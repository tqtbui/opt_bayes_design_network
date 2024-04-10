#######################################################
#This file contain functions for the gane model (i.e., the POW-DEG model in the paper)
#######################################################

library(MASS)
library(Matrix)
library(igraph)

#######################################################
#Make graph objects for the gane (POW-DEG) model so that the design criterion calculation can be sped up 
#######################################################

#' Extract information from an igraph object to calculate the design criterion for the gane (POW-DEG) model
#' 
#' @description This function extracts necessary information from an igraph object to make a new graph object. This help speed up the design criterion calculation for the gane (POW-DEG) model
#' @param g an igraph object (from the mc_funcs.RData file)
#' @return a graph object that contains graph information needed to calculate the design criterion for the gane (POw-DEG) model
gane_extract_graph <- function(g) {
  
  #adjacency matrix
  A <- as_adjacency_matrix(g)
  
  #degree
  deg <- degree(g)

  #number of nodes
  N <- length(deg)
  
  #return a list of info needed to calculate the design criterion for the gane (POW-DEG) model
  list(A = A, deg = deg, N = N)
}

load("funcs/mc_funcs.RData")
g_enron_gane <- gane_extract_graph(g_enron)
g_caltech_gane <- gane_extract_graph(g_fb_caltech)
g_umich_gane <- gane_extract_graph(g_fb_umich) 

#######################################################
#Functions to calculate the design criterion for the gane (POW-DEG) model
#######################################################

#' Power transformation
#' 
#' @description This function conducts a element-wise power transformation on a vector
#' @param lambda the power coefficient
#' @param vec the vector 
#' @param der the derivative of the power transformation x^lambda with respect to lambda. 0 is no derivative, 1 is first derivative and 2 is second derivative
#' @return the power-transformed vector
power_mod <- function(lambda, vec, der = 0) {
  pow <- sign(vec)*(abs(vec)^lambda)
  pow[is.nan(pow)] <- 0
  pow[is.na(pow)] <- 0
  pow[is.infinite(pow)] <- 0
  if (der == 0) {
    ret <- pow
  } else if (der == 1) {
    ret <- pow*log(vec) #derivative
    ret[is.nan(ret)] <- 0
    ret[is.na(ret)] <- 0
    ret[is.infinite(ret)] <- 0
  } else if (der == 2) {
    ret <- pow*log(vec)*log(vec) #second derivative
    ret[is.nan(ret)] <- 0
    ret[is.na(ret)] <- 0
    ret[is.infinite(ret)] <- 0
  }
  ret
}

#' Model matrix
#' 
#' @description This function calculates the model matrix for the gane (POW-DEG) model
#' @param g the graph object that contains relevant information for the calculation of the design criterion, for example, the g_enron_gane object
#' @param assignment the treatment assignment vector given to units on the graph
#' @return a model matrix for the gane (POW-DEG) model
gane_modmat <- function(g, assignment) {

  #general info
  N <- g$N #number of nodes
  A <- g$A #adjacency matrix
  deg <- g$deg #degree
  num_trt <- as.vector(A%*%assignment) #number of treated neighbors
  num_ctr <- deg - num_trt #number of controlled neighbors
  
  #model matrix
  m <- cbind(rep(1, N), 
             assignment, 
             num_trt, num_ctr) #no need data frame format because no fitting is needed
  
  colnames(m) <- c("m", "t", "g1", "g2")

  #return
  list(m = m, N = N, deg = deg)
}

#' Calculate the design criterion for the gane (POW-DEG) model
#' 
#' @description This function calculates the design criterion for the gane (POW-DEG) model
#' @param obj the model matrix object created using the gane_modmat() function
#' @param beta the beta parameter in the gane (POW-DEG) model
#' @param lambda the lambda parameter in the gane (POW-DEG) model
#' @param sigma the sigma parameter in the gane (POW-DEG) model
#' @return a numeric value of the design criterion
gane_crit <- function(obj, beta = c(0,0,0,0), 
                      lambda = 1, sigma = 1) {
  
  #input parameters
  omega <- sigma^2
  trans.func <- power_mod
  formals(trans.func)$lambda <- lambda
  
  #preparing input matrices
  m <- obj$m
  m[,3] <- trans.func(vec = obj$m[,3])
  m[,4] <- trans.func(vec = obj$m[,4])
  
  #derivative matrix
  m.der <- cbind(rep(0, obj$N), 
                 rep(0, obj$N), 
                 trans.func(vec = obj$m[,3], der = 1), 
                 trans.func(vec = obj$m[,4], der = 1))
  
  netfunc.der <- mean(trans.func(vec = obj$deg, der = 1))
  netfunc <- mean(trans.func(vec = obj$deg))
  
  #calculation is done in c++ for faster speed
  ret <- cpp_pow_crit(m, m.der, as.matrix(beta), netfunc, netfunc.der, obj$N, omega) 
  
  #return
  ret
}

#' Default prior for the gane (POW-DEG) model
#' 
#' @description This function randomly generates parameters for the gane (POW-DEG) model according to the default prior (as stated in the paper)
#' @param i any numeric value. This parameter is only used to comply with the Monte Carlo simulation function mc_gane()
#' @return a vector of values of beta, lambda, sigma drawn from the default prior distribution
gane_prior_default <- function(i) {
  mu <- rnorm(1, mean = 0, sd = 10) #beta[1]
  tau <- rnorm(1, mean = 0, sd = 10) #beta[2]
  gamma1 <- rnorm(1, mean = 0, sd = 10) #beta[3]
  gamma2 <- rnorm(1, mean = 0, sd = 10) #beta[4]
  lambda <- runif(1, min = 0, max = 1) #beta distribution = uniform
  sigma <- rgamma(1, shape = 1, rate = 1) #inverse gamma
  sigma^(-1/2)
  c(mu, tau, gamma1, gamma2, lambda, sigma)
}

#' Monte Carlo simulation to evaluate the Bayesian design criterion for the gane (POW-DEG) model
#' 
#' @description This function evaluates the Bayesian design criterion of a treatment assignment using Monte Carlo simulations for the gane (POW-DEG) model
#' @param g the graph object that contains relevant information for the calculation of the design criterion, for example, the g_enron_gane object
#' @param assignment the treatment assignment vector given to units on the graph
#' @param prior the prior function that return a vector of beta, lambda, sigma parameters. If NULL, the gane_prior_default() function will be used
#' @param nsim the number of parameter draws to use. Default is 5000.
#' @param timed whether the time elapsed should be printed
#' @param mc.cores the number of cores to distribute the simulation. Default 1 for Windows (other values do not work for Windows). But in the simulation run, we run the code on 8 cores on a Linux server
#' @param seed the seed value to be set for the simulation. The default is NULL, corresponding to "no seed is given"
#' @param trim the trim when calculating the Monte Carlo estimate, which is the mean of the design criterion values with model parameters generated from the prior
#' @param na.rm whether to remove the NA value when calculating the Monte Carlo (mean) estimate
#' @return a list of two objects: trace: vector of design criterion values for each parameter draw; value: the final numeric Monte Carlo estimate of the Bayesian design criterion
mc_gane <- function(g, assignment, 
                    prior = NULL, 
                    nsim = 5000, timed = FALSE, 
                    mc.cores = 1, seed = NULL, 
                    trim = 0, na.rm = TRUE) {

  #default prior
  if (is.null(prior)) {
    prior <- gane_prior_default
  }

  #initialization
  ret <- 0
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  #calculate the model matrix outside of the monte carlo simulation to save time
  modmat <- gane_modmat(g = g, assignment = assignment)
  
  #monte carlo integration
  start <- Sys.time()
  ret <- mclapply(1:nsim, function(i) {
    pars <- prior(i)
    gane_crit(obj = modmat, 
              beta = pars[1:4], lambda = pars[5], sigma = pars[6])
  }, mc.cores = mc.cores)
  end <- Sys.time()
  if (timed) {
    cat("Simulation ends in", hms_span(start, end), "\n")
  }
  
  #return
  ret <- unlist(ret)
  list(trace = ret, value = mean(ret, trim = trim, na.rm = na.rm))
}

#~~~~~~~~~~~~~~~~~~~~~~
#Save the workspace
#~~~~~~~~~~~~~~~~~~~~~~~

save(g_enron_gane, g_caltech_gane, g_umich_gane,
     power_mod, gane_extract_graph, gane_modmat, gane_crit, 
     gane_prior_default, mc_gane, 
     file = "funcs/gane_funcs.RData")

#######################################################
#Test if things work as expected
#######################################################

library(parallel)
library(Rcpp)
sourceCpp("rcpp_help_funcs.cpp")
#if Rcpp does not work, uncomment and run the command below
#source("rcpp_help_funcs.R")

load("funcs/mc_funcs.RData")
load("funcs/gane_funcs.RData")

set.seed(123456)
assignment <- assign_binomial(g_enron_gane$N)
check <- mc_gane(g_enron_gane, assignment = assignment, timed = TRUE)
check$value


