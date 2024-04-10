#######################################################
#This file contain functions for the basse model (i.e., the NS model in the paper)
#######################################################

library(Matrix)

#--------------------
#Make graph objects for the basse (NS) model
#--------------------

#These graph objects are useful to make the calculation of the Bayesian design criertion faster
load("funcs/gane_funcs.RData") #make use of the gane graph objects 

g_enron_basse <- list(N = g_enron_gane$N, #number of units
                      deg = g_enron_gane$deg, #degree, i.e., number of connections of each unit 
                      AA = t(g_enron_gane$A) %*% g_enron_gane$A) #A^T*A, where A is the adjacency matrix

g_caltech_basse <- list(N = g_caltech_gane$N,
                        deg = g_caltech_gane$deg,
                        AA = t(g_caltech_gane$A) %*% g_caltech_gane$A)

g_umich_basse <- list(N = g_umich_gane$N,
                      deg = g_umich_gane$deg,
                      AA = t(g_umich_gane$A) %*% g_umich_gane$A)

#--------------------
#calculate the design criteria
#--------------------

#' Calculate the design criterion for the basse (NS) model
#' 
#' @description This function calculates the design criterion for the basse (NS) model
#' @param g the graph object that contains relevant information for the calculation of the design criterion, for example, the g_enron_basse object
#' @param assignment the treatment assignment vector given to units on the graph
#' @param mu the mu parameter in the basse (NS) model
#' @param gamma the gamma parameter in the basse (NS) model
#' @param sigma the sigma parameter in the basse (NS) model
#' @return a numeric value of the design criterion
basse_crit <- function(g, assignment, 
                       mu = 1, gamma = 1, sigma = 2) {
  cpp_basse_crit(g$AA, g$deg, g$N, mu, gamma, sigma, assignment) #calculation is done in c++ for faster speed
}

#' Default prior for the basse (NS) model
#' 
#' @description This function randomly generates parameters for the basse (NS) model according to the default prior (as stated in the paper)
#' @param i any numeric value. This parameter is only used to comply with the Monte Carlo simulation function mc_basse()
#' @return a vector of values of mu, sigma, gamma drawn from the default prior distribution
basse_prior_default <- function(i) {
  mu <- rnorm(1, mean = 0, sd = 10)
  sigma <- rgamma(1, shape = 1, rate = 1) #inverse gamma
  sigma <- sigma^(-1/2)
  gamma <- rgamma(1, shape = 1, rate = 1) #inverse gamma
  gamma <- gamma^(-1/2)
  c(mu, sigma, gamma) 
}

#' Monte Carlo simulation to evaluate the Bayesian design criterion for the basse (NS) model
#' 
#' @description This function evaluates the Bayesian design criterion of a treatment assignment using Monte Carlo simulations for the basse (NS) model
#' @param g the graph object that contains relevant information for the calculation of the design criterion, for example, the g_enron_basse object
#' @param assignment the treatment assignment vector given to units on the graph
#' @param prior the prior function that return a vector of mu, sigma, gamma parameters. If NULL, the basse_prior_default() function will be used
#' @param nsim the number of parameter draws to use. Default is 5000
#' @param timed whether the time elapsed should be printed
#' @param mc.cores the number of cores to distribute the simulation. Default 1 for Windows (other values do not work for Windows). But in the simulation run, we run the code on 8 cores on a Linux server
#' @param seed the seed value to be set for the simulation. The default is NULL, corresponding to "no seed is given"
#' @param trim the trim when calculating the Monte Carlo estimate, which is the mean of the design criterion values with model parameters generated from the prior
#' @param na.rm whether to remove the NA value when calculating the Monte Carlo (mean) estimate
#' @return a list of two objects: trace: vector of design criterion values for each parameter draw; value: the final numeric Monte Carlo estimate of the Bayesian design criterion
mc_basse <- function(g, assignment, 
                     prior = NULL,  
                     nsim = 5000, timed = FALSE, 
                     mc.cores = 1, seed = NULL, 
                     trim = 0, na.rm = TRUE) {
  
  #changed assignment type
  assignment <- as.numeric(assignment)
  
  #default prior
  if (is.null(prior)) {
    prior <- basse_prior_default
  }

  #initialization
  ret <- 0
  if (!is.null(seed)) {
    set.seed(seed)
  }

  #monte carlo integration
  start <- Sys.time()

  ret <- mclapply(1:nsim, function(i) {
    pars <- prior(i)
    basse_crit(g, assignment = assignment, 
               mu = pars[1], gamma = pars[3], sigma = pars[2])
  }, mc.cores = mc.cores)
  
  end <- Sys.time()
  if (timed) {
    cat("Simulation ends in", hms_span(start, end), "\n")
  }
  
  #return
  if (is.list(ret)) {ret <- unlist(ret)}
  list(trace = ret, value = mean(ret, trim = trim, na.rm = na.rm))
}

#--------------------
#save the workspace
#--------------------

save(g_enron_basse, g_caltech_basse, g_umich_basse,
     basse_crit, basse_prior_default, mc_basse, 
     file = "funcs/basse_funcs.RData")

#######################################################
#Test if things work as expected
#######################################################

library(parallel)
library(Rcpp)
sourceCpp("rcpp_help_funcs.cpp")
#if Rcpp does not work, uncomment and run the command below
#source("rcpp_help_funcs.R")

load("funcs/mc_funcs.RData")
load("funcs/basse_funcs.RData")

set.seed(123456)
assignment <- assign_binomial(g_enron_basse$N)
check <- mc_basse(g_enron_basse, assignment = assignment, timed = TRUE)
check$value



