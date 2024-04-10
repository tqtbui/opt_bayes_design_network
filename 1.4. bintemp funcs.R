#######################################################
#This file contain functions for the bintemp model (i.e., the BNTAR model in the paper)
#######################################################

library(igraph)
library(MASS)
library(Matrix)

#######################################################
#Make graph objects for the bintemp (BNTAR) model so that the design criterion calculation can be sped up 
#######################################################

load("funcs/gane_funcs.RData")

g_enron_bintemp <- list(N = g_enron_gane$N,
                        A = g_enron_gane$A, 
                        deg = g_enron_gane$deg, 
                        gte = 0)

g_caltech_bintemp <- list(N = g_caltech_gane$N,
                          A = g_caltech_gane$A, 
                          deg = g_caltech_gane$deg, 
                          gte = 0)

g_umich_bintemp <- list(N = g_umich_gane$N,
                        A = g_umich_gane$A, 
                        deg = g_umich_gane$deg, 
                        gte = 0)

#######################################################
#Functions to calculate the design criterion for the bintemp (BNTAR) model
#######################################################

#' Response of the bintemp (BNTAR) model
#' 
#' @description This function generates responses for the bintemp (BNTAR) model
#' @param g the graph object that contains relevant information for the calculation of the design criterion, for example, the g_enron_bintemp object
#' @param assignment the treatment assignment vector given to units on the graph
#' @param tau the tau parameter in the bintemp (BNTAR) model
#' @param gamma the gamma parameter in the bintemp (BNTAR) model
#' @param mu the mu parameter in the bintemp (BNTAR) model. Default is -1.5 as in Eckles et al. (2016)
#' @param sigma the sigma parameter in the bintemp (BNTAR) model. Default is 1. 
#' @param time the number of time steps T to simulate the response. Default is 3. 
#' @return a vector of binary responses generated from the bintemp (BNTAR) model
response_time <- function(g, assignment, 
                          tau, gamma, 
                          mu = -1.5, sigma = 1, time = 3) {
  
  response_old <- rep(0, g$N)
  response_new <- rep(0, g$N)

  for (t in 1:time) {
    er <- rnorm(g$N, mean = 0, sd = sigma) #random error
    
    net <- as.vector(g$A%*%response_old)
    net <- net/g$deg
    net[is.infinite(net)] <- 0
    net[is.nan(net)] <- 0
    
    response_new <- mu + tau*assignment + gamma*net + er
    response_new <- (response_new > 0)
    response_old <- response_new
  }
  
  #return
  as.vector(response_new)
}

#' Default prior for the bintemp (BNTAR) model
#' 
#' @description This function randomly generates parameters for the bintemp (BNTAR) model according to the default prior (as stated in the paper)
#' @param i any numeric value. This parameter is only used to comply with the Monte Carlo simulation function mc_bintemp()
#' @return a vector of values of tau and gamma values drawn from the default prior distribution
bintemp_prior_default <- function(i) {
  tau <- runif(1, min = 0, max = 1)
  gamma <- runif(1, min = 0, max = 1)
  c(tau, gamma)
}

#' Monte Carlo simulation to evaluate the Bayesian design criterion for the bintemp (BNTAR) model
#' 
#' @description This function evaluates the Bayesian design criterion of a treatment assignment using Monte Carlo simulations for the bintemp (BNTAR) model
#' @param g the graph object that contains relevant information for the calculation of the design criterion, for example, the g_enron_bintemp object
#' @param assignment the treatment assignment vector given to units on the graph
#' @param prior the prior function that return a vector of tau and gamma parameters. If NULL, the bintemp_prior_default() function will be used
#' @param nsim the number of parameter draws to use. Default is 5000.
#' @param timed whether the time elapsed should be printed
#' @param mc.cores the number of cores to distribute the simulation. Default 1 for Windows (other values do not work for Windows). But in the simulation run, we run the code on 8 cores on a Linux server
#' @param seed the seed value to be set for the simulation. The default is NULL, corresponding to "no seed is given"
#' @param trim the trim when calculating the Monte Carlo estimate, which is the mean of the design criterion values with model parameters generated from the prior
#' @param na.rm whether to remove the NA value when calculating the Monte Carlo (mean) estimate
#' @param mse whether to calculate the mse from the true gte value stored in g. Default is true. If false, it will just calculate the mean response over different parameter values.
#' @return a list of two objects: trace: vector of design criterion values for each parameter draw; value: the final numeric Monte Carlo estimate of the Bayesian design criterion
mc_bintemp <- function(g, assignment, 
                       prior = NULL,
                       nsim = 5000, timed = FALSE,
                       mc.cores = 1, seed = NULL, 
                       trim = 0, na.rm = TRUE, 
                       mse = TRUE) {

  #default prior
  if (is.null(prior)) {
    prior <- bintemp_prior_default
  }

  #how to process the DiM estimator
  if (mse) {
    if ((sum(assignment == 0) == 0) | (sum(assignment == 1) == 0)) {
      gte_func <- function(response) {
        mean(response) #if all nodes are assigned to treatment/control, just calculate the mean
      }
    } else {
      gte_func <- function(response) {
        mean(response[assignment == 1]) - mean(response[assignment == 0]) #the difference in means estimator
      }
    }
    output_func <- function(gte) {
      (gte - g$gte)^2 
    }
    finalize_func <- function(vec) {
      sqrt(mean(vec, trim = trim, na.rm = na.rm)) #the root mean squared error
    }
  } else {
    gte_func <- function(response) {
      mean(response) #just use the mean response
    }
    output_func <- function(gte) {
      gte
    }
    finalize_func <- function(vec) {
      mean(vec, trim = trim, na.rm = na.rm) #mean over the different parameter values
    }
  }

  #initialization
  ret <- 0
  if (!is.null(seed)) {
    set.seed(seed)
  }

  #monte carlo integration
  start <- Sys.time()
  ret <- mclapply(1:nsim, function(i) {
    tmp <- prior(i)
    response <- response_time(g = g, assignment = assignment,
                              tau = tmp[1], gamma = tmp[2])
    gte <- gte_func(response)
    output_func(gte)
  }, mc.cores = mc.cores)
  end <- Sys.time()
  if (timed) {
    cat("Simulation ends in", hms_span(start, end), "\n")
  }

  #return
  ret <- unlist(ret)
  list(trace = ret, value = finalize_func(ret))
}

#######################################################
#Calculate the true GTE (averaged over prior of parameters) for the bintemp (BNTAR) model for each of the network and store them in the graph objects
#######################################################

#' Calculate the true GTE of the bintemp (BNTAR) model
#' 
#' @description This function calculates the true GTE (averaged over prior of parameters) for the bintemp (BNTAR) model given a network
#' @param g the graph object that contains relevant information for the calculation of the design criterion, for example, the g_enron_bintemp object
#' @param nsim the number of runs to generate parameters and response values. Default is 50000
#' @param seed the seed value to be set for the simulation. The default is 123456.
#' @return a numeric value of the true GTE
bintemp_true_gte <- function(g, nsim = 50000, 
                             seed = 123456) {
  a1 <- rep(1, g$N) #when everyone is treated
  a0 <- rep(0, g$N) #when everyone is controlled
  
  trace_1 <- mc_trace(mc.func = mc_bintemp, 
                      control.mc = list(g = g, assignment = a1, 
                                        nsim = 50000,  
                                        mse = FALSE)) #mean response when everyone is treated
  
  trace_0 <- mc_trace(mc.func = mc_bintemp, 
                      control.mc = list(g = g, assignment = a0, 
                                        nsim = 50000,  
                                        mse = FALSE)) #mean response when everyone is controlled
  
  mean(trace_1[50000,]) - mean(trace_0[50000,]) #GTE is the difference in mean response in the two cases
}

load("funcs/mc_funcs.RData")
g_enron_bintemp$gte <- bintemp_true_gte(g_enron_bintemp) #takes about 20 minutes
g_caltech_bintemp$gte <- bintemp_true_gte(g_caltech_bintemp) #takes about 80 minutes
g_umich_bintemp$gte <- bintemp_true_gte(g_umich_bintemp) #takes about 5 hours

#~~~~~~~~~~~~~~~~~~~~~~
#Save the workspace
#~~~~~~~~~~~~~~~~~~~~~~~

save(g_enron_bintemp, g_caltech_bintemp, g_umich_bintemp, 
     response_time, bintemp_prior_default, 
     mc_bintemp, file = "funcs/bintemp_funcs.RData")

#######################################################
#Test if things work as expected
#######################################################

library(parallel)
library(Rcpp)
sourceCpp("rcpp_help_funcs.cpp")
#if Rcpp does not work, uncomment and run the command below
#source("rcpp_help_funcs.R")

load("funcs/mc_funcs.RData")
load("funcs/bintemp_funcs.RData")

set.seed(123456)
assignment <- assign_binomial(g_enron_bintemp$N)
check <- mc_bintemp(g_enron_bintemp, assignment = assignment, timed = TRUE)
check$value


