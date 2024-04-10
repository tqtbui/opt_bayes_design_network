#######################################################
#This file contain functions for the locopt model (i.e., the CNAR model in the paper)
#######################################################

library(igraph)
library(MASS)
library(Matrix)

#######################################################
#Make graph objects for the locopt (CNAR) model so that the design criterion calculation can be sped up 
#######################################################

load("funcs/gane_funcs.RData")

g_enron_locopt <- list(N = g_enron_gane$N,
                       A = g_enron_gane$A, 
                       D = diag(g_enron_gane$deg))

g_caltech_locopt <- list(N = g_caltech_gane$N,
                         A = g_caltech_gane$A, 
                         D = diag(g_caltech_gane$deg))

g_umich_locopt <- list(N = g_umich_gane$N,
                       A = g_umich_gane$A, 
                       D = diag(g_umich_gane$deg))

#######################################################
#Functions to calculate the design criterion for the locopt (CNAR) model
#######################################################

#' Calculate the design criterion for the locopt (CNAR) model
#' 
#' @description This function calculates the design criterion for the locopt (CNAR) model
#' @param g the graph object that contains relevant information for the calculation of the design criterion, for example, the g_enron_locopt object
#' @param assignment the treatment assignment vector given to units on the graph
#' @param rho the rho parameter in the locopt (CNAR) model
#' @return a numeric value of the design criterion
locopt_crit <- function(g, assignment, rho) {
  x <- cbind(rep(1, g$N), assignment)
  cpp_locopt_crit(g$D, g$A, rho, x)
}

#' Default prior for the locopt (CNAR) model
#' 
#' @description This function randomly generates parameters for the locopt (CNAR) model according to the default prior (as stated in the paper)
#' @param i any numeric value. This parameter is only used to comply with the Monte Carlo simulation function mc_locopt()
#' @return a vector of values of rho drawn from the default prior distribution
locopt_prior_default <- function(i) {
  runif(1, min = 0, max = 0.99) 
}

#' Monte Carlo simulation to evaluate the Bayesian design criterion for the locopt (CNAR) model
#' 
#' @description This function evaluates the Bayesian design criterion of a treatment assignment using Monte Carlo simulations for the locopt (CNAR) model
#' @param g the graph object that contains relevant information for the calculation of the design criterion, for example, the g_enron_locopt object
#' @param assignment the treatment assignment vector given to units on the graph
#' @param prior the prior function that return a value of rho parameter. If NULL, the locopt_prior_default() function will be used
#' @param nsim the number of parameter draws to use. Default is 5000.
#' @param timed whether the time elapsed should be printed
#' @param mc.cores the number of cores to distribute the simulation. Default 1 for Windows (other values do not work for Windows). But in the simulation run, we run the code on 8 cores on a Linux server
#' @param seed the seed value to be set for the simulation. The default is NULL, corresponding to "no seed is given"
#' @param trim the trim when calculating the Monte Carlo estimate, which is the mean of the design criterion values with model parameters generated from the prior
#' @param na.rm whether to remove the NA value when calculating the Monte Carlo (mean) estimate
#' @return a list of two objects: trace: vector of design criterion values for each parameter draw; value: the final numeric Monte Carlo estimate of the Bayesian design criterion
mc_locopt <- function(g, assignment, 
                      prior = NULL,  
                      nsim = 5000, timed = FALSE, 
                      mc.cores = 1, seed = NULL, 
                      trim = 0, na.rm = TRUE) {
  
  #default prior
  if (is.null(prior)) {
    prior <- locopt_prior_default
  }

  #initialization
  ret <- 0
  if (!is.null(seed)) {
    set.seed(seed)
  }

  #monte carlo integration
  start <- Sys.time()
  ret <- mclapply(1:nsim, function(i) {
    rho <- prior(i)
    locopt_crit(g, assignment = assignment, rho = rho)
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

save(g_enron_locopt, g_caltech_locopt, g_umich_locopt,
     locopt_crit, locopt_prior_default, mc_locopt, 
     file = "funcs/locopt_funcs.RData")

#######################################################
#Test if things work as expected
#######################################################

library(parallel)
library(Rcpp)
sourceCpp("rcpp_help_funcs.cpp")
#if Rcpp does not work, uncomment and run the command below
#source("rcpp_help_funcs.R")

load("funcs/mc_funcs.RData")
load("funcs/locopt_funcs.RData")

set.seed(123456)
assignment <- assign_binomial(g_enron_locopt$N)
check <- mc_locopt(g_enron_locopt, assignment = assignment, timed = TRUE)
check$value


