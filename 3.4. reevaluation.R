#######################################################
#This file contains commands to re-evaluate designs found by different algorithms 
#######################################################

#-----------------------
#load the required packages and data
#-----------------------

library(igraph)
library(MASS)
library(Matrix)
library(parallel)
library(Rcpp)
sourceCpp("rcpp_help_funcs.cpp")
#if Rcpp does not work, uncomment and run the command below
#source("rcpp_help_funcs.R")

#------------------
#Re-evaluation function
#------------------

#' Results reevaluation
#' 
#' @description This function re-evaluate the designs found by different algorithms
#' @param g the graph object that contains relevant information for the calculation of the design criterion to be passed to function "mc_func"
#' @param obj the object containing design-finding results, e.g. output from function test_algo() or mass_assign() (these can be found in folder "results/designs")
#' @param mc_func the function to evaluate the design criteria
#' @param ... other parameters to be passed to the "mc_func" function
#' @return a new design object with design values re-evaluated
re_eval <- function(g, obj, mc_func, ...) {
  ret <- c() 
  n <- dim(obj$best_designs)[2]
  for (i in 1:n) {
    #evaluate the design according to a given design criteria
    tmp <- do.call(mc_func, args = list(g = g, 
                                        assignment = obj$best_designs[,i], 
                                        ...))
    ret <- c(ret, tmp$value)
  }
  obj$best_vals <- ret
  obj
}

#------------------
#Re-evaluating all designs
#------------------

#change this "model" parameter to either "basse" (NS), "bintemp" (BNTAR), "gane" (POW-DEG), or "locopt" (CNAR)
#to check the MC convergence for each model
model <- "gane" 

#change this "network" parameter to either "enron", "caltech", "umich"
network <- "enron"

#load the design criterion functions for the specified model
load(paste0("funcs/", model, "_funcs.RData"))
g <- get(paste0("g_", network, "_", model))

#load the design results
graph_based <- load(paste0("results/designs/", network, "_graph_based_uneval.RData"))
optimal_found <- load(paste0("results/designs/", network, "_", model, ".RData"))
list_names <- setdiff(c(graph_based, optimal_found), "name")
list_names #list of design objects' names
list_obj <- c(mget(list_names)) #list of design objects

#re-evaluating all designs in the list
for (i in 1:length(list_names)) {
  
  #re-evaluating the design over 50,000 parameter draws in MC eveluation to get an exact design criterion value
  eval(call("<-", as.name(list_names[i]), 
            re_eval(g, list_obj[[i]], mc_func = get(paste0("mc_", model2)), 
                    nsim = 50000, mc.cores = 8, trim = 0.05))) 
  
  cat(i, "/", length(list_names), "\n")
}

#save the re-evaluated results into a new object
save(list = list_names, 
     file = paste0("results/evaluation/", network, "_", model, "_reeval.RData"))
