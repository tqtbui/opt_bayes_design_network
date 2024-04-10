#######################################################
#This file contains commands to find designs by different algorithms multiple times 
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

#######################################################
#Graph based designs
#######################################################

#load the graph clustering functions and results from code file 3.2
load("results/clustering_compare.RData")
load("funcs/cluster_funcs.RData")

#------------------------------
#Assign treatment and control for graph clusterings
#------------------------------

#graph-cluster randomization
mass_assign(enron_louvain, enron_balanced_cluster, 
            caltech_louvain, caltech_balanced_cluster,  
            umich_louvain, umich_balanced_cluster,
            assign_func = assign_mixed, 
            control.assign = list(clustered = TRUE))

#balanced design (randomly assign half to treatment and half to control)
enron_balanced_unit <- enron_louvain #copy the louvain objects
caltech_balanced_unit <- caltech_louvain
umich_balanced_unit <- umich_louvain
mass_assign(enron_balanced_unit, caltech_balanced_unit, #assign
            umich_balanced_unit,
            assign_func = assign_mixed, 
            control.assign = list(clustered = FALSE))
enron_balanced_unit$clusters <- NULL #delete the clustering info
caltech_balanced_unit$clusters <- NULL
umich_balanced_unit$clusters <- NULL

#rename the algorithms for each object

enron_balanced_unit$algo <- "balanced randomization"
caltech_balanced_unit$algo <- "balanced randomization"
umich_balanced_unit$algo <- "balanced randomization"

enron_balanced_cluster$algo <- "balanced cluster randomization"
caltech_balanced_cluster$algo <- "balanced cluster randomization"
umich_balanced_cluster$algo <- "balanced cluster randomization"

enron_louvain$algo <- "imbalanced cluster randomization"
caltech_louvain$algo <- "imbalanced cluster randomization"
umich_louvain$algo <- "imbalanced cluster randomization"

#save the results
save(enron_balanced_unit, enron_balanced_cluster, enron_louvain, 
     file = "results/designs/enron_graph_based_uneval.RData")

save(caltech_balanced_unit, caltech_balanced_cluster, caltech_louvain, 
     file = "results/designs/caltech_graph_based_uneval.RData")

save(umich_balanced_unit, umich_balanced_cluster, umich_louvain, 
     file = "results/designs/umich_graph_based_uneval.RData")

#######################################################
#Finding optimal designs
#######################################################

#------------------------
#Load the required functions
#-----------------------

load("funcs/run_funcs.RData")
load("funcs/heuristics_funcs.RData")
load("funcs/neuralnet.RData")
load("funcs/bayesopt_funcs.RData")

#------------------------
#Load the functions for a specific model
#-----------------------

#change this "model" parameter to either "basse" (NS), "bintemp" (BNTAR), "gane" (POW-DEG), or "locopt" (CNAR)
#to check the MC convergence for each model
model <- "gane" 

#load the design criterion functions for the specified model
load(paste0("funcs/", model, "_funcs.RData"))

#------------------
#Run the design-finding algorithms for a specific network
#------------------

#change this "network" parameter to either "enron", "caltech", "umich"
network <- "enron"
g <- get(paste0("g_", network, "_", model)) #get the network

#NOTE: We only run genetic algorithm and tree-parzen estimator for umich network due to infeasibly long running time

#random search
test_algo(nrun = 30, name = "random search",
          name_obj = paste0(network, "_rs"), 
          file = paste0("results/designs/", network, "_", model, ".RData"),
          crit.func = get(paste0("mc_", model)),
          control.crit = list(g = g, trim = 0.05, mc.cores = 8),
          algo.func = random_search,
          control.algo = list(N = g$N)) #using default hyperparameters
cat("Done\n")

#tabu search
test_algo(nrun = 30, name = "tabu search",
          name_obj = paste0(network, "_ts"), 
          file = paste0("results/designs/", network, "_", model, ".RData"),
          crit.func = get(paste0("mc_", model)),
          control.crit = list(g = g, trim = 0.05, mc.cores = 8),
          algo.func = tabu_search,
          control.algo = list(N = g$N)) #using default hyperparameters
cat("Done\n")

#simulated annealing
test_algo(nrun = 30, name = "simulated annealing",
          name_obj = paste0(network, "_sa"), 
          file = paste0("results/designs/", network, "_", model, ".RData"),
          crit.func = get(paste0("mc_", model)),
          control.crit = list(g = g, trim = 0.05, mc.cores = 8),
          algo.func = simulated_annealing,
          control.algo = list(N = g$N)) #using default hyperparameters

#genetic algorithm
test_algo(nrun = 30, name = "genetic algorithm",
          name_obj = paste0(network, "_ga"), 
          file = paste0("results/designs/", network, "_", model, ".RData"),
          crit.func = get(paste0("mc_", model)),
          control.crit = list(g = g, trim = 0.05, mc.cores = 8),
          algo.func = genetic_algorithm,
          control.algo = list(N = g$N)) #using default hyperparameters

#tree-parzen estimator
test_algo(nrun = 30, name = "tree-parzen",
          name_obj = paste0(network, "_tree"), 
          file = paste0("results/designs/", network, "_", model, ".RData"),
          crit.func = get(paste0("mc_", model)),
          control.crit = list(g = g, trim = 0.05, mc.cores = 8),
          algo.func = tree_parzen,
          control.algo = list(N = g$N)) #using default hyperparameters

#local search using a surrogate model
test_algo(nrun = 30, name = "deep surrogate",
          name_obj = paste0(network, "_surrogate"), 
          file = paste0("results/designs/", network, "_", model, ".RData"),
          crit.func = get(paste0("mc_", model)),
          control.crit = list(g = g, trim = 0.05, mc.cores = 8),
          algo.func = bayesopt_surrogate,
          control.algo = list(N = g$N)) #using default hyperparameters

#bayesian optimization with reinforcement learning
test_algo(nrun = 30, name = "deep RL",
          name_obj = paste0(network, "_policy"), 
          file = paste0("results/designs/", network, "_", model, ".RData"),
          crit.func = get(paste0("mc_", model)),
          control.crit = list(g = g, trim = 0.05, mc.cores = 8),
          algo.func = bayesopt_surrogate_policy,
          control.algo = list(N = g$N)) #using default hyperparameters
