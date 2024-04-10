#######################################################
#This file contain commands perform small test run and reproduce selected (fast) results
#######################################################

#-----------------------
#Choose the model and the network
#-----------------------

#Specify the model and network. We will use the 
model <- "bintemp" #available options are "locopt" (CNAR), "basse" (NS), "gane" (POW-DEG), and "bintemp" (BNTAR) models
network <- "enron" #other options are "caltech" and "umich". We choose the smallest, i.e., the "enron" network, for fastest results

#-----------------------
#load the required packages and data sets
#-----------------------

library(igraph)
library(Matrix)
library(MASS)
library(parallel)
library(Rcpp)
sourceCpp("rcpp_help_funcs.cpp")
#if Rcpp does not work, uncomment and run the command below
#source("rcpp_help_funcs.R")

library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggh4x)
library(chron)

load("funcs/mc_funcs.RData")
load("funcs/run_funcs.RData")
load(paste0("funcs/", model, "_funcs.RData"))
load("funcs/heuristics_funcs.RData")
load("funcs/neuralnet.RData")
load("funcs/cluster_funcs.RData")
load("funcs/bayesopt_funcs.RData")

#-----------------------
#Calculate the design criterion 
#-----------------------

#get the relevant graph objects and function
g <- get(paste0("g_", network, "_", model))
mc_func <- get(paste0("mc_", model))

#generate a design, you can generate your own design, too
set.seed(123456)
assignment <- assign_binomial(g$N) 

#calculate the design criterion according to the default prior
check <- mc_func(g, assignment = assignment, timed = TRUE)
check$value #design criterion value

#you can also specify your own prior, as long as it return (randomly) a suitable set of parameters for the model you specified
prior <- function(i) {
  c(0.5, 0.5)
}
check <- mc_func(g, assignment = assignment, 
                 timed = TRUE, prior = prior, 
                 nsim = 5000) 
check$value #design criterion value

#-----------------------
#Reproduce a part of Figure S1
#-----------------------

#get the relevant graph objects and function
g <- get(paste0("g_", network, "_", model))
mc_func <- get(paste0("mc_", model))

#check the convergence of the monte carlo approximation
prior <- get(paste0(model, "_prior_default")) #default prior
tmp <- mc_trace(mc.func = mc_func, 
                control.mc = list(g = g, 
                                  assignment = assignment, 
                                  prior = prior, 
                                  nsim = 50000)) #about 10 minutes for bintemp model and enron network
tmp <- as.data.frame(tmp)
colnames(tmp) <- 1:ncol(tmp)
tmp <- reshape(tmp, varying = 1:ncol(tmp), v.names = "value", 
               timevar = "run", times = colnames(tmp), 
               direction = "long", new.row.names = 1:(ncol(tmp)*nrow(tmp)))
colnames(tmp) <- c("run", "value", "iteration")
ggline(data = tmp, x = "iteration", y = "value", 
       ylab = "Design Criterion \n", xlab = "\n Iteration", 
       plot_type = "l", size = 0.1, color = "run") +
  scale_color_manual(values = alpha(rep("black", 10), alpha = 0.2)) +
  theme(legend.position="none") +
  geom_vline(xintercept = 5000, colour = "black", linetype = "longdash")

#-----------------------
#Generate graph cluster randomized designs
#-----------------------

#load the corresponding network
if (network == "enron") {
  g <- get("g_enron")
} else {
  g <- get(paste0("g_fb_", network))
}
#choose one of these algorithms to perform graph clustering
tmp <- louvain_wrapper(g)     #Louvain clustering
tmp <- reLDG(g)               #restreaming linear deterministic greedy clustering
tmp <- social_hash(g)         #social hash clustering
tmp <- balanced_label_prop(g) #balanced label propagation clustering
#assign treatment and control using clustering info
assignment <- assign_mixed(tmp, clustered = TRUE)

#if you want to run a graph clustering algorithm multiple times and save the results
#look at code file 2.4 to see the documentation of the function and adjust the functions' parameters as you wish
tmp <- mass_clustering(g, 
                       cluster_func = reLDG, name = "reLDG", #choose the clustering algo
                       nrun = 30, #number of runs
                       file = NULL) #takes about 5 minutes
View(tmp$clusters) #clusters found
mean(times(tmp$times)) #average running time

#-----------------------
#Finding optimal designs using meta-heuristic/Bayesian optimization algorithms
#-----------------------

#get the relevant graph objects and function
g <- get(paste0("g_", network, "_", model)) #get the network
mc_func <- get(paste0("mc_", model)) #design criterion function with default prior

#choose one of these algorithms to find optimal designs
#look at code files 2.1 or 2.3 to see the documentation of these functions and adjust the functions' parameters as you wish
#these take a long time to run, especially when mc.cores = 1. For test run, recommand running meta-heuristic algos only.
tmp <- random_search(g$N, crit.func = mc_func, 
                     control.crit = list(g = g, trim = 0.05, mc.cores = 1), #set mc.cores > 1 for Mac or Linux
                     ndesign = 5000) 
tmp <- tabu_search(g$N, crit.func = mc_func, 
                   control.crit = list(g = g, trim = 0.05, mc.cores = 1), #set mc.cores > 1 for Mac or Linux
                   ndesign = 5000)
tmp <- simulated_annealing(g$N, crit.func = mc_func, 
                           control.crit = list(g = g, trim = 0.05, mc.cores = 1), #set mc.cores > 1 for Mac or Linux
                           ndesign = 5000)
tmp <- genetic_algorithm(g$N, crit.func = mc_func, 
                         control.crit = list(g = g, trim = 0.05, mc.cores = 1), #set mc.cores > 1 for Mac or Linux
                         ndesign = 5000)
tmp <- tree_parzen(g$N, crit.func = mc_func, 
                   control.crit = list(g = g, trim = 0.05, mc.cores = 1), #set mc.cores > 1 for Mac or Linux
                   ndesign = 5000)
tmp <- bayesopt_surrogate(g$N, crit.func = mc_func, 
                          control.crit = list(g = g, trim = 0.05, mc.cores = 1), #set mc.cores > 1 for Mac or Linux
                          ndesign = 5000)
tmp <- bayesopt_surrogate_policy(g$N, crit.func = mc_func, 
                                 control.crit = list(g = g, trim = 0.05, mc.cores = 1), #set mc.cores > 1 for Mac or Linux
                                 ndesign = 5000)

#summaries of results found
tmp$best_design
tmp$best_val
tmp$time

#plot the design criteria of designs visited in each iteration
plot(1:length(tmp$trace_val), tmp$trace_val, type = "l", 
     xlab = "iteration", ylab = "design criterion", 
     main = "results by iteration") 

#-----------------------
#Reproduce some selected results in Figure 1 and Figure 2
#-----------------------

#Here, we get some results from the imbalanced cluster randomization for faster results
#Use can get results from other algorithms in similar fashion as codes provided in 3.2 and 3.3

#get the relevant graph objects and function
if (network == "enron") {
  g <- get("g_enron")
} else {
  g <- get(paste0("g_fb_", network))
}
mc_func <- get(paste0("mc_", model)) 

#generate set of designs produced by imbalanced cluster randomization
tmp <- mass_clustering(g, 
                       cluster_func = louvain_wrapper, 
                       name = "imbalanced cluster randomization", 
                       nrun = 30, #number of runs
                       file = NULL) 
mass_assign(tmp, assign_func = assign_mixed, 
            control.assign = list(clustered = TRUE))

#re-evaluate the design criterion values of these designs
g <- get(paste0("g_", network, "_", model)) #get the network
ret <- c()
for (i in 1:dim(tmp$best_designs)[2]) { #takes about 15 minutes
  ret <- c(ret, do.call(mc_func, args = list(g = g, 
                                             assignment = tmp$best_designs[,i], 
                                             nsim = 50000, trim = 0.05))$value) 
}
tmp$best_vals <- ret

#----------------------->AVERAGE DESIGN CRITERION VALUES
mean(tmp$best_vals)

#----------------------->DESIGN CHARACTERISTICS
a <- tmp$best_designs[,1] #choose one design from the run, here I choose the first deisgn

#load the corresponding network
if (network == "enron") {
  g <- get("g_enron")
} else {
  g <- get(paste0("g_fb_", network))
}

A <- as_adjacency_matrix(g)
deg <- degree(g) 
N <- length(V(g))
between <- betweenness(g) 
close <- closeness(g)

#percentage of treated nodes
per_trt <- sum(a)/N
per_trt
  
#average difference in degree of treated vs controlled nodes
mean(deg[a == 1], na.rm = TRUE) - mean(deg[a == 0], na.rm = TRUE)

#average difference in betweenness of treated vs controlled nodes
mean(between[a == 1], na.rm = TRUE) - mean(between[a == 0], na.rm = TRUE)

#average difference in closeness of treated vs controlled nodes
mean(close[a == 1], na.rm = TRUE) - mean(close[a == 0], na.rm = TRUE)

#percentage of higher than expected proportion of neighbors sharing same treatment assignment
per_nei_trt <- as.vector(A%*%as.matrix(a))/deg 
per_nei_trt[is.nan(per_nei_trt) | is.infinite(per_nei_trt)] <- 0
per_nei_ctl <- 1 - per_nei_trt 
per_nei_ctl[deg == 0] <- 0
(sum(per_nei_trt[a==1] > per_trt) + sum(per_nei_ctl[a==0] > (1-per_trt)))/(sum(deg != 0))
