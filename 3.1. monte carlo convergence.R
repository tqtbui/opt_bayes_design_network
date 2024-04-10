#######################################################
#This file contains commands to investigate the iteration at which the Monte Carlo approximation converge
#######################################################

#-----------------------
#load the required packages and data
#-----------------------

library(igraph)
library(MASS)
library(Matrix)
library(parallel)
library(Rcpp)

library(ggplot2)
library(ggpubr)
library(dplyr)

load("funcs/mc_funcs.RData")
sourceCpp("rcpp_help_funcs.cpp")
#if Rcpp does not work, uncomment and run the command below
#source("rcpp_help_funcs.R")

#------------------------
#check the MC approximation
#------------------------

#change this "model" parameter to either "basse" (NS), "bintemp" (BNTAR), "gane" (POW-DEG), or "locopt" (CNAR)
#to check the MC convergence for each model
model <- "basse" 

#load the design criterion functions for the specified model
load(paste0("funcs/", model, "_funcs.RData"))

#check MC convergence of the specified model on the Enron network
#using 10 MC runs, each for 50,000 iterations
g <- get(paste0("g_enron_", model))
set.seed(123456)
a <- assign_binomial(g$N, 0.5) #fix a design on the Enron network
eval(call("<-", as.name(paste0("enron_", model, "_trace")),
          mc_trace(mc.func = get(paste0("mc_", model)),
                   control.mc = list(g = g, assignment = a,
                                     nsim = 50000, mc.cores = 1), 
                   trim = 0.05)))

#check MC convergence of the specified model on the Caltech network similarly
g <- get(paste0("g_caltech_", model))
set.seed(123456)
a <- assign_binomial(g$N, 0.5) #fix a design on the Caltech network
eval(call("<-", as.name(paste0("caltech_", model, "_trace")),
          mc_trace(mc.func = get(paste0("mc_", model)),
                   control.mc = list(g = g, assignment = a,
                                     nsim = 50000, mc.cores = 1), 
                   trim = 0.05)))

#run lines 33-50 for each of the 4 models ("model" parameter needs to be specified in line 30)
#then save the results from all models into file "mc_trace.RData" in a folder called "results"
save(enron_basse_trace, caltech_basse_trace,
     enron_bintemp_trace, caltech_bintemp_trace,
     enron_locopt_trace, caltech_locopt_trace,
     enron_gane_trace, caltech_gane_trace,
     file = "results/mc_trace.RData")

#-------------------------
#visualizing the traces
#-------------------------

library(igraph)
library(ggplot2)
library(ggpubr)
library(dplyr)

#' Flat the trace results
#' 
#' @description This function get the results from the "mc_trace.RData" file and merge them into a single data frame
#' @param list_obj the list of "mc_trace" objects
#' @param list_names the names of the "mc_trace" objects
#' @return a data frame that contains the results of the "mc_trace" simulations
flat_trace <- function(list_obj, list_names) {
  
  nobj <- length(list_obj)
  
  #reshape
  ret <- c()
  for (i in 1:nobj) {
    
    #extract the objects
    obj <- list_obj[[i]]
    nm <- list_names[i]
    
    #reshape
    obj <- as.data.frame(obj)
    colnames(obj) <- 1:ncol(obj)
    tmp <- reshape(obj, varying = 1:ncol(obj), v.names = "value", 
                   timevar = "run", times = colnames(obj), 
                   direction = "long", new.row.names = 1:(ncol(obj)*nrow(obj)))
    nms <- strsplit(nm, split = "_")
    tmp$network <- nms[[1]][1]
    tmp$model <- nms[[1]][2]
    
    #rbind the data frames
    ret <- rbind(ret, tmp)
  }
  
  #return
  ret
}

#load saved data
list_names <- load("results/mc_trace.RData") #YOUR results
# list_names <- load("saved results/mc_trace.RData") #SAVED results
list_obj <- mget(list_names)
traces <- flat_trace(list_obj, list_names)
colnames(traces)[3] <- "iteration" #change the column name from id -> iteration

#order networks and model names
traces$network <- factor(traces$network, 
                         levels = c("enron", "caltech"), 
                         labels = c("Enron", "Caltech"),
                         ordered = TRUE)
traces$model_names <- factor(traces$model,
                             levels = c("locopt", "basse", "gane", "bintemp"),
                             labels = c("CNAR", "NS", "POW-DEG", "BNTAR"),
                             ordered = TRUE)

#setting y-axis limits for each pannel
limits <- data.frame(model_names = levels(traces$model_names),
                     min = c(0,0,0,0),
                     max = c(0.0015,5000,5,0.03))
dataC <- inner_join(traces, limits) %>% filter(value > min, value < max)
dataC$model_names <- factor(dataC$model_names,
                            levels = c("CNAR", "NS", "POW-DEG", "BNTAR"),
                            labels = c("CNAR", "NS", "POW-DEG", "BNTAR"),
                            ordered = TRUE)

#plot the results ----------------> FIGURE S1
#wait about 2-3 minutes for the plot to appear
ggline(data = dataC, x = "iteration", y = "value", 
       ylab = "Design Criterion \n", xlab = "\n Iteration", 
       plot_type = "l", facet.by = c("model_names", "network"), 
       size = 0.1, color = "run") +
  scale_alpha_manual(values = 0.01) +
  geom_vline(xintercept = 5000, colour = "black", linetype = "longdash") + #the reference line at 5000
  facet_grid(model_names ~ network, scales="free") +
  theme(legend.position="none") +
  scale_color_manual(values = alpha(rep("black", 10), alpha = 0.2)) +
  theme(strip.text = element_text(size = 14), 
        axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12)) 
ggsave(filename = paste0("results/plots/figure_s1.png"),
       device = "png", dpi = 300,
       width = 9, height = 12, bg = "white") #save the plots


