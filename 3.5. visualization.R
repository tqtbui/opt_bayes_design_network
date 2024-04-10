#######################################################
#This file contains commands to visualize results found by different algorithms 
#######################################################

#-----------------------
#load the required packages and data
#-----------------------

library(igraph)
library(MASS)
library(Matrix)
library(chron)

library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggh4x)

#######################################################
#Running times
#######################################################

#-----------------------
#Function
#-----------------------

#' Summarize running times
#' 
#' @description This function summarizes the running times of design objects (for example, objects returned by test_algo() or the ones stored in results/evaluation)
#' @param list_obj the list of design objects
#' @return a data frame that contains the algorithm and running times 
designs_times <- function(list_obj) {
  
  nobj <- length(list_obj)
  ret <- NULL
  cnames <- c()
  for (i in 1:nobj) {
    #get the object
    obj <- list_obj[[i]]
    
    #get the algorithm
    algo <- gsub("_", " ", obj$algo)
    
    #calculate the running times as the average over 30 runs
    if (is.null(ret)) {
      ret <- data.frame(algorithm = algo, 
                        time = as.character(mean(times(obj$times))))
    } else {
      ret <- rbind.data.frame(ret, c(algo, as.character(mean(times(obj$times)))))
    }
  }
  
  #return
  ret
}

#-----------------------
#Running times - TABLE S2
#-----------------------

#change this "model" parameter to either "basse" (NS), "bintemp" (BNTAR), "gane" (POW-DEG), or "locopt" (CNAR)
#to check the MC convergence for each model
model <- "basse" 

#YOUR RESULTS
enron_reeval <- load(paste0("results/evaluation/enron_", model, "_reeval.RData"))
caltech_reeval <- load(paste0("results/evaluation/caltech_", model, "_reeval.RData"))

#SAVED RESULTS
# enron_reeval <- load(paste0("saved results/evaluation/enron_", model, "_reeval.RData"))
# caltech_reeval <- load(paste0("saved results/evaluation/caltech_", model, "_reeval.RData"))

enron_best_times <- designs_times(list_obj = mget(enron_reeval))
caltech_best_times <- designs_times(list_obj = mget(caltech_reeval))

colnames(enron_best_times) <- c("algorithm", "enron")
colnames(caltech_best_times) <- c("algorithm", "caltech")

best_times <- left_join(enron_best_times, caltech_best_times, by = "algorithm")
View(best_times)

#######################################################
#Summarizing design characteristics and design values
#######################################################

#-----------------------
#Summarizing functions
#-----------------------

#' Summarize a graph
#' 
#' @description This function summarizes the node characteristics of a graph
#' @param g an igraph object
#' @return a list containing the adjacency matrix, number of nodes, the degree, betweenness and closeness vectors
graph_summary <- function(g) {
  list(A = as_adjacency_matrix(g), N = length(V(g)), 
       deg = degree(g), between = betweenness(g), 
       close = closeness(g))
}

#' Summarize a design
#' 
#' @description This function summarizes the characteristics of a design
#' @param g a graph summary object returned by function graph_summary()
#' @param a a design
#' @return a list of design characteristics
design_summary <- function(g, a) {
  
  deg <- g$deg
  
  #percentage of treated nodes
  per_trt <- sum(a)/g$N
  
  #average difference in degree of treated vs controlled nodes
  deg_diff <- mean(deg[a == 1], na.rm = TRUE) - mean(deg[a == 0], na.rm = TRUE)
  
  #average difference in betweenness of treated vs controlled nodes
  bet_diff <- mean(g$between[a == 1], na.rm = TRUE) - mean(g$between[a == 0], na.rm = TRUE)
  
  #average difference in closeness of treated vs controlled nodes
  clo_diff <- mean(g$close[a == 1], na.rm = TRUE) - mean(g$close[a == 0], na.rm = TRUE)
  
  #percentage of higher than expected proportion of neighbors sharing same treatment assignment
  per_nei_trt <- as.vector(g$A%*%as.matrix(a))/deg
  per_nei_trt[is.nan(per_nei_trt) | is.infinite(per_nei_trt)] <- 0
  per_nei_ctl <- 1 - per_nei_trt
  per_nei_ctl[deg == 0] <- 0
  per_nei <- (sum(per_nei_trt[a==1] > per_trt) + sum(per_nei_ctl[a==0] > (1-per_trt)))/(sum(deg != 0))
  
  #return
  list(per_trt = per_trt,
       deg_diff = deg_diff,
       per_nei = per_nei, 
       bet_diff = bet_diff, 
       clo_diff = clo_diff)
}

#' Summarize many designs
#' 
#' @description This function summarizes the characteristics of designs contained in multiple design objects
#' @param g a graph summary object returned by function graph_summary()
#' @param list_obj a list of design objects, i.e., objects that contain designs returned from test_algo() or re_eval() functions (for example objects stored in folder results/evaluation)
#' @return a data frame that contain all information about the design objects together with corresponding design characteristics
designs_summaries <- function(g, list_obj) {
  
  #initialize a place holder
  nobj <- length(list_obj)
  ret <- data.frame(id = numeric(), #design number
                    algo = character(), #algorithm that produces the design
                    criteria = numeric(), #design criterion value
                    time = character(), #running time
                    chars = character(), #type of design characteristics
                    chars_value = numeric()) #value of design characteristics
  
  cnt <- 0
  for (i in 1:nobj) {
    #get the object
    obj <- list_obj[[i]]
    
    #algorithm
    algo <- gsub("_", " ", obj$algo)

    #get the summaries of designs inside each object
    for (j in 1:length(obj$best_vals)) {
      
      id <- j #design number
      crit <- obj$best_vals[j] #design criterion value
      time <- obj$times[j] #running times
      tmp <- design_summary(g, obj$best_designs[,j]) #design characteristics
      
      #add lines to the data frame
      cnt <- cnt+1
      ret[cnt,] <- c(id = id, 
                     algo = algo, 
                     criteria = crit, 
                     time = time, 
                     chars = "per_trt", 
                     chars_value = tmp$per_trt) #percentage of treated nodes
      
      cnt <- cnt+1
      ret[cnt,] <- c(id = id, 
                     algo = algo, 
                     criteria = crit, 
                     time = time, 
                     chars = "deg_diff", 
                     chars_value = tmp$deg_diff) #average difference in degree of treated vs controlled nodes
      
      cnt <- cnt+1
      ret[cnt,] <- c(id = id, 
                     algo = algo, 
                     criteria = crit, 
                     time = time, 
                     chars = "bet_diff", 
                     chars_value = tmp$bet_diff) #average difference in betweenness of treated vs controlled nodes
      
      cnt <- cnt+1
      ret[cnt,] <- c(id = id, 
                     algo = algo, 
                     criteria = crit, 
                     time = time, 
                     chars = "clo_diff", 
                     chars_value = tmp$clo_diff) #average difference in closeness of treated vs controlled nodes
      
      cnt <- cnt+1
      ret[cnt,] <- c(id = id, 
                     algo = algo, 
                     criteria = crit, 
                     time = time, 
                     chars = "per_nei", 
                     chars_value = tmp$per_nei) #percentage of higher than expected proportion of neighbors sharing same treatment assignment
    }
  }
  
  #return the data frame
  ret
}

#-----------------------
#Summarizing design objects
#-----------------------

#summarize the graphs
load("funcs/mc_funcs.RData")
g_enron_summary <- graph_summary(g_enron)
g_caltech_summary <- graph_summary(g_fb_caltech)
g_umich_summary <- graph_summary(g_fb_umich)

#change this "model" parameter to either "basse" (NS), "bintemp" (BNTAR), "gane" (POW-DEG), or "locopt" (CNAR)
#to check the MC convergence for each model
model <- "gane"

#get the design objects
enron_reeval <- load(paste0("results/evaluation/enron_", model, "_reeval.RData"))
caltech_reeval <- load(paste0("results/evaluation/caltech_", model, "_reeval.RData"))
umich_reeval <- load(paste0("results/evaluation/umich_", model, "_reeval.RData"))

#summarize the design objects
enron_designs <- designs_summaries(g_enron_summary, list_obj = mget(enron_reeval))
caltech_designs <- designs_summaries(g_caltech_summary, list_obj = mget(caltech_reeval))
umich_designs <- designs_summaries(g_umich_summary, list_obj = mget(umich_reeval))

#merge enron and caltech datasets
enron_designs$network <- "enron"
caltech_designs$network <- "caltech"
umich_designs$network <- "umich"
best_designs <- rbind(enron_designs, caltech_designs, umich_designs)

#order networks
best_designs$network <- factor(best_designs$network,
                               levels = c("enron", "caltech", "umich"),
                               ordered = TRUE)

#order algorithms
best_designs$algo[best_designs$algo == "deep surrogate"] <- "surrogate-based local search"
best_designs$algo[best_designs$algo == "deep RL"] <- "reinforcement learning"
best_designs$algo <- factor(best_designs$algo,
                            levels = c("random search", "tabu search", "simulated annealing", "genetic algorithm",
                                       "surrogate-based local search", "reinforcement learning", "tree-parzen",
                                       "balanced randomization", "balanced cluster randomization", "imbalanced cluster randomization"),
                            ordered = TRUE)

#types of algorithms
unique(best_designs$algo)
best_designs$type <- "bayesian-optimization"
best_designs$type[best_designs$algo %in% c("random search", "tabu search", "simulated annealing", "genetic algorithm")] <- "meta-heuristics"
best_designs$type[best_designs$algo %in% c("balanced randomization", "balanced cluster randomization", "imbalanced cluster randomization")] <- "graph-based"
best_designs$type <- factor(best_designs$type,
                            levels = c("meta-heuristics", "bayesian-optimization", "graph-based"),
                            ordered = TRUE)

#adjust types of variables
best_designs$criteria <- as.numeric(best_designs$criteria)
best_designs$chars_value <- as.numeric(best_designs$chars_value) 

#calculate eficiency
best_designs <- best_designs %>%
  group_by(algo, network) %>%
  mutate(mean_crit = mean(criteria, trim = 0.1)) %>%
  ungroup() %>%
  group_by(network) %>%
  mutate(ref = first(mean_crit[algo == "balanced randomization"])) %>%
  ungroup() %>%
  mutate(efficiency = (1-criteria/ref)*100) %>%
  select(-mean_crit, -ref)

#give names to characteristics
measure.vec <- c("per_trt",
                 "deg_diff",
                 "bet_diff", 
                 "clo_diff",
                 "per_nei")
measure.names <- c("Percentage of treated units",
                   "Difference in average degree \n between treated vs controlled units",
                   "Difference in average betweenness \n between treated vs controlled units",
                   "Difference in average closeness \n between treated vs controlled units",
                   "Percentage of units having \n a higher proportion of similarly assigned \n neighbors than expected")
best_designs$chars_names <- factor(best_designs$chars,
                                   levels = measure.vec,
                                   labels = measure.names,
                                   ordered = TRUE)

#reference lines for each design characteristics
best_designs$ref <- rep(0.5, nrow(best_designs))
best_designs$ref[best_designs$chars %in% c("deg_diff", "bet_diff", "clo_diff")] <- 0

#save the data frame into an object with names containing the models
eval(call("<-", as.name(paste0(model, "_best_designs")), best_designs))

#-----------------------
#Save results of all models into a .RData file
#-----------------------

save(basse_best_designs, bintemp_best_designs, gane_best_designs, locopt_best_designs,
     g_enron_summary, g_caltech_summary, g_umich_summary,
     file = "results/designs_summaries.RData")

#######################################################
#Plotting
#######################################################

library(igraph)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggh4x)

#load the design summaries
load(paste0("results/designs_summaries.RData")) #YOUR results
# load(paste0("saved results/designs_summaries.RData")) #SAVED results

#give model name to the design data frames
basse_best_designs$model <- "basse"
bintemp_best_designs$model <- "bintemp"
gane_best_designs$model <- "gane"
locopt_best_designs$model <- "locopt"

#merge the different design data frame into one
best_designs <- rbind.data.frame(basse_best_designs, bintemp_best_designs, 
                                 gane_best_designs, locopt_best_designs)
best_designs$model_names <- factor(best_designs$model,
                                   levels = c("locopt", "basse", "gane", "bintemp"),
                                   labels = c("CNAR", "NS", "POW-DEG", "BNTAR"),
                                   ordered = TRUE)
best_designs$network <- factor(best_designs$network, 
                               levels = c("enron", "caltech", "umich"), 
                               labels = c("Enron", "Caltech", "UMichigan"), 
                               ordered = TRUE)

#----------Performances of algorithms - FIGURE 1

obj_vals <- best_designs %>%
  distinct(algo, network, id, model, .keep_all = TRUE)

ggboxplot(data = obj_vals, x = "algo", y = "efficiency", 
          ylab = "\n Efficiency (%)", xlab = "Algorithm \n", 
          orientation = "horizontal",
          facet.by = c("model_names", "network")) +
  geom_hline(yintercept = 0, colour = "black", linetype = "longdash") +
  facet_grid(model_names ~ network, scales = "fixed") +
  coord_flip(ylim = c(-100, 100)) +
  theme(strip.text = element_text(size = 14), 
        axis.title = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12)) 
ggsave(filename = paste0("results/plots/figure_1.png"),
       device = "png", dpi = 300,
       width = 11, height = 12, bg = "white")

#----------Design characteristics - FIGURE 2

obj <- best_designs %>%
  filter(chars %in% c("deg_diff", "bet_diff", "per_nei"))

ggscatter(data = obj, 
          x = "chars_value", y = "efficiency", 
          xlab = "\n Value", ylab = "Efficiency (%) \n", 
          position = position_jitter(0.1), alpha = 0.1, 
          ylim = c(-100, 100)) + #lancet palette
  facet_nested(model_names + network ~ chars_names, scales = "free_x") +
  geom_vline(data = obj, aes(xintercept = ref), colour = "black", 
             linetype = "longdash") +
  theme(strip.text.y = element_text(size = 14), 
        strip.text.x = element_text(size = 11), 
        axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
ggsave(filename = paste0("results/plots/figure_2.png"),
       device = "png", dpi = 300,
       width = 11, height = 14, bg = "white")

#----------Other design characteristics - FIGURE S2
obj <- best_designs %>%
  filter(!chars %in% c("deg_diff", "bet_diff", "per_nei"))

ggscatter(data = obj, 
          x = "chars_value", y = "efficiency", 
          xlab = "\n Value", ylab = "Efficiency (%) \n", 
          position = position_jitter(0.1), alpha = 0.1) + 
  facet_nested(model_names + network ~ chars_names, scale = "free") +
  geom_vline(data = obj, aes(xintercept = ref), colour = "black", 
             linetype = "longdash") +
  theme(strip.text.y = element_text(size = 14), 
        strip.text.x = element_text(size = 11), 
        axis.title = element_text(size = 14, face = "bold"), 
        axis.text = element_text(size = 12), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
ggsave(filename = paste0("results/plots/figure_s2.png"),
       device = "png", dpi = 300,
       width = 9, height = 14, bg = "white")


