#######################################################
#This file contains commands to investigate the different clustering algorithms
#######################################################

#-----------------------
#load the required packages and data
#-----------------------

library(igraph)
library(MASS)
library(Matrix)

library(ggplot2)
library(ggpubr)
library(dplyr)

load("funcs/mc_funcs.RData")
load("funcs/cluster_funcs.RData")

#-----------------------
#Cluster graphs using different algorithms
#-----------------------

#~~~~~~~~~LOUVAIN CLUSTERING

enron_louvain <- mass_clustering(g_enron, 
                                 cluster_func = louvain_wrapper, 
                                 name = "louvain", 
                                 name_obj = "enron_louvain", 
                                 file = "results/clustering_compare.RData")

caltech_louvain <- mass_clustering(g_fb_caltech, 
                                  cluster_func = louvain_wrapper, 
                                  name = "louvain", 
                                  name_obj = "caltech_louvain", 
                                  file = "results/clustering_compare.RData")

umich_louvain <- mass_clustering(g_fb_umich, 
                                 cluster_func = louvain_wrapper, 
                                 name = "louvain", 
                                 name_obj = "umich_louvain", 
                                 file = "results/clustering_compare.RData")

#~~~~~~~~~BALANCED LABEL PROPAGATION

enron_balanced_cluster <- mass_clustering(g = g_enron,
                                          cluster_func = balanced_label_prop,
                                          name = "balanced label propagation",
                                          name_obj = "enron_balanced_cluster",
                                          file = "results/clustering_compare.RData",
                                          verbose = FALSE)
cat("Enron Done!\n")

caltech_balanced_cluster <- mass_clustering(g = g_fb_caltech,
                                            cluster_func = balanced_label_prop,
                                            name = "balanced label propagation",
                                            name_obj = "caltech_balanced_cluster",
                                            file = "results/clustering_compare.RData",
                                            verbose = FALSE)
cat("Caltech Done!\n")

umich_balanced_cluster <- mass_clustering(g = g_fb_umich,
                                          cluster_func = balanced_label_prop,
                                          name = "balanced label propagation",
                                          name_obj = "umich_balanced_cluster",
                                          file = "results/clustering_compare.RData",
                                          verbose = FALSE)
cat("UMich Done!\n")

#~~~~~~~~~reLDG

enron_reLDG <- mass_clustering(g = g_enron,
                               cluster_func = reLDG,
                               name = "reLDG",
                               name_obj = "enron_reLDG",
                               file = "results/clustering_compare.RData",
                               verbose = FALSE)
cat("Enron Done!\n")

caltech_reLDG <- mass_clustering(g = g_fb_caltech,
                                 cluster_func = reLDG,
                                 name = "reLDG",
                                 name_obj = "caltech_reLDG",
                                 file = "results/clustering_compare.RData",
                                 verbose = FALSE)
cat("Caltech Done!\n")

umich_reLDG <- mass_clustering(g = g_fb_umich,
                               cluster_func = reLDG,
                               name = "reLDG",
                               name_obj = "umich_reLDG",
                               file = "results/clustering_compare.RData",
                               verbose = FALSE)
cat("UMich Done!\n")

#~~~~~~~~~SOCIAL HASH

enron_social_hash <- mass_clustering(g = g_enron,
                                     cluster_func = social_hash,
                                     name = "social hash",
                                     name_obj = "enron_social_hash",
                                     file = "results/clustering_compare.RData",
                                     verbose = FALSE)
cat("Enron Done!\n")

caltech_social_hash <- mass_clustering(g = g_fb_caltech,
                                       cluster_func = social_hash,
                                       name = "social hash",
                                       name_obj = "caltech_social_hash",
                                       file = "results/clustering_compare.RData",
                                       verbose = FALSE)
cat("Caltech Done!\n")

umich_social_hash <- mass_clustering(g = g_fb_umich,
                                     cluster_func = social_hash,
                                     name = "social hash",
                                     name_obj = "umich_social_hash",
                                     file = "results/clustering_compare.RData",
                                     verbose = FALSE)
cat("UMich Done!\n")

#------------------------
#Comparing results
#------------------------

library(dplyr)
library(igraph)
library(chron)

#' Clustering summary
#' 
#' @description This function summarizes the modularity, range of cluster sizes, and the running times for a clustering object in the "cluster_compare.RData" file
#' @param nm the name of the R object to be summarized
#' @return a list that contains the average modularity, range of cluster sizes, and the running times over all the runs in the clustering the object
cluster_summary <- function(nm) {
  
  #get the clustering object and corresponding graph object
  obj <- get(paste0(nm))
  nm <- strsplit(nm, "_")[[1]]
  if (nm[1] == "enron") {
    g <- get(paste0("g_", nm[1]))
  } else {
    g <- get(paste0("g_fb_", nm[1]))
  }
  
  #place holder for results
  result <- data.frame(modularity = numeric(), 
                       range = numeric(), 
                       times = character())
  
  #calculate the modularity, cluster size range and running time for each run in the clustering object
  cnt <- 0
  for (i in 1:length(obj$times)) {
    cnt <- cnt+1
    tmp <- obj$clusters[,i]
    result[cnt, "modularity"] <- modularity(g,tmp)
    sizes <- sapply(1:length(unique(tmp)), function(i) {sum(tmp == i)})
    result[cnt, "range"] <- max(sizes) - min(sizes)
    result[cnt, "times"] <- obj$times[i]
  }
  
  #return
  list(network = nm, 
       algo = ifelse(length(nm) == 3, paste0(nm[2:3], collapse = " "), nm[2]), 
       modularity = mean(result$modularity), 
       range = mean(result$range), 
       times = as.character(mean(times(result$times))))
}

#' Clustering summaries
#' 
#' @description This function summarizes the modularity, range of cluster sizes, and the running times for a list of clustering objects in the "cluster_compare.RData" file
#' @param list_obj the list of names of R clustering objects to be summarized
#' @return a data frame that contains the summaries for each algorithm and network
cluster_summaries <- function(list_obj) {
  
  #place holder for results
  result <- data.frame(network = character(), 
                       algo = character(),
                       modularity = numeric(), 
                       range = numeric(), 
                       times = character())
  
  #summarize each object in the list
  cnt <- 0
  for (i in 1:length(list_obj)) {
    cnt <- cnt+1
    result[cnt,] <- cluster_summary(list_obj[i])
  }
  
  #return
  result
}

#summarize all objects into the "clustering_compare.RData" file
clusters <- load("results/clustering_compare.RData") #YOUR results
# clusters <- load("saved results/clustering_compare.RData") #SAVED results
load("funcs/mc_funcs.RData")

#------------------> TABLE S1
cluster_summaries(clusters)

