#######################################################
#This file contains functions for graph-cluster randomization
#######################################################

#-----------------------------
#Help functions
#-----------------------------

#' Resave an object into a file
#' 
#' @description This function save an R object into an existing .RData file
#' @param ... the objects to be saved
#' @param list (optional) list of objects to be saved
#' @param file the name of the .RData file
resave <- function(..., list = character(), file) {
  
  #load the file
  path <- paste0(getwd(),"/")
  previous  <- suppressWarnings(try(load(paste0(path, file)), silent = TRUE))
  
  #get the objects
  var.names <- c(list, as.character(substitute(list(...)))[-1L])
  for (var in var.names) assign(var, get(var, envir = parent.frame()))
  
  #save the the .RData file
  if ("try-error" %in% class(previous)) {
    save(list, file = file)
  } else {
    save(list = unique(c(previous, var.names)), file = file)
  }
}

#-----------------------------
#Graph clustering functions
#-----------------------------

#' Louvain clustering
#' 
#' @description This function is a wrapper for the built-in Louvain clustering function in igraph. It extracts the cluster label from the object returned by the built in function
#' @param g an igraph object
#' @return the cluster labels
louvain_wrapper <- function(g) {
  temp <- cluster_louvain(g)
  temp$membership
}

#' Balanced label propagation
#' 
#' @description This function performs balanced label propagation for balanced clustering on graphs
#' @param g an igraph object
#' @param ncom number of cluster (communities). Default is NULL, corresponding to the number of clusters returned by the Louvain algorithm
#' @param max_step maximum number of iterations. Default is 100
#' @param shuffle proportion of random cluster label shuffling in each iteration. Default is 0.05
#' @param tol threshold for modularity convergence
#' @param verbose whether the run progress should be printed
#' @return a vector of cluster labels
balanced_label_prop <- function(g, ncom = NULL, max_step = 100, 
                                shuffle = 0.05, tol = 1e-03, 
                                verbose = TRUE) {
  
  #initialization
  if (is.null(ncom)) {
    ncom <- length(unique(cluster_louvain(g)$membership)) #default number of clusters are decided by the louvain algorithm 
  }
  n <- length(V(g))
  
  #initialize cluster labels with balanced numbers of members in each cluster
  if (n%%ncom == 0) {
    member <- c(rep(1:ncom, n%/%ncom))
  } else {
    member <- c(rep(1:ncom, n%/%ncom), 1:(n%%ncom))
  }
  
  #membership and couples
  member_old <- member #old membership (cluster label)
  couples <- combn(1:n, 2) #pairs of nodes
  nc <- dim(couples)[2] #number of pairs of nodes
  shuffle <- round(shuffle*nc) #number of shuffling pairs
  
  #algorithm: alternate between two steps 
  #1. Label Propagation. Define gain as an increase of the number of edges between nodes of the same cluster. For every pair of nodes, switch their cluster labels in a greedy way until the gain cannot be increased any more.
  #2. Random Shuffling. Randomly select some pairs of nodes, and swap their labels.
  #Stop when label prop cannot increase gain 
  
  if (verbose) {
    cat("Start... ")
  }
  crit_out <- TRUE #break condition of outer loop
  k <- 0 #number of iteration
  
  start <- Sys.time()
  while (crit_out) {
    
    old_old_mod <- modularity(g, member_old) #outer-loop modularity
    k <- k+1
    
    #1. random shuffling
    if (k < max_step) {
      shuf <- sample(1:nc, shuffle, replace = FALSE)
      
      for (i in 1:shuffle) {
        node1 <- couples[1,shuf[i]]
        node2 <- couples[2,shuf[i]]
        
        #swap cluster label of two randomly selected nodes
        swap <- member[node1]
        member[node1] <- member[node2]
        member[node2] <- swap
      }
    }
    
    #2. greedy step
    crit_in <- TRUE #break condition in the inner loop
    while (crit_in) { 
      
      old_mod <- modularity(g, member) #inner-loop modularity
      
      #shuffle the orders of couples
      ord <- sample(1:nc, nc, replace = FALSE)
      
      for (i in 1:nc) {
        
        node1 <- couples[1,ord[i]] #first node
        node2 <- couples[2,ord[i]] #second node
        
        if (member[node1] != member[node2]) {
          
          #calculate the number of edges in a community
          temp <- which(member == member[node1])
          mod1 <- length(E(induced_subgraph(g, temp))) #number of edges in the cluster of node1
          temp <- setdiff(c(temp, node2), c(node1))
          mod2 <- length(E(induced_subgraph(g, temp))) #number of edges in the cluster of node1, but now node1 is replaced by node2
          
          temp <- which(member == member[node2])
          mod3 <- length(E(induced_subgraph(g, temp))) #number of edges in the cluster of node2
          temp <- setdiff(c(temp, node1), c(node2))
          mod4 <- length(E(induced_subgraph(g, temp))) #number of edges in the cluster of node2, but now node2 is replaced by node1
          
          if ((mod2+mod4) > (mod1+mod3)) {
            #swap two cluster labels
            swap = member[node1]
            member[node1] = member[node2]
            member[node2] = swap
          }
        }
      }
      
      new_mod <- modularity(g, member)
      crit_in <- abs(new_mod - old_mod) > tol #greedy until cannot get better
    }
    
    #update the cluster label after the iteration
    if (new_mod > old_old_mod) {
      member_old <- member
    } else {
      member <- member_old
      crit_out <- FALSE
    }
    
    #print updates
    if (verbose) {
      cat(round(new_mod,3), " (", round(abs(new_mod - old_old_mod),3), ") ... ")
    }
    
    crit_out <- (k <= max_step) & crit_out #criterion to stop outer loop
  }
  end <- Sys.time()
  if (verbose) {
    cat("Simulation ends in", hms_span(start, end), "\n") #print running time
  }
  
  #return
  member_old
}

#' Restreaming linear deterministic greedy clustering
#' 
#' @description This function clusters graphs using the restreaming linear deterministic greedy algorithm (Saveski et al. 2017)
#' @param g an igraph object
#' @param nrun number of iterations
#' @param max_size the maximum size of a cluster. Default is NULL, corresponding to maximum cluster size in balanced clustering
#' @param verbose whether the run progress should be printed
#' @return a vector of cluster labels
reLDG <- function(g, nrun = 20,
                  max_size = NULL, verbose = TRUE) {
  
  #initial clustering using the Louvain clustering algorithm
  member <- cluster_louvain(g)$membership
  n.block <- length(unique(member))
  n <- length(V(g))
  nodes <- 1:n
  
  #maximum cluster size
  if (is.null(max_size)) {
    max_size <- ceiling(n/n.block)
  }
  
  #starts
  if (verbose) {
    progress <- 0
    cat("Start... ")
  }
  start <- Sys.time()
  
  #iterate over multiple runs
  for (k in 1:nrun) {
    
    #randomize the orders of the nodes
    ord <- sample(nodes, n, replace = FALSE)
    
    #for each node
    for (u in ord) {
      #scores
      score <- rep(NA, n.block)
      nu <- neighbors(g, v = u)
      #calculate score of each cluster with respect to the node
      for (i in 1:n.block) { 
        ci <- nodes[member == i]
        score[i] <- length(intersect(ci, nu))*(1-length(ci)/max_size)
      }
      #assign a node to the cluster with highest score
      possible_cluster <- which(score == max(score))
      member[u] <- ifelse(length(possible_cluster) > 1, sample(possible_cluster, 1), possible_cluster) #change membership right away!
    }
    if (verbose) { #print progress
      if (k %/% (nrun %/% 10) > progress) {
        progress <- progress + 1
        cat(progress*10, "%... ")
      }
    }
  }
  end <- Sys.time()
  if (verbose) {
    cat("Simulation ends in", hms_span(start, end), "\n")
  }
  
  member #return
}

#' Social hash algorithm for clustering
#' 
#' @description This function adapts the social hash algorithm for bipartite graph (Shalita et al. 2016) to cluster undirected simple graphs
#' @param g an igraph object
#' @param n.block number of cluster. Default is NULL, corresponding to the number of clusters returned by the Louvain algorithm
#' @param nrun number of iterations
#' @param verbose whether the run progress should be printed
#' @return a vector of cluster labels
social_hash <- function(g, n.block = NULL, 
                        nrun = 20, verbose = TRUE) {
  
  #initialization
  if (is.null(n.block)) {
    n.block <- length(unique(cluster_louvain(g)$membership)) #default number of clusters are decided by the louvain algorithm 
  }
  n <- length(V(g))
  nodes <- 1:n
  
  #starts
  if (verbose) {
    progress <- 0
    cat("Start... ")
  }
  start <- Sys.time()
  
  #initial clustering - ensuring balance
  bucket <- nodes
  ord <- sample(nodes, n, replace = FALSE) #random order
  if (n%%n.block == 0) {
    bucket[ord] <- c(rep(1:n.block, n%/%n.block))
  } else {
    bucket[ord] <- c(rep(1:n.block, n%/%n.block), 1:(n%%n.block))
  }
  
  #local refinement
  for (k in 1:nrun) {
    
    #place holders for supporting matrices
    S <- matrix(0, nrow = n.block, ncol = n.block)
    prob <- matrix(0, nrow = n.block, ncol = n.block)
    target <- rep(NA, n)
    final_gain <- rep(NA, n)
    
    for (u in nodes) {
      #compute gain
      gain <- rep(NA, n.block)
      bucket_new <- bucket
      mod_old <- modularity(g, bucket)
      for (i in 1:n.block) {
        bucket_new[u] <- i
        gain[i] <- modularity(g, bucket_new) - mod_old
      }
      #find the best bucket
      possible_bucket <- which(gain == max(gain))
      target[u] <- ifelse(length(possible_bucket) > 1, sample(possible_bucket, 1), possible_bucket)
      final_gain[u] <- gain[target[u]]
      #update matrix
      if (gain[target[u]] > 0) {
        S[bucket[u], target[u]] <- S[bucket[u], target[u]] + 1
      }
    }
    
    #compute move probabilities
    for (i in 1:n.block) {
      for (j in 1:n.block) {
        prob[i,j] <- ifelse(S[i,j] == 0, 0, min(S[i,j], S[j,i])/S[i,j])
      }
    }
    
    #change buckets according to the move probabilities
    for (u in nodes) {
      if ((final_gain[u] > 0) & (runif(1) < prob[bucket[u],target[u]])) {
        bucket[u] <- target[u]
      }
    }
    
    #print the progress
    if (verbose) {
      if (k %/% (nrun %/% 10) > progress) {
        progress <- progress + 1
        cat(progress*10, "%... ")
      }
    }
  }
  end <- Sys.time()
  if (verbose) {
    cat("Simulation ends in", hms_span(start, end), "\n")
  }
  
  #return
  bucket
}

#' Clustering multiple times 
#' 
#' @description This function cluster a graph using one algorithm for multiple times and save the results into a .RData file
#' @param g an igraph object
#' @param cluster_func the function to perform clustering
#' @param name_obj the name of the R object to be saved
#' @param nrun number of run repetition
#' @param name name of the algorithm to be saved in the object
#' @param file the name of the .RData file to save the result object. If set to NULL, the results will not be saved in an .RData file
#' @param seed seed number given at the beginning of the runs
#' @param ... parameters to be passed to the clustering function cluster_func()
#' @return an object that contains the information about the runs, including the clustering labels, the running times and the name of the algorithm used
mass_clustering <- function(g, cluster_func, name_obj, 
                            nrun = 30, name = "balanced clustering",
                            file = "clustering.RData",
                            seed = 123456, ...) {
  
  #initialize place holders
  clusters <- c()
  times <- c()
  
  #set the randomization seed
  if (!is.null(seed)) {set.seed(seed)}
  
  for (i in 1:nrun) {
    start <- Sys.time()
    member <- do.call(cluster_func, args = list(g = g, ...)) #clustering
    if (is.list(member)) {
      member <- member$member
    }
    clusters <- cbind(clusters, member) #add in the results
    end <- Sys.time()
    times <- c(times, hms_span(start, end)) #add in the running times
    
    if (!is.null(file)) {
      eval(call("<-", as.name(name_obj), 
                list(clusters = clusters, times = times, algo = name)))
      resave(list = name_obj, file = file) #save into the file
    }
    
    cat(i, "/", nrun, "\n")
  }
  
  list(clusters = clusters, 
       times = times, 
       algo = name)
}

#------------------------------
#Treatment assignment
#------------------------------

#' Treatment assignment based on graph clustering
#' 
#' @description This function assigns treatments to units on a graph based on their cluster labels
#' @param member the vector of cluster label
#' @param p_assign the probability of assigning to treatment
#' @param clustered whether the randomization is done on the cluster-level or individual level. Default is TRUE, that is, randomly choose half of the clusters and assign all units in the selected clusters to treatment. If set to FALSE, randomly choose half the units and assign them to treatment
#' @return a vector of treatment assignment (1 for treatment and 0 for control) 
assign_mixed <- function(member, p_assign = 0.5, 
                         clustered = FALSE) {
  
  #initialization
  n <- length(member)
  assign <- rep(0, n)
  
  if (clustered) {
    #randomly choose half of the clusters to treatment/control
    n.block <- length(unique(member))
    samp <- sample(1:n.block, round(n.block*p_assign), replace = FALSE)
    for (i in 1:n.block) {
      if (i %in% samp) {
        assign[member == i] <- 1
      }
    }
  } else {
    #randomly choose half of the units to treatment/control
    samp <- sample(1:n, round(n*p_assign), replace = FALSE)
    assign[samp] <- 1
  }
  
  #return
  assign
}

#' Assigning treatments to multiple clustering
#' 
#' @description This function performs treatment assignments to an object returned by the mass_clustering() function, save the results to the same object
#' @param ... the clustering objects returned by the mass_clustering() function
#' @param assign_func the treatment assignment function
#' @param seed the seed for randomization given at the beginning of treatment assignment process
#' @param control.assign the list of parameters to be passed to the assign_func function
mass_assign <- function(..., assign_func, 
                        seed = 123456, 
                        control.assign = list()) {
  
  #get the names of all objects
  all_vars <- as.list(match.call())[-1]
  var_names <- sapply(all_vars, as.character)
  var_names <- var_names[!names(var_names) %in% c("assign_func", "seed", 
                                                  "control.assign")]
  nobj <- length(var_names)
  pf <- parent.frame()
  
  #start treatment assignment for each object
  for (i in 1:nobj) {
    
    #extract the objects
    nm <- var_names[[i]]
    obj <- get(nm) 
    
    #assigning 
    if (!is.null(seed)) {set.seed(seed)}
    ndesigns <- dim(obj$clusters)[2]
    best_designs <- c()
    
    for (j in 1:ndesigns) {
      assignment <- do.call(assign_func, args = c(list(member = obj$clusters[,j]),
                                                  control.assign))
      best_designs <- cbind(best_designs, assignment)
    }
    
    #assign the new results to current names
    eval(call("<-", as.name(nm), 
              list(best_designs = best_designs, clusters = obj$clusters,
                   times = obj$times, algo = obj$algo)), 
         envir = pf)
  }
}

save(louvain_wrapper, balanced_label_prop, reLDG, social_hash, 
     mass_clustering, resave,
     assign_mixed, mass_assign,
     file = "funcs/cluster_funcs.RData")

#######################################################
#Test if things work as expected
#######################################################

library(igraph)

load("funcs/mc_funcs.RData")
load("funcs/cluster_funcs.RData")

#different clustering algorithms
set.seed(123456)
tmp <- louvain_wrapper(g_enron)
modularity(g_enron, tmp) #check the modularity to see the quality of the clustering. higher modularity means better clustering. the function is built-in in the igraph package
tmp <- reLDG(g_enron) #about 6 secs
modularity(g_enron, tmp)
tmp <- social_hash(g_enron) #about 3 secs
modularity(g_enron, tmp)
tmp <- balanced_label_prop(g_enron) #about 5 minutes
modularity(g_enron, tmp)

#clustering multiple times
enron_louvain <- mass_clustering(g_enron, 
                                 cluster_func = louvain_wrapper, 
                                 name = "louvain", 
                                 name_obj = "enron_louvain", 
                                 file = "results/checks/checks_graph_cluster.RData")

rm(list = ls()) #clean the environment
load("results/checks/checks_graph_cluster.RData") #check the result file
View(enron_louvain) #object "enron_louvain" was created and saved in the "checks_graph_cluster.RData" file. The object contains information about 5 louvain clusterings

#mass treatment assignment
load("funcs/cluster_funcs.RData")
mass_assign(enron_louvain, assign_func = assign_mixed,
            control.assign = list(clustered = TRUE))
View(enron_louvain$clusters)
View(enron_louvain$best_designs)


