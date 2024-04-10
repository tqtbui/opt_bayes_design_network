#######################################################
#This file contains functions for the meta-heuristic optimization algorithms
#######################################################

#---------------------
#Meta-heuristic algorithms
#---------------------

#' Random search
#' 
#' @description This function finds a good design using random search
#' @param N the number of experimental units
#' @param crit.func the function to calculate design criterion
#' @param control.crit the list of parameters to be passed to the crit.func() function
#' @param timed_search whether the duration of the algorithm is printed. Default is TRUE
#' @param seed_search seed number given at the beginning of the algorithm. Default is NULL, which means "no seed is given"
#' @param ndesign number of designs to be considered in the search
#' @param lookahead when "ndesign" designs are considered, we can choose to continue the search by looking at "lookahead" more designs. If a better design is found, the algorithm will continue to search for another "lookahead" design, and so on
#' @return an object that contains the best design found by the algorithm, the design criterion at the best design, the running time, and other information about the search   
random_search <- function(N, crit.func, control.crit = list(), 
                          timed_search = TRUE, seed_search = NULL, 
                          ndesign = 5000, lookahead = 0) {
  
  #initialization
  best_design <- c()
  best_val <- 0
  trace_val <- c()
  
  #initialize counting
  cnt <- 0
  cnt_ahead <- 0
  is.continue <- TRUE
  
  #set seed if specified
  if (!is.null(seed_search)) {
    set.seed(seed_search)
  }
  
  cat("Searching... ")
  
  start <- Sys.time()
  #.....................
  #recursive search
  #.....................
  while (is.continue) {
    
    #generate a design
    a <- generate_random_design(N)
    
    #estimate the criteria over the prior
    val <- do.call(crit.func, c(list(assignment = a), control.crit))$value
    
    #update the trace
    trace_val <- c(trace_val, val)
    
    #conditions
    if (cnt < ndesign) {
      
      cnt <- cnt + 1

      #check for best of the batch
      if ((cnt == 1) | (val < best_val)) { 
        #smaller value of criterion is better
        best_design <- a
        best_val <- val
      }
      
      #update progress
      if (cnt %% (ndesign%/%10) == 0) {
        cat(cnt /(ndesign%/%10)*10, "%... ")
      }
      
      #stop searching and start intensifying
      if (cnt == ndesign) {
        cat("Intensifying... ")
      }
      
    } else {
      
      cnt_ahead <- cnt_ahead + 1
      
      if (val < best_val) {
        
        best_design <- a
        best_val <- val
        
        #refresh lookahead if there is new best design found
        cnt_ahead <- 0 
        cat("refreshed... ")
      }
    }
    
    #conditions to stop
    if ((cnt == ndesign) & (cnt_ahead == lookahead)) {
      is.continue <- FALSE
    }
  }
  end <- Sys.time()
  if (timed_search) {
    cat("Simulation ends in", hms_span(start, end), "\n")
  }
  
  #return
  list(best_val = best_val, best_design = best_design, 
       time = hms_span(start, end), 
       name = "random search", trace_val = trace_val)
}

#' Tabu search
#' 
#' @description This function finds a good design using tabu search
#' @param N the number of experimental units
#' @param crit.func the function to calculate design criterion
#' @param control.crit the list of parameters to be passed to the crit.func() function
#' @param timed_search whether the duration of the algorithm is printed. Default is TRUE
#' @param seed_search seed number given at the beginning of the algorithm. Default is NULL, which means "no seed is given"
#' @param ndesign number of designs to be considered in the search
#' @param lookahead when "ndesign" designs are considered, we can choose to continue the search by looking at "lookahead" more designs. If a better design is found, the algorithm will continue to search for another "lookahead" design, and so on
#' @param start_design the design to start with. The default is NULL, which means a design will be randomly generated to start the algorithm
#' @param ntabu the number of designs to be stored in the tabu list
#' @param nneighbor the number of neighbor designs to look at when finding the next design
#' @param aspiration_rate if a neighbor design is already contained in the tabu list, there is still "aspiration_rate" chance that this design will be selected as the next design 
#' @param shiftprob the "shiftprob" parameter to be passed to the generate_neighbor_deisgn() function to generate a neighbor design from a current design 
#' @return an object that contains the best design found by the algorithm, the design criterion at the best design, the running time, and other information about the search   
tabu_search <- function(N, crit.func, control.crit = list(), 
                        timed_search = TRUE, seed_search = NULL, 
                        ndesign = 5000, lookahead = 0, start_design = NULL,
                        ntabu = 100, nneighbor = 100, 
                        aspiration_rate = 0.1, shiftprob = 0.1) {
  
  #initialization
  best_design <- c()
  best_val <- 0
  tabu <- matrix(0, ncol = ntabu, nrow = N) #initialize a tabu list
  neighbors <- matrix(0, ncol = nneighbor, nrow = N) #initialize a list of neighbor designs
  val_neighbors <- rep(0, nneighbor) #initialize the design criterion values of neighbors
  
  #initialize counting
  cnt <- 0
  cnt_ahead <- 0
  is.continue <- TRUE
  progress <- 0
  
  #set seed if specified
  if (!is.null(seed_search)) {
    set.seed(seed_search)
  }
  
  #set up the first design
  cat("Searching... ")
  start <- Sys.time()
  if (is.null(start_design)) {
    start_design <- generate_random_design(N)
  }
  a_old <- start_design #set the current design
  tabu[,1] <- c(a_old) #add the design to the tabu list
  best_design <- a_old #set the currently best design
  best_val <- do.call(crit.func, c(list(assignment = a_old), control.crit))$value #approximate the design criterion of the first design using monte carlo
  trace_val <- best_val #add to the trace
  cnt <- cnt + 1
  
  #main algorithm
  while (is.continue) {
    
    #generate neighboring designs
    for (j in 1:nneighbor) {
      neighbors[,j] <- generate_neighbor_design(a_old, shiftprob = shiftprob)
      val_neighbors[j] <- do.call(crit.func, c(list(assignment = neighbors[,j]), control.crit))$value #approximate the deisgn criterion value using monte carlo
    }
    
    #order the neighboring designs according on values
    trace_val <- c(trace_val, val_neighbors)
    neighbors <- neighbors[,order(val_neighbors)]
    val_neighbors <- sort(val_neighbors)
    
    #choose the best neighbors according to tabu criteria
    is.next <- TRUE
    cnt_neighbors <- 0
    while (is.next) {
      cnt_neighbors <- cnt_neighbors + 1
      candidate <- neighbors[,cnt_neighbors]
      
      #check if the new design is in the tabu list 
      is.tabu <- sum(sapply(1:ntabu, function(k) {
        sum(candidate != tabu[,k]) == 0 #if new a is in the tabu list
      })) != 0
      
      #check if aspired or not
      is.aspire <- runif(1, min = 0, max = 1) < aspiration_rate

      if ((!is.tabu) | (is.tabu & is.aspire)) { #choose this design as the next design
        is.next <- FALSE
        val <- val_neighbors[cnt_neighbors] 
      } else if (cnt_neighbors == nneighbor) { #if no neighbor design satisfy the conditions, choose the next design as the current design
        is.next <- FALSE
        candidate <- a_old
        val <- trace_val[i]
      }
    }
    
    #update the new design
    a_old <- candidate
    
    #update the tabu list
    if (!is.tabu) {
      tabu <- cbind(candidate, tabu)
      tabu <- tabu[,1:ntabu]
    }
    
    #conditions
    if (cnt < ndesign) {
      
      cnt <- cnt + nneighbor
      
      #check for best of the batch
      if (val < best_val) { 
        #smaller value of criterion is better
        best_design <- candidate
        best_val <- val
      }
      
      #update progress
      if (cnt %/% (ndesign%/%10) > progress) {
        progress <- progress + 1
        cat(progress*10, "%... ")
      }
      
      #stop searching and start intensifying
      if (cnt >= ndesign) {
        cat("Intensifying... ")
      }
      
    } else {
      
      cnt_ahead <- cnt_ahead + nneighbor
      
      if (val < best_val) {
        
        best_design <- candidate
        best_val <- val
        
        #refresh lookahead if there is new best design found
        cnt_ahead <- 0
        cat("refreshed... ")
      }
    }
    
    #conditions to stop
    if ((cnt >= ndesign) & (cnt_ahead >= lookahead)) {
      is.continue <- FALSE
    }
  }
  end <- Sys.time()
  if (timed_search) {
    cat("Simulation ends in", hms_span(start, end), "\n")
  }
  
  #return
  list(best_val = best_val, best_design = best_design, 
       time = hms_span(start, end), 
       name = "tabu search", trace_val = trace_val)
}

#' Simulated annealing
#' 
#' @description This function finds a good design using simulated annealing
#' @param N the number of experimental units
#' @param crit.func the function to calculate design criterion
#' @param control.crit the list of parameters to be passed to the crit.func() function
#' @param timed_search whether the duration of the algorithm is printed. Default is TRUE
#' @param seed_search seed number given at the beginning of the algorithm. Default is NULL, which means "no seed is given"
#' @param ndesign number of designs to be considered in the search
#' @param lookahead when "ndesign" designs are considered, we can choose to continue the search by looking at "lookahead" more designs. If a better design is found, the algorithm will continue to search for another "lookahead" design, and so on
#' @param start_design the design to start with. The default is NULL, which means a design will be randomly generated to start the algorithm
#' @param cooling_step the temperature will be cooled down for every "cooling_step" iterations
#' @param cooling_coef the temeprature will be cooled down to "cooling_coef"*temperature
#' @param restart_step the number of iterations before the search is restarted, i.e., a new candidate design is chosen randomly. Default is NULL, which means the search will never be restarted.
#' @param max_rejection the maximum of design rejections before moving to the next iteration. Default is 100.
#' @param shiftprob the "shiftprob" parameter to be passed to the generate_neighbor_deisgn() function to generate a neighbor design from a current design 
#' @return an object that contains the best design found by the algorithm, the design criterion at the best design, the running time, and other information about the search   
simulated_annealing <- function(N, crit.func, control.crit = list(), 
                                timed_search = TRUE, seed_search = NULL, 
                                ndesign = 5000, lookahead = 0, start_design = NULL, 
                                cooling_step = 50, cooling_coef = 0.9, 
                                restart_step = NULL, max_rejection = 100, 
                                shiftprob = 0.1) {
  
  #initialization
  temp <- 1 #temperature
  cnt <- 0
  cnt_ahead <- 0
  is.continue <- TRUE
  progress <- 0
  
  #set seed if specified
  if (!is.null(seed_search)) {
    set.seed(seed_search)
  }
  
  #set the first design
  cat("Searching... ")
  start <- Sys.time()
  if (is.null(start_design)) {
    start_design <- generate_random_design(N)
  }
  a_old <- start_design #current design
  current_val <- do.call(crit.func, c(list(assignment = a_old), control.crit))$value #approximate the deisgn criterion value of the first design using monte carlo
  trace_val <- current_val #add the value to the trace
  best_design <- a_old #current best design
  best_val <- current_val #current best design criterion
  cnt <- cnt + 1 #number of designs evaluated by MC
  step <- 0 #count number of iterations
  
  #main algorithm
  while (is.continue) {
    
    step <- step + 1
    
    #restart if specified
    is.search <- TRUE #whether we should search for the next design in this iteration, or we should just restart everything
    if (!is.null(restart_step)) {
      if (step %% restart_step == 0) {
        a_old <- generate_random_design(N) #randomly choose a design to restart the search
        current_val <- do.call(crit.func, c(list(assignment = a_old), control.crit))$value #approximate the deisgn criterion value using monte carlo
        trace_val <- c(trace_val, current_val)
        cnt <- cnt + 1
        is.search <- FALSE #do not search for a next design in this iteration where the restart happens
      } 
    }
    
    #search for the next design
    if (is.search) {
      is.next <- TRUE  #whether we should consider the next design
      cnt_rejection <- 1 #number of designs evaluated in each acceptance step
      while (is.next) {
        
        #generate a neighbor of current design
        a_new <- generate_neighbor_design(a_old, shiftprob = shiftprob)
        
        #evaluate new design
        val <- do.call(crit.func, c(list(assignment = a_new), control.crit))$value 
        
        #add into the trace
        trace_val <- c(trace_val, val) #value of candidate design
        
        #accept or not
        if ((val < current_val) | (runif(1, min = 0, max = 1) < exp((current_val - val)/temp))) {
          is.next <- FALSE
          a_old <- a_new
          current_val <- val
        } else {
          cnt_rejection <- cnt_rejection + 1
        }
        
        #move to next iteration if the maximum of rejections is reached  
        if ((cnt_rejection == max_rejection) | ((cnt + cnt_rejection) >= ndesign)) {
          is.next <- FALSE
        }
      }
    }
    
    #conditions
    if (cnt < ndesign) {

      cnt <- cnt + cnt_rejection
      
      #check for best of the batch
      if (current_val < best_val) {
        best_design <- a_old
        best_val <- current_val
      }
      
      #update progress
      if (cnt %/% (ndesign%/%10) > progress) {
        progress <- progress + 1
        cat(progress*10, "%... ")
      }
      
      #stop searching and start intensifying
      if (cnt >= ndesign) {
        cat("Intensifying... ")
      }
      
    } else {
      
      cnt_ahead <- cnt_ahead + cnt_rejection
      
      if (current_val < best_val) {
        
        best_design <- a_old
        best_val <- current_val
        
        #refresh lookahead if there is new best design found
        cnt_ahead <- 0
        cat("refreshed... ")
      }
    }
    
    #conditions to stop
    if ((cnt >= ndesign) & (cnt_ahead >= lookahead)) { #stop when the number of designs considered is reached
      is.continue <- FALSE 
    }
    
    #cool the temperature
    if (step %% cooling_step == 0) {
      temp <- temp*cooling_coef
    }
  }
  end <- Sys.time()
  if (timed_search) {
    cat("Simulation ends in", hms_span(start, end), "\n")
  }
  
  #return
  list(best_val = best_val, best_design = best_design, 
       time = hms_span(start, end), 
       name = "simulated annealing", trace_val = trace_val)
}

#' Genetic algorithm
#' 
#' @description This function finds a good design using genetic algorithm
#' @param N the number of experimental units
#' @param crit.func the function to calculate design criterion
#' @param control.crit the list of parameters to be passed to the crit.func() function
#' @param timed_search whether the duration of the algorithm is printed. Default is TRUE
#' @param seed_search seed number given at the beginning of the algorithm. Default is NULL, which means "no seed is given"
#' @param ndesign number of designs to be considered in the search
#' @param lookahead when "ndesign" designs are considered, we can choose to continue the search by looking at "lookahead" more designs. If a better design is found, the algorithm will continue to search for another "lookahead" design, and so on
#' @param npop the number of designs in the population, which is maintained in each iteration
#' @param parentsprob the percentage of designs in the population eligible to become parents
#' @param eliteprob the percentage of top designs in the population who will be pass to the next iteration without being changed
#' @param mutateprob the percentage of units in the children designs whose treatment assignments will be randomly selected
#' @return an object that contains the best design found by the algorithm, the design criterion at the best design, the running time, and other information about the search   
genetic_algorithm <- function(N, crit.func, control.crit = list(), 
                              timed_search = TRUE, seed_search = NULL, 
                              ndesign = 5000, lookahead = 0, npop = 100, 
                              parentsprob = 0.5, eliteprob = 0.1, 
                              mutateprob = 0.1) {
  
  #initialization
  trace_val <- c()
  ntop <- round(npop*eliteprob) #number of designs in the population which will directly go to next iteration without being modified
  nparents <- round(npop*parentsprob) #number of eligible parents in the population
  nmutate <- round(npop*mutateprob) #number of units in children designs whose treatment assignments will be randomly generated
  
  #initialize counting
  cnt <- 0
  cnt_ahead <- 0
  is.continue <- TRUE
  progress <- 0

  #set seed if specified
  if (!is.null(seed_search)) {
    set.seed(seed_search)
  }
  
  #generate the initial population
  cat("Searching... ")
  start <- Sys.time()
  pop <- c()
  pop_val <- c()
  for (i in 1:npop) {
    a <- generate_random_design(N) #randomy generate a design 
    val <- do.call(crit.func, c(list(assignment = a), control.crit))$value #approximate the deisgn criterion value using monte carlo
    pop <- cbind(pop, a) #add to the population
    pop_val <- c(pop_val, val)
  }
  pop <- pop[,order(pop_val)] #order the population using the design criterion
  pop_val <- sort(pop_val)
  best_design <- pop[,1] #current best design
  best_val <- pop_val[1] #current best design criterion
  trace_val <- pop_val #trace
  cnt <- cnt + npop
  
  #main algorithm
  while (is.continue) {
    
    #choose the elite population to the new population
    new_pop <- pop[,1:ntop]
    new_val <- pop_val[1:ntop]
    
    #create children
    for (j in 1:(npop-ntop)) {
      
      #randomly choose parents
      idx_parents <- sample(1:nparents, 2, replace = FALSE)
      child <- rep(0, N)
      
      #choose the genes for mutation
      idx_mutate <- sample(1:N, nmutate, replace = FALSE)
      child[idx_mutate] <- generate_random_design(nmutate)
      idx_inherit <- setdiff(1:N, idx_mutate)
      
      #randomly choose genes from mom
      idx_mom <- sample(1:(N-nmutate),  round((N-nmutate)/2), replace = FALSE)
      idx_mom <- idx_inherit[idx_mom]
      child[idx_mom] <- pop[,idx_parents[1]][idx_mom]
      
      #rest of the genes are from dad
      idx_dad <- setdiff(1:N, union(idx_mutate, idx_mom))
      child[idx_dad] <- pop[,idx_parents[2]][idx_dad]
      
      #append the child to new population
      new_pop <- cbind(new_pop, child)
      
      #assess the design criterion of the child
      val <- do.call(crit.func, c(list(assignment = child), control.crit))$value
      new_val <- c(new_val, val)
    }
    
    pop <- new_pop
    pop_val <- new_val
    trace_val <- c(trace_val, pop_val)
    
    #order the population from smallest to largest design criterion
    pop <- pop[,order(pop_val)]
    pop_val <- sort(pop_val)
    
    #conditions
    if (cnt < ndesign) {
      
      cnt <- cnt + npop
      
      #check for best of the batch
      if (pop_val[1] < best_val) { 
        #smaller value of criterion is better
        best_design <- pop[,1]
        best_val <- pop_val[1]
      }
      
      #update progress
      if (cnt %/% (ndesign%/%10) > progress) {
        progress <- progress + 1
        cat(progress*10, "%... ")
      }
      
      #stop searching and start intensifying
      if (cnt >= ndesign) {
        cat("Intensifying... ")
      }
      
    } else {
      
      cnt_ahead <- cnt_ahead + npop
      
      if (pop_val[1] < best_val) { 

        best_design <- pop[,1]
        best_val <- pop_val[1]
        
        #refresh lookahead if there is new best design found
        cnt_ahead <- 0
        cat("refreshed... ")
      }
    }
    
    #conditions to stop
    if ((cnt >= ndesign) & (cnt_ahead >= lookahead)) {
      is.continue <- FALSE
    }
  }
  end <- Sys.time()
  if (timed_search) {
    cat("Simulation ends in", hms_span(start, end), "\n")
  }
  
  #return
  list(best_val = best_val, best_design = best_design, 
       time = hms_span(start, end), 
       name = "genetic algorithm", trace_val = trace_val)
}

#~~~~~~~~~~~~~~~~~~~~~~
#Save the workspace
#~~~~~~~~~~~~~~~~~~~~~~~

save(random_search, tabu_search, 
     simulated_annealing, genetic_algorithm, 
     file = "funcs/heuristics_funcs.RData")

#######################################################
#Test if things work as expected
#######################################################

library(igraph)
library(parallel)
library(MASS)
library(Matrix)
library(Rcpp)
sourceCpp("rcpp_help_funcs.cpp")
#if Rcpp does not work, uncomment and run the command below
#source("rcpp_help_funcs.R")

load("funcs/mc_funcs.RData")
load("funcs/gane_funcs.RData")
load("funcs/run_funcs.RData")
load("funcs/heuristics_funcs.RData")

#-----------------------
#setting the problem parameters
#-----------------------

g <- g_enron_gane
gane_prior_now <- function(i) { 
  c(0, 1, 0.5, 0.1, 0.5, 1) #set the prior to constant so that the test run is faster
}

#-----------------------
#running different algorithms
#-----------------------

#random search
check_rs <- random_search(N = g$N, crit.func = mc_gane, 
                          control.crit = list(g = g, 
                                              prior = gane_prior_now,
                                              nsim = 1)) #constant prior so only need 1 simulation to calculate the design criterion
plot(1:length(check_rs$trace_val), check_rs$trace_val, type = "l", 
     xlab = "iteration", ylab = "design criterion", 
     main = "random search results by iteration") #plot the design criteria of designs visited in each iteration

#tabu search
check_ts <- tabu_search(N = g$N, crit.func = mc_gane, 
                        control.crit = list(g = g, 
                                            prior = gane_prior_now,
                                            nsim = 1)) #constant prior so only need 1 simulation to calculate the design criterion
plot(1:length(check_ts$trace_val), check_ts$trace_val, type = "l", 
     xlab = "iteration", ylab = "design criterion", 
     main = "tabu search results by iteration") #plot the design criteria of designs visited in each iteration

#simulated annealing
check_sa <- simulated_annealing(N = g$N, crit.func = mc_gane, 
                                control.crit = list(g = g, 
                                                    prior = gane_prior_now,
                                                    nsim = 1)) #constant prior so only need 1 simulation to calculate the design criterion
plot(1:length(check_sa$trace_val), check_sa$trace_val, type = "l", 
     xlab = "iteration", ylab = "design criterion", 
     main = "simulated annealing results by iteration") #plot the design criteria of designs visited in each iteration

#genetic algorithm
check_ga <- genetic_algorithm(N = g$N, crit.func = mc_gane, 
                              control.crit = list(g = g, 
                                                  prior = gane_prior_now,
                                                  nsim = 1)) #constant prior so only need 1 simulation to calculate the design criterion
plot(1:length(check_ga$trace_val), check_ga$trace_val, type = "l", 
     xlab = "iteration", ylab = "design criterion", 
     main = "genetic algorithm results by iteration") #plot the design criteria of designs visited in each iteration

save(check_rs, check_ts, check_sa, check_ga, 
     file = "results/checks/checks_algo_runs.RData")

#-----------------------
#running an algorithm several times
#-----------------------

tmp <- test_algo(nrun = 5, crit.func = mc_gane, 
                 control.crit = list(g = g, 
                                     prior = gane_prior_now,
                                     nsim = 1), #constant prior so only need 1 simulation to calculate the design criterion
                 algo.func = genetic_algorithm, 
                 control.algo = list(N = g$N), 
                 file = "results/checks/try-meta-heuristics.RData") 

rm(list = ls()) #clean the environment
load("results/checks/try-meta-heuristics.RData") #check the result file
View(genetic_algorithm) #object "genetic_algorithm" was created and saved in the "try-meta-heuristics.RData" file. The object contains information about 5 genetic algorithm runs
