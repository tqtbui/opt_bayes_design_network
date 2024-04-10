#######################################################
#This file contains help functions to run optimization algorithms
#######################################################

#' Generate random design
#' 
#' @description This function generates a random design from the uniform distribution, where each node independently has 50% chance of getting treatment/control
#' @param n the number of experimental units
#' @return a vector containing treatment assignment (1 for treatment and 0 for control) for n units
generate_random_design <- function(n) {
  rbinom(n, 1, 0.5) #uniform distribution 
}

#' Generate neighbor design
#' 
#' @description This function generates a neighbor design from the current design by randomly choose a proportion of units to possibly change their treatment assignment from the current design
#' @param assignment_old the current design
#' @param shiftprob the proportion of units randomly chosen whose treatment assignment will be possibly change
#' @param random whether the treatment assignment change is random (independent bernoulli trial at 0.5 probability) or deterministic (who is previously treated will be assigned to control and vice versa)
#' @return a vector containing treatment assignment (1 for treatment and 0 for control) for n units
generate_neighbor_design <- function(assignment_old, 
                                     shiftprob = 0.1, 
                                     random = TRUE) {
  
  #randomly pick units whose treatment assignments will be possibly changed
  n <- length(assignment_old)
  nchange <- round(n*shiftprob)
  which.change <- sample(1:n, nchange, replace = FALSE)
  
  #change the treatment assignments of selected units
  ret <- assignment_old
  if (random) {
    ret[which.change] <- generate_random_design(nchange) #random change
  } else {
    ret[which.change] <- 1 - assignment_old[which.change] #deterministic change
  }
  
  #return  
  ret
}

#' Merge results of two objects
#' 
#' @description This function merges the results of two result objects (output of the test_algo() func)
#' @param obj1 the first result object
#' @param obj2 the second result object
#' @return a new result object which contains results from obj1 and obj2 
merge_results <- function(obj1, obj2) {
  best_designs <- cbind(obj1$best_designs, obj2$best_designs)
  best_vals <- c(obj1$best_vals, obj2$best_vals)
  times <- c(obj1$times, obj2$times)
  traces <- cbind(obj1$traces, obj2$traces)
  algo <- obj1$algo
  list(best_designs = best_designs, best_vals = best_vals,
       times = times, traces = traces, algo = algo)
}

#' Run an algorithm multiple times
#' 
#' @description This function runs an algorithm for multiple times and store the results in a .RData file
#' @param nrun the number of runs 
#' @param name the name of the algorithm to be run. Default is NULL, which means the name of the algorithm will be used
#' @param name_obj the name of the R object to be saved. Default is NULL, which means the name of the algorithm will be used
#' @param crit.func the function to calculate design criterion
#' @param control.crit the list of parameters to be passed to the crit.func() function
#' @param algo.func the design-finding algorithm function
#' @param control.algo the list of parameters to be passed to the algo.func() function
#' @param file the name of the .RData file to save the result object. Default is "", which means the results will not be saved into a .RData file
#' @return an object that contains the best designs, corresponding design criteria, runing times, etc. from the multiple runs 
test_algo <- function(nrun = 10, name = NULL, name_obj = NULL,
                      crit.func, control.crit = list(), 
                      algo.func, control.algo = list(), 
                      file = "") {
  
  #initialization
  if (is.null(name)) {
    name <- deparse(substitute(algo.func)) #name of the algorithm
  }
  if (is.null(name_obj)) {
    name_obj <- deparse(substitute(algo.func)) #name of the object to be saved in the .RData file
  }
  
  #initialize result holder
  best_vals <- c()
  best_designs <- c()
  times <- c()
  traces <- c()
  
  #run the algos several times
  for (i in 1:nrun) {
    
    tmp <- do.call(algo.func, c(list(crit.func = crit.func, control.crit = control.crit), 
                                control.algo)) #run algorithm

    best_designs <- cbind(best_designs, tmp$best_design) #best designs found
    best_vals <- c(best_vals, tmp$best_val) #design criterion values of best designs
    times <- c(times, tmp$time) #running time
    traces <- cbind(traces, tmp$trace_val) #traces
    
    if (file != "") { #save the object into a .RData file
      
      previous  <- suppressWarnings(try(load(file), silent = TRUE)) #load the file
      if ("try-error" %in% class(previous)) {
        save(name, file = file) #if the file is not existing, create a file with the specified name
      } 
      previous <- load(file) #load the file again
      
      #result object to be saved
      tmp_obj <- list(best_designs = tmp$best_design, best_vals = tmp$best_val, 
                      times = tmp$time, traces = tmp$trace_val, algo = name) #object to be saved
      
      if (name_obj %in% previous) {
        eval(call("<-", as.name(name_obj), 
                  merge_results(get(name_obj), tmp_obj))) #if an object of the same name exists in the .RData file, merge it with the current object 
      } else {
        eval(call("<-", as.name(name_obj), 
                  tmp_obj)) #if this object does not exist in the file, save the new object as it is
      }
      
      #save to the specified .RData file
      save(list = unique(c(previous, name_obj)), file = file)
    }
  }
  
  #return the result object
  list(best_designs = best_designs, best_vals = best_vals, 
       times = times, traces = traces,
       algo = name)
}

#~~~~~~~~~~~~~~~~~~~~~~
#Save the workspace
#~~~~~~~~~~~~~~~~~~~~~~~

save(generate_neighbor_design, generate_random_design, 
     test_algo, merge_results,
     file = "funcs/run_funcs.RData")
