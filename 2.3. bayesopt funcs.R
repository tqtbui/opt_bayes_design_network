#######################################################
#This file contains functions for the Bayesian optimization algorithms
#######################################################

#---------------------
#Bayesian optimization algorithms
#---------------------

#' Randomly generate data from multivariate independent Bernoulli distribution
#' 
#' @description This function randomly generates data from multivariate independent Bernoulli distribution
#' @param probvec The vector of success probability for each variable
#' @return a vector of values generated from the multivariate independent Bernoulli distribution
rbinom_mod <- function(probvec) {
  sapply(1:length(probvec), function(i) {rbinom(1, size = 1, prob = probvec[i])})
}

#' Likelihood of multivariate independent Bernoulli distribution
#' 
#' @description This function outputs the likelihood of an observation that follows a multivariate independent Bernoulli distribution  
#' @param x The vector of data (observation)
#' @param probvec The vector of success probability for each variable
#' @return the likelihood that the data is observed, given the data follow a multivariate independent Bernoulli distribution
pbinom_mod <- function(x, probvec) {
  prod(sapply(1:length(x), function(i) {ifelse(x[i] == 1, probvec[i], 1-probvec[i])}))
}

#' Tree-Parzen estimator
#' 
#' @description This function finds a good design using the Tree-Parzen estimator
#' @param N the number of experimental units
#' @param crit.func the function to calculate design criterion
#' @param control.crit the list of parameters to be passed to the crit.func() function
#' @param timed_search whether the duration of the algorithm is printed. Default is TRUE
#' @param seed_search seed number given at the beginning of the algorithm. Default is NULL, which means "no seed is given"
#' @param ndesign number of designs to be considered in the search
#' @param lookahead when "ndesign" designs are considered, we can choose to continue the search by looking at "lookahead" more designs. If a better design is found, the algorithm will continue to search for another "lookahead" design, and so on
#' @param npop the number of candidate designs to generate in each iteration
#' @param ninitial the number of randomly generated initial designs. These designs are used to estimate the distribution of design given that the design criterion is low/high 
#' @param threshold the quantile of design criterion used to classify whether a design is having low/high design criterion. This means that the smallest "threshold"*100% of observed design criterion values will be considered low. 
#' @return an object that contains the best design found by the algorithm, the design criterion at the best design, the running time, and other information about the search   
tree_parzen <- function(N, crit.func, control.crit = list(), 
                        timed_search = TRUE, seed_search = NULL, 
                        ndesign = 5000, lookahead = 0, npop = 100, 
                        ninitial = 1000, threshold = 0.15) {
  
  #initialize counting
  cnt <- 0
  cnt_ahead <- 0
  is.continue <- TRUE
  progress <- 0
  
  #initialization
  threshold_step <- ceiling(1/threshold) #number of iterations per which the threshold value is recalculated
  
  #set seed if specified
  if (!is.null(seed_search)) {
    set.seed(seed_search)
  }
  
  #generate an initial population
  cat("Searching... ")
  start <- Sys.time()
  designs <- c()
  designs_val <- c()
  for (i in 1:ninitial) {
    a <- generate_random_design(N) #randomly generate a design
    val <- do.call(crit.func, c(list(assignment = a), control.crit))$value #approximate the design criterion using monte carlo
    designs <- cbind(designs, a) #add the design to the design collection
    designs_val <- c(designs_val, val)
  }
  best_val <- min(designs_val) #current best design criterion
  best_design <- designs[,which.min(designs_val)] #current best design
  ustar <- quantile(designs_val, probs = threshold) #current threshold value
  cnt <- cnt + ninitial
  
  #main algorithm
  while (is.continue) {
    
    #calculate the probability
    probs_low <- rowMeans(designs[,designs_val <= ustar]) #probability for each unit to be assigned to treatment if a design has a low (lower than threshold) design criterion value
    probs_high <- rowMeans(designs[,designs_val > ustar]) #probability for each unit to be assigned to treatment if a design has a higher (higher than threshold) design criterion value
    
    #generate a pool of candidates
    pool_best_design <- c()
    pool_best_val <- 0.5 #we want to maximize P(z|low)/P(z|high), which is equivalent to maximize the expected improvement. We expect the value of P(z|low)/P(z|high) to be more than 1. we start with 0.5.
    for (j in 1:npop) {
      a <- rbinom_mod(probs_low) #generate candidates according to P(z|low)
      #evaluate the candidate design using P(z|low)/P(z|high)
      high_prob <- pbinom_mod(a, probs_high)
      low_prob <- pbinom_mod(a, probs_low)
      val <- ifelse(high_prob < 1e-6, 1000000, #for numerical stability in cases that P(z|high) is too low
                    low_prob/high_prob) 
      if ((j == 1) | (val > pool_best_val)) { #maximize
        pool_best_design <- a #update the best design so far
        pool_best_val <- val #update the best value of P(z|low)/P(z|high)
      }
    } #note that we avoided the heavy calculation of design criterion here. So in total we can consider more possible designs
    
    #add in the designs that maximize the expected improvement
    designs <- cbind(designs, pool_best_design)
    val <- do.call(crit.func, c(list(assignment = pool_best_design), control.crit))$value #approximate true design criterion value of the new design using monte carlo
    designs_val <- c(designs_val, val)
    
    #update ustar (threshold value)
    if (i %% threshold_step == 0) {
      ustar <- quantile(designs_val, probs = threshold) #current threshold
    }
    
    #conditions
    if (cnt < ndesign) {
      
      cnt <- cnt + 1
      
      #check for best of the batch
      if (val < best_val) { 
        #smaller value of criterion is better
        best_design <- pool_best_design
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
      
      cnt_ahead <- cnt_ahead + 1
      
      if (val < best_val) {
        
        best_design <- pool_best_design
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
  list(best_val = best_val, 
       best_design = best_design, 
       time = hms_span(start, end), 
       name = "tree-parzen", trace_val = designs_val)
}

#' Bayesian optimization via a local search with a surrogate model
#' 
#' @description This function finds a good design using bayesian optimization via a local search with a surrogate model
#' @param N the number of experimental units
#' @param crit.func the function to calculate design criterion
#' @param control.crit the list of parameters to be passed to the crit.func() function
#' @param timed_search whether the duration of the algorithm is printed. Default is TRUE
#' @param seed_search seed number given at the beginning of the algorithm. Default is NULL, which means "no seed is given"
#' @param ndesign number of designs to be considered in the search
#' @param lookahead when "ndesign" designs are considered, we can choose to continue the search by looking at "lookahead" more designs. If a better design is found, the algorithm will continue to search for another "lookahead" design, and so on
#' @param npop the number of candidate designs to generate in each iteration
#' @param ninitial the number of randomly generated initial designs. These designs are used to train the surrogate (network ensemble) model for the design criterion function 
#' @param nupdate the number of iterations per which the surrogate model is retrained
#' @param nensemble the number of feed forward neural networks in the surrogate (network ensemble)
#' @param predict_type the type of acquisition function to be passed to function predict.ensemble(). Possible values are "ucb", "thompson", "postmean". Default is "ucb"
#' @param surrogate_args the list of parameters to construct the ensemble of neural networks. These parameters will be passed to function ensemble_ffnn()
#' @param train_surrogate_args the list of parameters to train the ensemble. These will be passed to function train.ensemble()
#' @return an object that contains the best design found by the algorithm, the design criterion at the best design, the running time, and other information about the search   
bayesopt_surrogate <- function(N, crit.func, control.crit = list(), 
                               timed_search = TRUE, seed_search = NULL, 
                               ndesign = 5000, lookahead = 0, npop = 100, 
                               ninitial = 1000, nupdate = 500, 
                               nensemble = 10, predict_type = "ucb",
                               surrogate_args = list(layers = c(N, 4, 1), 
                                                     actFUN = c("identity", rep("sigmoid", 1), "identity"), 
                                                     lossFUN = "mse", learning_rates = rep(0.01, 2)), 
                               train_surrogate_args = list(epoch = 100, decay = 0.1, batchsize = 25)) {
  
  #initialize counting
  cnt <- 0
  cnt_ahead <- 0
  is.continue <- TRUE
  progress <- 0
  
  #set seed if specified
  if (!is.null(seed_search)) {
    set.seed(seed_search)
  }
  
  #generate an initial population
  cat("Searching... ")
  start <- Sys.time()
  designs <- c()
  designs_val <- c()
  for (i in 1:ninitial) {
    a <- generate_random_design(N) #randomly generate a design
    val <- do.call(crit.func, c(list(assignment = a), control.crit))$value #approximate the design criterion using monte carlo
    designs <- cbind(designs, a) #add to the design collections
    designs_val <- c(designs_val, val)
  }
  best_val <- min(designs_val) #current best design criterion
  best_design <- designs[,which.min(designs_val)] #current best design
  current_design <- best_design #current design
  cnt <- cnt + ninitial
  
  #train an ensemble of neural network to represent the blackbox (design criterion) function
  cnt_update <- 0
  surrogate <- ensemble_ffnn(ffnn_args = surrogate_args, rep = nensemble)
  surrogate <- do.call(train.ensemble, c(list(obj = surrogate, input = t(designs), 
                                              y = designs_val), train_surrogate_args))
  
  #main algorithm
  while (is.continue) {
    
    #generate a pool of candidates
    pool <- c()
    for (j in 1:npop) {
      a <- generate_neighbor_design(current_design) #generate neighbors of current designs
      pool <- rbind(pool, a)
    }
    pool_val <- predict.ensemble(surrogate, input = pool, type = predict_type) #use the surrogate the predict the design criterion values of these candidate designs
    current_design <- as.vector(pool[which.min(as.vector(pool_val)),]) #choose the candidate that has the smallest predicted value of design criterion
    
    #add in the designs
    designs <- cbind(designs, current_design)
    val <- do.call(crit.func, c(list(assignment = current_design), control.crit))$value #approximate true design criterion value of the new design using monte carlo
    designs_val <- c(designs_val, val)
    
    #retrain the surrogate when more data (designs, design criteria) are collected
    cnt_update <- cnt_update + 1
    if (cnt_update == nupdate) {
      cnt_update <- 0
      surrogate <- ensemble_ffnn(ffnn_args = surrogate_args, rep = nensemble)
      surrogate <- do.call(train.ensemble, c(list(obj = surrogate, input = t(designs), 
                                                  y = designs_val), train_surrogate_args))
    }
    
    #conditions
    if (cnt < ndesign) {
      
      cnt <- cnt + 1
      
      #check for best of the batch
      if (val < best_val) { 
        #smaller value of criterion is better
        best_design <- current_design
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
      
      cnt_ahead <- cnt_ahead + 1
      
      if (val < best_val) {
        
        best_design <- current_design
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
  list(best_val = best_val, 
       best_design = best_design, 
       time = hms_span(start, end), 
       name = "deep surrogate", trace_val = designs_val)
}

#' Bayesian optimization via reinforcement learning
#' 
#' @description This function finds a good design using bayesian optimization via a surrogate model and a policy network
#' @param N the number of experimental units
#' @param crit.func the function to calculate design criterion
#' @param control.crit the list of parameters to be passed to the crit.func() function
#' @param timed_search whether the duration of the algorithm is printed. Default is TRUE
#' @param seed_search seed number given at the beginning of the algorithm. Default is NULL, which means "no seed is given"
#' @param ndesign number of designs to be considered in the search
#' @param lookahead when "ndesign" designs are considered, we can choose to continue the search by looking at "lookahead" more designs. If a better design is found, the algorithm will continue to search for another "lookahead" design, and so on
#' @param npop the number of designs in the population, which is maintained in each iteration
#' @param ninitial the number of randomly generated initial designs. These designs are used to train the surrogate (network ensemble) model for the design criterion function 
#' @param nupdate the number of iterations per which the surrogate model is retrained
#' @param nensemble the number of feed forward neural networks in the surrogate (network ensemble)
#' @param predict_type the type of acquisition function to be passed to function predict.ensemble(). Possible values are "ucb", "thompson", "postmean". Default is "ucb"
#' @param surrogate_args the list of parameters to construct the ensemble of neural networks. These parameters will be passed to function ensemble_ffnn()
#' @param train_surrogate_args the list of parameters to train the ensemble. These will be passed to function train.ensemble()
#' @param epoch the number of epochs to train the policy network
#' @param policy_args the list of parameters to construct the policy neural network. These parameters will be passed to function ffnn()
#' @return an object that contains the best design found by the algorithm, the design criterion at the best design, the running time, and other information about the search   
bayesopt_surrogate_policy <- function(N, crit.func, control.crit = list(),
                                      timed_search = TRUE, seed_search = NULL,
                                      ndesign = 5000, lookahead = 0, npop = 100,
                                      ninitial = 1000, nupdate = 500,
                                      nensemble = 10, predict_type = "ucb", 
                                      surrogate_args = list(layers = c(N, 4, 1),
                                                            actFUN = c("identity", rep("sigmoid", 1), "identity"),
                                                            lossFUN = "mse", learning_rates = rep(0.01, 2)),
                                      train_surrogate_args = list(epoch = 30, decay = 0.1, batchsize = 25),
                                      epoch = 100, 
                                      policy_args = list(layers = c(N, 32, 8, N),
                                                         actFUN = c("identity", rep("sigmoid", 2), "sigmoid"),
                                                         lossFUN = "bernoulli",  learning_rates = rep(0.01, 3))) {
  
  #initialize counting
  cnt <- 0
  cnt_ahead <- 0
  is.continue <- TRUE
  progress <- 0
  
  #set seed if specified
  if (!is.null(seed_search)) {
    set.seed(seed_search)
  }
  
  #generate an initial population
  cat("Searching... ")
  start <- Sys.time()
  designs <- c()
  designs_val <- c()
  for (i in 1:ninitial) {
    a <- generate_random_design(N) #randomly generate a design 
    val <- do.call(crit.func, c(list(assignment = a), control.crit))$value #approximate true design criterion value using monte carlo
    designs <- cbind(designs, a) #add the design to the design collection
    designs_val <- c(designs_val, val)
  }
  ord <- order(designs_val)
  pop <- designs[,ord[1:npop]] #use the best npop designs as the current population
  pop_val <- designs_val[ord[1:npop]] #design criterion values of the current population
  best_design <- pop[,1] #current best design
  best_val <- pop_val[1] #current best design criterion
  cnt <- cnt + ninitial
  
  #train an ensemble of neural network to represent the blackbox design criterion function
  cnt_update <- 0
  surrogate <- ensemble_ffnn(ffnn_args = surrogate_args, rep = nensemble)
  surrogate <- do.call(train.ensemble, c(list(obj = surrogate, input = t(designs),
                                              y = designs_val), train_surrogate_args))
  
  #initiate a neural network for the policy
  policy <- do.call(ffnn, policy_args)
  
  #main algorithm
  while (is.continue) {
    
    #generate a new pool of candidates based on the policy network
    #and train the policy network at the same time
    pool <- pop
    pool_val <- as.vector(predict.ensemble(surrogate, input = t(as.matrix(pool)), type = predict_type)) #use the surrogate to predict the design criterion values of the pool
    for (t in 1:epoch) { #train the policy network for multiple epochs
      gradients <- as.list(rep(0, length(policy$weights))) #initialize gradients for the epoch
      pop_val <- c() #initialize vector of design criterion values estimated by the surrogate
      for (i in 1:npop) { #consider all designs in the population
        
        s <- pop[,i] #current design
        
        #sample the edit to the current design (action)
        ypred <- as.vector(predict.ffnn(policy, t(as.matrix(s))))
        changes <- rbinom_mod(ypred)
        idx <- which(changes == 1)
        snew <- s
        snew[idx] <- 1 - s[idx] #apply the edit to the design (these are the "actions") to get a new design
        
        #evaluate the reward
        rw_old <- as.vector(predict.ensemble(surrogate, input = t(as.matrix(s)), type = predict_type))
        rw_new <- as.vector(predict.ensemble(surrogate, input = t(as.matrix(snew)), type = predict_type))
        rw <- rw_old - rw_new #reward obtained by editing the old design into the new design
        
        gradients <- accumulate_weight.ffnn(gradients1 = gradients, 
                                            gradients2 = gradient.ffnn(policy, input = t(as.matrix(s)), y = snew)$gradients, 
                                            factor = rw/npop) #contribute gradients from this observation to the epoch's gradients
        #update population
        pop[,i] <- snew #add the newly created design to the population and the pool
        pool <- cbind(pool, snew)
        pop_val <- c(pop_val, rw_new) #add the estimated design criterion using surrogate
        pool_val <- c(pool_val, rw_new)
      }
      policy <- update.ffnn(policy, gradients = gradients) #update the weights of the policy network after each epoch
    }
    
    #choose the best designs (evaluated using the surrogate) from the whole pool (populations * epochs) into the next population
    ord <- order(pool_val)
    pop <- pool[,ord[1:npop]] #next population
    pop_val <- c() #initialize the design criterion values for the next population
    
    #evaluate the population
    for (i in 1:npop) {
      val <- do.call(crit.func, c(list(assignment = pop[,i]), control.crit))$value #approximate true design criterion value using monte carlo
      designs <- cbind(designs, pop[,i]) #add this population to the design collection 
      designs_val <- c(designs_val, val)
      pop_val <- c(pop_val, val) #add in the new population's design criterion values
    }
    
    #retrain the surrogate when more data (designs, design criteria) are collected
    cnt_update <- cnt_update + npop
    if (cnt_update >= nupdate) {
      cnt_update <- 0
      surrogate <- ensemble_ffnn(ffnn_args = surrogate_args, rep = nensemble)
      surrogate <- do.call(train.ensemble, c(list(obj = surrogate, input = t(designs),
                                                  y = designs_val), train_surrogate_args))
    }
    
    #conditions
    if (cnt < ndesign) {
      
      cnt <- cnt + npop
      
      #check for best of the population
      idx <- which.min(pop_val)
      if (pop_val[idx] < best_val) {
        #smaller value of criterion is better
        best_design <- pop[,idx]
        best_val <- pop_val[idx]
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
      
      idx <- which.min(pop_val)
      if (pop_val[idx] < best_val) {
        #smaller value of criterion is better
        best_design <- pop[,idx]
        best_val <- pop_val[idx]
        
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
  list(best_val = best_val,
       best_design = best_design,
       time = hms_span(start, end),
       name = "deep surrogate & policy",
       trace_val = designs_val)
}

#~~~~~~~~~~~~~~~~~~~~~~
#Save the workspace
#~~~~~~~~~~~~~~~~~~~~~~~

save(rbinom_mod, pbinom_mod, 
     tree_parzen, bayesopt_surrogate, 
     bayesopt_surrogate_policy, 
     file = "funcs/bayesopt_funcs.RData")

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
load("funcs/neuralnet.RData")
load("funcs/bayesopt_funcs.RData")

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

#tree-parzen
check_tree <- tree_parzen(N = g$N, crit.func = mc_gane, 
                          control.crit = list(g = g, 
                                              prior = gane_prior_now,
                                              nsim = 1)) #takes 20 minutes
plot(1:length(check_tree$trace_val), check_tree$trace_val, type = "l", 
     xlab = "iteration", ylab = "design criterion", 
     main = "tree-parzen results by iteration") #plot the design criteria of designs visited in each iteration

#local search with surrogate
check_surrogate <- bayesopt_surrogate(N = g$N, crit.func = mc_gane, 
                                      control.crit = list(g = g, 
                                                          prior = gane_prior_now,
                                                          nsim = 1)) #takes 12 hours
plot(1:length(check_surrogate$trace_val), check_surrogate$trace_val, type = "l", 
     xlab = "iteration", ylab = "design criterion", 
     main = "local search with surrogate results by iteration") #plot the design criteria of designs visited in each iteration

#reinforcement learning
check_reinforce <- bayesopt_surrogate_policy(N = g$N, crit.func = mc_gane, 
                                             control.crit = list(g = g, 
                                                                 prior = gane_prior_now,
                                                                 nsim = 1)) #takes 7 hours
plot(1:length(check_reinforce$trace_val), check_reinforce$trace_val, type = "l", 
     xlab = "iteration", ylab = "design criterion", 
     main = "reinforcement learning results by iteration") #plot the design criteria of designs visited in each iteration

load("results/checks/checks_algo_runs.RData")
save(check_rs, check_ts, check_sa, check_ga,
     check_tree, check_surrogate, check_reinforce, 
     file = "results/checks/checks_algo_runs.RData")

#NOTE: If you want to quickly plot the results without waiting for algorithms to run
#use the command below and run the plot commands above
load("saved results/checks/checks_algo_runs.RData")
