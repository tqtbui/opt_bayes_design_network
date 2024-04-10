#######################################################
#Functions for Monte Carlo evaluation
#######################################################

#' Print elapsed time
#' 
#' @description This function prints the elapsed time in hh:mm:ss format
#' @param start start (system) time
#' @param end end (system) time
#' @return a character string of the form hh:mm:ss
hms_span <- function(start, end) {
  dsec <- as.numeric(difftime(end, start, unit = "secs"))
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- dsec - 3600*hours - 60*minutes
  paste0(
    sapply(c(hours, minutes, seconds), function(x) {
      formatC(x, width = 2, format = "d", flag = "0")
    }), collapse = ":")
}

#' Binomial design
#' 
#' @description This function randomly generates a treatment assignment on n units
#' @param n number of experimental units
#' @param p_assign percentage of units assigned to treatment. p_assign = 0.5 corresponds to the design that randomly assigns half of the units to treatment and the other half to control
#' @return a vector containing treatment assignment (1 for treatment and 0 for control) for n units
assign_binomial <- function(n, p_assign = 0.5) {
  assign <- rep(0, n)
  samp <- sample(1:n, floor(n*p_assign), replace = FALSE)
  assign[samp] <- 1
  assign
}

#' Cumulative means
#' 
#' @description This function calculates the cumulative means for a vector
#' @param x the vector whose cumulative means are to be calculated
#' @param trim the fraction of observations to be trimmed on each end of x before the mean is computed. More info on base R's mean() function's documentation
#' @param ... arguments to be passed to the base R's mean() function 
#' @return a vector of cumulative means 
cummean <- function(x, trim = 0, ...) {
  sapply(1:length(x), function(i) {
    mean(x[1:i], trim = trim, ...)
  })
}

#' Trace of Monte Carlo simulation paths
#' 
#' @description This function runs Monte Carlo simulation and trace the Monte Carlo estimate as iteration increases
#' @param mc.func the function that generates Monte Carlo samples
#' @param control.mc the list of parameters that are passed to the mc.func function
#' @param paths the number of Monte Carlo simulations to be conducted
#' @param seed the seed number of the whole call of mc_trace function
#' @param timed whether the time elapsed should be printed
#' @param ... arguments to be passed to the cummean function that calculate the Monte Carlo estimate
#' @return a matrix, each column of which represents a path of Monte Carlo estimates
mc_trace <- function(mc.func, control.mc = list(), 
                     paths = 10, seed = NULL, timed = TRUE, ...) {
  
  #initialization
  ret <- c()
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  #monte carlo paths
  start <- Sys.time()
  for (i in 1:paths) {
    
    tmp <- do.call(mc.func, control.mc)
    ret <- cbind(ret, cummean(tmp$trace, ...))
    
    #update progress
    if (paths >= 10) {
      if (i%% (paths%/%10) == 0) {
        cat(i/(paths%/%10)*10, "%... ")
      }
    }
  }
  end <- Sys.time()
  if (timed) {
    cat("Simulation ends in", hms_span(start, end), "\n")
  }
  
  #return
  ret
}

#######################################################
#Data preprocessing
#######################################################

library(igraph)

#-----------------
#enron data
#-----------------

library(igraphdata) #enron network data taken from the igraph data package
data("enron")
g_enron <- as.undirected(enron) #make the network undirected
g_enron <- simplify(g_enron) #make the network simplified
rm(enron) #remove the "enron" object from the environement
length(V(g_enron)) #number of nodes
length(E(g_enron)) #number of edges

#-----------------
#Caltech network data 
#-----------------

#fb data Caltech
#the raw data file was downloaded from network repository https://networkrepository.com/socfb.php
#the current data file was modified from the raw data file by removing the first row
g_fb_caltech <- read_graph("data/socfb-Caltech36-mod.mtx", 
                           format = "edgelist") 
g_fb_caltech <- as.undirected(g_fb_caltech) #make the network undirected
g_fb_caltech <- simplify(g_fb_caltech) #make the network simplified
length(V(g_fb_caltech)) #number of nodes
length(E(g_fb_caltech)) #number of edges

#-----------------
#UMich network data 
#-----------------

#fb data UMich
#the raw data file was downloaded from network repository https://networkrepository.com/socfb.php
#the current data file was modified from the raw data file by removing the first row
g_fb_umich <- read_graph("data/socfb-Mich67-mod.mtx", 
                         format = "edgelist") #the file was downloaded from network repository https://networkrepository.com/socfb.php
g_fb_umich <- as.undirected(g_fb_umich)#make the network undirected
g_fb_umich <- simplify(g_fb_umich) #make the network simplified
length(V(g_fb_umich)) #number of nodes
length(E(g_fb_umich)) #number of edges

#~~~~~~~~~~~~~~~~~~~~~~
#Save the workspace
#~~~~~~~~~~~~~~~~~~~~~~~

save(g_enron, g_fb_caltech, g_fb_umich, 
     hms_span, cummean, assign_binomial,
     mc_trace, file = "funcs/mc_funcs.RData")
