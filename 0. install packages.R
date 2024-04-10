#######################################################
#Install packages for the code run
#######################################################

#-------------------
#Help functions
#-------------------

#' Install undownloaded packages 
#' 
#' @description This function checks if a package is installed. If not, then install it from cran
#' @param namevec a vector of names of the packages to be checked/installed
check.install <- function(namevec) {
  for (i in 1:length(namevec)) {
    name <- namevec[i]
    if (!require(name, character.only = TRUE)) {
      cat("Installing", name, "... \n")
      install.packages(paste0(name))
    } else {
      cat("Package", name, "has already been installed \n")
    }
  }
}

#-------------------
#Install packages needed for code run
#-------------------

# #install all packages at once
# check.install(c("igraph", "igraphdata", "MASS", "Matrix", 
#                 "Rcpp", "RcppEigen", "parallel", "neuralnet", 
#                 "ggplot2", "ggpubr", "dplyr", "chron", "ggh4x"))

#one-by-one installing in case some pop-up appears, for example, to requiring restart
check.install("igraph")
check.install("igraphdata")
check.install("MASS")
check.install("Matrix")
check.install("Rcpp")
check.install("RcppEigen")
check.install("parallel")
check.install("neuralnet")
check.install("ggplot2")
check.install("ggpubr")
check.install("dplyr")
check.install("chron")
check.install("ggh4x")

#-------------------
#Set working directory to the current file
#-------------------

path <- rstudioapi::getSourceEditorContext()$path
path <- substr(path, start = 1, 
               stop = gregexpr("/", path)[[1]][length(gregexpr("/", path)[[1]])])
setwd(path)

#-------------------
#Check if Rcpp is properly installed and ready to use
#-------------------

library(Rcpp)
sourceCpp("rcpp_help_funcs.cpp")

#if not work use this
source("rcpp_help_funcs.R")
