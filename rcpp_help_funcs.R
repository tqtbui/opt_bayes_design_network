#######################################################
#This file contains core functions to calculate target design criterion phi_0
#######################################################

#' Target design criterion for the gane (POW-DEG) model
#' 
#' @description This function calculate the target design criterion for the gane (POW-DEG) model
#' @param m the model matrix, denoted by F in the paper
#' @param m_der the derivative of the model matrix, denoted by F' in the paper
#' @param beta beta parameter
#' @param netfunc the third element of the d vector in the paper
#' @param netfun_der the fifth element of the d vector in the paper
#' @param N number of units
#' @param omega sigma squared
#' @return the value of the target design criterion for the gane (POW-DEG) model 
cpp_pow_crit <- function(m, m_der, beta, netfunc, netfunc_der, N, omega) {
  
  fish <- matrix(0, nrow = 6, ncol = 6)
  fish[1:4,1:4] <- (1/omega) * t(m) %*% m
  fish[6,6] <- N / (2*omega^2)
  fish[5,5] <- (1/omega) * t(as.matrix(beta)) %*% t(m_der) %*% m_der %*% as.matrix(beta)
  fish[1:4,5] <- (1/omega) * t(m) %*% m_der %*% as.matrix(beta)
  fish[5,1:4] <- fish[1:4,5]

  hess <- try(solve(fish), silent = TRUE)
  if ("try-error" %in% class(hess)) {
    hess <- ginv(hess) #use pseudo inverse when not invertible
  }
  
  dtotal <- rep(0, 5) 
  dtotal[2] <- 1
  dtotal[3] <- netfunc
  dtotal[4] <- -netfunc
  dtotal[5] <- (beta[3] - beta[4])*netfunc_der

  t(as.matrix(dtotal)) %*% hess[1:5, 1:5] %*% as.matrix(dtotal)
}

#' Target design criterion for the locopt (CNAR) model
#' 
#' @description This function calculate the target design criterion for the locopt (CNAR) model
#' @param D the diagonal matrix of units' degrees, denoted by K in the paper
#' @param W the adjacency matrix, denoted by A in the paper
#' @param rho rho parameter
#' @param X the model matrix, denoted by F in the paper
#' @return the value of the target design criterion for the locopt (CNAR) model 
cpp_locopt_crit <- function(D, W, rho, x) {
  
  fish <- t(x) %*% D %*% as.matrix(x) - rho * t(x) %*% W %*% x
  hess <- try(solve(fish), silent = TRUE)
  if ("try-error" %in% class(hess)) {
    hess <- ginv(hess) #use pseudo inverse when not invertible
  }
  
  hess[2,2]
}

#' Target design criterion for the basse (NS) model
#' 
#' @description This function calculate the target design criterion for the basse (NS) model
#' @param AA the adjacency matrix, denoted by A in the paper
#' @param d the vector of units' degree
#' @param N number of units
#' @param mu mu parameter
#' @param gamma gamma parameter
#' @param sigma sigma parameter
#' @param z1 the treatment assignment vector
#' @return the value of the target design criterion for the basse (NS) model 
cpp_basse_crit <- function(AA, d, N, mu, gamma, sigma, z1) {
  
  n1 <- sum(z1)
  z0 <- 1 - z1
  n0 <- sum(z0)
  AAz1 <- AA %*% as.matrix(z1)
  
  bias <- (1/n1) * sum(z1*d) - (1/n0) * sum(z0*d)
  var1 <- (1/n1) + (1/n0)
  var2 <- (1/n1^2) * sum(z1*AAz1) + (1/n0^2) * sum(t(as.matrix(z0))%*%AA%*%as.matrix(z0)) - (2/(n1*n0)) * sum(z0*AAz1)
  
  mu^2 * bias^2 + gamma^2 * var1 + sigma^2 * var2;
}
