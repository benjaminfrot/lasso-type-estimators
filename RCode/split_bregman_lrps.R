# An implementation of http://arxiv.org/abs/1110.3076

library(matrixcalc)

split.bregman.low.rank.plus.sparse.update.parameters <- function(ps, opts) {
  lp1 <- opts$lp1
  lp2 <- opts$lp2
  mu <- opts$mu
  C <- opts$Sigma
  A <- ps$A
  S <- ps$S
  L <- ps$L
  U <- ps$U

  # Update A  
  X1 <- mu * (S - L) - C - U
  X2 <- X1 %*% X1 + 4 * mu * diag(dim(X1)[1])
  X2 <- 0.5 * (X2 + t(X2))
  eig <- eigen(X2)
  sqrtX2 <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  A <- (X1 + sqrtX2) / (2 * mu)
  A <- 0.5 * (A + t(A))
  
  # Update S
  X1 <- A + L + (U / mu)
  X2 <- abs(X1) - lp1
  X2[X2 <= 0] <- 0
  S = X2 * sign(X1)
  S <- 0.5 * (S + t(S))
  
  # Update L
  X1 <- S - A - (U / mu)
  X1 <- 0.5 * (X1 + t(X1))
  eig <- eigen(X1)
  eigVal <- eig$values - lp2
  eigVal[eigVal < 0] <- 0
  L <- eig$vectors %*% diag(eigVal) %*% t(eig$vectors)
  L <- 0.5 * (L + t(L))
  
  #Update U
  U <- U + mu * (A - S + L)
  U <- 0.5 * (U + t(U))
  
  ps$A <- A
  ps$U <- U
  ps$S <- S
  ps$L <- L
  
  ps
}

split.bregman.low.rank.plus.sparse <- function(Sigma, Lambda1, Lambda2, init=NULL, maxiter=3000, mu=0.05, tol=1e-07) {
  
  # Sigma : Covariance Matrix
  # Lambda 1 : Penalty on S
  # Lambda 2 : Penalty on L
  
  options <- list()
  options$mu <- mu
  options$Sigma <- Sigma
  options$lp1 <- Lambda1 / mu
  options$lp2 <- Lambda2 / mu
  
  p <- dim(Sigma)[1]
  
  if (is.null(init)) {
    S <- diag(p)
    L <- S * 0.01
    A <- S - L + diag(rep(1, p))
    U <- diag(rep(1, p))
    parameters <- list(S=S, L=L, A=A, U=U)
  } else {
    parameters = init
  }
  
  diffs <- c()  
  for (i in 1:maxiter) {
    new_parameters <- split.bregman.low.rank.plus.sparse.update.parameters(parameters, options)
    
    # Check convergence of S, L
    if (frobenius.norm(parameters$L) > 0) {
      diff <- frobenius.norm(new_parameters$S - parameters$S) / frobenius.norm(parameters$S) + 
        frobenius.norm(new_parameters$L - parameters$L) / frobenius.norm(parameters$L)
    } else {
      diff <- frobenius.norm(new_parameters$S - parameters$S) / frobenius.norm(parameters$S) + 
        frobenius.norm(new_parameters$L - parameters$L)
    }
    diffs <- c(diffs, diff)
    parameters <- new_parameters
    if (diff < tol) {
      break()
    }
  }
  
  parameters$iter <-i 
  parameters$diffs <- diffs
  
  parameters
}