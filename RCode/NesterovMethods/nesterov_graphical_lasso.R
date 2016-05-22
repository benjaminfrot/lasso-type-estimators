require(matrixcalc)
#####
# This file contains an implementation of a proximal gradient
# descent method to solve (http://statweb.stanford.edu/~tibs/ftp/graph.pdf):
# argmin_{S} -logdet{S} + Trace{S.Sigma}
# l1 ||S||_1.
# with S > 0.
# We implement : 
# - Nesterov Y (2005) Smooth minimization of non-smooth functions. 
# Math Program 103: 127–152. doi: 10.1007/s10107-004-0552-5
# In order to avoid the typical "riples" that come with that method
# we also implement : 
# - http://statweb.stanford.edu/~candes/papers/adap_restart_paper.pdf

nesterov.graphical.lasso.compute.gradient <- function(Sigma, N, y) {
  # Compute the gradient at y given Sigma (the unscaled covmat) and
  # N (the sample size)
  
  S <- y$S
  p <- dim(Sigma)[1]
  
  S <- 0.5 * (S + t(S))
  
  # Check that S - L is positive semi-definite
  is.pos.sem.def <- is.positive.definite(S)
  if (!is.pos.sem.def) {
    return(NULL)
  }
  
  R <- chol(S)
  logdet <- 2 * sum(log(diag(R)))
  iSmL <- solve(S) # FIXME: Should use the chol dec. 
  
  dL <- N * iSmL - Sigma
  dL <- dL / N
  dL <- 0.5 * (dL + t(dL))
  dS <- -dL
  
  obj <- 0.5 * sum(diag((S) %*% Sigma)) - 0.5 * N * logdet
  obj <- obj / N
  
  result <- list()
  result$obj <- obj
  result$dS <- dS
  
  result
}

nesterov.graphical.lasso <- function(X, l1, init, tol = 1e-08, max.iter = 1000) {
  N <- dim(X)[1]
  p <- dim(X)[2]
  Sigma <- t(X) %*% X
  
  S <- init$S
  
  x <- list()
  y <- list()
  z <- list()
  x$S <- y$S <- z$S <- S
  
  l <- 1
  tk <- 0.66
  eta <- 1.5
  nTries <- 300
  # To store the value of the objective function at the prev. iter. :
  objs <- c()
  # Minimum number of iterations
  minIter <- 5
  
  for (i in 1:max.iter) {
    y$S <- (1 - tk) * x$S + tk * z$S
    
    y.r <- nesterov.graphical.lasso.compute.gradient(Sigma = Sigma, N = N, y = y)
    for (j in 1:nTries) {
      nz <- list()
      ny <- list()
      nz$S <- z$S - (1 / (l * tk)) * y.r$dS
      
      ny$S <- z$S - (1 / (l)) * y.r$dS
      
      # Compute the proximal operator of the l1-norm
      X1 <- 0.5 * (nz$S + t(nz$S))
      X2 <- abs(X1) - 2 * l1 / (l * tk)
      X2[X2 < 0] <- 0
      nz$S <- X2 * sign(X1)
      
      X1 <- 0.5 * (ny$S + t(ny$S))
      X2 <- abs(X1) - 2 * l1 / (l)
      X2[X2 < 0] <- 0
      ny$S <- X2 * sign(X1)
      
      # Compute the gradient at ny
      ny.r <- nesterov.graphical.lasso.compute.gradient(Sigma, N, ny)
      if (is.null(ny.r)) {
        isposdef <- FALSE
      } else {
        isposdef <- is.positive.definite(nz$S)
      }
      
      if (isposdef) {
        diffS <- ny$S - y$S
        scProd <- sum(diffS * y.r$dS)
        RHS <- y.r$obj + scProd + 0.5 * l * (sum((diffS)**2))
        if (ny.r$obj <= RHS + tol) {
          x <- ny
          z <- nz
          break()
        }  
      }
      l <- l * eta
    }
    
    r <- nesterov.graphical.lasso.compute.gradient(Sigma, N, x)
    pobj <- r$obj + l1 * sum(abs(x$S))
    tk <- (sqrt(tk**4 + 4 * tk ** 2) - tk**2) / 2
    objs <- c(objs, pobj)
    
    if(length(objs) > 1) {
      diff <- abs(objs[length(objs)] - objs[length(objs)-1])
      if(diff < tol) {
        if(i > minIter)
          break()
      }
#      if (pobj > objs[length(objs)-1]) {
#        tk <- 0.5 # Kill the momentum (see : http://statweb.stanford.edu/~candes/papers/adap_restart_paper.pdf)
#      }
    }
  }
  x$obj <- min(objs)
  x$objs <- objs
  x$diff <- diff
  x$iter <- i
  
  x
}
