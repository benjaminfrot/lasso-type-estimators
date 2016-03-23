require(matrixcalc)
#####
# This file contains an implementation of a proximal gradient
# descent method to solve 
# (http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003420)
# (http://www.jmlr.org/proceedings/papers/v28/wytock13.pdf) :
# We implement : 
# - Nesterov Y (2005) Smooth minimization of non-smooth functions. 
# Math Program 103: 127–152. doi: 10.1007/s10107-004-0552-5
# In order to avoid the typical "riples" that come with that method
# we also implement : 
# - http://statweb.stanford.edu/~candes/papers/adap_restart_paper.pdf

nesterov.scggm.compute.gradient <- function(SigmaZ, SigmaZX, SigmaX, N, y) {
  # Compute the gradient at y given Sigma (the unscaled covmat) and
  # N (the sample size)
  
  SX <- y$SX
  SZX <- y$SZX
  SX <- 0.5 * (SX + t(SX))
  
  # Check that SX is positive semi-definite
  is.pos.sem.def <- is.positive.definite(SX)
  if (!is.pos.sem.def) {
    return(NULL)
  }
  
  R <- chol(SX)
  logdet <- 2 * sum(log(diag(R)))
  iSX <- matrix.inverse(SX) # FIXME: Should use the chol dec. 
  txyityy <- SZX %*% iSX
  xtxth <- SigmaZ %*% txyityy
  txytxtxth <- t(SZX) %*% SigmaZ %*% txyityy
   
  l1 <- matrix.trace(SX %*% SigmaX)
  l2 <- matrix.trace(SigmaZX %*% t(SZX))
  l3 <- matrix.trace(txytxtxth)
  obj <- 0.5 * l1 + l2 + 0.5 * l3 - 0.5 * N * logdet
  obj <- obj / N
  
  dSZX <- (SigmaZX + xtxth) / N
  dSX <- 0.5 * (SigmaX - N * iSX - iSX %*% txytxtxth) / N
  dSX <- 0.5 * (dSX + t(dSX))
  
  result <- list()
  result$obj <- obj
  result$dSX <- dSX
  result$dSZX <- dSZX
  
  result
}

nesterov.scggm.lasso <- function(Z, X, l1, l2, init, tol = 1e-08, max.iter = 1000) {
  
  N <- dim(X)[1]
  p <- dim(X)[2]
  q <- dim(Z)[2]
  SigmaX <- t(X) %*% X
  SigmaZ <- t(Z) %*% Z
  SigmaZX <- t(Z) %*% X
  
  SX <- init$SX
  SZX <- init$SZX
  
  x <- list()
  y <- list()
  z <- list()
  x$SX <- y$SX <- z$SX <- SX
  x$SZX <- y$SZX <- z$SZX <- SZX
  
  l <- 1
  tk <- 0.66
  eta <- 1.5
  nTries <- 300
  # To store the value of the objective function at the prev. iter. :
  objs <- c()
  # Minimum number of iterations (to avoid being caught in the first "riples" of Nesterov's method)
  minIter <- 50
  
  for (i in 1:max.iter) {
    y$SX <- (1 - tk) * x$SX + tk * z$SX
    y$SZX <- (1 - tk) * x$SZX + tk * z$SZX
    
    y.r <- nesterov.scggm.compute.gradient(SigmaZ, SigmaZX, SigmaX, N = N, y = y)
    for (j in 1:nTries) {
      nz <- list()
      ny <- list()
      nz$SX <- z$SX - (1 / (l * tk)) * y.r$dSX
      nz$SZX <- z$SZX - (1 / (l * tk)) * y.r$dSZX
      
      ny$SX <- z$SX - (1 / (l)) * y.r$dSX
      ny$SZX <- z$SZX - (1 / (l)) * y.r$dSZX
      
      # Compute the proximal operator of the l1-norm
      X1 <- 0.5 * (nz$SX + t(nz$SX))
      X2 <- abs(X1) - 2 * l1 / (l * tk)
      X2[X2 < 0] <- 0
      nz$SX <- X2 * sign(X1)
      
      X1 <- 0.5 * (ny$SX + t(ny$SX))
      X2 <- abs(X1) - 2 * l1 / (l)
      X2[X2 < 0] <- 0
      ny$SX <- X2 * sign(X1)
      
      X1 <- nz$SZX
      X2 <- abs(X1) - 2 * l2 / (l * tk)
      X2[X2 < 0] <- 0
      nz$SZX <- X2 * sign(X1)
      
      X1 <- ny$SZX
      X2 <- abs(X1) - 2 * l2 / (l)
      X2[X2 < 0] <- 0
      ny$SZX <- X2 * sign(X1)
      
      
      # Compute the gradient at ny
      ny.r <- nesterov.scggm.compute.gradient(SigmaZ, SigmaZX, SigmaX, N, ny)
      if (is.null(ny.r)) {
        isposdef <- FALSE
      } else {
        isposdef <- is.positive.definite(nz$SX)
      }
      
      if (isposdef) {
        diffSX <- ny$SX - y$SX
        diffSZX <- ny$SZX - y$SZX
        scProd <- sum(diffSX * y.r$dSX) + sum(diffSZX * y.r$dSZX)
        RHS <- y.r$obj + scProd + 0.5 * l * (sum((diffSX)**2) + sum((diffSZX)**2))
        if (ny.r$obj <= RHS + tol) {
          x <- ny
          z <- nz
          break()
        }  
      }
      l <- l * eta
    }
    
    r <- nesterov.scggm.compute.gradient(SigmaZ, SigmaZX, SigmaX, N, x)
    pobj <- r$obj + l1 * sum(abs(x$SX)) + l2 * sum(abs(x$SZX))
    tk <- (sqrt(tk**4 + 4 * tk ** 2) - tk**2) / 2
    objs <- c(objs, pobj)
    
    if(length(objs) > 1) {
      diff <- abs(objs[length(objs)] - objs[length(objs)-1])
      if(diff < tol) {
        if(i > minIter)
          break()
      }
      if (pobj > objs[length(objs)-1]) {
        tk <- 0.5 # Kill the momentum (see : http://statweb.stanford.edu/~candes/papers/adap_restart_paper.pdf)
      }
    }
  }
  x$obj <- min(objs)
  x$objs <- objs
  x$diff <- diff
  x$iter <- i
  
  x
}
