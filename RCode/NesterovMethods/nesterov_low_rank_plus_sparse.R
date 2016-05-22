require(matrixcalc)

#####
# This file contains an implementation of a proximal gradient
# descent method to solve (http://projecteuclid.org/euclid.aos/1351602527):
# argmin_{S, L} -logdet{S - L} + Trace{(S -L)Sigma}
# l1 ||S||_1 + l2 ||L||_\ast.
# with S - L > 0, L >= 0.
# We implement : 
# - Nesterov Y (2005) Smooth minimization of non-smooth functions. 
# Math Program 103: 127–152. doi: 10.1007/s10107-004-0552-5
# In order to avoid the typical "riples" that come with that method
# we also implement : 
# - http://statweb.stanford.edu/~candes/papers/adap_restart_paper.pdf


nesterov.low.rank.plus.sparse.compute.gradient <- function(Sigma, N, y) {
  # Compute the gradient at y given Sigma (the unscaled covmat) and
  # N (the sample size)
  
  S <- y$S
  L <- y$L # Low rank matrix containing information about the latent vars.
  p <- dim(Sigma)[1]
  
  S <- 0.5 * (S + t(S))
  L <- 0.5 * (L + t(L))
  
  # Check that S - L is positive semi-definite
  is.pos.sem.def <- is.positive.definite(S - L)
  if (!is.pos.sem.def) {
    return(NULL)
  }
  
  R <- chol(S - L)
  logdet <- 2 * sum(log(diag(R)))
  iSmL <- solve(S - L) # FIXME: Should use the chol dec. 
  
  dL <- N * iSmL - Sigma
  dL <- dL / N
  dL <- 0.5 * (dL + t(dL))
  dS <- -dL
  
  obj <- 0.5 * sum(diag((S - L) %*% Sigma)) - 0.5 * N * logdet
  obj <- obj / N
  
  result <- list()
  result$obj <- obj
  result$dL <- dL
  result$dS <- dS
  
  result
}

nesterov.low.rank.plus.sparse <- function(X, l1, l2, init, tol = 1e-05, max.iter = 1000) {
  
  N <- dim(X)[1]
  p <- dim(X)[2]
  Sigma <- t(X) %*% X
  
  S <- init$S
  L <- init$L
  
  x <- list()
  y <- list()
  z <- list()
  x$S <- y$S <- z$S <- S
  x$L <- y$L <- z$L <- L
  
  l <- 1
  tk <- 0.5
  eta <- 1.5
  nTries <- 300
  # To store the value of the objective function at the prev. iter. :
  objs <- c()
	diffs <- c()
  # Minimum number of iterations
  minIter <- 100
	oX <- x
  
  for (i in 1:max.iter) {
    y$S <- (1 - tk) * x$S + tk * z$S
    y$L <- (1 - tk) * x$L + tk * z$L
    
    y.r <- nesterov.low.rank.plus.sparse.compute.gradient(Sigma = Sigma, N = N, y = y)
    for (j in 1:nTries) {
      nz <- list()
      ny <- list()
      nz$S <- z$S - (1 / (l * tk)) * y.r$dS
      nz$L <- z$L - (1 / (l * tk)) * y.r$dL
      
      ny$S <- z$S - (1 / (l)) * y.r$dS
      ny$L <- z$L - (1 / (l)) * y.r$dL
      
      # Compute the proximal operator of the l1-norm
      X1 <- 0.5 * (nz$S + t(nz$S))
      X2 <- abs(X1) - l1 / (l * tk)
      X2[X2 < 0] <- 0
      nz$S <- X2 * sign(X1)
      
      X1 <- 0.5 * (ny$S + t(ny$S))
      X2 <- abs(X1) - l1 / (l)
      X2[X2 < 0] <- 0
      ny$S <- X2 * sign(X1)
      
      # Now compute the proximal operator of the nuclear norm
      X1 <- 0.5 * (nz$L + t(nz$L))
      eig <- eigen(X1)
      eigVal <- eig$values - l2 / (l * tk)
      eigVal[eigVal < 0] <- 0
      nz$L <- eig$vectors %*% diag(eigVal) %*% t(eig$vectors)
      nz$L <- 0.5 * (nz$L + t(nz$L))
      
      X1 <- 0.5 * (ny$L + t(ny$L))
      eig <- eigen(X1)
      eigVal <- eig$values - l2 / l
      eigVal[eigVal < 0] <- 0
      ny$L <- eig$vectors %*% diag(eigVal) %*% t(eig$vectors)
      ny$L <- 0.5 * (ny$L + t(ny$L))
      
      # Compute the gradient at ny
      ny.r <- nesterov.low.rank.plus.sparse.compute.gradient(Sigma, N, ny)
      if (is.null(ny.r)) {
        isposdef <- FALSE
      } else {
        isposdef <- is.positive.definite(nz$S - nz$L)
        issemdef1 <- is.positive.semi.definite(ny$L)
        issemdef2 <- is.positive.semi.definite(nz$L)
      }

      if (isposdef && issemdef1 && issemdef2) {
        diffS <- ny$S - y$S
        diffL <- ny$L - y$L
        scProd <- sum(diffS * y.r$dS) + sum(diffL * y.r$dL)
        RHS <- y.r$obj + scProd + 0.5 * l * (sum((diffS)**2) + sum((diffL)**2))
        if (ny.r$obj <= RHS + tol) {
          x <- ny
          z <- nz
          break()
        }  
      }
      l <- l * eta
    }
    
    r <- nesterov.low.rank.plus.sparse.compute.gradient(Sigma, N, x)
    pobj <- r$obj + l1 * sum(abs(x$S)) + l2 * sum(diag(x$L))
    tk <- (sqrt(tk**4 + 4 * tk ** 2) - tk**2) / 2
    objs <- c(objs, pobj)

 		diff1 <- frobenius.norm(oX$S - x$S) / frobenius.norm(oX$S)	
		if (frobenius.norm(oX$L) > 0) {
		  diff2 <- frobenius.norm(oX$L - x$L) / frobenius.norm(oX$L)	
		} else {
		  diff2 <- frobenius.norm(oX$L - x$L)
		}
		diff <- diff1 + diff2
		diffs <- c(diffs, diff)
		oX <- x
   
    if(length(objs) > 1) {
      if(i > minIter) {
      	if(all(diffs[(length(diffs) - 100):length(diffs)] < tol)) {
          	break()
      	}
			}
      if (pobj > objs[length(objs)-1]) {
        tk <- 0.5 # Kill the momentum (see : http://statweb.stanford.edu/~candes/papers/adap_restart_paper.pdf)
      }
    }
  }
  x$obj <- min(objs)
  x$objs <- objs
  x$diff <- diff
	x$diffs <- diffs
  x$iter <- i
  
  x
}
