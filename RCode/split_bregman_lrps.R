.split.bregman.low.rank.plus.sparse.update.parameters <- function(ps, opts) {
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

#' Inverse covariance estimation via a low-rank plus sparse decomposition
#'
#' @export
#' @details Minimise the function \deqn{-logdet S + Tr(\Sigma S) + \lambda_1 ||S||_1 + \lambda_2 ||L||_\ast}
#' as per the estimator suggested in Chandrasekaran et. al (doi:10.1214/11-AOS949). We use Split-Bregman iterations
#' as described by Ye et al in https://arxiv.org/abs/1110.3076.
#' @param Sigma The sample covariance matrix
#' @param Lambda1 Penalty on the l1-norm of the sparse component S.
#' @param Lambda2 Penalty on the nuclear norm of the low-rank component L.
#' @param maxiter Maximum number of iterations.
#' @param mu Step size
#' @param tol Tolerance. The algorithm stops when ||S_k - S_k-1||_F / ||S_k-1||_F + ||L_k - L_k-1||_F / ||L_k-1||_F < tol.
#' @param init Initial parameters. Can be used for warm starts. See examples below.
#' @return A list L where L$S (resp. L$L) is the estimate of the sparse matrix (resp. low-rank matrix).
#' L$termcode is the termination code (0: good termination, -1: did not converge, -2: Lambda1 was so large that the estimated S is empty.)
#' L$iter is the number of iterations. L$diffs is the vector containing the sequence ||S_k - S_k-1||_F / ||S_k-1||_F + ||L_k - L_k-1||_F / ||L_k-1||_F
#' for k=1:L$iter, it shows the convergence speed of the algorithm.
#' @examples
#' # Generate some low-rank plus sparse data with 4 confounders
#' set.seed(0)
#' toy.data <- generate.conditional.data(n = 1000, 2, 5, 32)
#' Y <- cbind(toy.data$Z, toy.data$X)
#' Sigma = cov2cor(cov(Y))
#' # Fit the model with gamma = 0.2
#' gamma <- 0.2
#' lambda <- 0.1
#' fit <- fit.low.rank.plus.sparse(Sigma, gamma * lambda, (1-gamma)*lambda)
#' # Look at the convergence
#' plot(fit$diffs)
#' # Look at the sparse matrix
#' image(fit$S!=0)
#' # Now, fit for another value of lambda, but reuse the previous fit for a warm start
#' lambda <- 0.2
#' fit <- fit.low.rank.plus.sparse(Sigma, gamma * lambda, (1-gamma)*lambda, init = fit)
fit.low.rank.plus.sparse <- function(Sigma, Lambda1, Lambda2, init=NULL, maxiter=10000, mu=0.01, tol=1e-07) {
  require(matrixcalc)
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
  parameters$termcode <- -1
  for (i in 1:maxiter) {
    new_parameters <- .split.bregman.low.rank.plus.sparse.update.parameters(parameters, options)
    
    # Check convergence of S, L
    if(sum(new_parameters$S) == 0) {
      parameters$S <- diag(p)
      parameters$termcode <- -2
      break()
    }
    if ((matrixcalc::frobenius.norm(parameters$L) > 0)) {
      diff <- matrixcalc::frobenius.norm(new_parameters$S - parameters$S) / matrixcalc::frobenius.norm(parameters$S) +
        matrixcalc::frobenius.norm(new_parameters$L - parameters$L) / matrixcalc::frobenius.norm(parameters$L)
    } else {
      diff <- matrixcalc::frobenius.norm(new_parameters$S - parameters$S) / matrixcalc::frobenius.norm(parameters$S) +
        matrixcalc::frobenius.norm(new_parameters$L - parameters$L)
    }
    diffs <- c(diffs, diff)
    parameters <- new_parameters
    if (diff < tol) {
      parameters$termcode <- 0
      break()
    }
  }
  
  parameters$iter <-i
  parameters$diffs <- diffs
  
  parameters
}