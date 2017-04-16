rm(list=ls())
source("./split_bregman_lrps.R")
library(ggplot2)
library(R.matlab)

# Load a toy dataset.
# p = 32, n = 3000. 
toy.data <- readMat("./test_data.mat")
moramat <- toy.data$SX # True precision matrix
moramat <- 1 * (moramat !=0) # Here we look only the sparsity pattern.
moramat <- moramat - diag(diag(moramat)) # Remove the diagonal entries
Lmle <- toy.data$LXmle # This is an estimate of the latent piece obtained by computing the MLE with 5.10^5 samples.
# The true rank of L is 4.

# Compute the covariance matrix
Sigma <- cor(toy.data$X)

# Look at the precision recall curves for a range of values of gamma / lambda
gammas <- c(0.2, 0.5, 0.7)
# We compute the full path for multiple values of Gamma. A warm stat is used.
lambdas <- 1.1**(c(-50:10))
results <- matrix(NA, ncol=6)

fit = NULL

for (gamma in gammas) {
  print(paste("Now computing path for Gamma = ", gamma))
  for (lambda in lambdas) {
    l1 <- lambda * gamma
    l2 <- lambda * (1 - gamma)
    # Reuse the previous fit as initial value.
    fit <- fit.low.rank.plus.sparse(Sigma, l1, l2, init = fit, maxiter = 10000, mu=0.001, tol = 1e-7)
    if(fit$termcode == -2) {
      break()
    }
    rank <- sum(eigen(fit$L)$val > 1e-6)
    S <- fit$S - diag(diag(fit$S))
    S <- S != 0
    sparsity <- sum(S != 0) / (dim(S)[1]**2 - dim(S)[1])
    prec <- sum(S * moramat) / sum(S)
    rec <- sum(S * moramat) / sum(moramat)
    results <- rbind(results, c(lambda, gamma, prec, rec, rank, sparsity))
    print(paste("Lambda:", lambda, ", #Iterations:", fit$iter, ", Sparsity:", sparsity, ", Rank:", rank, "Precision", prec, "Recall", rec))
    if(sparsity == 0) {
      break()
    }
  }
}

df.results <- as.data.frame(results[2:nrow(results),])
names(df.results) <- c("Lambda", "Gamma", "Precision", "Recall", "Rank", "Sparsity")
idx <- is.nan(df.results$Precision)
df.results <- df.results[!idx,]
# Plot the Precision / Recall curves. One curve per value of Gamma
ggplot(df.results, aes(x=Recall, y=Precision)) + geom_point(aes(colour=as.factor(Gamma))) +
  facet_wrap(~Gamma)
# Plot the rank as function of Lambda, Gamma.
ggplot(df.results, aes(x=Lambda, y=Rank)) + geom_point(aes(colour=as.factor(Gamma))) +
  facet_wrap(~Gamma)
# Plot the sparsity as a function of Lambda, Gamma
ggplot(df.results, aes(x=Lambda, y=log(Sparsity))) + geom_point(aes(colour=as.factor(Gamma))) +
  facet_wrap(~Gamma)