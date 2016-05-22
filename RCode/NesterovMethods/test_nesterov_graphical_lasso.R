source("generate_data.R")
source("./nesterov_graphical_lasso.R")

set.seed(3)
N <- 1000
nO <- 40 # 40 Observed Variables
nL <- 3 # 3 Latent variables
sp_O_L <- 0 # NO EFECT OF THE LATENT VARIABLES
sp_O <- 0.01 # Sparsitiy of observed DAG
d <- generate_data(nO, nL, sp_O_L, sp_O, N)
X <- d$data[,d$obs]
# Standardise
for (i in 1:nO) {
  v <- X[,i]
  X[,i] <- (v - mean(v)) / sd(v)
}
init <- list()
init$S <- solve(cov(X)) # Initialise at the MLE since we have enough samples
# True moralised matrix
moramat = get.adjacency(d$ug, sparse=FALSE)[d$obs, d$obs]

# Look at the precision recall curves for a range of values of lambda
lambdas <- 1.1**(c(-100:10))
results <- matrix(NA, ncol=4)
for (lambda in lambdas) {
  l1 <- lambda
  fit <- nesterov.graphical.lasso(X, l1, init, max.iter = 10000)
  # Plot the value of the objective function as a function of the #iterations
  # plot(fit$objs)
  init <- list()
  init$S <- fit$S
  S <- fit$S - diag(diag(fit$S))
  S <- S != 0
  sparsity <- sum(S != 0) / (dim(S)[1]**2 - dim(S)[1])
  prec <- sum(S * moramat) / sum(S)
  rec <- sum(S * moramat) / sum(moramat)
  results <- rbind(results, c(1, lambda, prec, rec))
  print(paste("Lambda:", lambda, ", #Iterations:", fit$iter, ", Sparsity:", sparsity))
  if(sparsity == 0) {
    break()
  }
}

df.results <- as.data.frame(results[2:nrow(results),])
names(df.results) <- c("Method", "Lambda", "Precision", "Recall")
idx <- is.nan(df.results$Precision)
df.results <- df.results[!idx,]
plot(df.results$Recall, df.results$Precision)
title("Precision = f(Recall)")