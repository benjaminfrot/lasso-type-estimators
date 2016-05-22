source("./generate_data.R")
source("./NesterovMethods/nesterov_low_rank_plus_sparse.R")
require(ggplot2)

set.seed(3)
N <- 10000
nO <- 40 # 40 Observed Variables
nL <- 3 # 3 Latent variables
sp_O_L <- 0.7 # Sparsity of graph O -> Latent
sp_O <- 0.01 # Sparsitiy of observed DAG
d <- generate_data(nO, nL, sp_O_L, sp_O, N)
X <- d$data[,d$obs]
for (i in 1:nO) {
  v <- X[,]
  X[,i] <- (v - mean(v)) / sd(v)
}
init <- list()
init$S <- solve(cov(X))
init$L <- matrix(0, ncol=nO, nrow=nO)
# True moralised matrix
moramat = get.adjacency(d$ug, sparse=FALSE)[d$obs, d$obs]

# Look at the precision recall curves for a range of values of gamma / lambda
gammas <- c(0.2, 0.5, 0.7)
lambdas <- 1.1**(c(-50:10))
results <- matrix(NA, ncol=6)
for (gamma in gammas) {
  init <- list()
  init$S <- solve(cov(X))
  init$L <- matrix(0, ncol=nO, nrow=nO)
  print(paste("Now computing path for Gamma = ", gamma))
  for (lambda in lambdas) {
    l1 <- lambda * gamma
    l2 <- lambda * (1 - gamma)
    fit <- nesterov.low.rank.plus.sparse(X, l1, l2, init, max.iter = 10000)
    # Value of the objective function as a function of the iterations
    #plot(fit$objs)
    init <- list()
    init$S <- fit$S
    init$L <- fit$L
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
# Plot the Precision / Recall curve
ggplot(df.results, aes(x=Recall, y=Precision)) + geom_point(aes(colour=as.factor(Gamma))) +
  facet_wrap(~Gamma)
# Plot the rank as function of Gamma
ggplot(df.results, aes(x=Lambda, y=Rank)) + geom_point(aes(colour=as.factor(Gamma))) +
  facet_wrap(~Gamma)
ggplot(df.results, aes(x=Lambda, y=log(Sparsity))) + geom_point(aes(colour=as.factor(Gamma))) +
  facet_wrap(~Gamma)

# Both the graph structure and rank of the latent piece are correctly estimated for
# gamma = 0.5