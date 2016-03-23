source("./nesterov_sparse_conditional_gaussian_graphical_model.R")
source("generate_conditional_data.R")
set.seed(1)
N <- 10000
nZ <- 100 # Number of inputs
nX <- 40 # Number of outputs
nL <- 3 # Number of latent variables (confounders)
nP <- nZ # Number of latent variables mediating the effect of the inputs
sp_pph <- 0.05 # Sparsity of graph P -> Observed variables
sp_cph <- 0 # Sparsity of grpah Confounders -> Observed variables
sp_ph  <- 0.01 # Sparsity of DAG over X
sp_v <- 0.01  # Sparsity of graph : Z -> P
d <- generate_conditional_data(nZ, nX, nP, nL, sp_pph, sp_cph, sp_ph, sp_v, N)
full.dag <- d$dag
# Get the true matrix SZX
Zvars <- Xvars <- c()
for (n in d$names) {
  g1 <- grep("Ph", n)
  g2 <- grep("V", n)
  if (length(g1) > 0) {
    Xvars <- c(Xvars, n)
  }
  if (length(g2) > 0) {
    Zvars <- c(Zvars, n)
  }
}
true.SZX <- matrix(0,nrow=nZ, ncol=nX)
rownames(true.SZX) <- Zvars
colnames(true.SZX) <- Xvars
for(z in Zvars) {
  i <- which(Zvars == z)
  paths <- all_simple_paths(full.dag, from = z, to = Xvars)
  for (path in paths) {
    if(length(path) > 3) {
      next()
    }
    X_name <- path[3]$name
    j <- which(Xvars == X_name)
    true.SZX[i,j] <- 1
  }
}
# Get the true moralised graph SX
Zidx <- d$obs[1:nZ]
Xidx <- d$obs[(1 + nZ):(nZ + nX)]
true.SX <- get.adjacency(d$ug, sparse=FALSE)[Xidx, Xidx]

Z <- d$data[,Zidx]
X <- d$data[,Xidx]
# Standardise
for (i in 1:nZ) {
  v <- Z[,i]
  Z[,i] <- (v - mean(v)) / sd(v)
}
for (i in 1:nX) {
  v <- X[,i]
  X[,i] <- (v - mean(v)) / sd(v)
}

### Fit
init <- list()
init$SZX <- matrix(0, nrow=nZ, ncol=nX)
init$SX <- solve(cov(X))
lambda1 <- 0.02
lambda2 <- lambda1
fit <- nesterov.scggm.lasso(Z, X, lambda1, lambda2, init, max.iter = 10000)
plot(fit$objs)
SX <- fit$SX != 0
SX <- SX - diag(diag(SX))
SZX <- fit$SZX != 0
print(sum(SX * true.SX) / sum(true.SX)) # Recall for SX
print(sum(SX * true.SX) / sum(SX)) # Precision for SX
print(sum(SZX * true.SZX) / sum(true.SZX)) # Rec for SZX
print(sum(SZX * true.SZX) / sum(SZX)) # Prec for SZX
image(SZX)