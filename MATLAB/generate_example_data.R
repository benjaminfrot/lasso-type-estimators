library(MASS)
library(R.matlab)

args = commandArgs(trailingOnly=TRUE)
n <- as.numeric(args[[1]]) # Sample size
d_Z <- as.numeric(args[[2]]) # 2**d_Z = number of ``latent pathways'' = rank of KZX
d_H <- as.numeric(args[[3]]) # 2**d_H = number of confounders = rank of LX
fn <- (args[[4]]) # Filename: Where to save the data, e.g. toydata.mat

p <- 32 # Dimension of the problem

nn <- 100000
#### Common to all simulations
K_X <- matrix(0, ncol=p, nrow=p)
for(i in 1:p) {
  for(j in 1:p) {
    if(i == (j-1)) {
      if(i %% 5 > 0) {
        K_X[i,j] <- runif(n=1, min=0.3, max=1) * sign(rnorm(n=1))
      }
    }
  }
}
K_X <- K_X + t(K_X) + diag(runif(n=p, min=0.5, max=1.5))
Z <- matrix(rt(n=nn*p, df=4), ncol=p)
####

#d_Z <- 5 # The rank of the matrix K_{ZX} is 2^(d_Z)
#d_H <- 4 # The rank of the matrix K_{HX} is 2^(d_H)
rk <- 2**d_Z
nc <- p / rk
M_ZF <- matrix(0, ncol=rk, nrow = p)
for(j in 1:rk) {
  start <- (j - 1) * nc + 1
  stop <- j * nc
  M_ZF[start:stop,j] <- 1
}
M_ZF <- M_ZF[sample(p, size=p, replace=F),]
# Make the matrix of Fs as a function of the Zs
Fs <- matrix(NA, nrow=nn, ncol=rk)
for(j in 1:rk) {
  idx <- which(M_ZF[,j]!=0)
  A <- as.matrix(Z[,idx])
  Fs[,j] <- rowSums(A)
}

K_FX <- M_ZF
K_FX <- (K_FX[sample(p, size=p, replace=F),])
for(idx in which(K_FX!=0)) {
  K_FX[idx] <- sign(rnorm(n=1)) * runif(n=1, min=0.5, max=1.5)
}
K_FX <- t(K_FX)

rk <- 2**d_H
K_U <- diag(runif(n=rk, min=0.3, max=1))
nc <- p / rk
M_UX <- matrix(0, ncol=rk, nrow = p)
for(j in 1:rk) {
  start <- (j - 1) * nc + 1
  stop <- j * nc
  M_UX[start:stop,j] <- 1
}
M_UX <- M_UX[sample(p, size=p, replace=F),]
K_UX <- t(M_UX)
for(idx in which(K_UX != 0)) {
  K_UX[idx] <- sign(rnorm(n=1)) * runif(n=1, min=0.5, max=1.5)
}

####################
# Prevent it from being singular
K_XU <- t(K_UX)
K_XH <- cbind(rbind(K_X, K_UX), rbind(K_XU, K_U))
minEig <- min(eigen(K_XH)$values)
diag(K_XH) <- diag(K_XH) + 0.05 - minEig

# Generate the data
d_XH <- matrix(0, ncol=p + rk, nrow = nn)
for(i in 1:nn) {
  mu <- - solve(K_XH) %*% t(cbind(K_FX, K_FX[,1:rk]*0)) %*% matrix(Fs[i,], ncol=1)
  d_XH[i,] <- mvrnorm(n=1, mu=mu, Sigma=solve(K_XH))
}

X <- d_XH[,1:p]
H <- d_XH[,(p+1):(rk+p)]

for(i in 1:p) {
  v <- X[, i]
  X[, i] <- (v - mean(v)) / sd(v)
}

for(i in 1:ncol(Z)) {
  v <- Z[,i]
  Z[,i] <- (v - mean(v)) /sd(v)
}

for(i in 1:ncol(H)) {
  v <- H[,i]
  H[,i] <- (v - mean(v)) /sd(v)
}
nH <- ncol(H)
# Compute the MLEs
Sig_XH <- cor(cbind(X, H))
Sig_Z <- cor(Z)
Sig_Z_XH <- cor(Z,cbind(X,H))
KXHmle <- solve(Sig_XH - t(Sig_Z_XH) %*% solve(Sig_Z) %*% Sig_Z_XH)
KZXmle <- (-solve(Sig_Z) %*% Sig_Z_XH %*% KXHmle)[,1:p]
SXmle <- KXHmle[1:p, 1:p]
LXmle <- KXHmle[1:p,(p+1):(p+nH)] %*% solve(KXHmle[(p+1):(p+nH),(p+1):(p+nH)]) %*% t(KXHmle[1:p,(p+1):(p+nH)])

# Keep only the first n samples
X <- X[1:n,]
Z <- Z[1:n,]
for(i in 1:p) {
  v <- X[, i]
  X[, i] <- (v - mean(v)) / sd(v)
}

for(i in 1:ncol(Z)) {
  v <- Z[,i]
  Z[,i] <- (v - mean(v)) /sd(v)
}

ObsData <- cbind(Z,X)

writeMat(con=fn, SX=K_X, Z=Z, X=X, z=d_Z, u=d_H, ObsData=ObsData, Sig_XH=Sig_XH, 
         Sig_Z=Sig_Z, Sig_Z_XH=Sig_Z_XH, KXHmle=KXHmle, KZXmle=KZXmle, SXmle=SXmle, LXmle=LXmle)
