library(GpGp)
library(fields)
library(Rcpp)
library(RcppArmadillo)
library(numDeriv)

n1 <- 20
n2 <- 20
n <- n1*n2
locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
X <- cbind(rep(1,n),locs[,2])
covparms <- c(2, 0.3, 0.7, .1)
y <- fast_Gp_sim(covparms, "matern_isotropic", locs, 50 )
NNarray <- find_ordered_nn(locs,30)


#reorder so that first rank(X) rows are LI
r <- qr(X)$rank
i <- 2
j=nrow(X)
while(qr(X[1:r,])$rank != r){
  while(qr(X[1:i,])$rank != i){
    X <- X[c(1:(i-1),j, (i+1):(j-1), i) ,]
    locs <-locs[c(1:(i-1),j, (i+1):(j-1), i) ,]
    y <- y[c(1:(i-1),j, (i+1):(j-1), i)]
    j <- j-1
  }
  i <- i+1
}


sourceCpp("main.cpp")
vecchia_reml_loglik_grad_info(covparms,"matern_isotropic",y,X,locs,NNarray)
vecchia_profbeta_loglik_grad_info(covparms,"matern_isotropic",y,X,locs,NNarray)
