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
covparms <- c(2, 0.2, 0.75, 0)
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


sourceCpp("reml/reml/main.cpp")

#rl returns loglik, grad, and info as a function of covparms. rlCompute computes these
#in one pass through the data

rl = function(covps){
  rlCompute(covparms = covps, locs,NNarray,y,X,r)
}


rl(covparms)

#something is wrong with fisher information 



# rlg returns just the loglik
rlg = function(covps){
  rlCompute(covparms = covps,locs = locs,NNarray = NNarray,y=y,X=X,r=r)$loglik
}

# use finte difference derivative (very slow) to check if grad from rl is correct (requires library(numDeriv))
g = function(covps){
  return(grad(rlg,covps))
}

g(covparms)

# grad is correct, not sure how to check fisher information