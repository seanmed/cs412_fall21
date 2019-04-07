library('GpGp')
library('RColorBrewer')
library('fields')

# grid size for data locations
gsize <- 100
nvec <- c(gsize,gsize)
n <- prod(nvec)

x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))

# k is the number you want to group by for max-man-2. Should divide n.
k = 100
ngroups <- n/k
maxit <- 100
km <- kmeans(locs, ngroups, iter.max = maxit, algorithm = "Lloyd" )
if( km$iter == maxit + 1 ){
  cat("did not converge")
}


# loop over the groups
# each time assign indices to a list element
ind_group_list <- list()
for(j in 1:ngroups){
  ind_group_list[[j]] <- which(km$cluster == j )
}

# run max-min on the midpoints
ordm <- order_maxmin(km$centers)

# reorder ind_group_list according to ordm.
ord_ind_group_list = list()
ord_ind_group_list<-ind_group_list[ordm]

# Get the final ordering.
ord_full <- unlist(ord_ind_group_list)

# A plot of the grouping for max-min-2.
c <- n/k
color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
colvec <- sample(color,c,replace = T)
plot(locs[,1], locs[,2], col = colvec[km$cluster],pch=16,cex=1)

# Set covariance parameters and generate data.
covparms <- c(variance = 4, range = 0.1, smoothness = 1.0, nugget = 0)
y=fast_Gp_sim(covparms, "matern_isotropic",locs,40)

# An ordering using max-min.
ord <- order_maxmin(locs)
locsord =locs[ord,]
yord = y[ord]

# An ordering using max-min-2.
ord1 <- ord_full
locsord1 <- locs[ord1,]
yord1 <- y[ord1]

# find the ordered m nearest neighbors.
m <- 10
NNarray = find_ordered_nn(locsord,m)
NNarray1 = find_ordered_nn(locsord1,m)


# Comparisons: 
# Ungrouped
print( system.time( ll <- vecchia_loglik(covparms, "matern_isotropic", rep(0,n), locsord, NNarray)))
print( system.time( ll1 <- vecchia_loglik(covparms, "matern_isotropic", rep(0,n), locsord1, NNarray1)))

NNlist <- group_obs(NNarray, 2)
NNlist1 <- group_obs(NNarray1, 2)

# Grouped
print( system.time(llg <- vecchia_loglik_grouped(covparms, "matern_isotropic", rep(0,n), locsord, NNlist)))
print( system.time(llg1 <- vecchia_loglik_grouped(covparms, "matern_isotropic", rep(0,n), locsord1, NNlist1)))
