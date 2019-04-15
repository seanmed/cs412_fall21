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
k = 10
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
covparms <- c(variance = 4, range = 0.1, smoothness = 0.5, nugget = 0)
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
m <- 30
NNarray = find_ordered_nn(locsord,m)
NNarray1 = find_ordered_nn(locsord1,m)


# Comparisons: 
# Ungrouped
print( system.time( ll <- vecchia_loglik(covparms, "matern_isotropic", rep(0,n), locsord, NNarray)))
print( system.time( ll1 <- vecchia_loglik(covparms, "matern_isotropic", rep(0,n), locsord1, NNarray1)))
print( ll - ll1 )


NNlist <- group_obs(NNarray, 2)
NNlist1 <- group_obs(NNarray1, 2)

# Grouped
print( system.time(llg <- vecchia_loglik_grouped(covparms, "matern_isotropic", rep(0,n), locsord, NNlist)))
print( system.time(llg1 <- vecchia_loglik_grouped(covparms, "matern_isotropic", rep(0,n), locsord1, NNlist1)))


names(NNlist)

# list of all indices in U_1, U_2, ...
# listed out in order, one after the other
NNlist$all_inds

# information about the end of each U_i codified in last_ind_of_block
NNlist$last_ind_of_block[1]
# better name:
# NNlist$last_ind_of_U[12]


U12 <- NNlist$all_inds[ (NNlist$last_ind_of_block[11]+1):NNlist$last_ind_of_block[12] ]

# which of these are from B12

# all indicies in B1, B2, ...
# listed out in order, one after the other
NNlist$global_resp_inds

# this tells you where they end
NNlist$last_resp_of_block[12]
# better name
# last_ind_of_B[12]

B12 <- NNlist$global_resp_inds[ (NNlist$last_resp_of_block[11]+1):NNlist$last_resp_of_block[12] ]


B12 <- NNlist$global_resp_inds[ (NNlist$last_resp_of_block[11]+1):NNlist$last_resp_of_block[12] ]

sv <- NNlist$local_resp_inds[ (NNlist$last_resp_of_block[11]+1):NNlist$last_resp_of_block[12] ]
U12[sv]
B12


# construct global resp inds
# look in NNarray and grab union of all neighbors for each B
# assign to all_inds
# ...

