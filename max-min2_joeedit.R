library("GpGp")
library('RANN')

covparms <- c(variance = 4, range = 0.1, smoothness = 1.0, nugget = 0)

# grid size for data locations
gsize <- 100
nvec <- c(gsize,gsize)
n <- prod(nvec)

x1 <- (1:nvec[1])/nvec[1]
x2 <- (1:nvec[2])/nvec[2]
locs <- as.matrix(expand.grid(x1,x2))

# k is the number you want to group by,
k = 10
l= nn2(locs,locs,k)[[1]]
ind_group= matrix(unique(as.vector(t(l))),ncol = k, byrow = T)
ind_group_list = list()
for(i in 1:nrow(ind_group)){
  ind_group_list[[i]] = as.list(ind_group[i,])}

ngroups <- 1000
maxit <- 100
km <- kmeans(locs, ngroups, iter.max = maxit, algorithm = "Lloyd" )
if( km$iter == maxit + 1 ){ 
    cat("did not converge\n") 
}

plot( locs[,1], locs[,2], col = colvec[km$cluster] )
# loop over the 1000 (or however many) groups
# each time assign indices to a list element
ind_group_list = list()
for(j in 1:ngroups){
    ind_group_list[[j]] <- which( km$cluster == j )
}



# a plot of whats going on. red is first group, blue second, green third...
colorvec <- c("red","black","blue","green","magenta","orange")
colvec <- rep(colorvec,1000)
d = .1
a_1 = locs[1,]
plot(a_1[1], a_1[2], xlim=c(0, 1), ylim=c(0, 1), col='black', pch=20, xlab='x', ylab='y') 
#points(locs[,1], locs[,2], col='grey')  
#locs[c(ind_group[1,]),]
for(j in 1:length(ind_group_list)){
    points(locs[c(ind_group[j,]),][,1],locs[c(ind_group[j,]),][,2], col=colvec[j], pch=16)
}

points(locs[c(ind_group[1,]),][,1],locs[c(ind_group[1,]),][,2], col='red', pch=20)
points(locs[c(ind_group[2,]),][,1],locs[c(ind_group[2,]),][,2], col='blue', pch=20)
points(locs[c(ind_group[3,]),][,1],locs[c(ind_group[3,]),][,2], col='green', pch=20)

loc_group_mat = matrix(rep(NA,2*n),ncol = 2*k)
for(i in 1:(n/k)){
  for(j in seq(1,2*k,2)){
  loc_group_mat[i,seq(1,2*k,2)] = locs[ind_group[i,1:k],][,1]
  loc_group_mat[i,seq(2,2*k,2)] = locs[ind_group[i,1:k],][,2]
  }
}


mid = matrix(NA,nrow = nrow(ind_group),ncol = 2)
for(i in 1:nrow(loc_group_mat)){
  mid[i,1] = sum(loc_group_mat[i,seq(1,2*k,2)])/k
  mid[i,2] = sum(loc_group_mat[i,seq(2,2*k,2)])/k
}

# indices of max-min ordered points
ordm = order_maxmin(mid)
ordm = order_maxmin(km$centers)

#indices of all  points after max-min ordering
ord_full = rep(NA,n)
for(i in 1:k){
  ord_full[seq(i,n,k)] = ind_group[ordm,i]
}

