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


# a plot of whats going on. red is first group, blue second, green third...
d = .1
a_1 = locs[1,]
plot(a_1[1], a_1[2], xlim=c(0, a_1[1]+d), ylim=c(0, a_1[2]+d), col='black', pch=20, xlab='x', ylab='y') 
points(locs[,1], locs[,2], col='grey')  
locs[c(ind_group[1,]),]
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

#indices of all  points after max-min ordering
ord_full = rep(NA,n)
for(i in 1:k){
  ord_full[seq(i,n,k)] = ind_group[ordm,i]
}

