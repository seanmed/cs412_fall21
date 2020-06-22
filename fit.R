library(GpGp)
library(fields)
library(Rcpp)
library(RcppArmadillo)


n1 <- 20
n2 <- 20
n <- n1*n2
locs <- as.matrix( expand.grid( (1:n1)/n1, (1:n2)/n2 ) )
X <- cbind(rep(1,n),locs[,2])
covparms <- c(2, 0.2, 0.75, .1)
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


fit_model2 <- function(y, locs, X = NULL, covfun_name = "matern_isotropic",
                      NNarray = NULL, start_parms = NULL, reorder = TRUE, group = TRUE,
                      m_seq = c(10,30), max_iter = 40, fixed_parms = NULL,
                      silent = FALSE, st_scale = NULL, convtol = 1e-4){
  
  n <- length(y)
  
  # check that length of observation vector same as
  # number of locations
  if( nrow(locs) != n ){
    stop("length of observation vector y not equal
              to the number of locations (rows in locs)")
  }
  
  # check if design matrix is specified
  if( is.null(X) ){
    if(!silent) cat("Design matrix not specified, using constant mean \n")
    X <- rep(1,n)
  }
  X <- as.matrix(X)
  
  
  # detect and remove missing values
  not_missing <- apply( cbind(y,locs,X), 1,
                        function(x){
                          if( sum(is.na(x) | is.infinite(x)) > 0 ){
                            return(FALSE)
                          } else { return(TRUE) }
                        }
  )
  if( sum(not_missing) < n ){
    y <- y[not_missing]
    locs <- locs[not_missing,,drop=FALSE]
    X <- X[not_missing,,drop=FALSE]
    cat(paste0( n - sum(not_missing),
                " observations removed due to missingness or Inf\n"))
  }
  
  # redefine n
  n <- length(y)
  
  # check that start_parms is specified when fixed_parms is
  if( is.null(fixed_parms) ){
    if( is.null(start_parms) ){
      start <- get_start_parms(y,X,locs,covfun_name)
      start_parms <- start$start_parms
    } else {
      # check if start_parms has the right length
      start <- get_start_parms(y,X,locs,covfun_name)
      if(length(start_parms) != length(start$start_parms) ){
        stop(paste0("start_parms not correct length for ",covfun_name))
      }
    }
    # define the parameters we are not fixing
    active <- rep(TRUE, length(start_parms) )
  } else {
    if( is.null(start_parms) ){
      stop("start_parms must be specified whenever fixed_parms is")
    }
    # check if start_parms has the right length
    start <- get_start_parms(y,X,locs,covfun_name)
    if(length(start_parms) != length(start$start_parms) ){
      stop(paste0("start parms not correct length for ",covfun_name))
    }
    # check whether fixed_parms has appropriate values
    if( max( fixed_parms - floor(fixed_parms) ) > 0 ){
      stop("fixed_parms must contain indices of parms you want to fix")
    }
    if( min( fixed_parms < 1 ) || max(fixed_parms) > length(start_parms) ){
      stop("fixed_parms must be between 1 and number of parameters")
    }
    # define the parameters we are not fixing
    active <- rep(TRUE, length(start_parms) )
    active[fixed_parms] <- FALSE
  } 
  
  # get link functions
  linkfuns <- get_linkfun(covfun_name)
  link <- linkfuns$link
  dlink <- linkfuns$dlink
  invlink <- linkfuns$invlink
  invlink_startparms <- invlink(start_parms)
  lonlat <- linkfuns$lonlat
  if(lonlat){
    cat("Assuming columns 1 and 2 of locs are (longitude,latidue) in degrees\n")
  }
  space_time <- linkfuns$space_time
  
  penalty <- get_penalty(y,X,locs,covfun_name) 
  pen <- penalty$pen
  dpen <- penalty$dpen
  ddpen <- penalty$ddpen
  
  # get an ordering and reorder everything
  if(reorder){
    if(!silent) cat("Reordering...")
    if( n < 1e5 ){  # maxmin ordering if n < 100000
      ord <- order_maxmin(locs, lonlat = lonlat, space_time = space_time)
    } else {        # otherwise random order
      ord <- sample(n)
    }
    if(!silent) cat("Done \n")
  } else {
    ord <- 1:n
  }
  yord <- y[ord]
  locsord <- locs[ord,,drop=FALSE]
  Xord <- as.matrix( X[ord,,drop=FALSE] )
  
  # get neighbor array if not provided
  if( is.null(NNarray) ){
    if(!silent) cat("Finding nearest neighbors...")
    NNarray <- find_ordered_nn(locsord, m=max(m_seq), lonlat = lonlat,
                               st_scale = st_scale)
    if(!silent) cat("Done \n")
  }
  
  # refine the estimates using m in m_seq
  for(i in 1:length(m_seq)){
    m <- m_seq[i]
    if(group){
      
      NNlist <- group_obs(NNarray[,1:(m+1)])
      likfun <- function(logparms){
        
        lp <- rep(NA,length(start_parms))
        lp[active] <- logparms
        lp[!active] <- invlink_startparms[!active]
        
        likobj <- vecchia_grouped_profbeta_loglik_grad_info(
          link(lp),covfun_name,yord,Xord,locsord,NNlist)
        likobj$loglik <- -likobj$loglik - pen(link(lp))
        likobj$grad <- -c(likobj$grad)*dlink(lp) -
          dpen(link(lp))*dlink(lp)
        likobj$info <- likobj$info*outer(dlink(lp),dlink(lp)) -
          ddpen(link(lp))*outer(dlink(lp),dlink(lp))
        likobj$grad <- likobj$grad[active]
        likobj$info <- likobj$info[active,active]
        return(likobj)
        
      }
      
    } else {
      
      likfun <- function(logparms){
        
        lp <- rep(NA,length(start_parms))
        lp[active] <- logparms
        lp[!active] <- invlink_startparms[!active]
        
        likobj <- vecchia_reml_loglik_grad_info(
          link(lp),covfun_name,yord,Xord,locsord,NNarray[,1:(m+1)])
        likobj$loglik <- -likobj$loglik - pen(link(lp))
        likobj$grad <- -c(likobj$grad)*dlink(lp) -
          dpen(link(lp))*dlink(lp)
        likobj$info <- likobj$info*outer(dlink(lp),dlink(lp)) -
          ddpen(link(lp))*outer(dlink(lp),dlink(lp))
        likobj$grad <- likobj$grad[active]
        likobj$info <- likobj$info[active,active]
        return(likobj)
      }
    }
    fit <- fisher_scoring( likfun,invlink(start_parms)[active],
                           link,silent=silent, convtol = convtol, max_iter = max_iter )
    invlink_startparms[active] <- fit$logparms
    #start_parms[active] <- fit$covparms
    start_parms <- link(invlink_startparms)
    fit$loglik <- -fit$loglik - pen(start_parms)
    invlink_startparms <- invlink(start_parms)
  }
  
  # return fit and information used for predictions
  fit$covfun_name <- covfun_name
  #fit$covparms <- start_parms
  lp <- rep(NA,length(start_parms))
  lp[active] <- fit$logparms
  lp[!active] <- invlink_startparms[!active]
  fit$covparms <- link(lp)
  fit$y <- y
  fit$locs <- locs
  fit$X <- X
  class(fit) <- "GpGp_fit"
  return(fit)
}


