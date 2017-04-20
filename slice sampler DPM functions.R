#functions to implement dpm using slice sampler
require(mcclust)
require(MCMCpack)
require(mvtnorm)

determineKP <- function(min.u,
                        weights,
                        Kstar,
                        v,
                        alpha){
  counter <- Kstar
  if (sum(weights)<(1-min.u)){
    while(sum(weights)<(1-min.u)){
      v <- c(v,0)
      weights <- c(weights,0)
      v[counter+1] <- rbeta(1,1,alpha)
      weights[counter+1] <- v[counter+1]*prod(1-v[1:counter])
      counter <- counter+1
    }
  }
  #number of potential components
  KP <- counter-Kstar
  output <- list(v=v,weights=weights,KP=KP)
  return(output)
}

zupdate <- function(Kstar,
                    n,
                    y,
                    u,
                    weights,
                    mu,
                    Sigma){
  
  nk <- rep(0,Kstar)
  
  for (i in 1:n){
    zprob <- rep(0,Kstar)
    #zprob[weights>u[i]] <- sapply(which(weights>u[i]),function(k) -0.5*sum(log(eigen(Sigma[,,k])$values))-0.5*(y[,i]-mu[,k])%*%(solve(Sigma[,,k],y[,i]-mu[,k])))
    zprob[weights>u[i]] <- 
      sapply(which(weights>u[i]),
             function(k) dmvnorm(y[,i], mu[,k], Sigma[,,k], log = TRUE))
    
    zprob[weights>u[i]] <- exp(zprob[weights>u[i]]-max(zprob[weights>u[i]]))
    z[i] <- 1+sum(runif(1)>cumsum(zprob/sum(zprob)))
    nk[z[i]] <- nk[z[i]]+1
  }
  output <- list(nk=nk,z=z)
  return(output)
}

vupdate <- function(Kstar,nk,alpha,n){
  v <- rep(0,Kstar);
  weights <- rep(0,Kstar)
  v[1] <- rbeta(1,1+sum(z==1),alpha+sum(z>1))
  
  if (Kstar > 1){
    for (k in 2:Kstar){
      v[k] <- rbeta(1,1+nk[k],alpha+n-sum(nk[1:k]))
    }
  }
  
  return(v)
}

weightsupdate <- function(Kstar,v){
  weights <- rep(0,Kstar)
  weights[1] <- v[1]
  for (k in 2:Kstar){
    weights[k] <- v[k]*prod(1-v[1:k-1])
  }
  return(weights)
}

muSigmaupdate <- function(y,
                          z,
                          nk,
                          hypers,
                          P,
                          Kstar){
  
  Sigma <- array(0, c(P,P,Kstar))
  
  mu <- array(0, c(P,Kstar))
  
  for (k in 1:Kstar){ # can we vectorise?
    if(nk[k]>1){
      ybar <- apply(y[ , z==k],
                    1,
                    mean) # get the rowmeans (so what's the mean PC loading across all curves in this group?)
    } else {
      ybar <- y[ , z==k]
    }
    
    D <- as.matrix(y[ , z==k]-ybar) # centering
    
    Sigma[ , , k] <- 
      riwish(hypers$c0+nk[k],
             hypers$C0 + 
               with(hypers,
                    (nk[k]*N0/(nk[k]+N0))*(ybar-b0) %*% t(ybar-b0) + D%*%t(D) ))
    
    mu[ , k] <-
      rmvnorm(n=1,
              mean = (nk[k]*ybar + with(hypers,N0*b0))/(nk[k]+hypers$N0),
              sigma = Sigma[ , , k]/(hypers$N0+nk[k]))
  }
  
  output <- list(mu=mu,Sigma=Sigma)
  return(output)
}


#label switching
labelswitch <- function(Kstar,
                        z,
                        nk,
                        weights,
                        v,
                        alpha,
                        mu,
                        Sigma){
  move <- runif(1)
  
  if (move<=(1/2)){
    #choose two 'alive' components
    kchoose <- sample(1:Kstar, 2)
    logr <- (nk[kchoose[1]]-nk[kchoose[2]])*(log(weights[kchoose[2]])-log(weights[kchoose[1]]))
  } else {
    kchoose <- rep(0,2)
    #choose any two components from 1:Kstar-1 (ie. empties ok for this move)
    kchoose[1] <- sample(1:(Kstar-1),1)
    kchoose[2] <- kchoose[1]+1
    logr <- nk[kchoose[1]]*log(1-v[kchoose[2]])-nk[kchoose[2]]*log(1-v[kchoose[1]])
  }
  #if(move>(2/3)){
  #tmp <- sample(2:sum(nk>0), 1)
  #kchoose <- rep(0,2)
  #kchoose[1] <- which(nk>0)[tmp]
  #kchoose[2] <- which(nk>0)[tmp-1]
  #weightscomb <- weights[kchoose[1]]+weights[kchoose[2]]
  #R1 <- (1+alpha+nk[kchoose[2]]+sum(nk[(1:Kstar)>kchoose[2]]))/(alpha+nk[kchoose[2]]+sum(nk[(1:Kstar)>kchoose[2]]))
  #R2 <- (alpha+nk[kchoose[1]]+sum(nk[(1:Kstar)>kchoose[2]]))/(1+alpha+nk[kchoose[1]]+sum(nk[(1:Kstar)>kchoose[2]]))
  #weightsdash <- weights[kchoose[2]]*R1+weights[kchoose[1]]*R2
  #logr <- (nk[kchoose[1]]+nk[kchoose[2]])*(log(weightscomb/weightsdash))+nk[kchoose[2]]*log(R1)+nk[kchoose[1]]*log(R2)
  #}
  
  if (log(runif(1)) < min(0,logr)){
    ztmp <- z;
    vtmp <- v;
    nktmp <- nk;
    mutmp <- mu;
    sigmatmp <- Sigma
    
    z[ztmp==kchoose[1]] <- kchoose[2]
    z[ztmp==kchoose[2]] <- kchoose[1]
    nk[kchoose[1]] <- nktmp[kchoose[2]]
    nk[kchoose[2]] <- nktmp[kchoose[1]]
    mu[,kchoose[1]] <- mutmp[,kchoose[2]]
    mu[,kchoose[2]] <- mutmp[,kchoose[1]]
    Sigma[,,kchoose[1]] <- sigmatmp[,,kchoose[2]]
    Sigma[,,kchoose[2]] <- sigmatmp[,,kchoose[1]]
    if (move>(1/2)){
      v[kchoose[1]] <- vtmp[kchoose[2]]
      v[kchoose[2]] <- vtmp[kchoose[1]]
      weights <- weightsupdate(Kstar,v)
    }
    #if (move>(2/3)){
    #weightstmp <- weights;vtmp <- v
    #weights[kchoose[1]] <- weightstmp[kchoose[2]]*weightscomb*R1/weightsdash
    #weights[kchoose[2]] <- weightstmp[kchoose[1]]*weightscomb*R2/weightsdash
    #v[kchoose[2]] <- weights[kchoose[2]]/prod(1-v[1:(kchoose[2]-1)])
    #v[kchoose[1]] <- weights[kchoose[1]]/prod(1-v[1:(kchoose[1]-1)])
    #weights <- weightsupdate(Kstar,v)
    #}
  }
  output <- list(z=z,nk=nk,v=v,weights=weights,mu=mu,Sigma=Sigma)
  return(output)
}

updatealpha <- function(n, alpha, Kstar,
                        a=1, b=1){
  a <- 1
  b <- 1
  for (t in 1:50){
    phi1  <- rbinom(n = 1,
                    size = 1, 
                    prob = n/(n+alpha))
    
    phi   <- rbeta(n = 1,
                   shape1 = alpha+1,
                   shape2 = n)
    
    alpha <- rgamma(1,
                    shape = a+Kstar-phi1,
                    rate = (b-log(phi)))
  }
  return(alpha)
}


fit_dpm <- function(y,
                    MCMC = list(niter = 1000,
                                thin = 1,
                                counter = 1),
                    K = 10){
  
  r <- dim(y)[1]
  n <- dim(y)[2]
  
  #set hyperparameters
  hypers <- list(
    b0 = apply(y, 1, mean),
    N0 = 0.01,
    c0 = r + 1,
    C0 = 0.75 * cov(t(y)),
    alpha.a = 1,
    alpha.b = 1)
  
  #Initialisation
  #clustering structure 
  Kstar <- K
  tmp   <- kmeans(t(y), K) # naive k-means approach to initialise
  z     <- as.vector(tmp$cluster)
  nk    <- table(z)
  
  #cluster means (mu)
  mu <- t(tmp$centers)
  
  #cluster var-covar matrices (Sigma)
  Sigma <- array(0, c(r, r, K))
  
  for(k in 1:K){
    Sigma[ , , k] <- with(hypers, riwish(c0,C0))
  }
  
  #DP concentration parameter (alpha)
  alpha <- 1
  
  #stick breaking weights (v) and cluster weights (weights)
  v <- rep(0, K)
  v[1] <- rbeta(1, 1, alpha)
  weights <- v
  for(k in 2:K) {
    v[k] <- rbeta(1, 1, alpha)
    weights[k] <- v[k] * prod(1 - v[1:(k - 1)])
  }
  rm(tmp)
  
  
  
  #set up MCMC information
  
  MCMC.traces <-
    list(
      z = array(0, c(n, MCMC$niter / MCMC$thin)),
      alpha = rep(0, MCMC$niter / MCMC$thin),
      K = rep(0, MCMC$niter / MCMC$thin)
    )
  
  #set up progress bar
  pbc<-1;
  
  pb <-
    txtProgressBar(style = 3,
                   title = "MCMC progress",
                   #label = "0% done",
                   min = 0,
                   max = 1,
                   initial = 0)
  
  for (t in 1:MCMC$niter){
    #STEP 1: Determine current number of alive components
    Kstar <- max(z)
    nk <- nk[1:Kstar]
    
    #STEP 2: Update stick breaking weights for alive components
    v <- vupdate(Kstar, nk, alpha, n)
    weights <- weightsupdate(Kstar, v)
    
    #STEP 3: Update MVN parameters for alive components
    out<-muSigmaupdate(y,z,nk,hypers,r,Kstar)
    mu<-out$mu
    Sigma<-out$Sigma
    
    #STEP 4:label switching move
    out <- labelswitch(Kstar, z, nk, weights, v, alpha, mu, Sigma)
    z<-out$z
    nk<-out$nk
    v<-out$v
    weights<-out$weights
    mu<-out$mu
    Sigma<-out$Sigma
    
    #STEP 5:Update u
    u<-sapply(1:n,function(t) runif(1,0,weights[z[t]]))
    
    #STEP 6: Update alpha
    alpha<-updatealpha(n,alpha,Kstar,
                       hypers$alpha.a,hypers$alpha.b)
    
    #STEP 7: Determine if more components are required (potentials)
    out<-determineKP(min(u),weights,Kstar,v,alpha)
    KP<-out$KP
    Ktotal<-Kstar+KP
    if(KP>0){v<-out$v; weights<-out$weights}
    
    #STEP 8: Update theta for any potential components
    if (KP>0){
      mu<-cbind(mu,array(0,c(r,KP)))
      for (k in (Kstar+1):Ktotal){
        Sigma<-array(c(Sigma,matrix(0,r,r)),c(r,r,k))
        Sigma[,,k]<-riwish(hypers$c0,hypers$C0)
        mu[,k]<-rmvnorm(1,hypers$b0,Sigma[,,k]/hypers$N0)
      } #end of for loop
    }#end of if statement
    
    #STEP 9: Update z
    out<-zupdate(Ktotal,n,y,u,weights,mu,Sigma)
    nk<-out$nk
    z<-out$z
    
    #CLEAN UP: order by decreasing z
    ztmp<-z;counter<-1;
    tmp<-unique(ztmp)
    for(ii in tmp){
      z[ztmp==ii]<-counter
      counter<-counter+1
    }
    weights<-weights[tmp]
    v<-v[tmp]
    mu<-mu[,tmp]
    Sigma<-Sigma[,,tmp]
    nk<-nk[tmp]
    rm(tmp)
    
    #store traces
    if ((t/MCMC$thin)==MCMC$counter){
      MCMC.traces$z[,t]<-norm.label(z) #relabelled to 1:length(unique(z))
      MCMC.traces$alpha[t]<-alpha
      MCMC.traces$K[t]<-sum(nk>0) #number of occupied components
      MCMC$counter<-MCMC$counter+1
    }
    
    #update progress bar
    if ( t %% 10 == 0){
      #Sys.sleep(0.1)
      #info <- sprintf("%d%% done", round((t/MCMC$niter)*100))
      setTxtProgressBar(pb,t/MCMC$niter)
      #pbc<-pbc+1
    }
    
  } #end of t loop
  
  return(MCMC.traces)
  
}

spike.smoother <- function(x,
                           y,
                           xnew=x,
                           k=20){
  
  require(mgcv)
  
  dat <- data.frame(x = unlist(x),
                    y = c(unlist(y)))
  
  my.spline <- gam(y ~ s(x, bs="ps", k=k), data=dat)
  #my.spline <- loess(y ~ x, span=0.1)
  
  my.smooth <- predict(my.spline, 
                       data.frame(x=unlist(xnew)))
  
  return(my.smooth)
  
}