# slice sampler for DP mixture model.
# Nicole White, 2015
# email: nm.white@qut.edu.au
# Sam Clifford, 2017, minor edits
rm(list=ls())

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ggalt)

#load DPM functions
source("slice sampler DPM functions.R")

# #read in data; dim =rxn
# y <- read_csv("./data/waves2.csv",col_names = F) %>%
#   as.matrix
# 
# # from sam's simulation
y <- t(sim.pca$loadings)

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
K     <- 10
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
MCMC <- list(niter = 1000,
             thin = 1,
             counter = 1)

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
  for(i in tmp){
    z[ztmp==i]<-counter
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

#Post-processing: Find MAP clustering based on PEAR.
sim.mean <- comp.psm(t(MCMC.traces$z[,1:t-1]))

# order y by the first PC
o.y <- order(y[1,])
o.y <- sim.dat %>%
  group_by(curve, group) %>%
  distinct %>%
  ungroup %>%
  arrange(group) %>%
  select(curve) %>%
  .$curve

# take posterior similarity matrix and convert to tidy data frame
sim.image <- sim.mean[o.y, o.y] %>%
  data.frame(., x=1:nrow(.)) %>%
  gather(y, z, -x) %>%
  mutate(y=parse_number(y)) 

pdf(file = "sim.image.pdf")
ggplot(data=sim.image, aes(x=x, y=y)) +
  geom_raster(aes(fill=z)) +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_fill_gradientn(
    colours = magma(n=5),
    name="Posterior similarity")
dev.off()

zmap <- maxpear(sim.mean,
                t(MCMC.traces$z[,1:t-1]),
                method="draws")$cl 

zmap

yt <- data.frame(t(y),
                 z=zmap,
                 id=1:ncol(y)) 

yt %<>% bind_cols(sim.dat %>%
                    group_by(curve, group) %>%
                    distinct %>%
                    select(-c(curve, group)))

gather(yt, key, value, -c(z, id)) %>%
  mutate_each(funs(parse_number), key) %>%
  na.omit %>%
  ggplot(data=., aes(x=z, y=value)) +
  geom_point(alpha=0.1,
             position = position_jitter(h=0,
                                        w=0.1)) +
  facet_wrap( ~ key) +
  theme_bw() +
  xlab("Cluster") +
  ylab("Principle component loading") +
  ggtitle("Cluster properties")

pdf(file = "sim.pc.pdf", width = 7, height=5)
ggplot(data=yt, aes(x=PC1,
                    y=PC2,
                    fill=factor(zmap),
                    color=factor(zmap))) +
  #geom_point(alpha=0.5) +
  scale_color_brewer(palette="Set2",
                     name="Most likely clustering") +
  scale_fill_brewer(palette="Set2",
                     name="Most likely clustering") +
  theme_bw() +
  #coord_equal() +
  theme(legend.position = "bottom") +
  geom_encircle(s_shape=1, expand=0, alpha=0.25) +
  xlab("PC1") +
  ylab("PC2") +
  geom_point(aes(shape=factor(group))) +
  scale_shape(name="True grouping")
dev.off()
