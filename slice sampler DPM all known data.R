library(R.matlab)

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(viridis)
library(ggalt)

#load DPM functions
source("slice sampler DPM functions.R")

# data sourced from paper
# http://www.sciencedirect.com/science/article/pii/S0165027009004506
# data available http://www2.le.ac.uk/departments/engineering/research/bioengineering/neuroengineering-lab/software

data.files <- dir("./data/leacukcsn/",
                  pattern="mat",
                  full.names = T)

data.n <- length(data.files)

MCMC.all <- vector(mode = "list",
                   length = data.n)

for (i in 1:5){ # purrr approach would be lovely
  
  a <- readMat(data.files[i])
  
  # import EC recording as .Rdata
  
  #Extract sampling frequency (kHz)
  sampl.freq<-1/a$samplingInterval
  
  #Detrmine number of spikes and spike times
  N<-length(a$spike.times[[1]][[1]])
  
  #check lenght of recording in minutes (should be 2)
  T<-(length(a$data)/sampl.freq)/60000
  
  
  #spike times
  spike.times<-a$spike.times[[1]][[1]]
  
  #spike type (0,1 or 2)
  spike.type<-a$spike.class[[1]][[1]]
  
  #store individual spike trajectories (sup to 2.5ms after spike time; not aligned yet)
  spike.data <- sapply(1:64, function(t)
    a$data[1, (spike.times + t)])
  
  ## need to smooth the data
  
  spike.data %<>% data.frame
  spike.data$type <- t(spike.type)
  spike.data$ID <- 1:nrow(spike.data)
  
  spike.dat <- gather(spike.data,
                      time, 
                      value,
                      -c(ID, type)) %>%
    mutate(time = parse_number(time)) %>%
    split(.$ID) %>%
    map( ~ mutate(.x, 
                  value.smooth =
                    spike.smoother(x=.$time,
                                   y=.$value))) %>%
    bind_rows
  
  #locate maximum amplitude of each sampled spike
  # max.amp<-apply(spike.data[,1:64],
  #                1,
  #                function(x) which.max(x))
  
  max.amp <- spike.dat %>%
    group_by(ID) %>%
    summarise(ind = which.max(value.smooth)) %>%
    ungroup %>%
    select(ind) %>%
    unlist %>%
    c
  
  
  #overwrite spike.data to obtain spikes aligned at max. amp at time 20 (figures in paper)
  spike.data <- sapply(-20:44,
                       function(t) a$data[1,(spike.times+max.amp)+t]) %>% 
    data.frame %>%
    mutate(type = t(spike.type),
           ID = 1:nrow(spike.data))
  
  spike.dat <- gather(spike.data,
                      time, 
                      value,
                      -c(ID, type)) %>%
    mutate(time = parse_number(time)) %>%
    split(.$ID) %>%
    map( ~ mutate(.x, 
                  value.smooth =
                    spike.smoother(x=.$time,
                                   y=.$value))) %>%
    bind_rows
  
  # save the processed data for later
  saveRDS(spike.dat, 
          file = gsub(pattern = "mat",
                      replacement = "RDS",
                      x = data.files[i]))
  # 
  # windows()
  # ggplot(data=spike.dat, 
  #        aes(x=time, y=value.smooth)) +
  #   geom_line(alpha=0.1, aes(group=ID)) +
  #   facet_wrap( ~ type, ncol=1)
  
  # robust PCA, was originally 4 PCs now we leave it unrestricted
  spike.pca <- spike.dat %>%
    select(ID, y=value.smooth, x=time) %>%
    spread(key = ID, value = y) %>%
    select(-x) %>%
    data.matrix(.) %>%
    robpca
  
  #set.seed(1337)
  # 
  # yinds <- sample(1:nrow(spike.pca$loadings),
  #                 replace = F,
  #                 size=100)
  
  # can modify this with the above if you don't want to use ALL the data
  #yinds <- 1:nrow(spike.pca$loadings)
  
  #y <- t(spike.pca$loadings[yinds, ])
  
  #MCMC.all[[i]] <- fit_dpm(y, K=10)
  
}

saveRDS(MCMC.all, "MCMC.traces.RDS")

spike.mean.all <- MCMC.all %>% 
  map("z") %>%
  map(~comp.psm(t(.x[, 1:(t-1)])))

zmap.all <- spike.mean.all %>%
  map(~maxpear(.x))

# for storage of solutions 
yt <- vector("list", data.n)

for (i in 1:data.n){
  spike.dat <- readRDS(gsub(pattern = "mat",
                            replacement = "RDS",
                            x = data.files[i]))
  
  spike.pca <- spike.dat %>%
    select(ID, y=value.smooth, x=time) %>%
    spread(key = ID, value = y) %>%
    select(-x) %>%
    data.matrix(.) %>%
    robpca
  
  set.seed(1337)
  # 
  # yinds <- sample(1:nrow(spike.pca$loadings),
  #                 replace = F,
  #                 size=100)
  
  y <- t(spike.pca$loadings)
  
  yt[[i]] <- data.frame(t(y),
                   z=zmap[[i]],
                   id=1:ncol(y)) %>% 
    bind_cols(spike.dat %>%
                group_by(ID, type) %>%
                distinct %>%
                select(-c(ID, type)))
  
  
}
