## let's simulate some data


x <- seq(0, 90, by=1)

groups <- 6 # how many true groups are we simulating from?


make.data <- function(x, 
                      df,
                      bdeg,
                      n.groups,
                      n.obs,
                      seed=102){
  
  require(mvtnorm)
  require(tidyverse)
  require(splines)
  require(rospca)
  
  
  X <- bs(x = x, df = df, degree = bdeg)
  
  set.seed(seed)
  
  beta.0 <- rnorm(n=df)
  
  beta <- t(rmvnorm(n = groups, 
                    mean = beta.0,
                    sigma = diag(runif(n=df,
                                       max=0.15))))
  
  
  beta.df <- data.frame(beta) %>%
    gather(group, beta) %>%
    mutate(group = parse_number(group)) %>%
    group_by(group) %>%
    mutate(idx = 1:n())
  
  curves <- data.frame(
    group = sample(1:groups,
                   n.obs, 
                   replace=T)) %>%
    mutate(curve = 1:n())
  
  sim.dat <- 
    inner_join(curves, beta.df) %>%
    mutate(beta_ = rnorm(n=nrow(.),
                         mean = beta,
                         sd = 0.2)) %>%
    # maybe we want to tune that sd parameter
    split(.$curve) %>%
    map_df( ~ data.frame(y=X %*% .$beta_,
                         x=x) %>%
              mutate(y = y-mean(y)),
            .id="curve") %>%
    mutate(curve = parse_number(curve)) %>%
    inner_join(curves)
  
  return(sim.dat)
  
}

sim.dat <- make.data(x = x, df = 8, bdeg = 3,
          n.groups = 6, n.obs = 200, seed = 102)


ggplot(data=sim.dat, aes(x=x, y=y)) +
  geom_line(aes(group=curve,
                color=factor(group)))

# now we need to break down to the 4 principle components

sim.pca <- sim.dat %>%
  select(curve, y, x) %>%
  spread(key = curve, value = y) %>%
  select(-x) %>%
  data.matrix(.) %>%
  robpca(k=4)

sim.pca$scores

sim.pca$scores %>%
  data.frame %>%
  mutate(x=x) %>%
  gather(key, value, -x) %>%
  ggplot(data=., aes(x=x, y=value)) +
  geom_line(aes(group=key, color=factor(key)))

sim.pca$loadings # the data we want to use

ggplot(data=sim.pca$loadings %>%
         as.data.frame %>%
         mutate(group = curves$group ),
       aes(x=PC1, y=PC2)) +
  geom_point(aes(color=factor(group),
                 shape=factor(group))) +
  coord_equal()
