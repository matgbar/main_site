#Predicting Sunspots with Gaussian Process Models - More Time Series Applications: 
library(tidyverse)
library(rstan)
library(bayesplot)
library(rstanarm)

options(mc.cores=parallel::detectCores())
rstan_options(auto_write = TRUE)

data.folder<-'~/GitHub/main_site/DFs'
code.folder<-'~/GitHub/main_site/Code_base'

sun.dat<-read.table(
  paste0(
    data.folder,
    '/SN_m_tot_V2.0.txt'
  ),
  header = FALSE
  )

#Removing extraneous columns:
sun.dat<-sun.dat[,1:4]

#Naming Variables:
colnames(sun.dat)<-c(
  'Year', 
  'Month', 
  'Time',
  'Spots'
  )

#Exploratory Analyses: 
g1<-ggplot(
  data = sun.dat, 
  aes(
    x=Time, 
    y=Spots
    )
  )+
  geom_line()

g1

hist(sun.dat$Spots)

#Spectral/FFT approach to decomposing time series
del<-1/12 # sampling interval
x.spec <- spectrum(sun.dat$Spots,
                   log="no", 
                   span=10,
                   plot=FALSE
                   )

spx <- x.spec$freq/del
spy <- 2*x.spec$spec

plot(
  spy~spx,
  xlab="frequency",
  ylab="spectral density",
  type="l", 
  xlim = c(0,1)
  )

dom.freq<-spx[which.max(spy)]
dom.freq  #This is .1 or once every ten years - a good starting value for my model...

#Let's start with the simple frequency model. (may try to identify the second mini peak there)
#Three data points per cycle (assuming cycle is 10 years)
#Can think of this as a three Hz sampling rate 

Y<-sun.dat$Spots[seq(1, length(sun.dat[,1]), by=40)]
X<-sun.dat$Time[seq(1, length(sun.dat[,1]), by=40)]
Xpred<-seq(2016, 2017, by=1/12)
N1<-length(X)
N2<-length(Xpred)

stan.dat<-list(
  N1=N1, 
  N2=N2,
  X=X, 
  Y=Y,
  Xp=Xpred
  )

pars.to.monitor<-c(
  paste0('a', 1:2), 
  paste0('r', 1:3),
  'sigma_sq',
  'Ypred'
  )

set.seed(513)

init.list<-list(
  list(
    r3 = rnorm(1, 24, 1), 
    sigma_sq = rnorm(1,1500, 10)
  ),
  list(
    r3 = rnorm(1, 24, 1), 
    sigma_sq = rnorm(1,1500, 10)
  ),
  list(
    r3 = rnorm(1, 24, 1), 
    sigma_sq = rnorm(1,1500, 10)
  ),
  list(
    r3 = rnorm(1, 24, 1), 
    sigma_sq = rnorm(1,1500, 10)
  ),
  list(
    r3 = rnorm(1, 24, 1), 
    sigma_sq = rnorm(1,1500, 10)
  ),
  list(
    r3 = rnorm(1, 24, 1), 
    sigma_sq = rnorm(1,1500, 10)
  )
)

fit.stan10<-stan(
  file=paste0(code.folder, '/GP_sun.stan'),
  data = stan.dat, 
  warmup = 10000,
  iter = 20000,
  thin = 2,
  refresh=2000,
  chains = 6,
  init = init.list, 
  pars = pars.to.monitor,
  control = list(
    adapt_delta = .95, 
    max_treedepth = 15
    )
  )

traceplot(
  fit.stan10, 
  pars=c(paste0('a', 1:2), 
         paste0('r', 1:3), 
         'sigma_sq'
         )
  )

pairs(
  fit.stan10, 
  pars=c(paste0('a', 1:2), 
         paste0('r', 1:3), 
         'sigma_sq'
         )
  )

#Following examination of traceplots are still not great.
#Going to up the adapt delta

summary(fit.stan10)$summary
