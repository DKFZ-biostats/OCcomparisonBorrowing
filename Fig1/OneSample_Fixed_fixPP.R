########################################################################
########################################################################
########################################################################
### 1 sample, borrowing with fixed delta power prior, external data is fixed
###########################
rm(list = ls())
library(tidyverse)
library(parallel)

source("../AuxillaryFunctions/functions_1sample.R")

alpha <- 0.025

theta1 <- 0.5  #where power is evaluated, i.e. mean for current data
H0MU <- 0       #where type I error rate is evaluated

N <- 25
SIGMA <- 1   # assumed as known...

muE <- 0  ##   P(x>muE|.)> cE
cE <- .975

sigma <- 1
n <- 20

nsim <- 100

delta <- 0.5 # (fixed) power parameter

results_inc <- list()
set.seed(1, kind = "L'Ecuyer-CMRG")

i <- 0
for (thetaE in c(0.0, 0.5))  {
  #mean of external
  i <- i + 1
  ##start loop for generating historical data with this one parameter
  results_inc[[i]] <- as.data.frame(t(simplify2array(
    mclapply(1:nsim, 
             function(my.i) {
               histdata <- rnorm(n = n, mean = thetaE, sd = sigma)
               mean.histdata <- mean(histdata)
               
               # evaluate alphaB:
               alphaB <- PowerFix(
                 pprior = pPowerpriorFix,   # type I error rate when borrowing from this historical data set
                 current.mu = H0MU, # true mean in current data
                 current.n = N,   # sample size in current data
                 current.sd = SIGMA, # std dev. in current data (assumed known)
                 hist.mu = mean.histdata,  # mean in historical data
                 hist.n = n,  # sample size in historical data
                 hist.sd = sigma,  # std dev. in historical data (assumed known)
                 muE = muE,
                 cE = cE,
                 delta = delta ,
                 integrate = F,  # TRUE: use stats::integrate(); FALSE: custom integration function my.integrate()
                 subdivisions = 1e6,
                 min.q = 1e-9# only if integrate == FALSE, lower quantile of normal distribution based on current trial to specify width of grid for integration (1-min.q upper quantile), small for a wide grid. Defaults to 1e-6
                 
               )
               
               #evaluate power with borrowing
               powerwFix <- PowerFix(
                 pprior = pPowerpriorFix,  # power when borrowing from this historical data set
                 current.mu = theta1,  # true mean in current data
                 current.n = N, # sample size in current data
                 current.sd = SIGMA,  # std dev. in current data (assumed known)
                 hist.mu = mean.histdata,  # mean in historical data
                 hist.n = n, # sample size in historical data
                 hist.sd = sigma, # std dev. in historical data (assumed known)
                 muE = muE,
                 cE = cE,
                 delta = delta,
                 integrate = F, # TRUE: use stats::integrate(); FALSE: custom integration function my.integrate()
                 subdivisions = 1e6,
                 min.q = 1e-9
               )
               
               # frequentist power for test calibrated to alphaB
               powerwo_alphaB <- pnorm(sqrt(N) * (theta1 - H0MU) / SIGMA - qnorm(1 - alphaB))  
               
               #compare power for w and w/o borrowing calibrated to type I error rate with borrowing
               powerwFixdiff <- powerwFix - powerwo_alphaB
               
               return(as.numeric(
                 c(thetaE,
                   mean.histdata,
                   alphaB,
                   powerwo_alphaB,
                   powerwFix,
                   powerwFixdiff
                 )
               ))
               
             }, 
             mc.cores = 8)  # number of CPUs, set to 1 if only 1 core
  )))
}

results <- bind_rows(results_inc)
names(results) <- c("thetaE",
                   "mean.histdata",
                   "alphaB",
                   "powerwo_alphaB",
                   "powerwFix",
                   "powerwFixdiff")
write.csv(results, 
          file = "OneSampleFix_delta0.5_alpha0.025_nE20.csv")
