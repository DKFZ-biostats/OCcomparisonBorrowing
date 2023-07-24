########################################################################
########################################################################
########################################################################
### 1 sample, borrowing with Empirical Bayes power prior, fixed external data
###########################

rm(list = ls())
source("../AuxillaryFunctions/functions_1sample.R")

alpha = 0.025

theta1 = 0.5  #where I evaluate power, i.e. mean for current
H0MU = 0      #where I evaluate type I error rate

N = 25        #sample size current data (n in paper)
SIGMA = 1     #current sigma, assumed as known...

muE = 0       # P(x>muE|.)> cE
cE = .975

sigma = 1     #historical sigma, assumed as known...
n = 20        #sample size historical data (n_E in paper)


results_inc = rep(NA, 5)

# grid of mean d_E
for (dE in c(-300:400) / 100){
  # evaluate alphaB:
  alphaB = PowerFix(
    pprior = pPowerpriorEB,  # type I error when borrowing from this historical data set
    current.mu = H0MU, # true mean in current data
    current.n = N, # sample size in current data
    current.sd = SIGMA, # std dev. in current data (assumed known)
    hist.mu = dE, # mean in historical data
    hist.n = n, # sample size in historical data
    hist.sd = sigma, # std dev. in historical data
    muE = muE,
    cE = cE
  )
  
  #evaluate power with borrowing
  powerwEB =
    PowerFix(
      pprior = pPowerpriorEB, # power when borrowing from this historical data set
      current.mu = theta1, # true mean in current data
      current.n = N, # sample size in current data
      current.sd = SIGMA, # std dev. in current data (assumed known)
      hist.mu = dE, # mean in historical data
      hist.n = n, # sample size in historical data
      hist.sd = sigma, # std dev. in historical data
      muE = muE,
      cE = cE
    )
  
  #compare power for w and w/o borrowing, taking the adjusted typeI error
  powerwo_alphaB = pnorm(sqrt(N) * (theta1 - H0MU) / SIGMA - qnorm(1 - alphaB))  # frequentist power for adjusted alpha
  
  powerwEBdiff = powerwEB - powerwo_alphaB
  results_inc = rbind(results_inc, 
                      as.numeric(c(dE, alphaB, powerwo_alphaB, powerwEB, powerwEBdiff))
                      )
  
}

resultsFix = data.frame(results_inc)
resultsFix = resultsFix[-1, ]
names(resultsFix) = c("dE", "alphaB", "powerwo", "powerwEB", "powerwEBdiff")

# Power without borrowing at level alpha
resultsFix$powerwo_orig = Powerwo_1sample(current.mu = theta1,
                                          H0MU = H0MU, 
                                          current.sigma = SIGMA,
                                          current.n = N, 
                                          alpha = alpha)

write.csv(resultsFix, 
          file = "OneSampleFix_EB_nE20_dEgrid_alpha0.025.csv")

