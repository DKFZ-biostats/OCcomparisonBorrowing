########################################################################
########################################################################
########################################################################
### 1 sample, borrowing with fixed delta power prior, fixed external data
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

delta = 0.5   #(fixed) power parameter

grid <- c(-300:400) / 100  # grid of mean d_E

with = sapply(grid, 
              function (dE) {
                PowerFix(
                  pprior = pPowerpriorFix, # power when borrowing from this historical data set
                  current.mu = theta1, # true mean in current data
                  current.n = N, # sample size in current data
                  current.sd = SIGMA, # std dev. in current data (assumed known)
                  hist.mu = dE, # mean in historical data
                  hist.n = n, # sample size in historical data
                  hist.sd = sigma, # std dev. in historical data (assumed known)
                  muE = muE,
                  cE = cE,
                  delta = delta ,
                  integrate = F, # TRUE: use stats::integrate(); FALSE: custom integration function my.integrate()
                  subdivisions = 5e3,
                  min.q = 1e-5
                )
              })

# evaluate alphaB:
alphaB <- sapply(grid, function (dE) {
  PowerFix(
    pprior = pPowerpriorFix, # type I error when borrowing from this historical data set
    current.mu = H0MU, # true mean in current data
    current.n = N, # sample size in current data
    current.sd = SIGMA, # std dev. in current data (assumed known)
    hist.mu = dE, # mean in historical data
    hist.n = n, # sample size in historical data
    hist.sd = sigma, # std dev. in historical data (assumed known)
    muE = muE,
    cE = cE,
    delta = delta  ,
    integrate = F, # TRUE: use stats::integrate(); FALSE: custom integration function my.integrate()
    subdivisions = 1e4,
    min.q = 1e-7 # only if integrate == FALSE, lower quantile of normal distribution based on current trial to specify width of grid for integration (1-min.q upper quantile), small for a wide grid. Defaults to 1e-6
  )
})


# power w/o borrowing
powerwo <- sapply(1:length(grid), 
                  function (i) {
                    dE = grid[i]
                    #compare power for w and w/o borrowing, taking the adjusted type I error rate
                    powerwo_alphaB = pnorm(sqrt(N) * (theta1 - H0MU) / SIGMA - qnorm(1 - alphaB[i]))  # frequentist power for adjusted alpha
                    powerwo_alphaB
                  })

# Power without borrowing at level alpha
powerwo_orig <- Powerwo_1sample(current.mu = theta1,
                                H0MU = H0MU, 
                                current.sigma = SIGMA,
                                current.n = N, 
                                alpha = alpha)


resultsFix = data.frame(
  dE = grid,
  alphaB = alphaB,
  powerwo = powerwo,
  powerwFix = with,
  powerwFixdiff = with - powerwo,
  powerwo_orig = powerwo_orig
)






write.csv(resultsFix, 
          file = "OneSampleFix_fixPP_delta0.5_nE20_dEgrid_alpha0.025.csv")

