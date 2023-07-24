########################################################################
########################################################################
########################################################################
### 2 sample, hybrid control, borrowing with Empirical Bayes power prior, random external data
###########################

rm(list = ls())
source("../AuxillaryFunctions/functions_2sample.R")


CONTROL.MU = 0       #current control theta_c. We can keep this fixed and only vary mean d_E
CONTROL.N = 12       #current control sample size
CONTROL.SIGMA = 1    #current control sigma (assumed known)

TREAT.N = 15         #current sample size in treatment arm
TREAT.SIGMA = 1      #current sigma in treatment arm(assumed known)

H0Diff = 0           #theta_c - theta_t=0
H1Diff = 1           #where to evaluate power: theta_c - theta_t=1

HIST.N = 10          #external control sample size
HIST.SIGMA = 1       #external control sigma (assumed known)

muE = 0  # P(x>muE|.)> cE
cE = .975

nplotpoints = 21     #how many points are plotted
minplotpoints = -0.7 #min on x-axis
maxplotpoints = 1.3  #max on x-axis
bplotpoints = (maxplotpoints - minplotpoints) / (nplotpoints - 1)
aplotpoints = minplotpoints - bplotpoints

Delta = rep(NA, nplotpoints) #this is difference between observed mean historical and theta_c (E of current control)

alphaBEB = rep(NA, nplotpoints)
powerwEB = rep(NA, nplotpoints)

for (my.i in c(1:nplotpoints)) {
  Delta[my.i] = aplotpoints + my.i * bplotpoints
  
  #T1E rate with borrowing, alpha_B
  alphaBEB[my.i] <-           
    PowerRandomSim_2sample(
      pprior = pdiffPowerpriorEB, # function, e.g. pdiffPowerpriorEB, pdiffPowerpriorFix
      treat.mu = Delta[my.i] + H0Diff, # true mean in current active, here under H0
      treat.n = TREAT.N, # sample size in current active
      treat.sd = TREAT.SIGMA, # std dev. in current active (assumed known)
      control.mu = Delta[my.i], # true mean in current control
      control.n = CONTROL.N, # sample size in current control
      control.sd = CONTROL.SIGMA, # std dev. in current control (assumed known)
      hist.mu = 0, # mean in historical control, oBdA=0
      hist.n = HIST.N, # sample size in historical control
      hist.sd = HIST.SIGMA, # std dev. in historical control (assumed known)
      muE = 0,
      cE = cE, # P(treat-control>muE|.)> cE
      n.sim = 1e5,
      cpus = 1 # number of cpus for parallelization, 1 for no parallelization
    )
  
  
  #power with borrowing
  powerwEB[my.i] <-             
    PowerRandomSim_2sample(
      pprior = pdiffPowerpriorEB, # function, e.g. pdiffPowerpriorEB, pdiffPowerpriorFix
      treat.mu = Delta[my.i] + H1Diff, # true mean in current active, here under H0
      treat.n = TREAT.N, # sample size in current active
      treat.sd = TREAT.SIGMA, # std dev. in current active (assumed known)
      control.mu = Delta[my.i], # true mean in current control
      control.n = CONTROL.N, # sample size in current control
      control.sd = CONTROL.SIGMA, # std dev. in current control (assumed known)
      hist.mu = 0, # mean in historical control, oBdA=0
      hist.n = HIST.N, # sample size in historical control
      hist.sd = HIST.SIGMA, # std dev. in historical control (assumed known)
      muE = 0,
      cE = cE, # P(treat-control>muE|.)> cE
      n.sim = 1e5,
      cpus = 1 # number of cpus for parallelization, 1 for no parallelization
    )
}

######

#power w/o borrowing at original level alpha; depends only on H1Diff and nothing else:
powerwo = Powerwo_2sample(
  treat.mu = H1Diff,
  control.mu = 0,
  treat.sigma = TREAT.SIGMA,
  control.sigma = CONTROL.SIGMA,
  n = TREAT.N,
  alpha = 1 - cE
)

#T1E rate maximized over all possible theta_c
max_alphaBEB = max(alphaBEB)

#power w/o borrowing but level alphaB; depends only on H1Diff and nothing else
powerwo_alphaBEB = Powerwo_2sample(
  treat.mu = H1Diff,
  control.mu = 0,
  treat.sigma = TREAT.SIGMA,
  control.sigma = CONTROL.SIGMA,
  n = TREAT.N,
  alpha = max_alphaBEB
)

powerwo_vec = rep(powerwo, nplotpoints)
powerwo_alphaBEB_vec = rep(powerwo_alphaBEB, nplotpoints)
powerEBdiff = powerwEB - powerwo_alphaBEB_vec


write.table(
  cbind(
    Delta,
    alphaBEB,
    powerwEB,
    powerwo_vec,
    powerwo_alphaBEB_vec,
    powerEBdiff
  ),
  file = "TwoSample_EB_histRandom_alpha0.025.csv",
  sep = ",",
  append = TRUE,
  quote = FALSE,
  col.names = c(
    "Delta",
    "alphaBEB",
    "powerwEB",
    "powerwo",
    "powerwo_alphaBEB",
    "powerEBdiff"
  ),
  row.names = FALSE
)
