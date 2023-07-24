library(plyr)
library(doParallel)

# To install package ESS from biostats required for power prior calculations
# devtools::install_github("DKFZ-biostats/ESS")
library(ESS)


###############
# 2 dimensional numeric integration (solution using package cubature had been instable)
my.integrate2 <- function(f, # function to be integrated
                          grid.length = 50, # number of subintervals for integration
                          min.q = 1e-6, # lower quantile of normal distribution to specify width of grid for integration (1-min.q upper quantile), small for a wide grid
                          data.mean1, 
                          sigma1,
                          n1,
                          data.mean2, 
                          sigma2,
                          n2){
  nq1 = qnorm(min.q) * sigma1/sqrt(n1)
  grid1 = seq(data.mean1 + nq1, data.mean1 - nq1, length.out = grid.length)
  d.grid1 = diff(grid1[1:2])
  nq2 = qnorm(min.q) * sigma2/sqrt(n2)
  grid2 = seq(data.mean2 + nq2, data.mean2 - nq2, length.out = grid.length)
  d.grid2 = diff(grid2[1:2])
  g = expand.grid(grid1,grid2)
  sum(d.grid1 * d.grid2 * f(g[,1], g[,2]))
}



##############
# posterior probabilities for differences P(treat-control <=q|.)

  # empirical Bayes power prior
    pdiffPowerpriorEB <- Vectorize(
      function(
        treat.xbar,  # mean in active
        treat.n,     # sample size in active
        treat.sd,    # std dev. in active (assumed known)
        control.xbar, # mean in current control arm
        control.n, # sample size in current control
        control.sd,# std dev. in current control (assumed known)
        hist.mu,      # mean in historical control data
        hist.n,       # sample size in historical control data
        hist.sd,     # std dev. in historical control data
        q,              # quantile (threshold)
        lower.tail = TRUE #logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
      ) {
        require(RBesT)
        require(ESS)
        post.control= suppressMessages(
          postmix(
            priormix = as.powerprior(mixnorm(comp1 = c(1, hist.mu, hist.sd/sqrt(hist.n)), 
                                             sigma = control.sd)   # ??? or should that be hist.sd???
            ),
            m = control.xbar,
            se = control.sd/sqrt(control.n) 
          )
        )
        post.treat = mixnorm(comp1=c(1, treat.xbar, treat.sd/sqrt(treat.n)), 
                             sigma = treat.sd)    # assumes flat prior
        
        pmixdiff(
          mix1 = post.treat,
          mix2 = post.control,
          q = q, 
          lower.tail = lower.tail
        )
      }, 
      vectorize.args = c("treat.xbar", "control.xbar"))
  
  # power prior with fixed power parameter delta  
    pdiffPowerpriorFix <- Vectorize(
      function(
        treat.xbar,  # mean in active
        treat.n,     # sample size in active
        treat.sd,    # std dev. in active (assumed known)
        control.xbar, # mean in current control arm
        control.n,  # sample size in current control
        control.sd, # std dev. in current control (assumed known)
        hist.mu,      # mean in historical control data
        hist.n,       # sample size in historical control data
        hist.sd,     # std dev. in historical control data
        delta = 1, # fixed power parameter 0<=delta<=1, full borrwoing for delta=1, no borrowing delta=1e-100
        q,              # quantile (threshold)
        lower.tail = TRUE #logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
      ) {
        require(RBesT)
        hist.n <- hist.n * delta    # FIXED POWER PARAMETER -> SIMPLY multiply historical sample size by power parameter delta
        
        post.control= postmix(
          priormix = mixnorm(comp1=c(1, hist.mu, hist.sd/sqrt(hist.n)), 
                             sigma = control.sd),
          m = control.xbar,
          se = control.sd/sqrt(control.n) 
        )
        
        post.treat = mixnorm(comp1=c(1, treat.xbar, treat.sd/sqrt(treat.n)), 
                             sigma = treat.sd)    # assumes flat prior
        
        pmixdiff(
          mix1 = post.treat,
          mix2 = post.control,
          q = q, 
          lower.tail = lower.tail
        )
      }, 
      vectorize.args = c("treat.xbar", "control.xbar"))
  


###############
# Calculate Power 
    
# Without borrowing: Analytical power calculation assuming equal sample size in both groups
  Powerwo_2sample = function(treat.mu,
                             control.mu,
                             treat.sigma,
                             control.sigma,
                             n,
                             alpha) {
    
    pnorm(sqrt(n) * (treat.mu - control.mu) / sqrt(treat.sigma^2 + control.sigma^2) - qnorm(1 - alpha))
  }
    

# With borrowing        
  PowerFix_2sample <- Vectorize(
    function(
      treat.mu,  # true mean in active
      treat.n,     # sample size in active
      treat.sd,    # std dev. in active (assumed known)
      control.mu,  # true mean in control 
      control.n,     # sample size in control 
      control.sd,    # std dev. in control  (assumed known)
      hist.mu,      #  true mean in historical control
      hist.n,       # sample size in historical control
      hist.sd,     # std dev. in historical control
      muE=0, cE=.975,      # P(treat-control>muE|.)> cE
      pprior,   # function, e.g. pdiffPowerpriorEB, pdiffPowerpriorFix
      subdivisions = 150, # the maximum number of subintervals for integration
      min.q = 1e-6, # lower quantile of normal distribution to specify width of grid for integration (1-min.q upper quantile), small for a wide grid. 
      ...  # arguments passed to pprior, in particular to pRobustmixture() (ie. mixtureWeight etc) and delta in case of pPowerpriorFix
    ){
      CALL = as.list(match.call(expand.dots = F))
      dots <- CALL$...
      
      Fun <- function(Xbar.treat, Xbar.control) {
        dnorm(Xbar.treat,
              treat.mu,
              sd=treat.sd/sqrt(treat.n)) * 
          dnorm(Xbar.control,
                control.mu,
                sd=control.sd/sqrt(control.n)) * 
        (do.call(pprior,
                 c(list(
                   treat.xbar=Xbar.treat,  
                   treat.n=treat.n,  
                   treat.sd=treat.sd,   
                   control.xbar=Xbar.control,  
                   control.n=control.n,  
                   control.sd=control.sd,   
                   hist.mu=hist.mu,  
                   hist.n=hist.n, 
                   hist.sd=hist.sd,     
                   q=muE,
                   lower.tail=FALSE
                 ),
                 dots)
        ) > cE) 
      }
      
      my.integrate2(
        f = Fun,
        min.q = min.q,
        grid.length = subdivisions,
        data.mean1 = treat.mu,
        sigma1 = treat.sd,
        n1 = treat.n,
        data.mean2 = control.mu,
        sigma2 = control.sd,
        n2 = control.n
      )
      
    },
    vectorize.args = c("treat.mu", "control.mu", "hist.mu"))




  PowerRandomSim_2sample <- Vectorize(
    function(
      treat.mu,  # true mean in active
      treat.n,     # sample size in active
      treat.sd,    # std dev. in active (assumed known)
      control.mu,  # true mean in control 
      control.n,     # sample size in control 
      control.sd,    # std dev. in control  (assumed known)
      hist.mu,      #  true mean in historical control
      hist.n,       # sample size in historical control
      hist.sd,     # std dev. in historical control
      muE=0, cE=.975,      # P(treat-control>muE|.)> cE
      pprior,   # function, e.g. pdiffPowerpriorEB, pdiffPowerpriorFix
      n.sim = 1e4,
      cpus=8, # number of cpus for parallelization, 1 for no parallelization
      ...  # arguments passed to pprior, in particular to pRobustmixture() (ie. mixtureWeight etc) and delta in case of pPowerpriorFix
    ){
      CALL = as.list(match.call(expand.dots = F))
      dots <- CALL$...
      
      if (cpus>1) {
        registerDoParallel(cpus)
        on.exit(stopImplicitCluster())
      }
      Fun <- function(Xbar.treat, Xbar.control, Xbar.hist) {
        (do.call(pprior,
                 c(list(
                   treat.xbar=Xbar.treat,  
                   treat.n=treat.n,  
                   treat.sd=treat.sd,   
                   control.xbar=Xbar.control,  
                   control.n=control.n,  
                   control.sd=control.sd,   
                   hist.mu=Xbar.hist,  
                   hist.n=hist.n, 
                   hist.sd=hist.sd,     
                   q=muE,
                   lower.tail=FALSE
                 ),
                 dots)
        ) > cE) 
      }
      
      random.draws.control = rnorm(n = n.sim,
                                   mean = control.mu,
                                   sd = control.sd/sqrt(control.n))
      random.draws.hist = rnorm(n = n.sim,
                                mean = hist.mu,
                                sd = hist.sd/sqrt(hist.n))
      random.draws.treat = rnorm(n = n.sim,
                                 mean = treat.mu,
                                 sd = treat.sd/sqrt(treat.n))
      
      mean(unlist(
        mlply(
          data.frame(
            Xbar.treat = random.draws.treat, 
            Xbar.control = random.draws.control,
            Xbar.hist = random.draws.hist),
          Fun,
          .progress = ifelse(cpus>1, "none", "text"),
          .parallel = (cpus>1)
        )
      ))
    },
    vectorize.args = c("treat.mu", "control.mu", "hist.mu"))
  

