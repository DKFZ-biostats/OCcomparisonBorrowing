
#####################    
# Numerical integration 
# replaces stats::integrate() (optional) on window centered at data.mean with range based on standard error
my.integrate <- function(f, # function to be integrated
                         grid.length = 50, # number of subintervals for integration
                         min.q = 1e-6, # lower quantile of normal distribution to specify width of grid for integration (1-min.q upper quantile), small for a wide grid
                         data.mean, 
                         sigma,
                         n){
  nq = qnorm(min.q) * sigma/sqrt(n)
  grid = seq(data.mean + nq, data.mean - nq, length.out = grid.length)
  d.grid = diff(grid[1:2])
  sum(d.grid * f(grid))
}




###################
# posterior probabilities P(mu<=q|.)


# posterior probabilitiy with power prior 
pp_postprob <- Vectorize(function(Xbar,N,SIGMA, #current
                                  mu,n, sigma,  #historical
                                  muE, delta=NULL,
                                  lower.tail = TRUE){
  SE2 <- SIGMA^2/N 
  se2 <- sigma^2/n
  if (is.null(delta)) delta <- se2/(max((Xbar-mu)^2,
                                        SE2+se2)-SE2)
  
  postvar <- 1/( 1/(se2/delta) + 1/SE2 )
  mean <- postvar * ( mu/(se2/delta) + Xbar/SE2 )
  pnorm(muE,
        mean,
        sqrt(postvar) ,
        lower.tail = lower.tail)
  
},
vectorize.args = "Xbar")


# Convenience interfaces
  # empirical Bayes power prior
    pPowerpriorEB <- Vectorize(
      function(current.xbar,  # mean in current data
               current.n,     # sample size in current data
               current.sd,    # std dev. in current data (assumed known)
               hist.mu,       # mean in historical data
               hist.n,        # sample size in historical data
               hist.sd,       # std dev. in historical data
               q,             # quantile (threshold)
               lower.tail = TRUE #logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
               ) {
        pp_postprob(Xbar=current.xbar,N=current.n,SIGMA=current.sd, #current
                    mu=hist.mu,n=hist.n, sigma=hist.sd,  #historical
                    muE=q,
                    lower.tail= lower.tail)
      }, 
      vectorize.args = "current.xbar")

  # power prior with fixed power parameter delta  
    pPowerpriorFix <- Vectorize(
      function(current.xbar,  # mean in current data
               current.n,     # sample size in current data
               current.sd,    # std dev. in current data (assumed known)
               hist.mu,       # mean in historical data
               hist.n,        # sample size in historical data
               hist.sd,       # std dev. in historical data
               delta = 1,     # fixed power parameter 0<=delta<=1, full borrwoing for delta=1
               q,             # quantile (threshold)
               lower.tail = TRUE #logical; if TRUE (default), probabilities are P[X ≤ x] otherwise, P[X > x].
      ) {
         pp_postprob(Xbar = current.xbar, N = current.n, SIGMA = current.sd, #current
                    mu = hist.mu, n = hist.n, sigma = hist.sd,  #historical
                    delta = delta,
                    muE = q,
                    lower.tail = lower.tail)

      }, 
      vectorize.args = "current.xbar")
    

 
    
######################
# Calculate Power 
    
# Without borrowing       
  Powerwo_1sample = function(current.mu,
                             H0MU,
                             current.sigma,
                             current.n,
                             alpha) {
    pnorm(sqrt(current.n) * (current.mu - H0MU) / current.sigma - qnorm(1 - alpha))
  }
    
    
    
# With borrowing    
  # Power for fixed historical data, with P(x>muE|.)> cE
    PowerFix <- Vectorize(
      function(current.mu,  # true mean in current trial
               current.n,   # sample size in current trial
               current.sd,  # std dev. in current trial (assumed known)
               hist.mu,     #  "true" mean in historical trial
               hist.n,      # sample size in historical trial
               hist.sd,     # std dev. in historical trial
               muE,cE,      # P(x>muE|.)> cE
               pprior,      # function, pPowerpriorEB, pPowerpriorFix
               integrate = TRUE, # TRUE: use stats::integrate(); FALSE: custom integration function my.integrate()
               lower = -Inf, upper = Inf, # only if integrate == TRUE, the limits of integration
               subdivisions = 50, # the maximum number of subintervals for integration
               min.q = 1e-6, # only if integrate == FALSE, lower quantile of normal distribution based on current trial to specify width of grid for integration (1-min.q upper quantile), small for a wide grid. Defaults to 1e-6
               ...  # arguments passed to pprior, in particular delta in case of pPowerpriorFix()
      ){
        CALL = as.list(match.call(expand.dots = F))
        dots <- CALL$...
        
        Fun <- function(Xbar) {
          dnorm(Xbar,
                current.mu,
                sd=current.sd/sqrt(current.n)) * 
            (do.call(pprior,
                     c(list(current.xbar=Xbar,  
                            current.n=current.n,  
                            current.sd=current.sd,   
                            hist.mu=hist.mu,  
                            hist.n=hist.n, 
                            hist.sd=hist.sd,     
                            q=muE,
                            lower.tail=FALSE
                     ),
                     dots)
            ) > cE) 
        }
        
        if (integrate) {
          integrate(
            f = Fun,
            lower = lower,
            upper = upper,
            subdivisions = subdivisions
          )$value
        } else {
          my.integrate(
            f = Fun,
            min.q = min.q,
            grid.length = subdivisions,
            data.mean = current.mu,
            sigma = current.sd,
            n = current.n)
        }
      }, vectorize.args = "current.mu")
    
    
