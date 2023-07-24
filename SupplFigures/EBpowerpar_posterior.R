########################################################################
### 1 sample, borrowing with EB,
### generate plots with EB power parameter and posterior probability P(theta > 0 |current and historical data)
###########################

rm(list = ls())
library(RBesT)
# To install package ESS from biostats required for power prior calculations
# devtools::install_github("DKFZ-biostats/ESS")
library(ESS)
source("../AuxillaryFunctions/functions_1sample.R")

powerparameter2 = function(m.hist,
                           n.hist,
                           m.current = m.current,
                           n.current = n.current,
                           sigma = sd.current) {
  sapply(m.hist,
         function(m.hist.i)
           unname(
             powerparameter(
               postmix(flat,
                       m = m.hist.i,
                       n = n.hist),
               m = m.current,
               n = n.current,
               sigma = sd.current
             )
           ))
}

m.current = c(-3000:10000) / 10000 # current mean
n.current = 25                     # sample size current data
n.hist = 1000                      # sample size historical data
SIGMA = sigma = 1                  # standard deviation of current and of historical data
sd.hist = sd.current = 1           # standard deviation of current and of historical data

muE = 0                            # P(x>muE|.)> cE
cE = .975
H0MU = 0                           # where I evaluate type I error rate
theta1 = 0.5


for (my.j in c(0:10)) {
  setEPS()
  postscript(
    paste("Extreme_EB_delta_posterior_sameScale", my.j, ".eps"),
    width = 8,
    height = 8
  )
  
  par(mfcol = c(2, 2))
  for (my.k in c(0, 1)) {
    m.hist = (2 * my.j + my.k) / 100  # historical mean
    
    flat <- mixnorm(c(1, 0, Inf),
                    sigma = sd.hist)
    
    # "historical posterior"
    post.hist = postmix(flat,
                        m = m.hist,
                        n = n.hist)
    
    ppar = powerparameter2(
      m.hist = m.hist,
      n.hist = n.hist,
      m.current = m.current,
      n.current = n.current,
      sigma = sd.current
    )
    colnames(ppar) = m.hist
    rownames(ppar) = m.current
    
    plot(
      m.current,
      ppar,
      type = 'l',
      lwd = 2,
      las = 1,
      xlab = "",
      ylab = expression(paste("EB ", delta)),
      main = substitute(paste("EB ",
                              delta, ", ",
                              bar(d)[E], "=", m),
                        list(m = m.hist))
    )
    
    pprob = pPowerpriorEB(
      current.xbar = m.current, # mean in current data
      current.n = n.current, # sample size in current data
      current.sd = SIGMA, # std dev. in current data (assumed known)
      hist.mu = m.hist, # mean in historical data
      hist.n = n.hist, # sample size in historical data
      hist.sd = sigma, # std dev. in historical data  (assumed known)
      q = muE, # quantile (threshold)
      lower.tail = FALSE       #logical; if TRUE (default), probabilities are P[X <= x] otherwise, P[X > x].
    )
    
    #T1E rate with borrowing, alpha_B
    alphaB = PowerFix(
      pprior = pPowerpriorEB, # type I error when borrowing from this historical data set
      current.mu = H0MU, # true mean in current data
      current.n = n.current, # sample size in current data
      current.sd = SIGMA, # std dev. in current data (assumed known)
      hist.mu = m.hist, # mean in historical data
      hist.n = n.hist, # sample size in historical data
      hist.sd = sigma, # std dev. in historical data (assumed known)
      muE = muE,
      cE = cE
    )
    # power with borrowing
    powerwEB =
      PowerFix(
        pprior = pPowerpriorEB, # power when borrowing from this historical data set
        current.mu = theta1, # true mean in current data
        current.n = n.current, # sample size in current data
        current.sd = SIGMA, # std dev. in current data (assumed known)
        hist.mu = m.hist, # mean in historical data
        hist.n = n.hist, # sample size in historical data
        hist.sd = sigma, # std dev. in historical data(assumed known)
        muE = muE,
        cE = cE
      )
    
    
    plot(
      m.current,
      pprob,
      type = 'l',
      lwd = 2,
      las = 1,
      xlab = "current observed mean",
      ylab = expression(paste("P(", theta , " > 0 |current and historical data)")),
      main = substitute(
        paste("Posterior probability, ",
              bar(d)[E], "=", m),
        list(m = m.hist)
      ),
      ylim = c(0, 1)
    )
    abline(h = 0.975, col = "blue", lwd = 2)
    mtext(
      "0.975",
      side = 2,
      at = 0.93,
      col = "blue",
      las = 1,
      line = -3,
      cex = 0.8
    )
    text(
      0.3,
      0.4,
      adj = 0,
      substitute(paste(alpha[B](d[E]), "=", m),
                 list(m = round(alphaB, 3))),
      col = "red",
      cex = 1.2
    )
    text(
      0.3,
      0.3,
      adj = 0,
      paste("powerwEB = ", round(powerwEB, 3)),
      col = "green",
      cex = 1.2
    )
  }
  
  dev.off()
}
