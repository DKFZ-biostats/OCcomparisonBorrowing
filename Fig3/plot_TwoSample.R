## plots Fig 3

library(latex2exp)
rm(list = ls())

dta_histFix = read.csv(file = "TwoSample_PPFix_histFix_delta0.5_alpha0.025.csv")


setEPS()
postscript(
  "alphaB_powerwPPFix_2sample_histFix.eps",
  width = 8,
  height = 8,
  pointsize = 12
)
par(mfrow = c(2, 1),
    mar = c(4, 4, 1, 4))

plot(
  dta_histFix$Delta,
  dta_histFix$alphaBPPFix,
  xlim = c(min(dta_histFix$Delta), max(dta_histFix$Delta)),
  ylim = c(0, 1.03),
  las = 1,
  xlab = TeX("$\\frac{\\theta_c - \\bar{d_E} }{\\sigma}$"),
  ylab = expression(Probability ~ to ~ reject ~ H[0]),
  type = "n",
  xaxs = "i",
  yaxs = "i"
)

mtext(
  side = 2,
  at = 0.95,
  las = 1,
  line = 3,
  text = "a",
  cex = 1.5
)

lines(dta_histFix$Delta,
      dta_histFix$alphaBPPFix,
      col = "red",
      lwd = 2)
abline(
  h = dta_histFix$powerwo_alphaBPPFix[1],
  lty = 1,
  lwd = 2,
  col = "green"
)
lines(dta_histFix$Delta,
      dta_histFix$powerwPPFix,
      col = "blue",
      lwd = 2)

abline(v = 0)
abline(h = 0.025, lwd = 2, lty = 2)
abline(h = dta_histFix$powerwo[1], lty = 3) #power w/o borrowing

legend(
  x = 2.9,
  y = 0.5,
  col = c("black", "green", "blue", "red", "black"),
  lty = c(3, 1, 1, 1, 2),
  cex = 0.8,
  bg = "white",
  legend = c(
    TeX("power w/o borrowing, $\\alpha=0.025$"),
    "power calibrated to borrowing",
    "power with borrowing",
    TeX("$\\alpha_B(\\theta_c=\\theta_t; d_E)$"),
    TeX("$\\alpha=0.025$")
  )
)

plot(
  dta_histFix$Delta,
  dta_histFix$alphaBPPFix,
  xlim = c(min(dta_histFix$Delta), max(dta_histFix$Delta)),
  ylim = c(-0.4, 1.03),
  las = 1,
  xlab = TeX("$\\frac{\\theta_c - \\bar{d_E} }{\\sigma}$"),
  ylab = expression(Probability ~ to ~ reject ~ H[0]),
  type = "n",
  xaxs = "i",
  yaxs = "i"
)

mtext(
  side = 2,
  at = 0.95,
  las = 1,
  line = 3,
  text = "b",
  cex = 1.5
)

lines(dta_histFix$Delta,
      dta_histFix$alphaBPPFix,
      col = "red",
      lwd = 2)
lines(dta_histFix$Delta,
      dta_histFix$powerPPFixdiff,
      col = "blue",
      lwd = 2)

abline(v = 0)
abline(h = 0)
abline(h = 0.025, lwd = 2, lty = 2)

legend(
  x = 3,
  y = 0.5,
  col = c("blue", "red", "black"),
  lty = c(1, 1, 2),
  cex = 0.8,
  bg = "white",
  legend = c(
    "power difference",
    TeX("$\\alpha_B(\\theta_c=\\theta_t; d_E)$"),
    TeX("$\\alpha=0.025$")
  )
)

# Add new plot
par(new = TRUE) 
plot(
  dta_histFix$Delta,
  dta_histFix$alphaBPPFix,
  ylim = c(0, 1),
  type = "n",
  # Create second plot without axes
  axes = FALSE,
  xlab = "",
  ylab = ""
)
axis(side = 4, las = 1)  # Add second axis
mtext("Power difference",
      side = 4,
      line = 2.5,
      cex = 1)

dev.off()
