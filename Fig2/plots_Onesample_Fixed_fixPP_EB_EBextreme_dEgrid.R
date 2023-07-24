rm(list = ls())
library(latex2exp)

resultsPPFix = read.csv(file = "OneSampleFix_fixPP_delta0.5_nE20_dEgrid_alpha0.025.csv")
resultsEB = read.csv(file = "OneSampleFix_EB_nE20_dEgrid_alpha0.025.csv")
resultsEBextreme = read.csv(file = "OneSampleFix_EB_nE1000_dEgrid_alpha0.025.csv")


setEPS()
postscript(
  "alphaB_powerdiff_1sample_dEgrid_3plots.eps",
  width = 8,
  height = 8,
  pointsize = 12
)
par(mfrow = c(3, 1),
    mar = c(4, 4, 2, 4))

plot(
  resultsPPFix$alphaB,
  resultsPPFix$powerwFixdiff,
  ylim = c(0, 1),
  xlim = c(-3, 4),
  type = "n",
  las = 1,
  xaxs = "i",
  xlab = TeX("$\\bar{d_E}$"),
  ylab = expression(Probability ~ to ~ reject ~ H[0])
)

mtext(side = 3, text = "Fixed Power Prior")
mtext(
  side = 2,
  at = 1,
  las = 1,
  line = 3,
  text = "a",
  cex = 1.5
)
abline(h = 0)
abline(h = 0.025, lty = 2)
abline(h = resultsPPFix$powerwo_orig[1], lty = 3)   #power w/o borrowing

lines(resultsPPFix$dE,
      resultsPPFix$alphaB,
      col = "red",
      lwd = 2)
lines(resultsPPFix$dE,
      resultsPPFix$powerwFix,
      col = "green",
      lwd = 2)
lines(resultsPPFix$dE,
      resultsPPFix$powerwFixdiff,
      col = "blue",
      lwd = 2)

legend(
  x = -2.9,
  y = 1,
  col = c("black", "green", "red", "blue", "black"),
  lty = c(3, 1, 1, 1, 2),
  cex = 0.8,
  bg = "white",
  legend = c(
    TeX("power w/o borrowing, $\\alpha=0.025$"),
    "power w/ borrowing",
    TeX("$\\alpha_B(d_E)$"),
    "power difference",
    TeX("$\\alpha=0.025$")
  )
)


par(new = TRUE) # Add new plot
plot(
  resultsPPFix$alphaB,
  resultsPPFix$powerwFixdiff,
  ylim = c(0, 1),
  xlim = c(-3, 4),
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
      cex = 0.7)


plot(
  resultsEB$alphaB,
  resultsEB$powerwEBdiff,
  ylim = c(0, 1),
  xlim = c(-3, 4),
  type = "n",
  las = 1,
  xaxs = "i",
  xlab = TeX("$\\bar{d_E}$"),
  ylab = expression(Probability ~ to ~ reject ~ H[0])
)

mtext(side = 3,
      text = TeX("Empirical Bayes Power Prior, $n_E=20$"))
mtext(
  side = 2,
  at = 1,
  las = 1,
  line = 3,
  text = "b",
  cex = 1.5
)
abline(h = 0)
abline(h = 0.025, lty = 2)
abline(h = resultsEB$powerwo_orig[1], lty = 3)   #power w/o borrowing

lines(resultsEB$dE,
      resultsEB$alphaB,
      col = "red",
      lwd = 2)
lines(resultsEB$dE,
      resultsEB$powerwEB,
      col = "green",
      lwd = 2)
lines(resultsEB$dE,
      resultsEB$powerwEBdiff,
      col = "blue",
      lwd = 2)

par(new = TRUE) # Add new plot
plot(
  resultsPPFix$alphaB,
  resultsPPFix$powerwFixdiff,
  ylim = c(0, 1),
  xlim = c(-3, 4),
  type = "n",
  # Create second plot without axes
  axes = FALSE,
  xlab = "",
  ylab = ""
)
axis(side = 4, las = 1)      # Add second axis
mtext("Power difference",
      side = 4,
      line = 2.5,
      cex = 0.7)


xlimmin = min(resultsEBextreme$dE)
xlimmax = max(resultsEBextreme$dE)
plot(
  resultsEBextreme$alphaB,
  resultsEBextreme$powerwEBdiff,
  ylim = c(-0.4, 1),
  xlim = c(xlimmin, xlimmax),
  type = "n",
  las = 1,
  xaxs = "i",
  xlab = TeX("$\\bar{d_E}$"),
  ylab = expression(Probability ~ to ~ reject ~ H[0])
)

mtext(side = 3,
      text = TeX("Empirical Bayes Power Prior, $n_E=1000$"))
mtext(
  side = 2,
  at = 1,
  las = 1,
  line = 3,
  text = "c",
  cex = 1.5
)
abline(h = 0)
abline(h = 0.025, lty = 2)
abline(h = resultsEBextreme$powerwo_orig[1], lty = 3)   #power w/o borrowing

lines(resultsEBextreme$dE,
      resultsEBextreme$alphaB,
      col = "red",
      lwd = 2)
lines(resultsEBextreme$dE,
      resultsEBextreme$powerwEB,
      col = "green",
      lwd = 2)
lines(
  resultsEBextreme$dE,
  resultsEBextreme$powerwo,
  col = "orange",
  lwd = 2,
  lty = 2
)
lines(
  resultsEBextreme$dE,
  resultsEBextreme$powerwEBdiff,
  col = "blue",
  lwd = 2
)
legend(
  x = 0.5,
  y = 0.8,
  col = c("black", "green", "orange", "red", "blue", "black"),
  lty = c(3, 1, 2, 1, 1, 2),
  cex = 0.8,
  bg = "white",
  legend = c(
    TeX("power w/o borrowing, $\\alpha=0.025$"),
    "power w/ borrowing",
    "power calibrated to borrowing",
    TeX("$\\alpha_B(d_E)$"),
    "power difference",
    TeX("$\\alpha=0.025$")
  )
)


# Add new plot
par(new = TRUE)                            
plot(
  resultsPPFix$alphaB,
  resultsPPFix$powerwFixdiff,
  ylim = c(0, 1),
  xlim = c(-3, 4),
  type = "n",
  # Create second plot without axes
  axes = FALSE,
  xlab = "",
  ylab = ""
)
axis(side = 4, las = 1)      # Add second axis
mtext("Power difference",
      side = 4,
      line = 2.5,
      cex = 0.7)

dev.off()
