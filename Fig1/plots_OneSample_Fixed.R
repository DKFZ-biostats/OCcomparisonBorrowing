####Code for Figure 1 and descriptive statistics

rm(list = ls())
library(latex2exp)

resultsPPFix = read.csv(file = "OneSampleFix_delta0.5_alpha0.025_nE20.csv")
resultsEB = read.csv(file = "OneSampleEB_delta0.5_alpha0.025_nE20.csv")

ylimmax = max(abs(resultsPPFix$powerwFixdiff),
              abs(resultsEB$powerwEBdiff))
xlimmax = max(resultsPPFix$alphaB, resultsEB$alphaB)

resultsPPFix0 = subset(resultsPPFix, resultsPPFix$thetaE == 0)
resultsPPFix0.5 = subset(resultsPPFix, resultsPPFix$thetaE == 0.5)

# Figure 1
  jpeg(
    "alphaB_powerdiff_1sample.jpeg",
    width = 20,
    height = 20,
    units = "cm",
    res = 200,
    pointsize = 12
  )
  
  par(mfrow = c(2, 1),
      mar = c(4, 7, 2, 2))
  plot(
    resultsPPFix$alphaB,
    resultsPPFix$powerwFixdiff,
    ylim = c(-ylimmax, ylimmax),
    xlim = c(0, xlimmax),
    type = "n",
    las = 1,
    xlab = "T1E rate",
    ylab = "Power difference"
  )
  
  mtext(side = 3, line = -1, text = "Fixed Power Prior")
  mtext(
    side = 2,
    at = 0.0072,
    las = 1,
    line = 1,
    text = "a",
    cex = 1.5
  )
  abline(h = 0)
  abline(v = 0.025, lty = 2)
  abline(v = 0.0)
  points(resultsPPFix0$alphaB, 
         resultsPPFix0$powerwFixdiff, 
         col = "black")
  points(
    resultsPPFix0.5$alphaB,
    resultsPPFix0.5$powerwFixdiff,
    pch = 4,
    col = "red"
  )
  legend(
    "topright",
    pch = c(1, 4),
    col = c("black", "red"),
    legend = c(TeX("$\\theta_E=0"), TeX("$\\theta_E=0.5"))
  )
  axis(side = 1, at = 0.025, label = "")
  axis(
    side = 1,
    at = 0.025,
    line = 1.5,
    labels = TeX("\\alpha=0.025")
  )
  
  resultsEB0 = subset(resultsEB, resultsEB$thetaE == 0)
  resultsEB0.5 = subset(resultsEB, resultsEB$thetaE == 0.5)
  
  plot(
    resultsEB$alphaB,
    resultsEB$powerwEBdiff,
    ylim = c(-ylimmax, ylimmax),
    xlim = c(0, xlimmax),
    type = "n",
    las = 1,
    xlab = "T1E rate",
    ylab = "Power difference"
  )
  mtext(side = 3, line = -1, text = "Empirical Bayes Power Prior")
  mtext(
    side = 2,
    at = 0.0072,
    las = 1,
    line = 1,
    text = "b",
    cex = 1.5
  )
  abline(h = 0)
  abline(v = 0.025, lty = 2)
  abline(v = 0.0)
  points(resultsEB0$alphaB, resultsEB0$powerwEBdiff, col = "black")
  points(resultsEB0.5$alphaB,
         resultsEB0.5$powerwEBdiff,
         pch = 4,
         col = "red")
  legend(
    "topright",
    pch = c(1, 4),
    col = c("black", "red"),
    legend = c(TeX("$\\theta_E=0"), TeX("$\\theta_E=0.5"))
  )
  axis(side = 1, at = 0.025, label = "")
  axis(
    side = 1,
    at = 0.025,
    line = 1.5,
    labels = TeX("\\alpha=0.025")
  )
  
  dev.off()

# descriptive statistics for alphaB and power difference in both scenarios and both borrowing methods
  min(resultsPPFix0$alphaB)
  max(resultsPPFix0$alphaB)
  mean(resultsPPFix0$alphaB)
  median(resultsPPFix0$alphaB)
  
  min(resultsPPFix0$powerwFixdiff)
  max(resultsPPFix0$powerwFixdiff)
  mean(resultsPPFix0$powerwFixdiff)
  median(resultsPPFix0$powerwFixdiff)
  
  min(resultsPPFix0.5$alphaB)
  max(resultsPPFix0.5$alphaB)
  mean(resultsPPFix0.5$alphaB)
  median(resultsPPFix0.5$alphaB)
  
  min(resultsPPFix0.5$powerwFixdiff)
  max(resultsPPFix0.5$powerwFixdiff)
  mean(resultsPPFix0.5$powerwFixdiff)
  median(resultsPPFix0.5$powerwFixdiff)
  
  
  
  min(resultsEB0$alphaB)
  max(resultsEB0$alphaB)
  mean(resultsEB0$alphaB)
  median(resultsEB0$alphaB)
  
  min(resultsEB0$powerwEBdiff)
  max(resultsEB0$powerwEBdiff)
  mean(resultsEB0$powerwEBdiff)
  median(resultsEB0$powerwEBdiff)
  
  min(resultsEB0.5$alphaB)
  max(resultsEB0.5$alphaB)
  mean(resultsEB0.5$alphaB)
  median(resultsEB0.5$alphaB)
  
  min(resultsEB0.5$powerwEBdiff)
  max(resultsEB0.5$powerwEBdiff)
  mean(resultsEB0.5$powerwEBdiff)
  median(resultsEB0.5$powerwEBdiff)
