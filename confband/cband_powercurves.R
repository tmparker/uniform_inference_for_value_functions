# Plot the power curves from the simulation experiment.

alpha <- 0.05
mu1.plot <- -5:5

res.list <- c("cband_sim_size100_pvals.rda", "cband_sim_size500_pvals.rda", 
              "cband_sim_size1000_pvals.rda")
colvec <- c("navyblue", "firebrick", "forestgreen")

load(res.list[1])
pow <- apply(pvmat, 2, function(x) mean((x < alpha)))

pdf("cband_powercurves.pdf", width = 5, height = 4)
par(cex = 0.75)
plot(mu1.plot, pow, type = "l", ylim = c(0, 1), col = colvec[1], lty = 3,
      xlab = "local parameter", ylab = "empirical rejection probability",
      main = "Power curves, confidence band experiment")
      #main = expression(paste("Power curves: nominal ", alpha == 0.05, 
      # " rejection rate")))
abline(h = alpha, lty = 4, col = "gray")
abline(h = 0, lty = 3, col = "gray")

for (i in 2:3) {
  load(res.list[i])
  pow <- apply(pvmat, 2, function(x) mean((x < alpha)))
  lines(mu1.plot, pow, lty = 3 - i + 1, col = colvec[i])
}
legend("topleft", legend = paste("n =", c(100, 500, 1000)),
       lty = 3:1, col = colvec, bty = "n")

dev.off()
