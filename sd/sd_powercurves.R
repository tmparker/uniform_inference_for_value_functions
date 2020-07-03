# Plot the power curves from the simulation experiment.

alpha <- 0.05
mu_opt <- 1
mu1.plot <- seq(-3, 3, by = 0.5)

res.list <- c("sd_sim_size100_pvals.rda", "sd_sim_size500_pvals.rda", 
              "sd_sim_size1000_pvals.rda")
colvec <- c("navyblue", "firebrick", "forestgreen")

load(res.list[3])
pow <- apply(pval.mat, 2, function(x) mean((x < alpha), na.rm = TRUE))

pdf("sd_powercurves.pdf", width = 5, height = 4)
par(cex = 0.75)
plot(mu1.plot, pow, type = "n", ylim = range(pow),
      xlab = "local parameter", ylab = "empirical rejection probability",
      main = "Power curves, stochastic dominance experiment")
abline(h = alpha, lty = 4, col = "gray")
abline(h = 0, lty = 3, col = "gray")

for (i in 1:3) {
  load(res.list[i])
  pow <- apply(pval.mat, 2, function(x) mean((x < alpha)))
  lines(mu1.plot, pow, lty = 3 - i + 1, col = colvec[i])
}
legend("topright", legend = paste("n =", c(100, 500, 1000)),
       lty = 3:1, col = colvec, bty = "n")
dev.off()
