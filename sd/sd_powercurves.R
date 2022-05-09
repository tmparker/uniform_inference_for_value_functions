# Plot the power curves from the simulation experiment.

alpha <- 0.05
nvec <- c(100, 500, 1000)
aname <- c("-2", "-N^(-1/10)", "-N^(-1/6)", "-N^(-1/3)", "-N^(-1/2)",
          "-N^(-2/3)", "-N^(-1)", "0", "N^(-1)", "N^(-2/3)",
          "N^(-1/2)", "N^(-1/3)", "N^(-1/6)", "N^(-1/10)", "2")
hlname <- c("N^(-1/6)", "N^(-1/3)", "N^(-1/2)", "N^(-1)")
res_list <- c("sd_ss100_pvals.rda", "sd_ss500_pvals.rda", "sd_ss1000_pvals.rda")

# pdf("sd_powercurves.pdf", width = 12, height = 4)
# par(mfrow = c(1, length(res_list)), mar = c(3.1, 3.1, 1.1, 1.1),
#     mgp = 2:0, cex = 0.75)
# for (i in 1:length(res_list)) {
#   load(res_list[i])
#   pow <- lapply(pv_list, function(M) { rowMeans((M < alpha)) })
#   pmat <- do.call(rbind, pow)
#   matplot(1:length(aname), pmat, type = "l", ylim = c(0, 1),
#           col = c("red", rep("black", ncol(pmat) - 1)),
#           lwd = c(1.5, rep(1, ncol(pmat - 1))), xaxt = "n",
#           xlab = "Local alternative mean",
#           ylab = "empirical rejection probability",
#           main = paste("Sample size:", nvec[i]))
#   axis(side = 1, at = 1:length(aname), labels = aname)
#   abline(v = (length(aname) + 1) / 2, lty = 3, col = "gray50")
#   abline(h = alpha, lty = 4, col = "gray50")
#   abline(h = 0, lty = 3, col = "gray50")
#   legend("topright", legend = c("analytical", paste("eps. = ", hlname)),
#           lty = 1:5, lwd = c(1.5, rep(1, ncol(pmat) - 1)),
#           col = c("red", rep("black", ncol(pmat) - 1)), bty = "n", cex = 0.75)
# }
# dev.off()

library(Hmisc)
hln <- c("n^{-1/6}", "n^{-1/3}", "n^{-1/2}", "n^{-1}")
hl_ch <- paste0("$\\epsilon_n = ", hln, "$")
en <- c("-1/10", "-1/6", "-1/3", "-1/2", "-2/3", "-1")
an <- c("$-2$", paste0("$-n^{", en, "}$"), "$0$", rev(paste0("$n^{", en, "}$")), "$2$")

for (i in 1:length(res_list)) {
  load(res_list[i])
  cap_text <- paste0("Empirical rejection frequencies for a test that $L_A \\leq U_B$.  Analytic derivative estimates are described in the text, while Numeric derivative estimates use the method suggested by \\citet{HongLi18}.  Each row label gives the location parameter of one of the distributions, while the null assumes both parameters are zero.  Samples of size ", nvec[i], ", ", Rvec[i], " bootstrap repetitions per test, $1000$ simulation repetitions.")
  pow <- lapply(pv_list, function(M) { rowMeans((M < alpha)) })
  pmat <- do.call(rbind, pow)
  # latex(pmat, file = "", cgroup = c("Analytic", "Numeric"), n.cgroup = c(1, 4), colheads = ch)
  latex(pmat, file = paste0("../../fig/sd_tab_ss", nvec[i], ".tex"),
    cgroup = c("Analytic", "Numeric"), n.cgroup = c(1, 4), colheads = c("", hl_ch),
    rgroup = c("", "Null bndry.", ""), n.rgroup = c(7, 1, 7),
    rowlabel = "Location", rowname = an, caption = cap_text, 
    caption.loc = "bottom", label = paste0("sd_tab", nvec[i]))
}
