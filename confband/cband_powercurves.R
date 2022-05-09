# Plot the power curves from the simulation experiment.

library(colorspace)

alpha <- 0.05
nvec <- c(100, 500, 1000)
Rvec <- c(499, 999, 1999)
aname <- c("-2", "-N^(-1/10)", "-N^(-1/6)", "-N^(-1/3)", "-N^(-1/2)",
          "-N^(-2/3)", "-N^(-1)", "0", "N^(-1)", "N^(-2/3)",
          "N^(-1/2)", "N^(-1/3)", "N^(-1/6)", "N^(-1/10)", "2")
hlname <- c("N^(-1/6)", "N^(-1/3)", "N^(-1/2)", "N^(-1)")
res_list <- c("cb_ss100_pvals.rda", "cb_ss500_pvals.rda", "cb_ss1000_pvals.rda")
al <- length(aname)
hll <- length(hlname)

# pal1 <- sequential_hcl(hll + 1, palette = "purp")
# pal2 <- sequential_hcl(hll + 1, palette = "mint")

# pdf("cband_powercurves.pdf", width = 12, height = 4)
# par(mfrow = c(1, length(res_list)), mar = c(3.1, 3.1, 1.1, 1.1),
#     mgp = 2:0, cex = 0.75)
# for (i in 1:length(res_list)) {
#   load(res_list[i])
#   pow <- lapply(pv_list, function(M) { rowMeans((M < alpha)) })
#   pmat <- do.call(rbind, pow)
#   matplot(1:length(aname), pmat, type = "l", ylim = c(0, 1),
#           lty = c(1, rep(1:hll + 1, 2)),
#           col = c("red", pal1[1:hll], pal2[1:hll]),
#           lwd = c(2, rep(1, ncol(pmat - 1))), xaxt = "n",
#           xlab = "Local alternative mean",
#           ylab = "empirical rejection probability",
#           main = paste("Sample size:", nvec[i]))
#   axis(side = 1, at = 1:length(aname), labels = aname)
#   abline(v = (length(aname) + 1) / 2, lty = 3, col = "gray50")
#   abline(h = alpha, lty = 4, col = "gray50")
#   abline(h = 0, lty = 3, col = "gray50")
#   legend(1, 0.25, legend = c("analytical", rep(paste("eps. = ", hlname), 2)),
#           lty = c(1, rep(1:hll + 1, 2)), lwd = c(2, rep(1, ncol(pmat) - 1)),
#           col = c("red", pal1[1:hll], pal2[1:hll]),
#           bty = "n")#, cex = 0.75)
# }
# dev.off()

library(Hmisc)
hln <- c("n^{-1/6}", "n^{-1/3}", "n^{-1/2}", "n^{-1}")
hl_ch <- paste0("$\\epsilon_n = ", hln, "$")
en <- c("-1/10", "-1/6", "-1/3", "-1/2", "-2/3", "-1")
an <- c("$-2$", paste0("$-n^{", en, "}$"), "$0$", rev(paste0("$n^{", en, "}$")), "$2$")

for (i in 1:length(res_list)) {
  load(res_list[i])
  cap_text <- paste0("Empirical rejection frequencies for a test that $L_{X+Y} \\equiv L_0$ for known $L_0$.  Analytic derivative estimates are described in the text, while Numeric derivative estimates use the method suggested by \\citet{HongLi18}.  Each row label gives the location parameter of one of the distributions, while the null assumes both parameters are zero.  Samples of size ", nvec[i], ", ", Rvec[i], " bootstrap repetitions per test, $1000$ simulation repetitions.")
  pow <- lapply(pv_list, function(M) { rowMeans((M < alpha)) })
  pmat <- do.call(rbind, pow)
  latex(pmat, file = paste0("../../fig/cband_tab_ss", nvec[i], ".tex"),
    size = "footnotesize", cgroup = c("Analytic", "Numeric (null known)",
    "Numeric (null estimated)"), n.cgroup = c(1, 4, 4), colheads = c("", rep(hl_ch, 2)),
    rgroup = c("", "Null is true", ""), n.rgroup = c(7, 1, 7),
    rowlabel = "Location", rowname = an, caption = cap_text,
    caption.loc = "bottom", label = paste0("cband_tab", nvec[i]),
    landscape = TRUE)
}

