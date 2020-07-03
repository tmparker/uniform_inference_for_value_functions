# Put things back together.

N <- c(100, 500, 1000)
for (ss in 1:length(N)) {
  load(paste0("./parts/size", N[ss], "run1.rda"))
  pval.mat <- pval
  for (run in 2:100) {
    load(paste0("./parts/size", N[ss], "run", run, ".rda"))
    pval.mat <- rbind(pval.mat, pval)
  }
  save(pval.mat, file = paste0("sd_sim_size", N[ss], "_pvals.rda"))
}

