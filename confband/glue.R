# Put things back together.

N <- c(100, 500, 1000)
for (ss in 1:length(N)) {
  load(paste0("./parts/size", N[ss], "run1.rda"))
  pvmat <- pval
  for (run in 2:100) {
    load(paste0("./parts/size", N[ss], "run", run, ".rda"))
    pvmat <- rbind(pvmat, pval)
  }
  save(pvmat, file = paste0("cband_sim_size", N[ss], "_pvals.rda"))
}

