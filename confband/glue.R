# Put things back together.

N <- 1000
nparts <- 500

load(paste0("./parts/ss", N, "part001.rda"))
pv_list <- pv
nd <- length(pv) # number of alternative designs
altnames <- names(pv)

for (j in 2:nparts) {
  jform <- sprintf("%03d", j)
  load(paste0("./parts/ss", N, "part", jform, ".rda"))
  pv_list <- lapply(1:nd, function(i) cbind(pv_list[[i]], pv[[i]]))
}
names(pv_list) <- altnames
save(pv_list, file = paste0("cb_ss", N, "_pvals.rda"))

