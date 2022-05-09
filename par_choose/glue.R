# Put things back together.

N <- 40
nparts <- 3
avec <- 1:3# seq(0.5, 5, by = 0.5)
bvec <- 1:4#seq(0.5, 5, by = 0.5)

load(paste0("./parts/ss", N, "part1.rda"))
pv_list <- pvals

for (j in 2:nparts) {
  load(paste0("./parts/ss", N, "part", j, ".rda"))
  pv_list <- lapply(1:length(avec), function(i) rbind(pv_list[[i]], pvals[[i]]))
}
names(pv_list) <- paste("an =", avec)
save(pv_list, file = paste0("par_ss", N, "_pvals.rda"))

