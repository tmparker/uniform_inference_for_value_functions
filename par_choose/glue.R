# Glue things back together

N <- c(100, 500, 1000)

for (n in 1:length(N)) {
  load(paste0("./parts/ssize", N[n], "part1.rda"))
  # loads a list of two matrices giving numbers of rejections for the tests run 
  # in that particular part.
  rej.list <- anslist
  for (k in 1:100) {
    load(paste0("./parts/ssize", N[n], "part", k, ".rda"))
    rej.list <- mapply("+", rej.list, anslist, SIMPLIFY = FALSE)
    pval.list <- lapply(rej.list, function(x) x / 1000)
    save(pval.list, file = paste0("pvalues_ssize", N[n], ".rda"))
  }
}
