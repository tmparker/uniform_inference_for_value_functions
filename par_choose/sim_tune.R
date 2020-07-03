# This file is for choosing tuning parameters in the two-sided testing case.
# Very similar to the confidence band experiment, just with varying tuning
# parameter for epsilon-maximization.

library(Rcpp)
sourceCpp("./bound_localpower.cpp", cacheDir="./lib_cpp")
library(rlecuyer)
nstream <- 100 # number of parallel random number streams

.lec.SetPackageSeed(6:1)
stream.names <- paste0("st", 1:nstream)
.lec.CreateStream(stream.names)
# To see all of them, they are stored in a global variable called 
# .lec.Random.seed.table.

N <- c(100, 500, 1000)
R <- c(499, 999, 1999)
gfine <- 0.05
reps <- 10
alpha <- 0.05
mu.big <- -5:5
Kvec <- 1:4 / 10

# Construct population lower-bound function for 2 normal distributions
L <- function(x, mu1, mu0) {ifelse(x >= mu1 - mu0, 
                         2 * pnorm((x - mu1 + mu0) / 2) - 1, 0)}

cband.pval <- function(N, mu, R, gr.fine, K) {
  co <- rnorm(N)
  tr <- rnorm(N, mean = mu)
  ends <- c(min(tr) - max(co), max(tr) - min(co))
  grd <- seq(ends[1] - 0.1, ends[2] + 0.1, by = gr.fine)
  bnds <- lohi(grd, tr, co)
  lpro <- sqrt(2 * N) * (bnds[1, ] - L(grd, 0, 0))
  stats.lo <- c(max(abs(lpro)), l2norm(grd, lpro))
  kap <- K * log(log(2 * N)) / sqrt(2 * N)
  b.cband <- lapply(boot_cb(grd, tr, co, R, kap), function(x) x * sqrt(2 * N))
  pval.ks <- mean(b.cband$KS > stats.lo[1])
  pval.cvm <- mean(b.cband$CvM > stats.lo[2])
  pv.vec <- c(pval.ks, pval.cvm)
  return(pv.vec)
}

runname <- as.numeric(commandArgs())
runname <- runname[length(runname)] # for some reason it is a vector.
.lec.CurrentStream(paste0("st", runname))

for (n in 1:length(N)) {
  parr.ks <- parr.cvm <- array(dim = c(length(mu.big), length(Kvec), reps))
  mu.local <- mu.big / sqrt(N[n])
  # dimensions: 1 = mu.big, 2 = Kvec, 3 = reps
  for (i in 1:length(mu.big)) {
    for (j in 1:length(Kvec)) {
      tmp <- replicate(reps, cband.pval(N[n], mu.local[i], R[n], gfine, Kvec[j]))
      parr.ks[i, j, ] <- tmp[1, ]
      parr.cvm[i, j, ] <- tmp[2, ]
    }
  }
  rej.ks <- apply(parr.ks, 1:2, function(x) sum(x < alpha))
  rej.cvm <- apply(parr.cvm, 1:2, function(x) sum(x < alpha))
  dimnames(rej.ks) <- dimnames(rej.cvm) <- list(mu.big, Kvec)
  anslist <- list(rej.ks = rej.ks, rej.cvm = rej.cvm)
  save(anslist, file = paste0("./parts/ssize", N[n], "part", runname, ".rda"))
}

.lec.CurrentStreamEnd()

