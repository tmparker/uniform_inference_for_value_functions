# This file runs a test using a normal location-shift experiment to verify that
# uniform confidence bands have the correct coverage probability.
# This focuses on the lower bound function.

library(Rcpp)
sourceCpp("../bound.cpp", cacheDir="../lib_cpp")
#set.seed(123)

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

L <- function(x, mu1, mu0) {ifelse(x >= mu1 - mu0, 
                         2 * pnorm((x - mu1 + mu0) / 2) - 1, 0)}

# This function generates one bootstrap p-value
cband.pval <- function(N, mu, R, gr.fine) {
  co <- rnorm(N)
  tr <- rnorm(N, mean = mu)
  ends <- c(min(tr) - max(co), max(tr) - min(co))
  grd <- seq(ends[1] - 0.1, ends[2] + 0.1, by = gr.fine)
  bnds <- lohi(grd, tr, co)
  lpro <- sqrt(2 * N) * (bnds[1, ] - L(grd, 0, 0))
  ks.lo <- max(abs(lpro))
  kap <- 0.2 * log(log(2 * N)) / sqrt(2 * N)
  b.cband <- sqrt(2 * N) * boot_cb(grd, tr, co, R, kap)$simL
  pval.lo <- mean(b.cband > ks.lo)
  return(pval.lo)
}

runname <- as.numeric(commandArgs())
runname <- runname[length(runname)] # for some reason it is a vector.
.lec.CurrentStream(paste0("st", runname))

pval <- matrix(0, reps, length(mu.big))
for (ss in 1:length(N)) {
  mu.local <- mu.big / sqrt(N[ss])
  for (i in 1:length(mu.big)) {
    pval[, i] <- replicate(reps, cband.pval(N[ss], mu.local[i], R[ss], gfine))
  }
  save(pval, file = paste0("./parts/size", N[ss], "run", runname, ".rda"))
}

.lec.CurrentStreamEnd()

