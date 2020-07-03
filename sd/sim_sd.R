# A simulation experiment to show the breakdown frontier test example.

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
reps <- 10

alpha <- 0.05
mu_bndry <- 1
mu_grid <- seq(-2, 1, by = 0.25) 
gfine <- 0.01

# This function generates one bootstrap p-value
sd.pval <- function(N, mu, R, gr.fine) {
  co <- runif(N)
  A <- runif(N, min = mu, max = mu + 1)
  B <- runif(N)
  grd <- seq(-1, mu + 1, by = gr.fine)
  bndA <- lohi(grd, A, co)
  bndB <- lohi(grd, B, co)
  diffproc <- bndA[2, ] - bndB[1, ]
  sd.stat <- sqrt(2 * N) * l2norm(grd, pmax(diffproc, 0))
  kap <- 0.2 * log(log(2 * N)) / sqrt(2 * N)
  con <- (abs(diffproc) <= 3 * log(log(2 * N)) / sqrt(2 * N))
  b.one <- sqrt(2 * N) * boot_sd(grd, co, A, B, R, kap)
  bstats <- apply(pmax(b.one, 0) * con, 2, l2norm, eval = grd)
  pv <- mean(bstats > sd.stat) 
  return(pv)
}

runname <- as.numeric(commandArgs())
runname <- runname[length(runname)] # for some reason it is a vector.
.lec.CurrentStream(paste0("st", runname))

pval <- matrix(0, reps, length(mu_grid))
for (ss in 1:length(N)) {
  mu_local <- mu_bndry + mu_grid / sqrt(2 * N[ss])
  for (i in 1:length(mu_grid)) {
    pval[, i] <- replicate(reps, sd.pval(N[ss], mu_local[i], R[ss], gfine))
  }
  save(pval, file = paste0("./parts/size", N[ss], "run", runname, ".rda"))
}

.lec.CurrentStreamEnd()
