# This file is used for tuning parameter selection.
# Tests using sums of two normal random variables, focusing on the lower bound.
# This file focuses on setting the epsilon-maximizer- and contact set estimation
# parameters right under the null hypothesis (in this implementation, the null
# is that both distributions are standard normal).

library(rlecuyer)
nstream <- 100 # number of parallel random number streams
.lec.SetPackageSeed(6:1)
stream.names <- paste0("st", 1:nstream)
.lec.CreateStream(stream.names)
# To see all of them, they are stored in a global variable called
# .lec.Random.seed.table.

# Fix sample size, loop over tuning parameters
N <- 40 # sample size
R <- 9 # bootstrap repetitions
avec <- seq(0.5, 5, by = 0.5)
bvec <- seq(0.5, 5, by = 0.5)
reps <- 5 # simulation repetitions
gfine <- 0.05 # grid size for pointwise-supremum-finding
alpha <- 0.05 # nominal test size
mu <- 0 # this simulation is about size

# Just calculate a lower bound function from two samples
# ev is a vector of evaluation points
# sam1, sam0 are samples
L_sum <- function(ev, sam1, sam0) {
  ne <- length(ev)
  F1 <- ecdf(sam1)
  F1_eval <- F1(ev)
  L <- double(ne)
  for (i in 1:ne) {
    shift0 <- sam0 - ev[i]
    F0 <- ecdf(shift0)
    dproL <- F1_eval + F0(-ev) - 1
    L[i] <- max(dproL)
  }
  pmax(0, L)
}

# Like L_sum but also keeping track of epsilon-maximizers and contact sets
# an = parameter for epsilon-maximizer set estimation
# bn = parameter for contact set estimation
eps_sum_lb <- function(ev, sam1, sam0, an, bn) {
  ne <- length(ev)
  F1 <- ecdf(sam1)
  F1_eval <- F1(ev)
  epsmax <- matrix(0, ne, ne)
  Lv <- double(ne)
  for (i in 1:ne) {
    shift0 <- sam0 - ev[i]
    F0 <- ecdf(shift0)
    dproL <- F1_eval + F0(-ev) - 1
    Lv[i] <- max(max(dproL), 0)
    epsmax[, i] <- (dproL >= Lv[i] - an)
  }
  is0 <- (abs(Lv) <= bn)
  isbig <- (Lv > bn)
  list(epsmax = epsmax, is0 = is0, isbig = isbig)
}

# Bootstrap the lower bound function for confidence bands
# R is number of bootstrap repetitions, all other args the same
# Calls eps_sum_lb before going to bootstrap loop
boot_cb_lb <- function(ev, sam1, sam0, R, an, bn) {
  ne <- length(ev)
  F1 <- ecdf(sam1)
  F1_eval <- F1(ev)
  elist <- eps_sum_lb(ev, sam1, sam0, an, bn)
  max_on_U <- double(ne)
  blp <- double(R)
  max0 <- maxp <- double(R)
  for (r in 1:R) {
    s0 <- sample(sam0, replace = TRUE)
    s1 <- sample(sam1, replace = TRUE)
    bF1 <- ecdf(s1)
    bF1_eval <- bF1(ev)
    for (i in 1:ne) {
      shift0 <- sam0 - ev[i]
      bshift0 <- s0 - ev[i]
      F0 <- ecdf(shift0)
      bF0 <- ecdf(bshift0)
      dpro <- bF1_eval - F1_eval - bF0(-ev) + F0(-ev)
      max_on_U[i] <- max(dpro * elist$epsmax[, i])
    }
    max0[r] <- max(max(max_on_U[elist$is0]), 0)
    maxp[r] <- max(abs(max_on_U[elist$isbig]))
  }
  pmax(max0, maxp)
}

# Bound for two normal distributions, one standard, one with mean muB.
L <- function(x, muB) pmax(2 * pnorm((x - muB) / 2) - 1, 0)

# This function generates one bootstrap p-value
# gr_fine is how finely to set the evaluation grid
# aconst = constant in epsilon-maximum set estimation
# bconst = constant in contact set estimation
cband_pval <- function(N, mu, R, gr_fine, aconst, bconst) {
  co <- rnorm(N)
  tr <- rnorm(N, mean = mu)
  ends <- c(min(tr) - max(co), max(tr) - min(co))
  grd <- seq(ends[1] - 0.1, ends[2] + 0.1, by = gr_fine)
  bnd_lo <- L_sum(grd, tr, co)
  lpro <- sqrt(N) * (bnd_lo - L(grd, mu))
  ks_lo <- max(abs(lpro))
  an <- aconst * sqrt(log(log(N)) / N)
  bn <- bconst * log(log(N)) / sqrt(N)
  b_cband <- boot_cb_lb(grd, tr, co, R, an, bn)
  pval_lo <- mean(b_cband > ks_lo)
  return(pval_lo)
}

pname <- commandArgs(trailingOnly = TRUE)
.lec.CurrentStream(paste0("st", pname))

tfun <- function(A, B) {
  replicate(reps, cband_pval(N, 0, R, gfine, A, B))
}

pvals <- lapply(avec, function(A) sapply(bvec, tfun, A=A))
names(pvals) <- paste("an =", avec)
for (i in 1:length(avec)) colnames(pvals[[i]]) <- paste("bn =", bvec)
save(pvals, file = paste0("./parts/", "ss", N, "part", pname, ".rda"))

.lec.CurrentStreamEnd()

