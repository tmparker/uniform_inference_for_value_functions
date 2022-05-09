# This file runs a test using a normal location-shift experiment to check the
# power of uniform tests for the lower bound of a CDF example.
# Alternative sequences defined as in Hong and Li (2018), and simulated
# derivatives as in H&L (2018) as well, although not entirely appropriate for
# the problem.

# rlecuyer library used for generating random numbers in parallel.
library(rlecuyer)

Nvec <- c(100, 500, 1000)
Rvec <- c(499, 999, 1999)

# Use nproc processors, each making 1000/nproc runs for an N/R pair.
# Call one N/R pair a "design".  There are 3 here.
# Job array index values go from 0 to nproc-1, nproc to 2nproc-1, etc.
# So passing the entire array to SLURM uses array=0-(3*nproc-1).

nproc <- 500
reps <- 1000 / nproc
arg <- as.numeric(commandArgs(trailingOnly = TRUE))
dnum <- floor(arg / nproc) + 1 # design number

N <- Nvec[dnum]
R <- Rvec[dnum]
gfine <- 0.05
alpha <- 0.05

# H&L alternatives
alt <- c(-2, -(2 * N)^(-1/10), -(2 * N)^(-1/6), -(2 * N)^(-1/3), -(2 * N)^(-1/2),
          -(2 * N)^(-2/3), -(2 * N)^(-1), 0, (2 * N)^(-1), (2 * N)^(-2/3),
          (2 * N)^(-1/2), (2 * N)^(-1/3), (2 * N)^(-1/6), (2 * N)^(-1/10), 2)
hl_eps <- c((2 * N)^(-1/6), (2 * N)^(-1/3), (2 * N)^(-1/2), (2 * N)^(-1))

aname <- c("-2", "-N^(-1/10)", "-N^(-1/6)", "-N^(-1/3)", "-N^(-1/2)",
          "-N^(-2/3)", "-N^(-1)", "0", "N^(-1)", "N^(-2/3)",
          "N^(-1/2)", "N^(-1/3)", "N^(-1/6)", "N^(-1/10)", "2")
hlname <- c("N^(-1/6)", "N^(-1/3)", "N^(-1/2)", "N^(-1)")

# Compute the lower bound function of the CDF of a sum of RVs
L_sum <- function(ev, sam1, sam0) {
  ne <- length(ev)
  F1 <- ecdf(sam1)
  F1_eval <- F1(ev)
  L <- double(ne)
  for (i in 1:ne) {
    F0 <- ecdf(-sam0 + ev[i])
    dproL <- F1_eval - F0(ev)
    L[i] <- max(max(dproL), 0)
  }
  pmax(0, L)
}

# Estimate level set estimates used in the analytical derivative estimates in
# boot_cb_combined() below.
eps_sum_lb <- function(ev, sam1, sam0, an, bn) {
  ne <- length(ev)
  F1 <- ecdf(sam1)
  F1_eval <- F1(ev)
  epsmax <- matrix(0, ne, ne)
  Lv <- double(ne)
  for (i in 1:ne) {
    F0 <- ecdf(-sam0 + ev[i])
    dproL <- F1_eval - F0(ev)
    Lv[i] <- max(max(dproL), 0)
    epsmax[, i] <- (dproL >= Lv[i] - an)
  }
  is0 <- (abs(Lv) <= bn)
  isbig <- (Lv > bn)
  list(epsmax = epsmax, is0 = is0, isbig = isbig)
}

# Bound for two normal distributions, one standard, one with mean mu.
L <- function(x, mu) pmax(2 * pnorm((x - mu) / 2) - 1, 0)

# Bootstrap reference distributions using analytical approximation of the
# derivative and Hong and Li's (2018) estimator.
# These options use an L_0 since it doesn't really work without it.
boot_cb_combined <- function(ev, sam1, sam0, R, an, bn, epsn) {
  ne <- length(ev)
  F1 <- ecdf(sam1)
  F1_eval <- F1(ev)
  elist <- eps_sum_lb(ev, sam1, sam0, an, bn) # for estimated one
  Lhat <- L_sum(ev, sam1, sam0)
  Lev <- L(ev, 0) # null hypothesis
  Lest <- max(abs(Lhat - Lev))
  max_on_U <- double(ne)
  analytical <- double(R)
  L_eps <- matrix(0, nrow = ne, ncol = length(epsn))
  sim1 <- sim2 <- matrix(0, nrow = R, ncol = length(epsn))
  for (r in 1:R) {
    # Generate samples
    s0 <- sample(sam0, replace = TRUE)
    s1 <- sample(sam1, replace = TRUE)
    bF1 <- ecdf(s1)
    bF1_eval <- bF1(ev)
    for (i in 1:ne) {
      F0 <- ecdf(-sam0 + ev[i])
      bF0 <- ecdf(-s0 + ev[i])
      pro <- bF1_eval - bF0(ev) - F1_eval + F0(ev)
      max_on_U[i] <- max(pro * elist$epsmax[, i])
      pro_eps <- sapply(epsn, function(e) {F1_eval - F0(ev) +
                          e * sqrt(2 * N) * pro})
      L_eps[i, ] <- apply(pro_eps, 2, function(x) max(x, 0))
    }
    if (sum(elist$is0) > 0) { max0 <- max(max(max_on_U[elist$is0]), 0)
    } else { max0 <- 0 }
    if (sum(elist$isbig) > 0) { maxp <- max(abs(max_on_U[elist$isbig]))
    } else { maxp <- 0 }
    analytical[r] <- max(max0, maxp)
    sim1[r, ] <- apply(L_eps, 2, function(x) {max(abs(x - Lev)) - Lest}) / epsn
    sim2[r, ] <- apply(L_eps, 2, function(x) {max(abs(x - Lhat))}) / epsn
  }
  list(analytical = analytical, sim1 = sim1, sim2 = sim2)
}

# Generate data and conduct tests, return p-values.
cband_pval <- function(N, mu, R, gr_fine, HL_eps) {
  X <- rnorm(N)
  Y <- rnorm(N, mean = mu)
  qgrd <- seq(0, 1, by = 0.01)
  qX <- quantile(X, probs = qgrd, type = 1)
  qYr <- quantile(-Y, probs = qgrd, type = 1)
  gmin <- min(qX - qYr) - 2 * gr_fine
  gmax <- max(X) + max(Y) + 2 * gr_fine
  grd <- seq(gmin, gmax, by = gr_fine)
  bnd_lo <- L_sum(grd, Y, X)
  ks_lo <- max(abs(sqrt(2 * N) * (bnd_lo - L(grd, 0))))
  an <- 0.5 * sqrt(log(log(2 * N)) / (2 * N))
  bn <- 3.5 * log(log(2 * N)) / sqrt(2 * N)
  boots <- boot_cb_combined(grd, Y, X, R, an, bn, HL_eps)
  pval_an <- mean(sqrt(2 * N) * boots$analytical > ks_lo)
  pval_si <- apply(boots$sim1, 2, function(x) mean(x + 1e-06 > ks_lo))
  pval_si2 <- apply(boots$sim2, 2, function(x) mean(x + 1e-06 > ks_lo))
  c(pval_an, pval_si, pval_si2)
}

# Generate random numbers in parallel here
nstream <- 3 * nproc # parallel RNG streams is #designs * #processors
.lec.SetPackageSeed(6:1)
stream.names <- paste0("st", 0:(nstream - 1))
.lec.CreateStream(stream.names)
.lec.CurrentStream(paste0("st", arg))

inner <- function(M) {replicate(reps, cband_pval(N, M, R, gfine, hl_eps))}
pv <- lapply(alt, inner)
names(pv) <- paste("mu =", aname)
for (i in 1:length(pv)) {
  rownames(pv[[i]]) <- c("analytic", rep(paste("eps =", hlname), 2))
}

part_print <- sprintf("%03d", arg %% nproc + 1) # print the part number
save(pv, file = paste0("./parts/ss", N, "part", part_print, ".rda"))

options(warn = -1)
.lec.CurrentStreamEnd() # The warning seems to be about other RNGs.
