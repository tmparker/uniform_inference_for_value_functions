# This file runs simulations for the stochastic dominance example, example B in
# the text.

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
gfine <- 0.02
alpha <- 0.05

# H&L alternatives - for SD experiment, -1 is the natural center
alt <- -1 + c(-2, -(2 * N)^(-1/10), -(2 * N)^(-1/6), -(2 * N)^(-1/3), -(2 * N)^(-1/2),
          -(2 * N)^(-2/3), -(2 * N)^(-1), 0, (2 * N)^(-1), (2 * N)^(-2/3), (2 * N)^(-1/2),
          (2 * N)^(-1/3), (2 * N)^(-1/6), (2 * N)^(-1/10), 2)
hl_eps <- c((2 * N)^(-1/6), (2 * N)^(-1/3), (2 * N)^(-1/2), (2 * N)^(-1))

aname <- c("-2", "-N^(-1/10)", "-N^(-1/6)", "-N^(-1/3)", "-N^(-1/2)",
          "-N^(-2/3)", "-N^(-1)", "0", "N^(-1)", "N^(-2/3)",
          "N^(-1/2)", "N^(-1/3)", "N^(-1/6)", "N^(-1/10)", "2")
hlname <- c("N^(-1/6)", "N^(-1/3)", "N^(-1/2)", "N^(-1)")

# Compute lower and upper bounds for CDF of difference X_1 - X_0
# This function isn't used anymore but the other functions are patterned after
# it.
lohi <- function(ev, sam1, sam0) {
  ne <- length(ev)
  F1 <- ecdf(sam1)
  F1_eval <- F1(ev)
  L <- U <- double(ne)
  for (i in 1:ne) {
    F0 <- ecdf(sam0 + ev[i])
    dpro <- F1_eval - F0(ev)
    L[i] <- max(max(dpro), 0)
    U[i] <- 1 + min(min(dpro), 0)
  }
  rbind(L, U)
}

# Calculate the L2 norm of the difference of two bounds
# gr_fine is used for trapezoid-rule integration
Lambda <- function(ev, sam0, samA, samB, gr_fine) {
  ne <- length(ev)
  FA <- ecdf(samA)
  FA_eval <- FA(ev)
  FB <- ecdf(samB)
  FB_eval <- FB(ev)
  LA <- UB <- double(ne)
  for (i in 1:ne) {
    F0 <- ecdf(sam0 + ev[i])
    F0_eval <- F0(ev)
    dproA <- FA_eval - F0_eval
    dproB <- FB_eval - F0_eval
    LA[i] <- max(max(dproA), 0)
    UB[i] <- 1 + min(min(dproB), 0)
  }
  d2pos <- pmax(LA - UB, 0)^2
  sqrt(0.5 * gr_fine * sum(d2pos[-ne] + d2pos[-1]))
}

# Estimate level set estimates used in the analytical derivative estimates in
# boot_sd_combined() below.
eps_diff <- function(ev, sam0, samA, samB, an, bn) {
  ne <- length(ev)
  FA <- ecdf(samA)
  FA_eval <- FA(ev)
  FB <- ecdf(samB)
  FB_eval <- FB(ev)
  emA <- matrix(0, ne, ne)
  emB <- matrix(0, ne, ne)
  LA <- UB <- double(ne)
  for (i in 1:ne) {
    F0 <- ecdf(sam0 + ev[i])
    F0_eval <- F0(ev)
    dproA <- FA_eval - F0_eval
    dproB <- FB_eval - F0_eval
    LA[i] <- max(dproA) # ev should be wide enough that 0 is included
    UB[i] <- min(dproB) # 1 isn't necessary for comparison below
    emA[, i] <- (dproA >= LA[i] - an)
    emB[, i] <- (dproB <= UB[i] + an)
  }
  X0_hat <- (abs(LA - UB) <= bn)
  if (sum(X0_hat) == 0) {X0_hat = rep(TRUE, ne)}
  list(emLA = emA, emUB = emB, X0_hat = X0_hat)
}

# Generate a bootstrap reference distribution
boot_sd_combined <- function(ev, sam0, samA, samB, R, an, bn, epsn, gr_fine) {
  ne <- length(ev)
  FA <- ecdf(samA)
  FA_eval <- FA(ev)
  FB <- ecdf(samB)
  FB_eval <- FB(ev)
  elist <- eps_diff(ev, sam0, samA, samB, an, bn)
  Lam <- Lambda(ev, sam0, samA, samB, gr_fine)
  LpA <- UpB <- double(ne)
  analytical <- double(R)
  LAeps <- UBeps <- matrix(0, ne, length(epsn))
  simulated <- matrix(0, R, length(epsn))
  for (r in 1:R) {
    s0 <- sample(sam0, replace = TRUE)
    sA <- sample(samA, replace = TRUE)
    sB <- sample(samB, replace = TRUE)
    bFA <- ecdf(sA)
    bFA_eval <- bFA(ev)
    bFB <- ecdf(sB)
    bFB_eval <- bFB(ev)
    for (i in 1:ne) {
      F0 <- ecdf(sam0 + ev[i])
      bF0 <- ecdf(s0 + ev[i])
      F0_eval <- F0(ev)
      bF0_eval <- bF0(ev)
      dproA <- bFA_eval - bF0_eval - FA_eval + F0_eval
      dproB <- bFB_eval - bF0_eval - FB_eval + F0_eval
      LpA[i] <- max(dproA * elist$emLA[, i])
      UpB[i] <- min(dproB * elist$emUB[, i])
      LAeps[i, ] <- sapply(epsn, function(e) {
        max(FA_eval - F0_eval + e * sqrt(2 * N) * dproA)})
      UBeps[i, ] <- 1 + sapply(epsn, function(e) {
        min(FB_eval - F0_eval + e * sqrt(2 * N) * dproB)})
    }
    dpr2 <- pmax((LpA - UpB) * elist$X0_hat, 0)^2
    analytical[r] <- sqrt(0.5 * gr_fine * sum(dpr2[-ne] + dpr2[-1]))
    dsim2 <- pmax(LAeps - UBeps, 0)^2
    L2sims <- sqrt(0.5 * gr_fine * colSums(dsim2[-ne, ] + dsim2[-1, ]))
    simulated[r, ] <- (L2sims - Lam) / epsn
  }
  list(analytical = analytical, simulated = simulated)
}

# Generate data and conduct tests, return p-values.
sd_pval <- function(N, mu, R, gr_fine, HL_eps) {
  co <- runif(N)
  A <- runif(N, min = mu, max = mu + 1)
  B <- runif(N)
  qgrd <- seq(0, 1, by = 0.01)
  q0 <- quantile(co, probs = qgrd, type = 1)
  qA <- quantile(A, probs = qgrd, type = 1)
  qB <- quantile(B, probs = qgrd, type = 1)
  gminA <- min(qA - q0) - 2 * gr_fine
  gmaxA <- max(A) - min(co) + 2 * gr_fine
  gminB <- min(B) - max(co) - 2 * gr_fine
  gmaxB <- max(qB - q0) + 2 * gr_fine
  grd <- seq(min(gminA, gminB), max(gmaxA, gmaxB), by = gr_fine)
  sd_stat <- sqrt(2 * N) * Lambda(grd, co, A, B, gr_fine)
  an <- 0.5 * sqrt(log(log(2 * N))) / sqrt(2 * N)
  bn <- 3.5 * log(log(2 * N)) / sqrt(2 * N)
  boots <- boot_sd_combined(grd, co, A, B, R, an, bn, HL_eps, gr_fine)
  pval_an <- mean(sqrt(2 * N) * boots$analytical + 1e-06 > sd_stat) # if all 0
  pval_si <- apply(boots$simulated, 2, function(x) mean(x + 1e-06 > sd_stat))
  c(pval_an, pval_si)
}

# Generate random numbers in parallel here
nstream <- 3 * nproc # parallel RNG streams is #designs * #processors
.lec.SetPackageSeed(6:1)
stream.names <- paste0("st", 0:(nstream - 1))
.lec.CreateStream(stream.names)
.lec.CurrentStream(paste0("st", arg))

inner <- function(M) {replicate(reps, sd_pval(N, M, R, gfine, hl_eps))}
pv <- lapply(alt, inner)
names(pv) <- paste("mu =", aname)
for (i in 1:length(pv)) {
  rownames(pv[[i]]) <- c("analytic", paste("eps =", hlname))
}

part_print <- sprintf("%03d", arg %% nproc + 1) # print the part number
save(pv, file = paste0("./parts/ss", N, "part", part_print, ".rda"))

options(warn = -1)
.lec.CurrentStreamEnd() # The warning seems to be about other RNGs.
