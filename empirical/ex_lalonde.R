# Confidence bands for the experimental Lalonde data and fun things to do with
# them.  Produces four figures: 
# - lalonde_cdfs.pdf is the empirical CDFs from the control and treatment
# distributions.
# - cdf_bounds_lalonde.pdf shows uniform confidence bands around each bound
# function (lower and upper) using the Lalonde data.
# - cdf_confband_lalonde.pdf shows the upper bound of the upper confidence band
# and lower bound of the lower confidence band to make a uniform confidence band
# for the CDF of the treatment effect, along with a few other estiamtes obtained
# using simple identification restrictions.
# - bounds_compare.pdf compares a collection of pointwise confidence intervals
# with the uniform confidence bands.

library(Rcpp)
sourceCpp("../bound.cpp", cacheDir="../lib_cpp")
load("lalonde.RData")
set.seed(8675309)

# alpha for CI of coverage rate 1 - alpha
# grd.diam is diameter of the grid used to calculate bound functions
alpha <- 0.05
grd.diam <- 100

# Read in data and make a few estimates
co <- lalonde.exp[lalonde.exp$treat == 0,]$re78
tr <- lalonde.exp[lalonde.exp$treat == 1,]$re78
pool <- deval(co, tr)
Fco <- ecdf(co)
Ftr <- ecdf(tr)
nco <- length(co)
ntr <- length(tr)
n <- nco + ntr

# Plot of the two empirical distribution functions
cvec <- c("firebrick", "navyblue")
pdf("lalonde_cdfs.pdf", width = 6, height = 4)
par(mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
plot(Ftr, do.points = FALSE, col = cvec[2], main = "", xlab = "outcome (dollars)", 
      ylab = "CDF", lty = 4)
lines(Fco, do.points = FALSE, col = cvec[1])
legend("bottomright", legend = c("control", "treatment", ""), col = cvec,
       lty = c(1, 4, 0), bty = "n")
abline(v = 0, lty = 3, col = "grey")
dev.off()

# For bootstrap/subsampling
kap <- 0.2 * log(log(n)) / sqrt(n)
B <- 1999
q <- 0.95
nlen <- 15
nbco <- floor(q^(1:nlen) * nco)
nbtr <- floor(q^(1:nlen) * ntr)

# Calculate on a grid
ends <- c(min(tr) - max(co), max(tr) - min(co))
grd <- seq(ends[1] - 1000, ends[2] + 1000, by = grd.diam)
G <- length(grd)
bnds <- lohi(grd, tr, co)

# Plot of the two empirical distribution functions
#par(mfrow = c(2, 1), mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
#plot(Fco, do.points = FALSE)
#lines(Ftr, col = "red", lty = 2, do.points = FALSE)

# Pointwise confidence intervals for each bound function.
# Using Fan & Park (2010) notation:
M <- m <- double(G)
argmax <- argmin <- double(G)
for (i in 1:G) {
  fd <- Ftr(pool) - Fco(pool - grd[i])
  M[i] <- max(fd)
  m[i] <- min(fd)
  argmax[i] <- pool[which.max(fd)]
  argmin[i] <- pool[which.min(fd)]
}

# Estimated asymptotic variance functions
sigL <- (n / ntr) * Ftr(argmax) * (1 - Ftr(argmax)) + 
        (n / nco) * Fco(argmax - grd) * (1 - Fco(argmax - grd))
sigU <- (n / ntr) * Ftr(argmin) * (1 - Ftr(argmin)) + 
        (n / nco) * Fco(argmin - grd) * (1 - Fco(argmin - grd))

# Pointwise confidence intervals using Fan & Park Thm. 3.2
Llo <- Lhi <- double(G)
Llo[M > 0] <- bnds[1, M > 0] - qnorm(1 - alpha) * sqrt(sigL[M > 0] / n) 
Lhi[M > 0] <- bnds[1, M > 0] + qnorm(1 - alpha) * sqrt(sigL[M > 0] / n)
Llo[M == 0] <- bnds[1, M == 0]
Lhi[M == 0] <- bnds[1, M == 0] + qnorm(((1 - alpha) + 1) / 2) * 
                sqrt(sigL[M == 0] / n)
Llo <- pmax(0, Llo)
Lhi <- pmin(1, Lhi)

Ulo <- Uhi <- rep(1, G)
Ulo[m < 0] <- bnds[2, m < 0] - qnorm(1 - alpha) * sqrt(sigU[m < 0] / n)
Uhi[m < 0] <- bnds[2, m < 0] + qnorm(1 - alpha) * sqrt(sigU[m < 0] / n)
Ulo[m == 0] <- bnds[2, m == 0] + qnorm(alpha / 2) * sqrt(sigU[m == 0] / n)
Uhi[m == 0] <- bnds[2, m == 0] 
Ulo <- pmax(0, Ulo)
Uhi <- pmin(1, Uhi)

#plot(grd, bnds[1,], type = "s", main = "Pointwise asymptotic confidence intervals",
#      xlab = "effect", ylab = "bounds")# xlim = c(-40000, 40000))
#lines(grd, bnds[2,])
#
#lines(grd, Llo, col = "red", lty = 2)
#lines(grd, Lhi, col = "red", lty = 2)
#lines(grd, Ulo, col = "blue", lty = 2)
#lines(grd, Uhi, col = "blue", lty = 2)

#####
# Pointwise bootstrap confidence intervals as described by Fan & Park (2010)
loboot <- array(0, dim = c(B, G, nlen))
hiboot <- array(0, dim = c(B, G, nlen))
for (k in 1:nlen) {
  lo <- matrix(0, B, G)
  hi <- matrix(0, B, G)
  for (b in 1:B) {
    nb <- nbco[k] + nbtr[k]
    costar <- sample(co, nbco[k], replace = TRUE)
    trstar <- sample(tr, nbtr[k], replace = TRUE)
    Fcost <- ecdf(costar)
    Ftrst <- ecdf(trstar)
    bndstar <- lohi(grd, trstar, costar)
    maxst <- minst <- double(G)
    for (i in 1:G) {
      fd <- Ftrst(pool) - Fcost(pool - grd[i])
      maxst[i] <- pool[which.max(fd)]
      minst[i] <- pool[which.min(fd)]
    }
    sigLst <- (nb / nbtr[k]) * Ftrst(maxst) * (1 - Ftrst(maxst)) + 
              (nb / nbco[k]) * Fcost(maxst - grd) * (1 - Fcost(maxst - grd)) 
    sigUst <- (nb / nbtr[k]) * Ftrst(minst) * (1 - Ftrst(minst)) +
              (nb / nbco[k]) * Fcost(minst - grd) * (1 - Fcost(minst - grd))
    TstarL <- sqrt(nb) * (bndstar[1, ] - bnds[1, ]) / sqrt(sigLst)
    TstarU <- sqrt(nb) * (bndstar[2, ] - bnds[2, ]) / sqrt(sigUst)
    lo[b, ] <- TstarL
    hi[b, ] <- TstarU
  }
  loboot[, , k] <- apply(lo, 2, sort, na.last = TRUE)
  hiboot[, , k] <- apply(hi, 2, sort, na.last = TRUE)
}

ldiff <- hdiff <- matrix(0, nlen - 1, G)
for (k in 1:(nlen - 1)) {
  for (j in 1:G) {
    if (any(!is.na(loboot[, j, k])) && any(!is.na(loboot[, j, k + 1]))) {
      ldiff[k, j] <- max(abs(ecdf(loboot[, j, k])(grd) - ecdf(loboot[, j, k + 1])(grd)))
    } else {ldiff[k, j] <- 2}
    if (any(!is.na(hiboot[, j, k])) && any(!is.na(hiboot[, j, k + 1]))) {
      hdiff[k, j] <- max(abs(ecdf(hiboot[, j, k])(grd) - ecdf(hiboot[, j, k + 1])(grd)))
    } else {hdiff[k, j] <- 2}
  }
}

lstar.lo <- apply(ldiff, 2, which.min)
lstar.hi <- apply(hdiff, 2, which.min)
lo.c1 <- lo.c2 <- double(G)
hi.c1 <- hi.c2 <- double(G)
for (j in 1:G) {
    lo.c1[j] <- quantile(loboot[, j, lstar.lo[j]], 1 - alpha, na.rm = TRUE)
    lo.c2[j] <- quantile(loboot[, j, lstar.lo[j]], alpha, na.rm = TRUE)
    hi.c1[j] <- quantile(hiboot[, j, lstar.hi[j]], 1 - alpha, na.rm = TRUE)
    hi.c2[j] <- quantile(hiboot[, j, lstar.hi[j]], alpha, na.rm = TRUE)
}

Llo.b <- bnds[1, ] - lo.c1 * sqrt(sigL / n)
Lhi.b <- bnds[1, ] - lo.c2 * sqrt(sigL / n)
Ulo.b <- bnds[2, ] - hi.c1 * sqrt(sigU / n)
Uhi.b <- bnds[2, ] - hi.c2 * sqrt(sigU / n)

Llo.b[is.finite(Llo.b)] <- pmax(Llo.b[is.finite(Llo.b)], 0)
Lhi.b[is.finite(Lhi.b)] <- pmin(Lhi.b[is.finite(Lhi.b)], 1)
Ulo.b[is.finite(Ulo.b)] <- pmax(Ulo.b[is.finite(Ulo.b)], 0)
Uhi.b[is.finite(Uhi.b)] <- pmin(Uhi.b[is.finite(Uhi.b)], 1)

#####
# Uniform bounds using KS statistics
uboot <- boot_cb(grd, tr, co, B, kap)
uq_lo <- quantile(uboot[[1]], 1 - alpha)
uq_hi <- quantile(uboot[[2]], 1 - alpha)
Llo.unif <- pmax(bnds[1, ] - uq_lo, 0)
Lhi.unif <- pmin(bnds[1, ] + uq_lo, 1)
Ulo.unif <- pmax(bnds[2, ] - uq_hi, 0)
Uhi.unif <- pmin(bnds[2, ] + uq_hi, 1)

# Uniform confidence band for the true CDF needs alpha/2 limits.
unifL <- pmax(bnds[1, ] - quantile(uboot[[1]], 1 - alpha / 2), 0)
unifU <- pmin(bnds[2, ] + quantile(uboot[[2]], 1 - alpha / 2), 1)

# Under rank invariance, QTEs give you the CDF of treatment effects.
tgrd <- seq(0.01, 0.99, by = 0.005)
Qco <- function(tau) {quantile(co, tau, type = 1)}
Qtr <- function(tau) {quantile(tr, tau, type = 1)}
QTE.ri <- sort(Qtr(tgrd) - Qco(tgrd))
FTE.ri <- stepfun(QTE.ri, c(tgrd, 1))

save.image(file = "lalonde_cdf_data.rda")

load("lalonde_cdf_data.rda")

#####
# Plots
cvec <- c("red", "blue")
pdf("cdf_bounds_lalonde.pdf", width = 6, height = 4)
par(mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
plot(grd, bnds[1,], type = "s", main = "",
      xlab = "effect", ylab = "CDF/bounds", xlim = c(-45000, 45000), 
      col = cvec[1])
lines(grd, bnds[2,], type = "s", col = cvec[2])
lines(grd, Llo.unif, type = "s", lty = 3, col = cvec[1])
lines(grd, Lhi.unif, type = "s", lty = 3, col = cvec[1])
lines(grd, Ulo.unif, type = "s", lty = 4, col = cvec[2])
lines(grd, Uhi.unif, type = "s", lty = 4, col = cvec[2])
#lines(grd, FTE.ri(grd), type = "s", col = "purple")
legend("topleft", legend = c("lower", "upper"), col = cvec, 
       lty = c(3, 4), bty = "n")
dev.off()

cvec <- c("black", "purple", "forestgreen", "firebrick")
pdf("cdf_confband_lalonde.pdf", width = 6, height = 4)
par(mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
plot(grd, unifL, type = "s", main = "",
      xlab = "effect", ylab = "CDF", xlim = c(-25000, 25000), 
      ylim = 0:1, col = cvec[1], lty = 1)
lines(grd, unifU, type = "s", lty = 1, col = cvec[1])
lines(grd, FTE.ri(grd), type = "s", lty = 2, col = cvec[2])
abline(v = 0, lty = 3, col = cvec[3])
abline(v = mean(tr) - mean(co), lty = 4, col = cvec[4])
legend("topleft", legend = c("confidence band", "comonotonic est.", 
        "degenerate at zero", "est. ATE"), col = cvec, lty = 1:4, bty = "n")
dev.off()

cvec <- c("red", "blue", "firebrick", "navyblue")
pdf("bounds_compare.pdf", width = 6, height = 4)
par(mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
plot(grd, bnds[1,], type = "s", main = "", 
      xlab = "effect", ylab = "bounds", xlim = c(-45000, 45000))
lines(grd, bnds[2,], type = "s")
lines(grd, Llo.b, type = "s", col = cvec[3], lty = 3)
lines(grd, Lhi.b, type = "s", col = cvec[3], lty = 3)
lines(grd, Ulo.b, type = "s", col = cvec[4], lty = 3)
lines(grd, Uhi.b, type = "s", col = cvec[4], lty = 3)
lines(grd, Llo.unif, type = "s", lty = 2, col = cvec[1])
lines(grd, Lhi.unif, type = "s", lty = 2, col = cvec[1])
lines(grd, Ulo.unif, type = "s", lty = 2, col = cvec[2])
lines(grd, Uhi.unif, type = "s", lty = 2, col = cvec[2])
legend("topleft", legend = c("lower conf. bands", "upper conf. bands", 
        "lower conf. intervals", "upper conf. intervals"), col = cvec, 
       lty = c(2, 2, 3, 3), bty = "n")
dev.off()

