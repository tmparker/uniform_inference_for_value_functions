# Get empirical rejection probabilities from the simulations

nreps <- 1000
avec <- seq(0.5, 5, by = 0.5)
bvec <- seq(0.5, 5, by = 0.5)

load("par_ss100_pvals.rda") # pv_list has matrices of p-values
pv100 <- pv_list
load("par_ss500_pvals.rda") # pv_list has matrices of p-values
pv500 <- pv_list
load("par_ss1000_pvals.rda") # pv_list has matrices of p-values
pv1000 <- pv_list

emp_check <- function(x, al) {colMeans(x < al) - al}
r100_05 <- sapply(pv100, emp_check, al = 0.05)
r100_10 <- sapply(pv100, emp_check, al = 0.10)
r500_05 <- sapply(pv500, emp_check, al = 0.05)
r500_10 <- sapply(pv500, emp_check, al = 0.10)
r1000_05 <- sapply(pv1000, emp_check, al = 0.05)
r1000_10 <- sapply(pv1000, emp_check, al = 0.10)

# # What has rejection probability closest to 5% or 10%?
# r05 <- sapply(pv_list, function(x) colMeans(x < 0.05) - 0.05)
# r10 <- sapply(pv_list, function(x) colMeans(x < 0.1) - 0.1)

# par(mfrow = c(1, 2), mar = c(0.1, 0.5, 0.1, 0.1), mgp = 2:0)
# pl <- persp(avec, bvec, r05, theta = 120, col = "lightblue", ticktype = "detailed", main = "5%")
# lines(trans3d(x = avec[1], y = bvec, z = 0, pl))
# lines(trans3d(x = avec[length(avec)], y = bvec, z = 0, pl))
# lines(trans3d(x = avec, y = bvec[1], z = 0, pl))
# lines(trans3d(x = avec, y = bvec[length(bvec)], z = 0, pl))

# pl <- persp(avec, bvec, r10, theta = 120, col = "lightblue", ticktype = "detailed", main = "10%")
# lines(trans3d(x = avec[1], y = bvec, z = 0, pl))
# lines(trans3d(x = avec[length(avec)], y = bvec, z = 0, pl))
# lines(trans3d(x = avec, y = bvec[1], z = 0, pl))
# lines(trans3d(x = avec, y = bvec[length(bvec)], z = 0, pl))

# # What p-value distribution is closest to uniform?
# unif_dist <- function(x) {
#   L <- length(x)
#   dev <- sort(x) - 1:L / L
#   max(abs(dev))
# }

# # uniftests <- sapply(pv_list, function(x) apply(x, 2, unif_dist))
# # par(mar = c(0.1, 0.5, 0.1, 0.1), mgp = 2:0)
# # pl <- persp(avec, bvec, uniftests, theta = -90, col = "lightblue", ticktype = "detailed")

# If it seems like a = 0.5 is clear, then check p-value plots for b?
# First element of pv_list corresponds to a = 0.5.

load("par_ss100_pvals.rda") # pv_list has matrices of p-values
pv100 <- pv_list[[1]]
load("par_ss500_pvals.rda") # pv_list has matrices of p-values
pv500 <- pv_list[[1]]
load("par_ss1000_pvals.rda") # pv_list has matrices of p-values
pv1000 <- pv_list[[1]]

library(colorspace)
pal <- qualitative_hcl(10, palette = "Dark3")
pvp_fun <- function(vec) {sort(vec) - 1:length(vec) / length(vec)}
par(mfrow = c(4, 1), mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
matplot(apply(pv100, 2, pvp_fun), col = pal, type = "l", lty = 1, cex = 0.25, xlim = c(0, 200), ylim = c(-0.03, 0.01))
abline(h = 0, lty = 2, lwd = 3)
matplot(apply(pv500, 2, pvp_fun), col = pal, type = "l", lty = 1, cex = 0.25, xlim = c(0, 200), ylim = c(-0.03, 0.01))
abline(h = 0, lty = 2, lwd = 3)
matplot(apply(pv1000, 2, pvp_fun), col = pal, type = "l", lty = 1, cex = 0.25, xlim = c(0, 200), ylim = c(-0.03, 0.01))
abline(h = 0, lty = 2, lwd = 3)
plot(1:10, rep(0, 10), col = pal, pch = 20, cex = 3)
text(1:10, rep(-0.5, 10), labels = seq(0.5, 5, by = 0.5))
