# Make a table of empirical rejection probabilities from the different tuning
# parameter choices.

alpha <- 0.05
nvec <- c(100, 500, 1000)
Rvec <- c(499, 999, 1999)

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
r500_05 <- sapply(pv500, emp_check, al = 0.05)
r1000_05 <- sapply(pv1000, emp_check, al = 0.05)
r100_10 <- sapply(pv100, emp_check, al = 0.10)
r500_10 <- sapply(pv500, emp_check, al = 0.10)
r1000_10 <- sapply(pv1000, emp_check, al = 0.10)
tab_05 <- rbind(r100_05, r500_05, r1000_05)
tab_10 <- rbind(r100_10, r500_10, r1000_10)
dimnames(tab_05) <- list(paste("$c_b$ =", rep(avec, 3)), paste("$c_a$ =", bvec))
dimnames(tab_10) <- list(paste("$c_b$ =", rep(avec, 3)), paste("$c_a$ =", bvec))

library(Hmisc)

cap_05 <- paste0("Empirical minus nominal 5\\% rejection frequency for a test that $L_{X+Y} \\equiv L_0$ for known $L_0$.  Numbers in the table closer to zero are better.  Rows describe tuning parameters for contact set estimation and columns for epsilon-maximizer estimation.")
cap_10 <- paste0("Empirical minus nominal 10\\% rejection frequency for a test that $L_{X+Y} \\equiv L_0$ for known $L_0$.  Numbers in the table closer to zero are better.  Rows describe tuning parameters for contact set estimation and columns for epsilon-maximizer estimation.")
tab_p05 <- round(tab_05 * 100, digits = 2)
tab_p10 <- round(tab_10 * 100, digits = 2)

latex(tab_p05, file = paste0("../../fig/tab_05.tex"),
  size = "footnotesize", rgroup = paste("s. size", nvec),
  caption = cap_05, caption.loc = "bottom", label = "size_05_tab")

latex(tab_p10, file = paste0("../../fig/tab_10.tex"),
  size = "footnotesize", rgroup = paste("s. size", nvec),
  caption = cap_10, caption.loc = "bottom", label = "size_10_tab")

