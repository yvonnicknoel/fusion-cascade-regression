## Wall-clock cost of a per-design Monte Carlo calibration (second revision, AE-15).
source("calibrate_c_fast.R"); source("calibrate_c_reg.R")
LOG <- "time_calibration.log"; cat("", file = LOG)
logm <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
logm("Single core, B = 8000 null draws, calibrate_cvec_mc_fast():")
for (np in list(c(200, 20), c(1000, 50), c(3536, 135), c(10000, 500))) {
  n <- np[1]; p <- np[2]; set.seed(1); X <- matrix(rnorm(n * p), n, p)
  tt <- system.time(calibrate_cvec_mc_fast(X, 0.05, B = 8000, seed = 7))[["elapsed"]]
  logm("  n = %5d, p = %3d : %7.1f s", n, p, tt)
}
n <- 1000; p <- 50; set.seed(1); X <- matrix(rnorm(n * p), n, p)
tt <- system.time(calibrate_cvec_mc(X, 0.05, B = 2000))[["elapsed"]]
logm("Reference lm()-based calibrate_cvec_mc(), n = 1000, p = 50, B = 2000: %.1f s  (x4 for B = 8000: ~%.0f s)", tt, 4 * tt)
logm("Done.")
