## ============================================================
## Fast, vectorized Monte Carlo calibration of c_vec(m), m = 1..p
## ------------------------------------------------------------
## Mathematically IDENTICAL to calibrate_cvec_mc() in
## calibrate_c_reg.R (full-model t-statistics under the global
## null; random size-m subset max-|t| via permutation + cummax),
## but computed with matrix algebra against a single precomputed
## QR of the fixed standardized design, so that even n = 10000,
## p = 500 is cheap.
##
## Per replicate the old code did one lm(y ~ .) fit: O(n p^2).
## Here the QR is factorized ONCE per design; each null draw is
## then O(n p) (one matvec Q'y) + O(p^2) back-substitution, and
## all B draws are batched.
## ============================================================

## Vectorized full-model |t|-statistics for B null responses.
##  D    : n x (p+1) intercept-augmented STANDARDIZED design
##  Returns p x B matrix of |t| for the p predictors (intercept dropped).
.abs_t_null_batch <- function(Q, R, diagcov_pred, dfree, Y) {
  ## Q: n x (p+1) orthonormal; R: (p+1) x (p+1) upper triangular
  ## Y: n x b matrix of N(0,1) null responses
  Z    <- crossprod(Q, Y)                 # (p+1) x b   = Q'Y
  RSS  <- colSums(Y * Y) - colSums(Z * Z) # length b    = ||(I-H)y||^2
  Bhat <- backsolve(R, Z)                 # (p+1) x b   solves R Bhat = Z
  s2   <- RSS / dfree                      # length b    sigma_hat^2
  ## SE_jb = sqrt(s2_b * diagcov_j); t = Bhat / SE  (predictors only)
  Bp   <- Bhat[-1, , drop = FALSE]         # p x b
  denom <- sqrt(outer(diagcov_pred, s2))   # p x b
  abs(Bp / denom)
}

calibrate_cvec_mc_fast <- function(X, alpha = 0.05, B = 5000,
                                   chunk = 1000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
  if (p < 1) stop("Need p >= 1.")

  ## Standardize EXACTLY as fusion_cascade_reg / calibrate_cvec_mc:
  col_means <- colMeans(X)
  col_sds   <- apply(X, 2, sd)
  if (any(col_sds == 0)) stop("Zero-variance predictor(s) in calibration.")
  Xs <- sweep(sweep(X, 2, col_means, "-"), 2, col_sds, "/")

  ## Intercept-augmented design + one QR
  D  <- cbind(1, Xs)                       # n x (p+1)
  qr <- qr(D)
  Q  <- qr.Q(qr)                           # n x (p+1)
  R  <- qr.R(qr)                           # (p+1) x (p+1)
  ## diag of (D'D)^{-1} = Rinv Rinv';  Rinv = R^{-1}
  Rinv     <- backsolve(R, diag(ncol(R)))
  diagcov  <- rowSums(Rinv * Rinv)
  diagcov_pred <- diagcov[-1]              # predictors only, length p
  dfree    <- n - (p + 1L)                 # residual df, matches lm

  ## M[m, b] = max|t| over a random size-m subset for replicate b
  M <- matrix(NA_real_, nrow = p, ncol = B)
  done <- 0L
  while (done < B) {
    b  <- min(chunk, B - done)
    Y  <- matrix(rnorm(n * b), nrow = n, ncol = b)
    aT <- .abs_t_null_batch(Q, R, diagcov_pred, dfree, Y)   # p x b
    ## one independent permutation per replicate, then cummax
    for (j in seq_len(b)) {
      M[, done + j] <- cummax(aT[sample.int(p), j])
    }
    done <- done + b
  }

  c_vec <- apply(M, 1, stats::quantile, probs = 1 - alpha, type = 8)
  cummax(as.numeric(c_vec))                # enforce monotonicity in m
}


## ============================================================
## Self-test: fast method must reproduce lm()-based t-statistics
## (statistic-for-statistic, same y) up to floating point.
## Run:  Rscript calibrate_c_fast.R
## ============================================================
if (sys.nframe() == 0L) {
  set.seed(1)
  n <- 60L; p <- 8L
  X  <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("X", 1:p)
  Xs <- sweep(sweep(X, 2, colMeans(X), "-"), 2, apply(X, 2, sd), "/")
  D  <- cbind(1, Xs); qr <- qr(D); Q <- qr.Q(qr); R <- qr.R(qr)
  Rinv <- backsolve(R, diag(ncol(R))); dc <- rowSums(Rinv * Rinv)[-1]

  max_diff <- 0
  for (rep in 1:200) {
    y <- rnorm(n)
    t_lm <- abs(coef(summary(lm(y ~ Xs)))[-1, "t value"])
    t_fast <- as.numeric(.abs_t_null_batch(Q, R, dc, n - (p + 1L),
                                           matrix(y, ncol = 1)))
    max_diff <- max(max_diff, max(abs(t_lm - t_fast)))
  }
  cat(sprintf("Self-test max |t_lm - t_fast| over 200 draws: %.3e\n", max_diff))
  stopifnot(max_diff < 1e-8)
  cat("OK: fast t-statistics reproduce lm() exactly.\n")

  ## Timing sanity at the largest target design
  gen <- function(n, p) matrix(rnorm(n * p), n, p)
  for (cfg in list(c(1000, 60), c(10000, 500))) {
    nn <- cfg[1]; pp <- cfg[2]
    t0 <- proc.time()[3]
    cv <- calibrate_cvec_mc_fast(gen(nn, pp), B = 5000, chunk = 500, seed = 7)
    cat(sprintf("n=%5d p=%3d : %.1fs  (c(1)=%.3f c(%d)=%.3f)\n",
                nn, pp, proc.time()[3] - t0, cv[1], pp, cv[pp]))
  }
}
