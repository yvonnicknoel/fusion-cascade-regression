## ============================================================
## Negatively correlated predictors (second revision, AE-9).
##
## The cascade works on |t|, and reversing the sign of a column X_j
## negates t_j and every correlation involving X_j while leaving all
## |t| unchanged. Any correlation structure that can be turned into a
## positive manifold by sign reversals therefore yields IDENTICAL
## critical values; AR(1) with -rho is such a case (reverse every
## second column). This script
##   (1) verifies that identity numerically on one AR(1) design, and
##       shows that independent +rho / -rho designs give the same
##       global-null selection rate;
##   (2) tests structures that are NOT sign-reducible: blocks of four
##       predictors mutually correlated at -0.30 (the most negative
##       equicorrelation a 4x4 block admits is -1/3), a dense 3-factor
##       model with mixed-sign loadings, and weak negative
##       equicorrelation across all 20 predictors (-0.05; the bound is
##       -1/(p-1)). For each: the realized global-null selection rate
##       of the cascade with the closed-form (rho = 0) thresholds, and
##       the signed relative difference between those thresholds and a
##       Monte Carlo recalibration on the design itself.
## Outputs: negative_correlation_check.csv, negative_correlation_check.log
## Env: B (10000 null draws per design), BCAL (8000 calibration draws), CORES.
## ============================================================
suppressMessages({ library(MASS); library(parallel) })
source("cascade_reg.R"); source("get_critical_value.R"); source("method_wrappers.R"); source("calibrate_c_fast.R")
B <- as.integer(Sys.getenv("B", "10000")); BCAL <- as.integer(Sys.getenv("BCAL", "8000"))
CORES <- as.integer(Sys.getenv("CORES", "40"))
n <- 200; p <- 20; alpha <- 0.05
LOG <- "negative_correlation_check.log"; cat("", file = LOG)
logm <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
t0 <- Sys.time(); el <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

ar1   <- function(p, r) r^abs(outer(1:p, 1:p, "-"))
block <- function(p, bs, r) { S <- diag(p); for (s in seq(1, p, by = bs)) { i <- s:min(s + bs - 1, p); S[i, i] <- r }; diag(S) <- 1; S }
equi  <- function(p, r) { S <- matrix(r, p, p); diag(S) <- 1; S }
mixed_factor <- function(p, q = 3, seed = 11) {  # dense, mixed-sign loadings, full rank
  set.seed(seed); L <- matrix(rnorm(p * q, 0, 0.55), p, q)
  h <- rowSums(L^2); L <- L * sqrt(pmin(1, 0.85 / h)); L %*% t(L) + diag(1 - rowSums(L^2)) }
mkX <- function(S, seed) { set.seed(seed); X <- MASS::mvrnorm(n, rep(0, nrow(S)), S); colnames(X) <- paste0("X", 1:ncol(X)); X }
sel <- function(expr) length(tryCatch(expr, error = function(e) integer(0))) > 0

## share of predictor pairs whose correlation sign cannot be made positive by any column sign reversal
## Balanced-signed-graph test: assign column signs by breadth-first search over the
## nonzero correlations; a sign-reducible matrix yields no violated pair, and the share of
## violated pairs is reported otherwise (0 = reducible to a positive manifold).
irreducibility <- function(S) { p <- nrow(S); A <- sign(S); diag(A) <- 0; s <- rep(NA_real_, p)
  for (root in 1:p) { if (!is.na(s[root])) next; s[root] <- 1; queue <- root
    while (length(queue)) { i <- queue[1]; queue <- queue[-1]
      for (j in which(A[i, ] != 0)) if (is.na(s[j])) { s[j] <- s[i] * A[i, j]; queue <- c(queue, j) } } }
  pr <- which(upper.tri(S) & A != 0, arr.ind = TRUE)
  if (nrow(pr) == 0) return(0); mean(A[pr] != s[pr[, 1]] * s[pr[, 2]]) }
neg_share <- function(S) mean(S[upper.tri(S)] < 0)

null_rate <- function(X) {
  cc <- get_critical_value(1:p, n, p)
  mean(unlist(mclapply(1:B, function(b) { set.seed(5000 + b); sel(select_fusion(X, rnorm(n), alpha, ccrit = cc)) }, mc.cores = CORES)))
}
recal_err <- function(X) {               # signed relative error of the closed form vs MC on this design
  cc <- get_critical_value(1:p, n, p); cm <- calibrate_cvec_mc_fast(X, alpha, B = BCAL, seed = 99)
  e <- (cm - cc) / cm; c(median = median(e), max_abs = max(abs(e)))
}
SE <- sqrt(alpha * (1 - alpha) / B); tol <- alpha + 3 * SE
logm("n=%d p=%d, B=%d null draws per design, BCAL=%d, %d cores; tol = alpha + 3SE = %.3f", n, p, B, BCAL, CORES, tol)

## ---- (1) the sign identity on AR(1) ----
Xp <- mkX(ar1(p, 0.7), 42)
Xn <- Xp %*% diag((-1)^(1:p))             # same design, every second column reversed: corr = (-0.7)^|j-k|
colnames(Xn) <- colnames(Xp)
cp <- calibrate_cvec_mc_fast(Xp, alpha, B = BCAL, seed = 99); cn <- calibrate_cvec_mc_fast(Xn, alpha, B = BCAL, seed = 99)
logm("\n(1) Sign identity. Max |c_+(k) - c_-(k)| after reversing every second column of an AR(1) rho=0.7 design: %.2e", max(abs(cp - cn)))
logm("    realized corr(X1,X2): %.3f (original) vs %.3f (reversed)", cor(Xp)[1, 2], cor(Xn)[1, 2])

## ---- (2) designs ----
designs <- list(
  "AR(1) rho=+0.7"              = list(S = ar1(p, 0.7),           seed = 42),
  "AR(1) rho=-0.7"              = list(S = ar1(p, -0.7),          seed = 43),
  "AR(1) rho=+0.9"              = list(S = ar1(p, 0.9),           seed = 50),
  "AR(1) rho=-0.9"              = list(S = ar1(p, -0.9),          seed = 44),
  "Blocks(5x4) at -0.30"        = list(S = block(p, 4, -0.30),    seed = 45),
  "Blocks(5x4) at +0.30"        = list(S = block(p, 4, 0.30),     seed = 46),
  "3-factor, mixed-sign loadings" = list(S = mixed_factor(p),     seed = 47),
  "Equicorr -0.05 (all pairs)"  = list(S = equi(p, -0.05),        seed = 48),
  "Independent (reference)"     = list(S = diag(p),               seed = 49))
rows <- list()
for (nm in names(designs)) {
  S <- designs[[nm]]$S; X <- mkX(S, designs[[nm]]$seed)
  fw <- null_rate(X); re <- recal_err(X)
  rows[[nm]] <- data.frame(design = nm, neg_share = neg_share(S), irreducible_share = irreducibility(S),
                           null_rate = fw, closed_form_rel_err_median = re["median"], closed_form_rel_err_maxabs = re["max_abs"])
  logm("  %-30s neg pairs %4.0f%%  irreducible %4.0f%%  null rate %.3f %s  closed-form vs MC: median %+.2f%%, max |.| %.2f%%   (%s)",
       nm, 100 * neg_share(S), 100 * irreducibility(S), fw, if (fw <= tol) "OK " else "***",
       100 * re["median"], 100 * re["max_abs"], el())
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, "negative_correlation_check.csv", row.names = FALSE)
logm("\nDone in %s.", el())
