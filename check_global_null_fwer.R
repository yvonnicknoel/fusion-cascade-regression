## ============================================================
## Regression check for the CENTRAL guarantee of the Fusion Cascade:
## under the global null (all beta = 0), P(select any predictor) <= alpha.
## Runs the exact production path -- select_fusion() with the closed-form
## thresholds get_critical_value(1:p, n, p) -- so it guards the claim the
## manuscript makes. Fails loudly if any cell exceeds alpha by > 3 MC SE.
## Run:  Rscript check_global_null_fwer.R
## ============================================================
source("cascade_reg.R"); source("get_critical_value.R")
source("generate_data.R"); source("method_wrappers.R")

alpha <- 0.05
B     <- as.integer(Sys.getenv("B", "2000"))
grid  <- expand.grid(n = c(200, 1000), p = c(10, 20), rho = c(0, 0.7))

fwer_null <- function(sel_fun, n, p, rho, B) {
  X <- gen_X_ar1(n, p, rho); cc <- get_critical_value(1:p, n, p, alpha)
  hit <- 0L
  for (b in seq_len(B)) {
    y   <- rnorm(n)                                   # pure global null
    sel <- sel_fun(X, y, cc)
    if (length(sel) > 0) hit <- hit + 1L
  }
  hit / B
}
fusion_sel <- function(X, y, cc) suppressMessages(select_fusion(X, y, alpha = alpha, ccrit = cc))
holm_sel   <- function(X, y, cc) select_holm(X, y, alpha = alpha)

set.seed(20260707)
res <- data.frame()
for (i in seq_len(nrow(grid))) {
  n <- grid$n[i]; p <- grid$p[i]; rho <- grid$rho[i]
  fF <- fwer_null(fusion_sel, n, p, rho, B)
  fH <- fwer_null(holm_sel,   n, p, rho, B)
  se <- sqrt(alpha * (1 - alpha) / B)
  res <- rbind(res, data.frame(n, p, rho, Fusion = fF, Holm = fH,
                               pass = fF <= alpha + 3 * se))
  cat(sprintf("n=%4d p=%2d rho=%.1f | Fusion FWER=%.3f  Holm=%.3f  %s\n",
              n, p, rho, fF, fH, if (fF <= alpha + 3 * se) "OK" else "*** EXCEEDS alpha ***"))
}
cat(sprintf("\nMC SE at B=%d is %.4f;  tolerance = alpha + 3 SE = %.4f\n",
            B, sqrt(alpha*(1-alpha)/B), alpha + 3*sqrt(alpha*(1-alpha)/B)))
if (all(res$pass)) cat("ALL CELLS PASS: global-null FWER controlled at alpha.\n") else
  stop("REGRESSION: at least one cell exceeds the global-null FWER tolerance.")
