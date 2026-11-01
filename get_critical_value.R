## ============================================================
## FUNCTION TO OBTAIN CRITICAL VALUES c(m; n, p)
## Model I4x (extended):
##   log c ~ [ logn + logp + logn*logp + p/n + (p/n)^2 + (p/n)^3
##             + 1/(n-p-1) ]   (7 terms + intercept = 8 functions of n,p)
##           x { 1, logm, logm^2, logm^3 }        = 32 parameters
##
## The residual-degrees-of-freedom term 1/(n-p-1) is what makes the
## approximation hold in the p ~ n regime: when few residual df remain
## the t-statistics have heavy tails and c must inflate sharply, which a
## smooth polynomial in (logn, logp, p/n) cannot reproduce.
## ============================================================

#' Get critical values for the Fusion Cascade in regression
#'
#' @param m Vector of family sizes (1 to p)
#' @param n Sample size
#' @param p Number of predictors
#' @param alpha Nominal level (default 0.05)
#' @return Vector of critical values c(m; n, p) for m = 1, ..., length(m)
#'
#' @details
#' Parametric approximation of Monte Carlo critical values, calibrated
#' under the global null (beta = 0) on an AR(1), rho = 0 design over the
#' grid 50 <= n <= 10000, 5 <= p <= 500 (p < n), B = 8000 replications.
#'
#' Model performance (rho = 0, 8025 calibration points):
#'   R^2 = 0.99999
#'   RMSE = 0.0043 (log scale)
#'   relative error on c: median 0.24%, 90th pct 0.68%, max 3.2%
#'
#' Robustness to AR(1) correlation: the rho = 0 model tracks calibrated
#' values for rho in {0.3, 0.7} to within ~3% (roughly centered, i.e. NOT
#' strictly conservative on the extended grid -- a small two-sided
#' deviation). If strict conservativeness is required, recalibrate for the
#' target rho or add a safety margin.
#'
#' Validity domain:
#'   50 <= n <= 10000
#'   5  <= p <= 500
#'   p < n
#'
#' Outside this domain, recalibrate via Monte Carlo with
#' calibrate_cvec_mc() (or the fast calibrate_cvec_mc_fast()).
#'
#' @examples
#' c_vec <- get_critical_value(m = 1:20, n = 100, p = 20)
#' c_10  <- get_critical_value(m = 10,   n = 100, p = 20)
#'
#' @export
get_critical_value <- function(m, n, p, alpha = 0.05) {

  ## ---- checks ----
  if (any(m < 1) || any(m > p)) stop("m must be between 1 and p")
  if (p >= n) warning("p >= n: model unreliable. Use Monte Carlo recalibration.")
  if (n < 50 || n > 10000 || p < 5 || p > 500)
    warning("Configuration outside calibration domain (50<=n<=10000, 5<=p<=500). ",
            "Interpret with caution or recalibrate.")
  if (alpha != 0.05)
    warning("Model calibrated for alpha = 0.05; recalibrate for other levels.")

  ## ---- 32 coefficients, canonical order ----
  ## For g in {1, logn, logp, logn*logp, p/n, (p/n)^2, (p/n)^3, 1/(n-p-1)},
  ## for logm power i in 0:3.
  coefs_I4x <- c(
    0.7054284794, 0.1437927441, 0.001498086222, -0.001410593589,
    -0.00549368934, 0.009345604222, -0.004126503409, 0.0004903252236,
    -0.004499392673, 0.009598396503, -0.004328679397, 0.0005403726442,
    0.001034283422, -0.001842456591, 0.0008179671675, -0.0001008237671,
    -0.06740363138, 0.09356878396, -0.03784321301, 0.003879507247,
    0.1655793457, -0.2572207047, 0.1186475656, -0.01366102945,
    -0.1158347286, 0.2007287284, -0.1005315503, 0.01216643169,
    1.156370102, 0.3633966919, -0.01017409274, -0.003331390545
  )

  ## ---- (n,p) basis (length 8), scalar n,p ----
  logn <- log(n); logp <- log(p); pn <- p / n
  g <- c(1, logn, logp, logn * logp, pn, pn^2, pn^3, 1 / (n - p - 1))

  ## coefficient matrix 8 x 4 (row = basis function, col = logm power)
  cf <- matrix(coefs_I4x, nrow = 8, ncol = 4, byrow = TRUE)

  ## effective coefficients on {1, logm, logm^2, logm^3}
  gc <- as.numeric(g %*% cf)                       # length 4

  ## evaluate over m
  lm1  <- log(m)
  logc <- gc[1] + gc[2] * lm1 + gc[3] * lm1^2 + gc[4] * lm1^3
  exp(as.numeric(logc))
}

## ============================================================
## USAGE EXAMPLES
## ============================================================
if (FALSE) {
  # Standard configuration
  c_vec <- get_critical_value(m = 1:20, n = 100, p = 20)
  print(head(c_vec, 10))

  # Extended domain (previously out of range)
  c_big <- get_critical_value(m = 1:500, n = 10000, p = 500)

  # Use in fusion_cascade_reg
  source("cascade_reg.R")
  set.seed(123); n <- 100; p <- 20
  X <- matrix(rnorm(n * p), n, p); y <- rnorm(n)
  c_vec  <- get_critical_value(m = 1:p, n = n, p = p)
  result <- fusion_cascade_reg(X, y, alpha = 0.05, c_vec = c_vec)
  print(result$active_predictors)
}
