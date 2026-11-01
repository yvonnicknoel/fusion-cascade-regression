## ============================================================
## Extended Monte Carlo tabulation of c(m; n, p, rho)
## ------------------------------------------------------------
## Extends the calibration box to n up to 10000 and p up to 500,
## using the fast vectorized calibrator (calibrate_c_fast.R).
##
## Runs incrementally by BATCH so the full compute need not be
## committed up front:
##   BATCH=extend_n  new n levels {2000,5000,10000} at old p (<=60)
##   BATCH=extend_p  new p levels {100,150,250,350,500} at old n (<=1000)
##   BATCH=corner    new n AND new p (the joint high corners)
##   BATCH=new       everything outside the original box (default)
##   BATCH=full      recompute the ENTIRE union grid (clean, one methodology)
##
## Env vars: BATCH, B (default 5000), CORES (default 90), SEED (default 20260706)
## Output  : c_table_ext_<BATCH>.csv   (plain CSV: alpha,n,p,rho,m,c,batch)
##           tabulate_c_extended_<BATCH>.log
## ============================================================

source("calibrate_c_fast.R")
suppressMessages(library(parallel))

## ---- configuration ----------------------------------------
BATCH <- Sys.getenv("BATCH", "new")
B     <- as.integer(Sys.getenv("B",     "5000"))
CORES <- as.integer(Sys.getenv("CORES", "90"))
SEED  <- as.integer(Sys.getenv("SEED",  "20260706"))
OUT   <- sprintf("c_table_ext_%s.csv", BATCH)
LOG   <- sprintf("tabulate_c_extended_%s.log", BATCH)

logmsg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", file = LOG, append = TRUE, sep = "")
}

## ---- grid definition --------------------------------------
n_old <- c(50, 100, 150, 200, 500, 1000)
p_old <- c(5, 10, 20, 30, 40, 60)
n_new <- c(2000, 5000, 10000)
p_new <- c(100, 150, 250, 350, 500)
n_all <- sort(union(n_old, n_new))
p_all <- sort(union(p_old, p_new))
rho_vec <- c(0.0, 0.3, 0.7)
alpha   <- 0.05

grid <- expand.grid(alpha = alpha, n = n_all, p = p_all, rho = rho_vec,
                    KEEP.OUT.ATTRS = FALSE)
grid <- subset(grid, p < n)
## classify each cell
is_old <- with(grid, n %in% n_old & p %in% p_old)
new_n  <- with(grid, n %in% n_new & p %in% p_old)
new_p  <- with(grid, p %in% p_new & n %in% n_old)
corner <- with(grid, n %in% n_new & p %in% p_new)
grid$batch <- ifelse(is_old, "old",
               ifelse(new_n, "extend_n",
               ifelse(new_p, "extend_p", "corner")))

sel <- switch(BATCH,
  extend_n = grid$batch == "extend_n",
  extend_p = grid$batch == "extend_p",
  corner   = grid$batch == "corner",
  new      = grid$batch != "old",
  full     = rep(TRUE, nrow(grid)),
  stop("Unknown BATCH: ", BATCH))
grid <- grid[sel, , drop = FALSE]
grid <- grid[order(grid$p, grid$n, grid$rho), ]   # cheap-first-ish
rownames(grid) <- NULL
n_scen <- nrow(grid)

## ---- AR(1) design generator (Cholesky, fast for large p) ---
gen_X_ar1 <- function(n, p, rho) {
  Z <- matrix(rnorm(n * p), n, p)
  if (rho != 0) {
    idx <- seq_len(p)
    L   <- chol(rho ^ abs(outer(idx, idx, "-")))  # upper chol
    Z   <- Z %*% L
  }
  colnames(Z) <- paste0("X", seq_len(p))
  Z
}

## ---- per-scenario worker ----------------------------------
run_one <- function(s_id) {
  a <- grid$alpha[s_id]; n <- grid$n[s_id]; p <- grid$p[s_id]; rho <- grid$rho[s_id]
  set.seed(SEED + 1000L * s_id)            # deterministic regardless of scheduling
  t0 <- proc.time()[3]
  X  <- gen_X_ar1(n, p, rho)
  cv <- calibrate_cvec_mc_fast(X, alpha = a, B = B, chunk = 500L)
  dt <- proc.time()[3] - t0
  logmsg("[%3d/%3d] n=%5d p=%4d rho=%.1f  %6.1fs  c(1)=%.3f c(%d)=%.3f",
         s_id, n_scen, n, p, rho, dt, cv[1], p, cv[p])
  data.frame(alpha = a, n = n, p = p, rho = rho, m = seq_len(p),
             c = cv, batch = grid$batch[s_id], stringsAsFactors = FALSE)
}

## ---- launch -----------------------------------------------
cat("", file = LOG)                        # truncate log
logmsg("BATCH=%s  scenarios=%d  B=%d  cores=%d  seed=%d",
       BATCH, n_scen, B, CORES, SEED)
logmsg("Grid: n in {%s}, p in {%s}, rho in {%s}",
       paste(n_all, collapse=","), paste(p_all, collapse=","),
       paste(rho_vec, collapse=","))
t_start <- proc.time()[3]

res <- mclapply(seq_len(n_scen), run_one,
                mc.cores = CORES, mc.preschedule = FALSE)

## surface any worker failures instead of silently dropping them
ok <- vapply(res, is.data.frame, logical(1))
if (any(!ok)) {
  logmsg("WARNING: %d scenario(s) failed:", sum(!ok))
  for (i in which(!ok)) logmsg("  scenario %d (n=%d p=%d rho=%.1f): %s",
                               i, grid$n[i], grid$p[i], grid$rho[i],
                               conditionMessage(attr(res[[i]], "condition")))
}
tab <- do.call(rbind, res[ok])
write.csv(tab, OUT, row.names = FALSE)
logmsg("Done in %.1f min. Wrote %d rows -> %s",
       (proc.time()[3] - t_start) / 60, nrow(tab), OUT)
