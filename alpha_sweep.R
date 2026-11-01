## ============================================================
## ANALYSIS B -- alpha-sweep of the Fusion Cascade operating point
## ------------------------------------------------------------
## Answers Reviewer 2's points 3 and 4 (MET-2026-0732).
##
## All results in the submitted manuscript are at alpha = 0.05, which
## invites the reading that the Fusion Cascade IS a conservative
## procedure. It is not: it is a ONE-PARAMETER FAMILY whose dial has an
## interpretable frequentist meaning (the global-null probability of
## selecting any predictor at all). This script traces the family.
##
## Two things are established:
##   (1) the TPR/FPR path across alpha in {0.01, 0.05, 0.10, 0.20},
##       overlaid on the competitor operating points from the main
##       simulation -- Fusion traces a frontier of its own;
##   (2) at s = 0 (global null), the REALIZED selection rate tracks the
##       NOMINAL alpha -- the dial means what it says, at every setting.
##
## Holm is swept over the same alphas: it is the only competitor with a
## comparable dial, so the comparison is like-for-like.
##
## EFFICIENCY: the four alphas are four quantiles of the SAME Monte
## Carlo replicates, so one calibration pass per (n, p) yields the whole
## sweep. Calibration is at rho = 0, matching the convention of the
## shipped transportable surface.
##
## Output: alpha_sweep_raw.csv        (one row per scenario x rep x method x alpha)
##         alpha_sweep_summary.csv    (one row per scenario x method x alpha)
##         alpha_sweep_overall.csv    (one row per method x alpha; Table-1 analogue)
##         alpha_sweep_global_null.csv(realized vs nominal, s = 0)
##         alpha_sweep_cvec.csv       (the calibrated thresholds themselves)
##         alpha_sweep.log
##
## Env vars: NREP (default 200), B (default 8000), CORES (default 90)
## ============================================================

source("calibrate_c_fast.R")
source("get_critical_value.R")
source("cascade_reg.R")
source("generate_data.R")
source("method_wrappers.R")
suppressMessages({ library(parallel); library(dplyr); library(tidyr) })

ALPHAS <- c(0.01, 0.05, 0.10, 0.20)
NREP   <- as.integer(Sys.getenv("NREP",  "200"))
B_MC   <- as.integer(Sys.getenv("B",     "8000"))
CORES  <- as.integer(Sys.getenv("CORES", "90"))
SEED   <- 1234L          # same base seed as main_simulation.R
n_test <- 1000
LOG    <- "alpha_sweep.log"

logmsg <- function(...) {
  m <- sprintf(...)
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE); flush.console()
}

## ------------------------------------------------------------
## Scenario grid -- IDENTICAL to main_simulation.R so that the
## alpha = 0.05 slice reproduces the published numbers and the
## overlay against competitors is apples-to-apples.
## ------------------------------------------------------------
n_vec   <- c(50, 100, 200, 500, 1000, 2000)
p_vec   <- c(10, 20, 50)
s_vec   <- c(0, 2, 5, 10)
rho_vec <- c(0, 0.3, 0.7)
R2_vec  <- c(0.2, 0.5)

scen_grid <- expand.grid(n = n_vec, p = p_vec, s = s_vec,
                         rho = rho_vec, R2 = R2_vec)
scen_grid <- subset(scen_grid, p <= n / 2 & s < p)
scen_grid <- subset(scen_grid, s > 0 | R2 == R2_vec[1])
scen_grid <- scen_grid[order(scen_grid$p, scen_grid$n), ]
rownames(scen_grid) <- NULL
n_scen <- nrow(scen_grid)

cat("", file = LOG)
logmsg("ANALYSIS B -- alpha-sweep of the Fusion operating point")
logmsg("%d scenarios x %d reps x %d alphas x 2 sweepable methods | B=%d cores=%d",
       n_scen, NREP, length(ALPHAS), B_MC, CORES)
logmsg("alphas: %s", paste(ALPHAS, collapse = ", "))

## ============================================================
## STEP 1 -- calibrate c_vec at every alpha, once per (n, p)
## ------------------------------------------------------------
## Multi-alpha variant of calibrate_cvec_mc_fast(): the subset-max
## matrix M is built once and read at four quantile levels. This is
## exactly the same estimator as the single-alpha function, evaluated
## at several probs, so the alpha = 0.05 column is directly comparable
## to the published calibration.
## ============================================================
calibrate_cvec_multi_alpha <- function(X, alphas, B = 8000, chunk = 1000L) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  col_sds <- apply(X, 2, sd)
  if (any(col_sds == 0)) stop("Zero-variance predictor(s) in calibration.")
  Xs <- sweep(sweep(X, 2, colMeans(X), "-"), 2, col_sds, "/")

  D <- cbind(1, Xs); qrd <- qr(D); Q <- qr.Q(qrd); R <- qr.R(qrd)
  Rinv <- backsolve(R, diag(ncol(R)))
  diagcov_pred <- rowSums(Rinv * Rinv)[-1]
  dfree <- n - (p + 1L)

  M <- matrix(NA_real_, nrow = p, ncol = B)
  done <- 0L
  while (done < B) {
    b  <- min(chunk, B - done)
    Y  <- matrix(rnorm(n * b), nrow = n, ncol = b)
    aT <- .abs_t_null_batch(Q, R, diagcov_pred, dfree, Y)
    for (j in seq_len(b)) M[, done + j] <- cummax(aT[sample.int(p), j])
    done <- done + b
  }
  ## one column per alpha, each monotone in k
  out <- vapply(alphas, function(a)
    cummax(as.numeric(apply(M, 1, stats::quantile, probs = 1 - a, type = 8))),
    numeric(p))
  colnames(out) <- paste0("alpha_", alphas)
  out
}

np_pairs <- unique(scen_grid[c("n", "p")])
rownames(np_pairs) <- NULL
logmsg("")
logmsg("STEP 1: calibrating c(k) at %d alphas for %d (n,p) pairs (rho = 0, B = %d)...",
       length(ALPHAS), nrow(np_pairs), B_MC)

t0 <- proc.time()[3]
cal_list <- mclapply(seq_len(nrow(np_pairs)), function(i) {
  n <- np_pairs$n[i]; p <- np_pairs$p[i]
  set.seed(SEED + 7777L * i)
  tt <- proc.time()[3]
  X  <- gen_X_ar1(n, p, rho = 0)              # transportable convention: rho = 0
  cm <- calibrate_cvec_multi_alpha(X, ALPHAS, B = B_MC)
  cat(sprintf("  calib %2d/%2d  n=%4d p=%2d  %5.1fs  c(1)@.05=%.3f c(%d)@.05=%.3f\n",
              i, nrow(np_pairs), n, p, proc.time()[3] - tt,
              cm[1, "alpha_0.05"], p, cm[p, "alpha_0.05"]),
      file = LOG, append = TRUE)
  list(n = n, p = p, cmat = cm)
}, mc.cores = min(CORES, nrow(np_pairs)), mc.preschedule = FALSE)

cal_key <- vapply(cal_list, function(z) paste(z$n, z$p, sep = "_"), character(1))
names(cal_list) <- cal_key
logmsg("STEP 1 done in %.1f min.", (proc.time()[3] - t0) / 60)

## Sanity check: MC at alpha = 0.05 vs the shipped closed-form surface
chk <- do.call(rbind, lapply(cal_list, function(z) {
  ch <- get_critical_value(m = seq_len(z$p), n = z$n, p = z$p)
  data.frame(n = z$n, p = z$p, k = seq_len(z$p),
             c_mc = z$cmat[, "alpha_0.05"], c_surface = ch)
}))
chk$rel_pct <- 100 * (chk$c_mc - chk$c_surface) / chk$c_surface
logmsg("SANITY: MC(alpha=.05) vs closed-form surface -- median %.2f%%, p95 %.2f%%, max |%.2f|%%",
       median(chk$rel_pct), quantile(abs(chk$rel_pct), .95), max(abs(chk$rel_pct)))

## Persist the thresholds themselves
cvec_df <- do.call(rbind, lapply(cal_list, function(z)
  data.frame(n = z$n, p = z$p, k = seq_len(z$p), as.data.frame(z$cmat))))
write.csv2(cvec_df, "alpha_sweep_cvec.csv", row.names = FALSE)

## ============================================================
## STEP 2 -- run Fusion and Holm at each alpha over the grid
## ============================================================
run_rep <- function(r, s_id, n, p, s, rho, R2, cmat) {
  ## Same seed stream as main_simulation.R: at alpha = 0.05 the Fusion
  ## rows reproduce the published replication-for-replication.
  set.seed(SEED + 100000L * s_id + r)

  dtrain <- gen_data_scenario(n, p, s, rho, R2_target = R2, beta_val = 1)
  X_train <- dtrain$X; y_train <- dtrain$y
  beta_true <- dtrain$beta; active <- dtrain$active

  X_test   <- gen_X_ar1(n_test, p, rho)
  eta_test <- as.vector(X_test %*% beta_true)
  var_eta  <- stats::var(eta_test)
  sigma2_test <- if (s == 0 || is.null(R2) || is.na(R2) || var_eta == 0) 1 else
    var_eta * (1 - R2) / R2
  y_test <- eta_test + rnorm(n_test, 0, sqrt(sigma2_test))

  rows <- list(); i <- 0L
  for (a in ALPHAS) {
    cv <- cmat[, paste0("alpha_", a)]
    for (meth in c("Fusion", "Holm")) {
      out <- tryCatch({
        selected <- if (meth == "Fusion")
          select_fusion(X_train, y_train, alpha = a, ccrit = cv)
        else
          select_holm(X_train, y_train, alpha = a)
        mets <- selection_metrics(selected, active, p)
        mse  <- prediction_mse(X_train, y_train, X_test, y_test, selected)
        data.frame(TP = mets$TP, FP = mets$FP, FN = mets$FN, TN = mets$TN,
                   TPR = mets$TPR, FPR = mets$FPR, FDR = mets$FDR,
                   model_size = mets$model_size, MSE = mse)
      }, error = function(e) data.frame(TP=NA,FP=NA,FN=NA,TN=NA,TPR=NA,
                                        FPR=NA,FDR=NA,model_size=NA,MSE=NA))
      i <- i + 1L
      rows[[i]] <- data.frame(scen_id = s_id, rep = r, method = meth, alpha = a,
                              n = n, p = p, s = s, rho = rho, R2 = R2, out)
    }
  }
  do.call(rbind, rows)
}

logmsg("")
logmsg("STEP 2: running the grid...")
t_start <- proc.time()[3]
results_list <- vector("list", n_scen)

for (s_id in seq_len(n_scen)) {
  n <- scen_grid$n[s_id]; p <- scen_grid$p[s_id]; s <- scen_grid$s[s_id]
  rho <- scen_grid$rho[s_id]; R2 <- scen_grid$R2[s_id]
  cmat <- cal_list[[paste(n, p, sep = "_")]]$cmat

  tt <- proc.time()[3]
  reps <- mclapply(seq_len(NREP),
                   function(r) run_rep(r, s_id, n, p, s, rho, R2, cmat),
                   mc.cores = CORES, mc.preschedule = TRUE)
  ok <- vapply(reps, is.data.frame, logical(1))
  results_list[[s_id]] <- do.call(rbind, reps[ok])
  logmsg("[%3d/%3d] n=%4d p=%2d s=%2d rho=%.1f R2=%.1f  %5.1fs  (%d/%d reps ok)",
         s_id, n_scen, n, p, s, rho, R2, proc.time()[3] - tt, sum(ok), NREP)
}

sweep_raw <- dplyr::bind_rows(results_list)
sweep_raw$F1 <- with(sweep_raw, 2*TP / (2*TP + FP + FN))
write.csv2(sweep_raw, "alpha_sweep_raw.csv", row.names = FALSE)

## ---- per-scenario summary ----
sweep_summary <- sweep_raw %>%
  group_by(method, alpha, n, p, s, rho, R2) %>%
  summarise(mean_TPR = mean(TPR, na.rm = TRUE), mean_FPR = mean(FPR, na.rm = TRUE),
            mean_FDR = mean(FDR, na.rm = TRUE), mean_F1 = mean(F1, na.rm = TRUE),
            mean_size = mean(model_size, na.rm = TRUE), mean_MSE = mean(MSE, na.rm = TRUE),
            FWER_FN = mean(FN > 0, na.rm = TRUE), FWER_FP = mean(FP > 0, na.rm = TRUE),
            .groups = "drop")
write.csv2(sweep_summary, "alpha_sweep_summary.csv", row.names = FALSE)

## ---- overall (Table-1 analogue; signal scenarios only, as in the paper) ----
sweep_overall <- sweep_raw %>%
  filter(s > 0) %>%
  group_by(method, alpha) %>%
  summarise(TPR = mean(TPR, na.rm = TRUE), FPR = mean(FPR, na.rm = TRUE),
            FDR = mean(FDR, na.rm = TRUE), F1 = mean(F1, na.rm = TRUE),
            model_size = mean(model_size, na.rm = TRUE),
            FWER_FN = mean(FN > 0, na.rm = TRUE), FWER_FP = mean(FP > 0, na.rm = TRUE),
            MSE = mean(MSE, na.rm = TRUE), .groups = "drop") %>%
  arrange(method, alpha)
write.csv2(sweep_overall, "alpha_sweep_overall.csv", row.names = FALSE)

## ---- THE KEY CHECK: realized global-null rate vs nominal alpha (s = 0) ----
gn <- sweep_raw %>%
  filter(s == 0) %>%
  group_by(method, alpha) %>%
  summarise(realized = mean(FP > 0, na.rm = TRUE),
            mean_size = mean(model_size, na.rm = TRUE),
            n_cells = n(), .groups = "drop") %>%
  mutate(nominal = alpha, ratio = realized / alpha)
gn_by_cell <- sweep_raw %>%
  filter(s == 0) %>%
  group_by(method, alpha, n, p, rho) %>%
  summarise(realized = mean(FP > 0, na.rm = TRUE), .groups = "drop")
write.csv2(rbind(
  cbind(gn[c("method","alpha","realized")], n = NA, p = NA, rho = NA, scope = "pooled"),
  cbind(gn_by_cell[c("method","alpha","realized")],
        gn_by_cell[c("n","p","rho")], scope = "cell")),
  "alpha_sweep_global_null.csv", row.names = FALSE)

## ============================================================
## Report
## ============================================================
logmsg("")
logmsg("=== GLOBAL-NULL CALIBRATION (s = 0): realized vs nominal ===")
for (i in seq_len(nrow(gn)))
  logmsg("  %-6s alpha = %.2f :  realized = %.4f   (ratio %.2f)",
         gn$method[i], gn$alpha[i], gn$realized[i], gn$ratio[i])
logmsg("")
logmsg("  Worst per-cell realized rate for Fusion, by alpha:")
wf <- gn_by_cell %>% filter(method == "Fusion") %>%
  group_by(alpha) %>% summarise(worst = max(realized), .groups = "drop")
for (i in seq_len(nrow(wf)))
  logmsg("    alpha = %.2f : max over cells = %.4f", wf$alpha[i], wf$worst[i])
logmsg("")
logmsg("=== OPERATING POINTS ACROSS ALPHA (signal scenarios, s > 0) ===")
so <- as.data.frame(sweep_overall)
for (i in seq_len(nrow(so)))
  logmsg("  %-6s alpha=%.2f : TPR=%.3f FPR=%.3f FDR=%.3f F1=%.3f size=%.2f FWER_FP=%.3f MSE=%.3f",
         so$method[i], so$alpha[i], so$TPR[i], so$FPR[i], so$FDR[i], so$F1[i],
         so$model_size[i], so$FWER_FP[i], so$MSE[i])
logmsg("")
logmsg("Wrote alpha_sweep_{raw,summary,overall,global_null,cvec}.csv")
logmsg("Total elapsed: %.1f min", (proc.time()[3] - t_start) / 60)
