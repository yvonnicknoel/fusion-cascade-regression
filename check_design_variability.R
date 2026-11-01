## ============================================================
## ANALYSIS A -- Design-realization variability of c(k)
## ------------------------------------------------------------
## Answers the Associate Editor's question (MET-2026-0732) about
## fixed vs random predictors.
##
## The Fusion Cascade's guarantee is CONDITIONAL on the realized
## design matrix X: the calibration draws one X, holds it fixed,
## and averages over Y only. The *transportable* surface
## c(k; n, p) is however a function of (n, p) alone, so it
## necessarily marginalizes over the design distribution.
##
## The error committed by an analyst who uses the transportable
## surface instead of calibrating on their own X is therefore the
## DESIGN-TO-DESIGN VARIABILITY of c(k). This script measures it:
## for each (n, p, rho) cell we draw R independent design matrices,
## calibrate c(k) on each, and report the across-draw spread.
##
## Output: design_variability_raw.csv   (one row per cell x draw x k)
##         design_variability_cells.csv (one row per cell)
##         design_variability.log
##
## Env vars: R_DRAWS (default 20), B (default 8000), CORES (default 90)
## ============================================================

source("calibrate_c_fast.R")
source("get_critical_value.R")
source("generate_data.R")
suppressMessages({ library(parallel); library(dplyr) })

R_DRAWS <- as.integer(Sys.getenv("R_DRAWS", "20"))
B_MC    <- as.integer(Sys.getenv("B",       "8000"))
CORES   <- as.integer(Sys.getenv("CORES",   "90"))
SEED    <- 20260821L
LOG     <- "design_variability.log"

logmsg <- function(...) {
  m <- sprintf(...)
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE); flush.console()
}

## ------------------------------------------------------------
## Cells: the (n, p) pairs used in the paper's simulation and
## calibration grids, plus two deliberately hard p ~ n/2 corners
## where design variability should be largest (small residual df).
## ------------------------------------------------------------
base_grid <- expand.grid(n = c(50, 100, 200, 500, 1000),
                         p = c(5, 10, 20, 50))
base_grid <- subset(base_grid, p <= n / 2)
corner    <- data.frame(n = c(50, 100), p = c(40, 60))   # df = 9, 39
np_grid   <- rbind(base_grid[c("n", "p")], corner)

cells <- merge(np_grid, data.frame(rho = c(0, 0.3, 0.7)))
cells <- cells[order(cells$p, cells$n, cells$rho), ]
rownames(cells) <- NULL
n_cells <- nrow(cells)

## ------------------------------------------------------------
## Job list. TWO arms per cell, so that the two sources of spread
## can be separated:
##
##   arm = "design" : R draws, each a FRESH X and a fresh MC stream.
##                    Spread = design variability + Monte Carlo error.
##   arm = "mc"     : R draws, all on ONE FIXED X, fresh MC stream only.
##                    Spread = Monte Carlo error alone.
##
## Since the two sources are independent,
##     Var_design = max(0, Var_total - Var_mc),
## which is what an analyst using the transportable surface instead of
## calibrating on their own X actually incurs. Reporting the raw
## across-draw spread alone would overstate it, because a large part of
## it is simply the Monte Carlo error of a (1-alpha) quantile at B
## replicates -- an error the analyst would also incur by calibrating.
## ------------------------------------------------------------
jobs <- do.call(rbind, lapply(seq_len(n_cells), function(i)
  rbind(
    data.frame(cell_id = i, arm = "design", draw = seq_len(R_DRAWS),
               n = cells$n[i], p = cells$p[i], rho = cells$rho[i]),
    data.frame(cell_id = i, arm = "mc", draw = seq_len(R_DRAWS),
               n = cells$n[i], p = cells$p[i], rho = cells$rho[i]))))
n_jobs <- nrow(jobs)

cat("", file = LOG)
logmsg("ANALYSIS A -- design-realization variability of c(k)")
logmsg("%d cells (n,p,rho) x %d draws x 2 arms (design | mc-only) = %d calibrations | B=%d cores=%d",
       n_cells, R_DRAWS, n_jobs, B_MC, CORES)

## ------------------------------------------------------------
## One calibration.
##  arm "design": fresh X each draw.
##  arm "mc"    : the cell's single reference X every draw (seeded on
##                cell_id alone), only the null responses differ.
## ------------------------------------------------------------
run_job <- function(j) {
  n <- jobs$n[j]; p <- jobs$p[j]; rho <- jobs$rho[j]
  cid <- jobs$cell_id[j]; arm <- as.character(jobs$arm[j]); dr <- jobs$draw[j]

  if (arm == "design") {
    set.seed(SEED + 1000L * cid + dr)           # fresh design realization
    X <- gen_X_ar1(n, p, rho)
  } else {
    set.seed(SEED + 1000L * cid)                # the cell's fixed reference X
    X <- gen_X_ar1(n, p, rho)
  }
  ## independent MC stream in both arms
  set.seed(SEED + 500000L + 1000L * cid + dr + (arm == "mc") * 250000L)

  t0 <- proc.time()[3]
  cv <- calibrate_cvec_mc_fast(X, alpha = 0.05, B = B_MC, chunk = 1000L)
  el <- proc.time()[3] - t0
  cat(sprintf("  job %4d/%4d [%-6s] n=%5d p=%3d rho=%.1f draw=%2d  %5.1fs\n",
              j, n_jobs, arm, n, p, rho, dr, el),
      file = LOG, append = TRUE)
  data.frame(cell_id = cid, arm = arm, draw = dr,
             n = n, p = p, rho = rho, k = seq_len(p), c = cv)
}

t_start <- proc.time()[3]
res <- mclapply(seq_len(n_jobs), run_job, mc.cores = CORES, mc.preschedule = FALSE)
ok  <- vapply(res, is.data.frame, logical(1))
if (any(!ok)) logmsg("WARNING: %d job(s) failed.", sum(!ok))
raw <- dplyr::bind_rows(res[ok])
write.csv2(raw, "design_variability_raw.csv", row.names = FALSE)
logmsg("Calibrations done in %.1f min (%d rows).",
       (proc.time()[3] - t_start) / 60, nrow(raw))

## ------------------------------------------------------------
## Across-draw spread of c(k), per (cell, arm, k), then decomposed
## ------------------------------------------------------------
by_arm <- raw %>%
  group_by(cell_id, n, p, rho, arm, k) %>%
  summarise(mean_c = mean(c), sd_c = sd(c),
            min_c = min(c), max_c = max(c), .groups = "drop") %>%
  mutate(cv_pct    = 100 * sd_c / mean_c,
         range_pct = 100 * (max_c - min_c) / mean_c)

## Var_design = Var_total - Var_mc  (independent sources; floored at 0)
by_k <- by_arm %>%
  select(cell_id, n, p, rho, k, arm, mean_c, sd_c, cv_pct, range_pct) %>%
  tidyr::pivot_wider(names_from = arm,
                     values_from = c(mean_c, sd_c, cv_pct, range_pct)) %>%
  mutate(
    ## total (design arm) and Monte-Carlo-only components
    cv_total_pct  = cv_pct_design,
    cv_mc_pct     = cv_pct_mc,
    ## net design-to-design variability, on the CV scale
    cv_design_pct = 100 * sqrt(pmax(0, sd_c_design^2 - sd_c_mc^2)) / mean_c_design,
    ## share of the observed spread attributable to Monte Carlo error
    mc_share_pct  = 100 * pmin(1, (sd_c_mc^2) / pmax(sd_c_design^2, .Machine$double.eps))
  )

## Per cell: worst case over family sizes k
by_cell <- by_k %>%
  group_by(cell_id, n, p, rho) %>%
  summarise(df_res         = first(n) - first(p) - 1,
            med_cv_total   = median(cv_total_pct),
            med_cv_mc      = median(cv_mc_pct),
            med_cv_design  = median(cv_design_pct),
            max_cv_design  = max(cv_design_pct),
            max_range_pct  = max(range_pct_design),
            .groups = "drop") %>%
  arrange(desc(max_cv_design))
write.csv2(by_cell, "design_variability_cells.csv", row.names = FALSE)

## ------------------------------------------------------------
## Sanity check: does the mean MC c(k) match the closed-form
## surface at rho = 0? (the surface is fitted at rho = 0)
## ------------------------------------------------------------
surf <- by_k %>%
  filter(rho == 0) %>%
  rowwise() %>%
  mutate(c_hat = get_critical_value(m = k, n = n, p = p)) %>%
  ungroup() %>%
  mutate(rel_err_pct = 100 * (mean_c_design - c_hat) / c_hat)

## ------------------------------------------------------------
## Report
## ------------------------------------------------------------
logmsg("")
logmsg("=== VARIANCE DECOMPOSITION OF THE OBSERVED SPREAD IN c(k) ===")
logmsg("  (pooled over all cells and family sizes k; CV = SD / mean, in %%)")
logmsg("  total   (fresh X each draw) : median %.3f  p95 %.3f  max %.3f",
       median(by_k$cv_total_pct), quantile(by_k$cv_total_pct, .95), max(by_k$cv_total_pct))
logmsg("  MC only (X held fixed)      : median %.3f  p95 %.3f  max %.3f",
       median(by_k$cv_mc_pct), quantile(by_k$cv_mc_pct, .95), max(by_k$cv_mc_pct))
logmsg("  DESIGN net of MC            : median %.3f  p95 %.3f  max %.3f",
       median(by_k$cv_design_pct), quantile(by_k$cv_design_pct, .95), max(by_k$cv_design_pct))
logmsg("  -> Monte Carlo error accounts for a median of %.1f%% of the observed variance.",
       median(by_k$mc_share_pct))
logmsg("")
logmsg("=== BY SAMPLE SIZE (net design CV%% of c(k)) ===")
bn <- by_k %>% group_by(n) %>%
  summarise(med = median(cv_design_pct), mx = max(cv_design_pct), .groups = "drop")
for (i in seq_len(nrow(bn)))
  logmsg("  n = %5d :  median design CV = %.3f%%   max = %.3f%%", bn$n[i], bn$med[i], bn$mx[i])
logmsg("")
logmsg("=== BY CORRELATION (net design CV%% of c(k)) ===")
br <- by_k %>% group_by(rho) %>%
  summarise(med = median(cv_design_pct), mx = max(cv_design_pct), .groups = "drop")
for (i in seq_len(nrow(br)))
  logmsg("  rho = %.1f :  median design CV = %.3f%%   max = %.3f%%", br$rho[i], br$med[i], br$mx[i])
logmsg("")
logmsg("=== TEN MOST DESIGN-SENSITIVE CELLS ===")
h <- head(as.data.frame(by_cell), 10)
for (i in seq_len(nrow(h)))
  logmsg("  n=%5d p=%3d rho=%.1f (df=%4d) : design CV median %.3f%% max %.3f%%  (total %.3f%%, MC %.3f%%)",
         h$n[i], h$p[i], h$rho[i], h$df_res[i],
         h$med_cv_design[i], h$max_cv_design[i], h$med_cv_total[i], h$med_cv_mc[i])
logmsg("")
logmsg("=== CLOSED-FORM SURFACE vs MEAN MC VALUE (rho = 0) ===")
logmsg("  relative error (%%):  median %.3f  p95 %.3f  max %.3f",
       median(surf$rel_err_pct), quantile(abs(surf$rel_err_pct), .95), max(abs(surf$rel_err_pct)))
logmsg("")
logmsg("Wrote design_variability_raw.csv and design_variability_cells.csv")
logmsg("Total elapsed: %.1f min", (proc.time()[3] - t_start) / 60)
