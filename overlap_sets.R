## ============================================================
## Overlap of the selected predictor sets across methods (R2, AE-13).
##
## Question (Associate Editor): do the other methods merely *extend* the
## Fusion set, or do they select different predictors?
##
## Two parts:
##  (A) SAPA application: the 50 train/test splits of spi_cv.R are re-run
##      with the SAME seeds (so the sets are those behind Figure 5), and the
##      selected sets are kept. For each competitor M we report
##        containment  |F ∩ M| / |F|   share of Fusion's items that M also selects
##        coverage     |F ∩ M| / |M|   share of M's items that Fusion also selects
##        Jaccard      |F ∩ M| / |F ∪ M|
##      plus how often Holm's set is nested in Fusion's.
##  (B) Simulation slice with known truth: n = 200, p = 20, s = 5, R² = 0.5,
##      ρ ∈ {0, 0.7}, NREP replications, all eight methods. Same statistics,
##      and the predictors M selects *beyond* Fusion are split into true and
##      false positives (likewise the predictors Fusion keeps that M drops).
##
## Outputs: overlap_spi.csv, overlap_spi_sets.rds, overlap_sim.csv,
##          overlap_sim_raw.csv; progress in overlap_sets.log.
## Env: CORES (50), NREP (200), PART ("spi,sim").
## ============================================================
suppressMessages({ library(psychTools); library(parallel); library(dplyr) })
source("get_critical_value.R"); source("cascade_reg.R")
source("generate_data.R"); source("method_wrappers.R")

CORES <- as.integer(Sys.getenv("CORES", "50"))
NREP  <- as.integer(Sys.getenv("NREP",  "200"))
PART  <- strsplit(Sys.getenv("PART", "spi,sim"), ",")[[1]]
LOG <- "overlap_sets.log"; cat("", file = LOG)
logm <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
t0 <- Sys.time()
elapsed <- function() sprintf("%.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))

methods <- c("Fusion","Holm","AIC","BIC","Lasso","ElasticNet","SCAD","RelaxedLasso")
others  <- setdiff(methods, "Fusion")

overlap_stats <- function(F, M) {
  i <- length(intersect(F, M)); u <- length(union(F, M))
  c(containment = if (length(F) > 0) i / length(F) else NA_real_,
    coverage    = if (length(M) > 0) i / length(M) else NA_real_,
    jaccard     = if (u > 0) i / u else 1,
    nested_in_F = as.numeric(all(M %in% F)),
    F_nested    = as.numeric(all(F %in% M)))
}

## ------------------------------------------------------------
## (A) SAPA
## ------------------------------------------------------------
if ("spi" %in% PART) {
  data(spi)
  demo  <- c("age","sex","health","education","p1edu","p2edu","wellness","exer","smoke","ER")
  items <- setdiff(colnames(spi), demo)
  d <- spi[!is.na(spi$health), ]
  X <- as.matrix(d[, items]); y <- as.numeric(d$health)
  n <- nrow(X); p <- ncol(X)
  R <- 50L
  logm("(A) SAPA: %d items, n = %d, %d splits, %d cores", p, n, R, CORES)

  run_split <- function(r) {                     # identical to spi_cv.R::run_rep
    set.seed(1000 + r)
    tr <- sample(n, floor(0.7 * n))
    Xtr <- X[tr, ]; ytr <- y[tr]
    ntr <- nrow(Xtr)
    cc <- get_critical_value(1:p, ntr, p)
    pick <- function(m) switch(m,
      Fusion = select_fusion(Xtr, ytr, alpha = 0.05, ccrit = cc), Holm = select_holm(Xtr, ytr, 0.05),
      AIC = select_aic(Xtr, ytr), BIC = select_bic(Xtr, ytr), Lasso = select_lasso(Xtr, ytr),
      ElasticNet = select_enet(Xtr, ytr), SCAD = select_scad(Xtr, ytr),
      RelaxedLasso = select_relaxed_lasso(Xtr, ytr))
    sets <- list()
    for (m in methods) { s <- tryCatch(pick(m), error = function(e) integer(0)); sets[[m]] <- s[!is.na(s)] }
    cat(sprintf("  split %2d done (%s)\n", r, elapsed()), file = LOG, append = TRUE)
    sets
  }
  spi_sets <- mclapply(seq_len(R), run_split, mc.cores = min(CORES, R), mc.preschedule = FALSE)
  ok <- vapply(spi_sets, is.list, logical(1)); spi_sets <- spi_sets[ok]
  saveRDS(list(items = items, sets = spi_sets), "overlap_spi_sets.rds")

  rows <- lapply(seq_along(spi_sets), function(r) {
    F <- spi_sets[[r]]$Fusion
    do.call(rbind, lapply(others, function(m) {
      st <- overlap_stats(F, spi_sets[[r]][[m]])
      data.frame(split = r, method = m, size_F = length(F), size_M = length(spi_sets[[r]][[m]]),
                 t(st))
    }))
  })
  raw <- bind_rows(rows)
  spi_summary <- raw %>% group_by(method) %>%
    summarise(size_F = mean(size_F), size_M = mean(size_M),
              containment = mean(containment, na.rm = TRUE),
              coverage = mean(coverage, na.rm = TRUE),
              jaccard = mean(jaccard),
              F_nested_pct = 100 * mean(F_nested), M_nested_in_F_pct = 100 * mean(nested_in_F),
              .groups = "drop") %>% arrange(desc(containment))
  write.csv(spi_summary, "overlap_spi.csv", row.names = FALSE)
  logm("\n=== (A) SAPA: overlap with Fusion over %d splits ===", length(spi_sets))
  logm("%-12s %6s %6s %12s %9s %8s %10s %10s", "method", "|F|", "|M|", "containment", "coverage", "Jaccard", "F⊂M %", "M⊂F %")
  for (i in seq_len(nrow(spi_summary))) with(spi_summary[i, ],
    logm("%-12s %6.1f %6.1f %12.3f %9.3f %8.3f %10.0f %10.0f",
         method, size_F, size_M, containment, coverage, jaccard, F_nested_pct, M_nested_in_F_pct))
}

## ------------------------------------------------------------
## (B) Simulation slice with known truth
## ------------------------------------------------------------
if ("sim" %in% PART) {
  n <- 200; p <- 20; s <- 5; R2 <- 0.5; rho_vec <- c(0, 0.7)
  cc <- get_critical_value(1:p, n, p)
  logm("\n(B) Simulation slice: n=%d p=%d s=%d R2=%.1f rho in {%s}, %d reps each",
       n, p, s, R2, paste(rho_vec, collapse = ","), NREP)

  run_rep <- function(idx) {
    rho <- rho_vec[(idx - 1) %/% NREP + 1]; r <- (idx - 1) %% NREP + 1
    set.seed(777000 + 1000 * which(rho_vec == rho) + r)
    dd <- gen_data_scenario(n, p, s, rho, R2_target = R2, beta_val = 1)
    Xtr <- dd$X; ytr <- dd$y; act <- dd$active
    pick <- function(m) switch(m,
      Fusion = select_fusion(Xtr, ytr, alpha = 0.05, ccrit = cc), Holm = select_holm(Xtr, ytr, 0.05),
      AIC = select_aic(Xtr, ytr), BIC = select_bic(Xtr, ytr), Lasso = select_lasso(Xtr, ytr),
      ElasticNet = select_enet(Xtr, ytr), SCAD = select_scad(Xtr, ytr),
      RelaxedLasso = select_relaxed_lasso(Xtr, ytr))
    sets <- list()
    for (m in methods) { z <- tryCatch(pick(m), error = function(e) integer(0)); sets[[m]] <- z[!is.na(z)] }
    F <- sets$Fusion
    out <- do.call(rbind, lapply(others, function(m) {
      M <- sets[[m]]; st <- overlap_stats(F, M)
      extra <- setdiff(M, F); missing <- setdiff(F, M)
      data.frame(rho = rho, rep = r, method = m, size_F = length(F), size_M = length(M), t(st),
                 extra_TP = sum(extra %in% act), extra_FP = sum(!(extra %in% act)),
                 missing_TP = sum(missing %in% act), missing_FP = sum(!(missing %in% act)))
    }))
    if (r %% 20 == 0) cat(sprintf("  rho=%.1f rep %3d done (%s)\n", rho, r, elapsed()), file = LOG, append = TRUE)
    out
  }
  sim_raw <- bind_rows(mclapply(seq_len(NREP * length(rho_vec)), run_rep, mc.cores = CORES, mc.preschedule = FALSE))
  write.csv(sim_raw, "overlap_sim_raw.csv", row.names = FALSE)
  sim_summary <- sim_raw %>% group_by(rho, method) %>%
    summarise(size_F = mean(size_F), size_M = mean(size_M),
              containment = mean(containment, na.rm = TRUE), coverage = mean(coverage, na.rm = TRUE),
              jaccard = mean(jaccard),
              F_nested_pct = 100 * mean(F_nested), M_nested_in_F_pct = 100 * mean(nested_in_F),
              ## share computed from the raw counts BEFORE the means below shadow the column names
              extra_FP_share = sum(extra_FP) / pmax(1, sum(extra_TP + extra_FP)),
              missing_TP_share = sum(missing_TP) / pmax(1, sum(missing_TP + missing_FP)),
              extra_TP = mean(extra_TP), extra_FP = mean(extra_FP),
              missing_TP = mean(missing_TP), missing_FP = mean(missing_FP),
              .groups = "drop") %>% arrange(rho, desc(containment))
  write.csv(sim_summary, "overlap_sim.csv", row.names = FALSE)
  for (rh in rho_vec) {
    logm("\n=== (B) rho = %.1f: overlap with Fusion (mean over %d reps) ===", rh, NREP)
    logm("%-12s %5s %5s %11s %8s %7s %6s %6s %8s %8s %8s %9s", "method", "|F|", "|M|", "containment", "coverage",
         "Jacc", "F⊂M%", "M⊂F%", "extraTP", "extraFP", "FPshare", "missTP/FP")
    sub <- sim_summary[sim_summary$rho == rh, ]
    for (i in seq_len(nrow(sub))) with(sub[i, ],
      logm("%-12s %5.2f %5.2f %11.3f %8.3f %7.3f %6.0f %6.0f %8.2f %8.2f %8.2f %4.2f/%4.2f",
           method, size_F, size_M, containment, coverage, jaccard, F_nested_pct, M_nested_in_F_pct,
           extra_TP, extra_FP, extra_FP_share, missing_TP, missing_FP))
  }
}
logm("\nDone in %s.", elapsed())
