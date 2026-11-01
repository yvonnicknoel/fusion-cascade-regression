## ============================================================
## run_all.R -- integrated reproduction script
##
## "A Fusion Cascade for Predictor Selection in Linear Regression"
##
## Reproduces every analysis, table and figure in the manuscript,
## end to end, from this directory.
##
##   Rscript run_all.R                 # everything
##   STAGES=sim,figures Rscript run_all.R
##   CORES=8 NREP=50 Rscript run_all.R # smaller/faster pass
##   Rscript run_all.R --list          # show stages and exit
##   Rscript run_all.R --dry-run       # show what would run
##
## Each stage runs in its OWN R subprocess. The individual scripts
## were written to be run that way, they share global variable names,
## and isolating them keeps this orchestrator from changing any
## result: every script here is byte-identical to the one that
## produced the published numbers.
##
## Environment variables (all optional):
##   STAGES   comma-separated stage names (default: all enabled ones)
##   CORES    parallel workers                        (default 90)
##   NREP     replications per simulation scenario    (default 200)
##   B        Monte Carlo replicates for calibration  (default 8000)
##   R_DRAWS  design draws for the variability check  (default 20)
##   FORCE    "1" to also run stages marked optional  (default off)
##
## A full default run is heavy (the calibration grid and the main
## simulation dominate). The two stages that are optional by default
## -- 'ctable' and 'fitmodel' -- regenerate the critical-value
## surface from scratch; their outputs ship with this folder
## (c_table_ext_full.csv, and the coefficients hardcoded in
## get_critical_value.R), so they are skipped unless FORCE=1.
## ============================================================

CORES   <- Sys.getenv("CORES",   "90")
NREP    <- Sys.getenv("NREP",    "200")
B_MC    <- Sys.getenv("B",       "8000")
R_DRAWS <- Sys.getenv("R_DRAWS", "20")
FORCE   <- Sys.getenv("FORCE",   "0") == "1"

LOG <- "run_all.log"

## ------------------------------------------------------------
## Stage table
##   name      short selector
##   script    file to run
##   optional  skipped unless FORCE=1 (outputs already shipped)
##   outputs   files the stage is expected to produce
## ------------------------------------------------------------
STAGE_TABLE <- list(
  list(name = "ctable", script = "tabulate_c_extended.R", optional = TRUE,
       what = "Monte Carlo calibration of c(k) over the extended (n,p,rho) grid",
       outputs = "c_table_ext_full.csv"),

  list(name = "fitmodel", script = "fit_final_model.R", optional = TRUE,
       what = "Fit the 32-parameter closed-form surface (Model I4x) to the grid",
       outputs = c("model_I4x_final.rds", "coefficients_model_I4x.csv", "coef_block_I4x.R")),

  list(name = "oog", script = "validate_surface_oog.R", optional = FALSE,
       what = "Out-of-grid validation: surface vs fresh MC at held-out designs",
       outputs = "validation_oog.rds"),

  list(name = "designvar", script = "check_design_variability.R", optional = FALSE,
       what = "REVISION -- design-realization variability of c(k), MC error partialled out",
       outputs = c("design_variability_raw.csv", "design_variability_cells.csv")),

  list(name = "sim", script = "main_simulation.R", optional = FALSE,
       what = "Main simulation study: 8 methods over the full scenario grid (Table 2)",
       outputs = c("simulation_results_reg_raw_aic.csv", "simulation_overall_with_aic.csv",
                   "simulation_by_n_rho_aic.csv", "simulation_global_null_fwer.csv")),

  list(name = "globalnull", script = "check_global_null_fwer.R", optional = FALSE,
       what = "Global-null FWER check of the cascade decision rule",
       outputs = character(0)),

  list(name = "corrrobust", script = "check_correlation_robustness.R", optional = FALSE,
       what = "Robustness beyond AR(1): block, factor and equicorrelated designs",
       outputs = character(0)),

  list(name = "negcorr", script = "check_negative_correlation.R", optional = FALSE,
       what = "REVISION 2 -- negatively correlated predictors: sign identity and non-reducible designs",
       outputs = c("negative_correlation_check.csv", "negative_correlation_check.log")),

  list(name = "timing", script = "time_calibration.R", optional = FALSE,
       what = "REVISION 2 -- wall-clock cost of a per-design Monte Carlo calibration",
       outputs = "time_calibration.log"),

  list(name = "alphasweep", script = "alpha_sweep.R", optional = FALSE,
       what = "REVISION -- alpha-sweep of the Fusion (and Holm) operating point",
       outputs = c("alpha_sweep_overall.csv", "alpha_sweep_global_null.csv",
                   "alpha_sweep_summary.csv", "alpha_sweep_cvec.csv")),

  list(name = "grouped", script = "grouped_demo.R", optional = FALSE,
       what = "Grouped / categorical predictors (referenced in the Discussion)",
       outputs = character(0)),

  list(name = "spi", script = "spi_cv.R", optional = FALSE,
       what = "SAPA personality-health application: 50 random 70/30 splits",
       outputs = "spi_cv_results.rds"),

  list(name = "overlap", script = "overlap_sets.R", optional = FALSE,
       what = "REVISION 2 -- overlap/nesting of the selected sets across methods (SAPA splits + simulation cell)",
       outputs = c("overlap_spi.csv", "overlap_sim.csv", "overlap_sim_raw.csv", "overlap_spi_sets.rds")),

  list(name = "fig1", script = "plot_validation_I4x.R", optional = FALSE,
       what = "Figure 1 -- approximation accuracy vs Monte Carlo",
       outputs = "plot_validation_I4x.pdf"),

  list(name = "fig2", script = "plot_decision_mechanics.R", optional = FALSE,
       what = "Figure 2 -- effective |t| keep-threshold of each selection rule",
       outputs = "plot_decision_mechanics.pdf"),

  list(name = "fig3", script = "plot_sim_vs_n_rho07.R", optional = FALSE,
       what = "Figure 3 -- behaviour under collinearity (rho = 0.7) vs n",
       outputs = "plot_sim_vs_n_rho07.pdf"),

  list(name = "fig4", script = "plot_alpha_frontier.R", optional = FALSE,
       what = "REVISION -- Figure 4: operating space and calibration of the alpha dial",
       outputs = "plot_alpha_frontier.pdf"),

  list(name = "fig5", script = "plot_spi_operating_point.R", optional = FALSE,
       what = "Figure 5 -- SAPA operating point (model size vs out-of-sample R^2)",
       outputs = "plot_spi_operating_point.pdf")
)

## Convenience groups
GROUPS <- list(
  calibration = c("ctable", "fitmodel", "oog", "designvar"),
  sim         = c("sim", "globalnull", "corrrobust", "negcorr", "alphasweep", "grouped"),
  application = c("spi", "overlap"),
  figures     = c("fig1", "fig2", "fig3", "fig4", "fig5"),
  revision    = c("designvar", "alphasweep", "fig4"),  # new in the first revision
  revision2   = c("overlap", "negcorr", "timing")      # new in the second revision
)

stage_names <- vapply(STAGE_TABLE, `[[`, character(1), "name")

## ------------------------------------------------------------
## CLI
## ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if ("--list" %in% args) {
  cat("\nStages (in run order):\n\n")
  for (s in STAGE_TABLE)
    cat(sprintf("  %-11s %-9s %s\n", s$name,
                if (s$optional) "[opt]" else "", s$what))
  cat("\nGroups: ", paste(names(GROUPS), collapse = ", "), "\n")
  cat("Example: STAGES=revision Rscript run_all.R\n\n")
  quit(status = 0)
}
DRY <- "--dry-run" %in% args

sel <- Sys.getenv("STAGES", "")
if (nzchar(sel)) {
  want <- trimws(strsplit(sel, ",")[[1]])
  expanded <- unlist(lapply(want, function(w)
    if (!is.null(GROUPS[[w]])) GROUPS[[w]] else w), use.names = FALSE)
  bad <- setdiff(expanded, stage_names)
  if (length(bad))
    stop("Unknown stage(s): ", paste(bad, collapse = ", "),
         "\nRun 'Rscript run_all.R --list' to see the available stages.")
  ## keep canonical run order
  chosen <- stage_names[stage_names %in% expanded]
  explicit <- TRUE
} else {
  chosen <- stage_names
  explicit <- FALSE
}

## ------------------------------------------------------------
## Logging
## ------------------------------------------------------------
logmsg <- function(...) {
  m <- sprintf(...)
  cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE); flush.console()
}
hr <- function() logmsg(strrep("-", 72))

## ------------------------------------------------------------
## Package check (fail early, with an actionable message)
## ------------------------------------------------------------
REQUIRED <- c("MASS", "glmnet", "ncvreg", "grpreg", "dplyr", "tidyr",
              "parallel", "ggplot2", "ggrepel", "patchwork", "psychTools")
missing <- REQUIRED[!vapply(REQUIRED, requireNamespace, logical(1), quietly = TRUE)]

cat("", file = LOG)
logmsg("============================================================")
logmsg("Fusion Cascade for regression -- full reproduction")
logmsg("started %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
logmsg("R %s | CORES=%s NREP=%s B=%s R_DRAWS=%s FORCE=%s",
       getRversion(), CORES, NREP, B_MC, R_DRAWS, FORCE)
logmsg("============================================================")

if (length(missing)) {
  logmsg("")
  logmsg("MISSING PACKAGES: %s", paste(missing, collapse = ", "))
  logmsg('Install with: install.packages(c(%s))',
         paste0('"', missing, '"', collapse = ", "))
  if (!DRY) stop("Install the packages above, then re-run.")
}

## ------------------------------------------------------------
## Run
## ------------------------------------------------------------
child_env <- c(paste0("CORES=",   CORES),
               paste0("NREP=",    NREP),
               paste0("B=",       B_MC),
               paste0("R_DRAWS=", R_DRAWS))

results <- list()
t_all <- proc.time()[3]

for (nm in chosen) {
  s <- STAGE_TABLE[[match(nm, stage_names)]]

  ## optional stages run only when asked for explicitly, or with FORCE=1
  if (s$optional && !FORCE && !explicit) {
    logmsg("")
    logmsg("SKIP  %-11s %s", s$name, s$what)
    logmsg("      (outputs ship with this folder; set FORCE=1 to regenerate)")
    results[[nm]] <- list(status = "skipped", secs = 0)
    next
  }
  if (!file.exists(s$script)) {
    logmsg("")
    logmsg("ERROR %-11s script not found: %s", s$name, s$script)
    results[[nm]] <- list(status = "missing-script", secs = 0)
    next
  }

  hr()
  logmsg("STAGE %-11s %s", s$name, s$what)
  logmsg("      script: %s", s$script)
  if (DRY) { results[[nm]] <- list(status = "dry-run", secs = 0); next }

  t0 <- proc.time()[3]
  code <- system2("Rscript", shQuote(s$script), env = child_env,
                  stdout = "", stderr = "")
  el <- proc.time()[3] - t0

  produced <- if (length(s$outputs)) file.exists(s$outputs) else logical(0)
  ok <- (code == 0L) && all(produced)

  logmsg("      -> %s in %.1f min (exit %d)",
         if (ok) "OK" else "FAILED", el / 60, code)
  if (length(s$outputs) && !all(produced))
    logmsg("      MISSING OUTPUT: %s", paste(s$outputs[!produced], collapse = ", "))

  results[[nm]] <- list(status = if (ok) "ok" else "failed", secs = el)
}

## ------------------------------------------------------------
## Summary
## ------------------------------------------------------------
hr()
logmsg("")
logmsg("=== SUMMARY ===")
for (nm in names(results)) {
  r <- results[[nm]]
  logmsg("  %-11s %-14s %6.1f min", nm, r$status, r$secs / 60)
}
nfail <- sum(vapply(results, function(r) r$status %in% c("failed", "missing-script"),
                    logical(1)))
logmsg("")
logmsg("Total wall clock: %.1f min", (proc.time()[3] - t_all) / 60)
logmsg("finished %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
if (nfail) {
  logmsg("")
  logmsg("%d stage(s) did not complete cleanly -- see %s and the stage output above.",
         nfail, LOG)
  quit(status = 1)
}
logmsg("All requested stages completed.")
