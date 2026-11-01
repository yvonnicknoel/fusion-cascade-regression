# Fusion Cascade for Regression

**A sequentially rejective procedure for predictor selection in linear regression, with Monte-Carlo-calibrated family-wise control of spurious selection and a closed-form approximation of its critical values.**

This repository accompanies the article:

> Noel, Y. (2026). *A Fusion Cascade for predictor selection in linear regression.* **Psychological Methods.** (Accepted 1 November 2026.)

It contains the reference R implementation (`cascade_reg.R`, `get_critical_value.R`), the Monte Carlo
calibration table and the fitted 32-parameter approximation of the critical values, and every script
needed to reproduce the simulations, the robustness checks, the SAPA application and the figures
reported in the paper. The companion repository for the ANOVA version of the method is
[`fusion-cascade`](https://github.com/yvonnicknoel/fusion-cascade).

## Using the method on your own data

```r
source("get_critical_value.R"); source("cascade_reg.R")
c_vec  <- get_critical_value(m = 1:ncol(X), n = nrow(X), p = ncol(X), alpha = 0.05)
result <- fusion_cascade_reg(X, y, alpha = 0.05, c_vec = c_vec)
result$active_predictors
```

`X` is the numeric predictor matrix (observations in rows), `y` the response. The closed-form
critical values are valid for 50 ≤ n ≤ 10 000 and 5 ≤ p ≤ 500 with p < n at α = 0.05; outside that range,
or at another level, calibrate on your own design with `calibrate_cvec_mc_fast()` from
`calibrate_c_fast.R` (seconds for most designs, about a minute for n = 10 000, p = 500).

## Licence

MIT (see `LICENSE`).

---

# Reproduction package

All scripts are plain R and use paths relative to this folder, so run them with this
directory as the working directory. Each script's header documents its inputs and
outputs.

## Quick start

```bash
Rscript run_all.R --list      # show every stage
Rscript run_all.R             # reproduce everything
STAGES=figures Rscript run_all.R
STAGES=revision Rscript run_all.R   # what was new in the first revision
STAGES=revision2 Rscript run_all.R  # what is new in the second revision
```

`run_all.R` is the integrated entry point. It runs each stage in its own R
subprocess — the individual scripts share global variable names and were written to be
run that way, so isolating them guarantees the orchestrator changes no result. Every
script here is byte-identical to the one that produced the reported numbers.

Settings are environment variables: `CORES` (default 90), `NREP` (200), `B` (8000),
`R_DRAWS` (20), `FORCE` (0). A smaller pass:

```bash
CORES=8 NREP=50 B=2000 Rscript run_all.R
```

Two stages are **optional by default** (`ctable`, `fitmodel`): they regenerate the
critical-value surface from scratch, which is the single most expensive part of the
project. Their outputs ship with this folder (`c_table_ext_full.csv`, and the fitted
coefficients hardcoded in `get_critical_value.R`), so they are skipped unless you set
`FORCE=1` or name them explicitly.

## Which script supports which claim

Every file here is on the reproduction path: each is either sourced by another script or wired
into a `run_all.R` stage. If you only want to check one result, this is where to look.

| Claim in the paper | Script | Stage |
|---|---|---|
| The algorithm itself | `cascade_reg.R` | — (library) |
| Applying the method to new data | `get_critical_value.R` | — (library) |
| Figure 1 — approximation vs Monte Carlo | `plot_validation_I4x.R` | `fig1` |
| Approximation accuracy (R² = 0.99999, max 3.2%) | `fit_final_model.R` | `fitmodel` |
| Out-of-grid validation (median 0.32%, max 2.0%) | `validate_surface_oog.R` | `oog` |
| Figure 2 — keep-thresholds of each rule | `plot_decision_mechanics.R` | `fig2` |
| Sensitivity to the realized design (93.9% is MC error) | `check_design_variability.R` | `designvar` |
| Table 2 — overall performance | `main_simulation.R` | `sim` |
| Global-null control under the cascade rule | `check_global_null_fwer.R` | `globalnull` |
| Robustness beyond AR(1); the equicorrelation limit | `check_correlation_robustness.R` | `corrrobust` |
| Negatively correlated predictors (sign identity; non-reducible designs) | `check_negative_correlation.R` | `negcorr` |
| Cost of a per-design calibration (Discussion timings) | `time_calibration.R` | `timing` |
| Figure 3 — behaviour under collinearity | `plot_sim_vs_n_rho07.R` | `fig3` |
| Figure 4 and the α-sweep | `alpha_sweep.R`, `plot_alpha_frontier.R` | `alphasweep`, `fig4` |
| Factor-level results cited in the Discussion | `grouped_demo.R` | `grouped` |
| Figure 5 and the SAPA application | `spi_cv.R`, `plot_spi_operating_point.R` | `spi`, `fig5` |
| Nesting of the selected sets across methods (application and simulation cell) | `overlap_sets.R` | `overlap` |

To apply the method without reproducing anything, three files suffice:
`get_critical_value.R` (self-contained — its coefficients are hardcoded, and the `source` call in
its file is inside a demo block that does not run on load), `cascade_reg.R`, and
`calibrate_c_reg.R`, which `cascade_reg.R` sources at load time.

## Requirements

R (≥ 4.0) with: `MASS`, `glmnet`, `ncvreg`, `grpreg`, `dplyr`, `tidyr`, `parallel`,
`ggplot2`, `ggrepel`, `patchwork`. The SAPA example additionally uses `psychTools`
(which ships the `spi` data set).

```r
install.packages(c("MASS","glmnet","ncvreg","grpreg","dplyr","tidyr",
                   "ggplot2","ggrepel","patchwork","psychTools"))
```

`run_all.R` checks for these up front and stops with an actionable message if any are
absent.

## 1. Method: algorithm and transportable critical values

| File | Purpose |
|------|---------|
| `cascade_reg.R` | The Fusion Cascade algorithm: `fusion_cascade_reg(X, y, alpha, c_vec / ccrit)`. Sequential backward elimination of the smallest \|t\|, halting when the smallest surviving \|t\| exceeds its step-down critical value. |
| `get_critical_value.R` | **Production entry point.** `get_critical_value(m, n, p, alpha = 0.05)` evaluates the closed-form 32-parameter approximation (Model I4x) of `c(k; n, p)`. Self-contained (coefficients hardcoded); the single source of truth for applying the method. |
| `calibrate_c_reg.R` | Reference Monte Carlo calibration of `c_vec` under the global null (`calibrate_cvec_mc`). |
| `calibrate_c_fast.R` | Vectorised, QR-based rewrite (`calibrate_cvec_mc_fast`) — mathematically identical to the reference (reproduces `lm` t-statistics to 1e-15), fast enough for n up to 10 000 and p up to 500. Run it directly to execute its self-test. |
| `tabulate_c_extended.R` | Drives the Monte Carlo calibration over the extended design grid (parallel) to produce the calibration table the approximation is fitted on. |
| `fit_final_model.R` | Fits Model I4x (8 basis functions in (n, p), including the residual-df term 1/(n−p−1), × 4 powers of log k) by OLS on the calibration grid. |
| `coef_block_I4x.R` | The 32 fitted coefficients as an R literal (auto-generated); provenance for the values embedded in `get_critical_value.R`. |

Minimal use on new data:

```r
source("get_critical_value.R"); source("cascade_reg.R")
c_vec  <- get_critical_value(m = 1:ncol(X), n = nrow(X), p = ncol(X), alpha = 0.05)
result <- fusion_cascade_reg(X, y, alpha = 0.05, c_vec = c_vec)
result$active_predictors
```

## 2. Validation and robustness checks

| File | Purpose |
|------|---------|
| `validate_surface_oog.R` | Out-of-grid validation: fresh Monte Carlo at twelve held-out designs lying strictly between grid nodes, confirming the surface interpolates rather than overfits. |
| `check_global_null_fwer.R` | Direct check that the cascade's step-down decision rule holds the global-null family-wise selection rate at α. |
| `check_correlation_robustness.R` | Behaviour beyond AR(1): block (subscale), dense factor-model, and equicorrelated designs — including the extreme equicorrelation regime reported in the Limitations. |
| `check_negative_correlation.R` | **New in the second revision.** Negatively correlated predictors: verifies that reversing column signs leaves the calibrated thresholds exactly unchanged (so AR(1) with −ρ equals +ρ), and tests designs that no sign reversal can make positive (blocks at −0.30, a mixed-sign 3-factor model, weak negative equicorrelation): realized global-null selection rate with the closed-form thresholds, and closed form vs per-design Monte Carlo. |
| `time_calibration.R` | **New in the second revision.** Single-core wall-clock time of `calibrate_cvec_mc_fast()` with B = 8000 for n up to 10 000 and p up to 500; the timings quoted in the Discussion. |
| `check_design_variability.R` | **New in this revision.** Design-realization variability of `c(k)`, with Monte Carlo error partialled out by a second arm holding the design fixed. Answers the Associate Editor's fixed-vs-random predictors question. |

## 3. Simulation study

| File | Purpose |
|------|---------|
| `generate_data.R` | Data generation: AR(1) design, sparse coefficient vectors, target-R² scenarios, and selection metrics (TPR, FPR, FDR, F1, …). |
| `method_wrappers.R` | Uniform wrappers for every competitor: Holm, AIC- and BIC-backward selection, LASSO, elastic net, relaxed LASSO, SCAD, plus out-of-sample MSE. |
| `main_simulation.R` | The main study (parallel): all eight methods across the sample-size / sparsity / signal / correlation grid, including the `s = 0` global-null cells. Produces the results behind Table 2. |
| `alpha_sweep.R` | **New in this revision.** Recalibrates `c(k)` at α ∈ {0.01, 0.05, 0.10, 0.20} and re-runs the full grid for Fusion and Holm, tracing the operating-point path and checking that the realized global-null rate tracks the nominal level. The four levels are four quantiles of the same Monte Carlo replicates, so one calibration pass yields the whole sweep. |
| `grouped_demo.R` | Grouped and categorical predictors: group-level cascade against grouped Holm, BIC, and group LASSO. `RHO` toggles correlated factors. Referenced in the Discussion. |

## 4. Application

| File | Purpose |
|------|---------|
| `spi_cv.R` | The reported personality–health application: 50 repeated 70/30 train–test splits for all methods (model size and out-of-sample R²). |
| `overlap_sets.R` | **New in the second revision.** Overlap of the selected predictor sets across methods: re-runs the 50 SAPA splits with the same seeds keeping each method's set (containment of Fusion's items in each competitor's set, and conversely), and a simulation cell (n = 200, p = 20, s = 5, R² = 0.5, ρ ∈ {0, 0.7}) where the predictors each competitor adds beyond the Fusion set are split into true and false positives. |

## 5. Figures

| Script | Manuscript figure | Reads |
|--------|-------------------|-------|
| `plot_validation_I4x.R` | Figure 1 — approximation accuracy vs Monte Carlo | `c_table_ext_full.csv` |
| `plot_decision_mechanics.R` | Figure 2 — effective \|t\| keep-threshold of each rule | `get_critical_value.R` (closed form only) |
| `plot_alpha_frontier.R` | **Figure 4 (new)** — operating space with the α-indexed Fusion path, and calibration of the dial | `simulation_overall_with_aic.csv`, `alpha_sweep_overall.csv`, `alpha_sweep_global_null.csv` |
| `plot_sim_vs_n_rho07.R` | Figure 3 — behaviour under collinearity (ρ = 0.7) vs n | `simulation_by_n_rho_aic.csv` |
| `plot_spi_operating_point.R` | Figure 5 — SAPA operating point (model size vs out-of-sample R²) | `spi_cv_results.rds` |

All figure scripts are black-and-white safe: methods are separable by line type and point
shape, never by colour alone.

## 6. Shipped results

The aggregated outputs the figures depend on are included, so every figure can be
regenerated in seconds without rerunning the simulations:

| File | Contents |
|------|----------|
| `c_table_ext_full.csv` | Monte Carlo calibration table over the extended (n, p, ρ) grid |
| `simulation_overall_with_aic.csv` | Table 2 — overall means by method |
| `simulation_results_reg_with_aic.csv` | Per-scenario summaries |
| `simulation_by_n_aic.csv`, `simulation_by_n_rho_aic.csv` | Slices by sample size and correlation |
| `simulation_global_null_fwer.csv` | Global-null (s = 0) selection rates by method |
| `alpha_sweep_overall.csv`, `alpha_sweep_summary.csv` | α-sweep operating points |
| `alpha_sweep_global_null.csv` | Realized vs nominal α at the global null |
| `alpha_sweep_cvec.csv` | The calibrated thresholds at each α |
| `design_variability_cells.csv` | Per-cell design-variability decomposition |
| `spi_cv_results.rds` | SAPA application results (50 splits) |
| `overlap_spi.csv`, `overlap_sim.csv` | Nesting of the selected sets across methods (application; simulation cell). `overlap_sim_raw.csv` holds the per-replication counts and `overlap_spi_sets.rds` the selected item sets per split |

The full per-replication raw files are not included (the α-sweep raw output alone is
~40 MB); rerun the corresponding stage to regenerate them.

## Reproduction order

`run_all.R` encodes the order below and handles it automatically.

1. `ctable` → `fitmodel` — regenerate the critical-value surface *(optional; outputs shipped)*
2. `oog`, `designvar` — validate the surface
3. `sim`, `globalnull`, `corrrobust`, `negcorr`, `alphasweep`, `grouped` — simulation studies; `timing`
4. `spi`, `overlap` — application and set-overlap analysis
5. `fig1`…`fig5` — figures
