source("calibrate_c_reg.R")
source("get_critical_value.R")
source("cascade_reg.R")
source("generate_data.R")
source("method_wrappers.R")
suppressMessages({ library(parallel); library(dplyr); library(tidyr) })

## ============================================================
## Grille de scénarios
## ============================================================
n_vec   <- c(50, 100, 200, 500, 1000, 2000)
p_vec   <- c(10, 20, 50)
s_vec   <- c(0, 2, 5, 10)   # nb de prédicteurs non nuls (s=0 = null global : contrôle de garantie FWER)
rho_vec <- c(0, 0.3, 0.7)   # corrélations
R2_vec  <- c(0.2, 0.5)      # intensité du signal

scen_grid <- expand.grid(n = n_vec, p = p_vec, s = s_vec,
                         rho = rho_vec, R2 = R2_vec)
## garder p <= n/2 (les compétiteurs basés sur le modèle complet, Holm/AIC/BIC,
## exigent n nettement > p) et s < p (il faut des prédicteurs nuls)
scen_grid <- subset(scen_grid, p <= n / 2 & s < p)
## s = 0 (null global) : R2 est sans objet -> ne garder qu'une seule valeur de R2 (évite les doublons)
scen_grid <- subset(scen_grid, s > 0 | R2 == R2_vec[1])
scen_grid <- scen_grid[order(scen_grid$p, scen_grid$n), ]
rownames(scen_grid) <- NULL
n_scen <- nrow(scen_grid)

n_rep  <- as.integer(Sys.getenv("NREP",  "200"))   # réplications par scénario
CORES  <- as.integer(Sys.getenv("CORES", "90"))
n_test <- 1000                                       # échantillon test pour la MSE
SEED   <- 1234
LOG    <- "sim_run.log"

logmsg <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }

## Les 8 méthodes comparées (RelaxedLasso ajouté EN DERNIER pour ne pas
## perturber le flux aléatoire -- donc les nombres des 7 autres sont inchangés)
method_names <- c("Fusion", "Holm", "AIC", "BIC", "Lasso", "ElasticNet", "SCAD", "RelaxedLasso")

## ============================================================
## Une réplication : renvoie un data.frame (une ligne par méthode)
## ============================================================
run_rep <- function(r, s_id, n, p, s, rho, R2, ccrit_s) {
  set.seed(SEED + 100000L * s_id + r)

  ## 1) Données train
  dtrain    <- gen_data_scenario(n, p, s, rho, R2_target = R2, beta_val = 1)
  X_train   <- dtrain$X; y_train <- dtrain$y
  beta_true <- dtrain$beta; active <- dtrain$active

  ## 2) Données test (même β, même structure de corrélation)
  X_test   <- gen_X_ar1(n_test, p, rho)
  eta_test <- as.vector(X_test %*% beta_true)
  var_eta  <- stats::var(eta_test)
  sigma2_test <- if (s == 0 || is.null(R2) || is.na(R2) || var_eta == 0) 1 else var_eta * (1 - R2) / R2
  y_test   <- eta_test + rnorm(n_test, 0, sqrt(sigma2_test))

  pick <- function(meth) switch(meth,
    "Fusion"     = select_fusion(X_train, y_train, alpha = 0.05, ccrit = ccrit_s),
    "Holm"       = select_holm(X_train, y_train, alpha = 0.05),
    "AIC"        = select_aic(X_train, y_train),
    "BIC"        = select_bic(X_train, y_train),
    "Lasso"      = select_lasso(X_train, y_train),
    "ElasticNet" = select_enet(X_train, y_train),
    "SCAD"       = select_scad(X_train, y_train),
    "RelaxedLasso" = select_relaxed_lasso(X_train, y_train))

  rows <- lapply(method_names, function(meth) {
    out <- tryCatch({
      selected <- pick(meth)
      mets <- selection_metrics(selected, active, p)
      mse  <- prediction_mse(X_train, y_train, X_test, y_test, selected)
      data.frame(TP = mets$TP, FP = mets$FP, FN = mets$FN, TN = mets$TN,
                 TPR = mets$TPR, FPR = mets$FPR, FDR = mets$FDR,
                 model_size = mets$model_size, MSE = mse)
    }, error = function(e) data.frame(TP=NA,FP=NA,FN=NA,TN=NA,TPR=NA,FPR=NA,
                                      FDR=NA,model_size=NA,MSE=NA))
    data.frame(scen_id = s_id, rep = r, method = meth,
               n = n, p = p, s = s, rho = rho, R2 = R2, out)
  })
  do.call(rbind, rows)
}

## ============================================================
## Boucle sur les scénarios ; réplications en parallèle
## ============================================================
cat("", file = LOG)
logmsg("Simulation: %d scénarios x %d réplications x %d méthodes | cores=%d",
       n_scen, n_rep, length(method_names), CORES)
t_start <- proc.time()[3]
results_list <- vector("list", n_scen)

for (s_id in seq_len(n_scen)) {
  n <- scen_grid$n[s_id]; p <- scen_grid$p[s_id]; s <- scen_grid$s[s_id]
  rho <- scen_grid$rho[s_id]; R2 <- scen_grid$R2[s_id]
  ccrit_s <- get_critical_value(m = 1:p, n = n, p = p)   # seuils Fusion (Model I4x)

  t0 <- proc.time()[3]
  reps <- mclapply(seq_len(n_rep),
                   function(r) run_rep(r, s_id, n, p, s, rho, R2, ccrit_s),
                   mc.cores = CORES, mc.preschedule = TRUE)
  ok <- vapply(reps, is.data.frame, logical(1))
  results_list[[s_id]] <- do.call(rbind, reps[ok])
  logmsg("[%2d/%2d] n=%3d p=%2d s=%d rho=%.1f R2=%.1f  %5.1fs  (%d/%d reps ok)",
         s_id, n_scen, n, p, s, rho, R2, proc.time()[3] - t0, sum(ok), n_rep)
}

sim_results <- dplyr::bind_rows(results_list)
## F1 = 2TP / (2TP + FP + FN); precision = TP/(TP+FP) = 1 - FDR
sim_results$F1  <- with(sim_results, 2*TP / (2*TP + FP + FN))
write.csv2(sim_results, file = "simulation_results_reg_raw_aic.csv", row.names = FALSE)

## ---- résumé par scénario ----
summary_results <- sim_results %>%
  group_by(method, n, p, s, rho, R2) %>%
  summarise(mean_TPR = mean(TPR, na.rm = TRUE), mean_FPR = mean(FPR, na.rm = TRUE),
            mean_FDR = mean(FDR, na.rm = TRUE), mean_F1 = mean(F1, na.rm = TRUE),
            mean_size = mean(model_size, na.rm = TRUE), mean_MSE = mean(MSE, na.rm = TRUE),
            FWER_FN = mean(FN > 0, na.rm = TRUE), FWER_FP = mean(FP > 0, na.rm = TRUE),
            .groups = "drop")
write.csv2(summary_results, file = "simulation_results_reg_with_aic.csv", row.names = FALSE)

## ---- moyennes globales (analogue Table 1) ----
overall <- sim_results %>%
  filter(s > 0) %>%                     # tables/figures de comparaison : scénarios avec signal
  group_by(method) %>%
  summarise(TPR = mean(TPR, na.rm = TRUE), FPR = mean(FPR, na.rm = TRUE),
            FDR = mean(FDR, na.rm = TRUE), F1 = mean(F1, na.rm = TRUE),
            model_size = mean(model_size, na.rm = TRUE),
            FWER_FN = mean(FN > 0, na.rm = TRUE), FWER_FP = mean(FP > 0, na.rm = TRUE),
            MSE = mean(MSE, na.rm = TRUE), .groups = "drop") %>%
  arrange(method)
write.csv2(overall, file = "simulation_overall_with_aic.csv", row.names = FALSE)

## ---- tranche par (méthode, n) : effet de la taille d'échantillon ----
by_n <- sim_results %>%
  filter(s > 0) %>%
  group_by(method, n) %>%
  summarise(TPR = mean(TPR, na.rm = TRUE), FPR = mean(FPR, na.rm = TRUE),
            FDR = mean(FDR, na.rm = TRUE), model_size = mean(model_size, na.rm = TRUE),
            FWER_FP = mean(FP > 0, na.rm = TRUE), MSE = mean(MSE, na.rm = TRUE),
            .groups = "drop")
write.csv2(by_n, file = "simulation_by_n_aic.csv", row.names = FALSE)

## ---- tranche par (méthode, n, rho) : effet de la colinéarité selon n ----
by_n_rho <- sim_results %>%
  filter(s > 0) %>%
  group_by(method, n, rho) %>%
  summarise(TPR = mean(TPR, na.rm = TRUE), FPR = mean(FPR, na.rm = TRUE),
            FDR = mean(FDR, na.rm = TRUE), model_size = mean(model_size, na.rm = TRUE),
            FWER_FP = mean(FP > 0, na.rm = TRUE), MSE = mean(MSE, na.rm = TRUE),
            .groups = "drop")
write.csv2(by_n_rho, file = "simulation_by_n_rho_aic.csv", row.names = FALSE)

## ---- CONTRÔLE DE GARANTIE : null global (s = 0) ----
## FWER = P(sélectionner au moins un prédicteur | tous les beta = 0). Doit être <= alpha
## pour Fusion (et Holm) ; documente aussi l'absence de contrôle des méthodes pénalisées.
global_null <- sim_results %>%
  filter(s == 0) %>%
  group_by(method, n, p, rho) %>%
  summarise(FWER = mean(FP > 0, na.rm = TRUE), mean_size = mean(model_size, na.rm = TRUE),
            .groups = "drop")
write.csv2(global_null, file = "simulation_global_null_fwer.csv", row.names = FALSE)
gnf <- global_null$FWER[global_null$method == "Fusion"]
logmsg("GLOBAL-NULL (s=0) Fusion FWER : min=%.3f max=%.3f mean=%.3f  (cible <= 0.05)",
       min(gnf), max(gnf), mean(gnf))

logmsg("Done in %.1f min. Rows: %d. Wrote raw / with_aic / overall / by_n / by_n_rho / global_null.",
       (proc.time()[3] - t_start) / 60, nrow(sim_results))
cat("\n=== Overall means across all scenarios ===\n")
print(as.data.frame(overall), digits = 3)
cat("\n=== FPR by sample size (effect of n) ===\n")
print(as.data.frame(tidyr::pivot_wider(by_n[c("method","n","FPR")],
      names_from = n, values_from = FPR)), digits = 2)
cat("\n=== AIC under collinearity (rho=0.7): FPR by n ===\n")
print(as.data.frame(subset(by_n_rho, rho == 0.7 & method %in% c("AIC","BIC","Fusion"))[
      c("method","n","TPR","FPR","FWER_FP")]), digits = 2)
