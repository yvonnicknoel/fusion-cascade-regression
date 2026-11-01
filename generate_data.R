## Packages nécessaires
library(MASS)     # mvrnorm
library(glmnet)   # lasso, elastic net
library(ncvreg)   # SCAD
library(dplyr)    # pour le résumé (optionnel, mais pratique)

##-------------------------------------------------
## Génération X ~ N(0, Sigma) avec structure AR(1)
##-------------------------------------------------
gen_X_ar1 <- function(n, p, rho) {
  if (rho == 0) {
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  } else {
    idx <- 1:p
    Sigma <- rho ^ abs(outer(idx, idx, "-"))
    X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  }
  colnames(X) <- paste0("X", 1:p)
  X
}

##-------------------------------------------------
## Génération du vecteur β vrai
## - p : nb total de prédicteurs
## - s : nb de prédicteurs non nuls
## - beta_val : amplitude (fixe ici pour simplifier)
##-------------------------------------------------
gen_beta <- function(p, s, beta_val = 1) {
  beta <- rep(0, p)
  if (s > 0) {
    act <- sample.int(p, size = s, replace = FALSE)
    signs <- sample(c(-1, 1), size = s, replace = TRUE)
    beta[act] <- beta_val * signs
  } else {
    act <- integer(0)
  }
  list(beta = beta, active = act)
}

##-------------------------------------------------
## Génération (X, y) pour un scénario donné
## - R2_target : R² désiré (NA ou NULL si pas de contrainte)
##-------------------------------------------------
gen_data_scenario <- function(n, p, s, rho, R2_target = 0.5, beta_val = 1) {
  X <- gen_X_ar1(n, p, rho)
  binfo <- gen_beta(p, s, beta_val = beta_val)
  beta <- binfo$beta
  active <- binfo$active

  eta <- as.vector(X %*% beta)
  var_eta <- stats::var(eta)

  if (s == 0 || is.null(R2_target) || is.na(R2_target) || var_eta == 0) {
    sigma2 <- 1
  } else {
    sigma2 <- var_eta * (1 - R2_target) / R2_target
  }

  y <- eta + rnorm(n, mean = 0, sd = sqrt(sigma2))

  list(
    X = X, y = y,
    beta = beta,
    active = active,
    sigma2 = sigma2
  )
}

##-------------------------------------------------
## Métriques de sélection
## selected : indices sélectionnés (par ex. 1:p)
## active   : indices réellement non nuls
## p        : nb total de prédicteurs
##-------------------------------------------------
selection_metrics <- function(selected, active, p) {
  selected <- sort(unique(selected))
  active   <- sort(unique(active))

  TP <- length(intersect(selected, active))
  FP <- length(setdiff(selected, active))
  FN <- length(setdiff(active, selected))
  TN <- p - TP - FP - FN

  TPR <- if (length(active) > 0) TP / length(active) else NA_real_
  FDR <- if (length(selected) > 0) FP / length(selected) else NA_real_
  FPR <- if (p - length(active) > 0) FP / (p - length(active)) else NA_real_

  list(
    TP = TP, FP = FP, FN = FN, TN = TN,
    TPR = TPR, FPR = FPR, FDR = FDR,
    model_size = length(selected)
  )
}
