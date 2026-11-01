##-------------------------------------------------
##  Méthode 1 : Fusion Cascade
##-------------------------------------------------
select_fusion <- function(X, y, alpha = 0.05, ccrit, ...) {
  fc <- fusion_cascade_reg(
    X, y,
    alpha        = alpha,
    ccrit        = ccrit,
    standardize  = TRUE,
    verbose      = FALSE
  )
  # indices 1:p des prédicteurs actifs
  fc$active_indices
}

##-------------------------------------------------
##  Méthode 2 : Holm sur les coefficients du modèle complet
##-------------------------------------------------
select_holm <- function(X, y, alpha = 0.05) {
  df <- data.frame(y = y, X)
  fit <- lm(y ~ ., data = df)
  sfit <- summary(fit)
  # lignes 2:(p+1) dans summary$coef (sans intercept)
  pvals <- sfit$coefficients[-1, 4]
  padj  <- p.adjust(pvals, method = "holm")
  which(padj < alpha)  # indices 1..p
}

##-------------------------------------------------
##  Méthode 3 : Backward BIC
##-------------------------------------------------
select_bic <- function(X, y) {
  df <- data.frame(y = y, X)
  full <- lm(y ~ ., data = df)
  # BIC = AIC avec k = log(n)
  bic_model <- step(full, direction = "backward",
                    k = log(nrow(X)), trace = 0)
  coefs <- names(coef(bic_model))
  # les colonnes Xj sont dans la formule sous ces noms-là
  selected_names <- setdiff(coefs, "(Intercept)")
  match(selected_names, colnames(X))  # indices 1..p
}

##-------------------------------------------------
##  Méthode 3bis : Backward AIC
##-------------------------------------------------
select_aic <- function(X, y) {
  df <- data.frame(y = y, X)
  full <- lm(y ~ ., data = df)
  # AIC : k = 2 (pénalité par défaut de step)
  aic_model <- step(full, direction = "backward",
                    k = 2, trace = 0)
  coefs <- names(coef(aic_model))
  selected_names <- setdiff(coefs, "(Intercept)")
  match(selected_names, colnames(X))  # indices 1..p
}

##-------------------------------------------------
##  Méthode 4 : Lasso (λ.1SE)
##-------------------------------------------------
select_lasso <- function(X, y, alpha_lasso = 1) {
  cvfit <- cv.glmnet(X, y, alpha = alpha_lasso,
                     standardize = TRUE,
                     family = "gaussian")
  beta_hat <- coef(cvfit, s = "lambda.1se")
  # beta_hat : intercept + p coefficients
  idx <- which(as.vector(beta_hat)[-1] != 0)  # indices 1..p
  idx
}

##-------------------------------------------------
##  Méthode 5 : Elastic Net (α = 0.5, λ.1SE)
##-------------------------------------------------
select_enet <- function(X, y, alpha_enet = 0.5) {
  cvfit <- cv.glmnet(X, y, alpha = alpha_enet,
                     standardize = TRUE,
                     family = "gaussian")
  beta_hat <- coef(cvfit, s = "lambda.1se")
  idx <- which(as.vector(beta_hat)[-1] != 0)
  idx
}

##-------------------------------------------------
##  Méthode 5bis : Relaxed Lasso (cv.glmnet, relax = TRUE, λ.1SE)
##  Chemin lasso puis réajustement non pénalisé des variables actives ;
##  gamma choisi par validation croisée (Hastie, Tibshirani & Tibshirani, 2020).
##-------------------------------------------------
select_relaxed_lasso <- function(X, y) {
  cvfit <- cv.glmnet(X, y, relax = TRUE,
                     standardize = TRUE, family = "gaussian")
  beta_hat <- coef(cvfit, s = "lambda.1se")   # gamma sélectionné par CV
  idx <- which(as.vector(beta_hat)[-1] != 0)
  idx
}

##-------------------------------------------------
##  Méthode 6 : SCAD (cv.ncvreg, λ.min)
##-------------------------------------------------
select_scad <- function(X, y) {
  cvfit <- cv.ncvreg(X, y, family = "gaussian", penalty = "SCAD")
  beta_hat <- coef(cvfit, lambda = cvfit$lambda.min)
  # beta_hat[1] = intercept, les suivants = prédicteurs
  idx <- which(beta_hat[-1] != 0)
  idx
}

##-------------------------------------------------
## MSE de prédiction pour un ensemble de variables sélectionnées
##-------------------------------------------------
prediction_mse <- function(X_train, y_train, X_test, y_test, selected) {
  if (length(selected) == 0) {
    # pas de prédicteur : prédire la moyenne
    y_hat <- rep(mean(y_train), length(y_test))
  } else {
    df_train <- data.frame(y = y_train, X_train[, selected, drop = FALSE])
    df_test  <- data.frame(X_test[, selected, drop = FALSE])
    colnames(df_test) <- colnames(df_train)[-1]
    fit <- lm(y ~ ., data = df_train)
    y_hat <- predict(fit, newdata = df_test)
  }
  mean((y_test - y_hat)^2)
}
