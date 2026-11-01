##-------------------------------------------------
## Calibration Monte Carlo du vecteur c_vec(m), m = 1..p
## pour une matrice de design X (sans intercept).
## - X : n x p (sans intercept)
## - alpha : niveau global
## - B : nombre de répétitions Monte Carlo
## Retour : numeric de longueur p, c_vec[m] = seuil pour une famille de taille m.
##-------------------------------------------------

calibrate_cvec_mc <- function(X, alpha = 0.05, B = 5000) {

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  if (p < 1) stop("Il faut au moins un prédicteur (p >= 1).")

  ## Standardisation EXACTEMENT comme dans fusion_cascade_reg :
  ## - centrage par la moyenne de colonne
  ## - division par l'écart-type (sd de R, n-1)
  col_means <- colMeans(X)
  col_sds   <- apply(X, 2, sd)

  if (any(col_sds == 0)) {
    idx0 <- which(col_sds == 0)
    stop("Prédicteur(s) à variance nulle pendant la calibration : ",
         paste(colnames(X)[idx0], collapse = ", "))
  }

  X_std <- sweep(X, 2, col_means, FUN = "-")
  X_std <- sweep(X_std, 2, col_sds,   FUN = "/")

  ## Data frame standardisé pour utiliser lm(y ~ .)
  dfX <- as.data.frame(X_std)
  if (is.null(colnames(dfX))) {
    colnames(dfX) <- paste0("X", seq_len(p))
  }

  ## Matrice B x p : pour chaque b, M[b, m] = max |T_j|
  ## sur un sous-ensemble aléatoire de taille m (via une permutation)
  M <- matrix(NA_real_, nrow = B, ncol = p)

  for (b in seq_len(B)) {
    ## y simulé sous le global null
    y_null <- rnorm(n, mean = 0, sd = 1)

    ## Modèle linéaire complet : y ~ X1 + X2 + ... + Xp (standardisés)
    df_null  <- cbind.data.frame(y = y_null, dfX)
    fit_null <- lm(y ~ ., data = df_null)

    ## Coefficients + SE + t-values via coef(summary())
    coef_tab <- coef(summary(fit_null))
    ## coef_tab est une matrice (p+1) x 4 : Estimate, Std. Error, t value, Pr(>|t|)

    ## On enlève l'intercept : indices 2..(p+1)
    t_pred <- coef_tab[-1, "t value"]
    absT   <- abs(t_pred)

    ## Subsets aléatoires : permutation des p prédicteurs, puis cummax
    perm      <- sample.int(p, size = p, replace = FALSE)
    absT_perm <- absT[perm]

    ## M[b, m] = max_{j dans subset de taille m} |T_j|
    M[b, ] <- cummax(absT_perm)
  }

  ## Quantiles (1 - alpha) pour chaque m = 1..p
  c_vec <- apply(M, 2, stats::quantile, probs = 1 - alpha, type = 8)

  ## Sécurisation : c(m) non décroissante
  c_vec <- cummax(c_vec)

  c_vec
}


##-------------------------------------------------
## Calibration Monte Carlo du vecteur d_vec(m), m = 1..p
## Version "Holm-like" : seuils multiplicatifs pour p-values
## Règle de décision : p-value ≤ alpha * d(m) / m
##
## Cette formulation maximise la ressemblance avec l'ANOVA Fusion Cascade
## où on comparait p-value ≤ alpha * c(k) / k
##
## - X : n x p (sans intercept)
## - alpha : niveau global
## - B : nombre de répétitions Monte Carlo
## Retour : numeric de longueur p, d_vec[m] = multiplicateur pour une famille de taille m
##-------------------------------------------------

calibrate_dvec_mc <- function(X, alpha = 0.05, B = 5000) {

  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  if (p < 1) stop("Il faut au moins un prédicteur (p >= 1).")

  ## Standardisation EXACTEMENT comme dans fusion_cascade_reg :
  ## - centrage par la moyenne de colonne
  ## - division par l'écart-type (sd de R, n-1)
  col_means <- colMeans(X)
  col_sds   <- apply(X, 2, sd)

  if (any(col_sds == 0)) {
    idx0 <- which(col_sds == 0)
    stop("Prédicteur(s) à variance nulle pendant la calibration : ",
         paste(colnames(X)[idx0], collapse = ", "))
  }

  X_std <- sweep(X, 2, col_means, FUN = "-")
  X_std <- sweep(X_std, 2, col_sds,   FUN = "/")

  ## Data frame standardisé pour utiliser lm(y ~ .)
  dfX <- as.data.frame(X_std)
  if (is.null(colnames(dfX))) {
    colnames(dfX) <- paste0("X", seq_len(p))
  }

  ## Matrice B x p : pour chaque b, M_pval[b, m] = min(p-value)
  ## sur un sous-ensemble aléatoire de taille m (via une permutation)
  M_pval <- matrix(NA_real_, nrow = B, ncol = p)

  for (b in seq_len(B)) {
    ## y simulé sous le global null
    y_null <- rnorm(n, mean = 0, sd = 1)

    ## Modèle linéaire complet : y ~ X1 + X2 + ... + Xp (standardisés)
    df_null  <- cbind.data.frame(y = y_null, dfX)
    fit_null <- lm(y ~ ., data = df_null)

    ## Coefficients + SE + t-values + p-values via coef(summary())
    ## coef_tab est une matrice (p+1) x 4 : Estimate, Std. Error, t value, Pr(>|t|)
    coef_tab <- coef(summary(fit_null))

    ## On enlève l'intercept : indices 2..(p+1)
    pvals <- coef_tab[-1, "Pr(>|t|)"]

    ## Subsets aléatoires : permutation des p prédicteurs, puis cummin
    ## (car on cherche le MIN des p-values dans chaque subset)
    perm       <- sample.int(p, size = p, replace = FALSE)
    pvals_perm <- pvals[perm]

    ## M_pval[b, m] = min_{j dans subset de taille m} p-value_j
    M_pval[b, ] <- cummin(pvals_perm)
  }

  ## Pour chaque taille de famille m, on veut d(m) tel que :
  ## P(min p-value ≤ alpha * d(m) / m | H0) ≈ alpha
  ##
  ## Équivalent : on cherche le quantile alpha de la distribution
  ## de (m/alpha) * min(p-value), qui donne d(m)

  d_vec <- numeric(p)
  for (m in 1:p) {
    # Ratio : combien de fois plus libéral que Holm standard ?
    # Si d(m) = 1, on retrouve exactement Holm
    # Si d(m) > 1, la procédure est plus libérale
    ratio <- (m / alpha) * M_pval[, m]
    d_vec[m] <- quantile(ratio, probs = alpha, type = 8)
  }

  ## Sécurisation : d(m) non décroissante (cohérence avec la structure séquentielle)
  d_vec <- cummax(d_vec)

  return(d_vec)
}
