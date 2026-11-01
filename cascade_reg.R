source("calibrate_c_reg.R")

## ------------------------------------------------------------
## Fusion Cascade pour la régression linéaire avec c(m)
## - X : matrice (ou data.frame) de prédicteurs, sans intercept (n x p)
## - y : vecteur réponse (longueur n)
## - alpha : niveau global (utile si c_vec et ccrit sont NULL)
## - c_vec : vecteur de seuils c(m), m = 1..p (prioritaire si non NULL)
## - ccrit : seuil critique global scalaire (optionnel).
##           S'il est fourni et c_vec est NULL, on met c_vec = rep(ccrit, p).
## - B : nb de simulations pour calibrer c_vec si c_vec et ccrit sont NULL
## - standardize : TRUE = standardisation des prédicteurs
## - verbose : messages d'avancement
## ------------------------------------------------------------

fusion_cascade_reg <- function(X, y,
                               alpha = 0.05,
                               c_vec = NULL,
                               ccrit = NULL,
                               B = 2000,
                               standardize = TRUE,
                               verbose = TRUE) {
  X <- as.matrix(X)
  y <- as.numeric(y)

  n <- nrow(X)
  p <- ncol(X)

  if (length(y) != n) {
    stop("X et y n'ont pas le même nombre de lignes.")
  }
  if (p < 1) {
    stop("Il faut au moins un prédicteur (p >= 1).")
  }

  predictor_names <- colnames(X)
  if (is.null(predictor_names)) {
    predictor_names <- paste0("X", seq_len(p))
    colnames(X) <- predictor_names
  }

  ## 1. Standardisation des prédicteurs (colonnes de X)
  col_means <- rep(0, p)
  col_sds   <- rep(1, p)

  if (standardize) {
    col_means <- colMeans(X)
    col_sds   <- apply(X, 2, sd)

    if (any(col_sds == 0)) {
      idx0 <- which(col_sds == 0)
      stop("Prédicteur(s) à variance nulle : ",
           paste(predictor_names[idx0], collapse = ", "))
    }

    X <- sweep(X, 2, col_means, FUN = "-")
    X <- sweep(X, 2, col_sds,   FUN = "/")
  }

  ## Matrice de régression avec intercept
  Xd <- cbind("(Intercept)" = 1, X)

  ## ------------------------------------------------------------
  ## 2. Détermination de c_vec(m)
  ## ------------------------------------------------------------

  B_used <- NA_integer_

  if (!is.null(c_vec)) {
    # c_vec fourni : vérifier la taille et la monotonie
    if (length(c_vec) < p) {
      stop("c_vec doit être de longueur au moins p.")
    }
    c_vec <- c_vec[1:p]
    if (verbose) {
      message("Vecteur c_vec fourni par l'utilisateur.")
    }
  } else if (!is.null(ccrit)) {
    if (length(ccrit) > 1) {
      # ccrit vectoriel : traité comme c_vec (seuils par taille de famille c(m))
      if (length(ccrit) < p) stop("ccrit vectoriel doit être de longueur au moins p.")
      c_vec <- ccrit[1:p]
      if (verbose) message("ccrit vectoriel fourni : traité comme c_vec (seuils c(m)).")
    } else {
      # ccrit scalaire fourni : seuil constant
      c_vec <- rep(ccrit, p)
      if (verbose) {
        message(sprintf("Seuil scalaire fourni (ccrit = %.3f). Utilisation de c_vec constant.", ccrit))
      }
    }
  } else {
    # Ni c_vec ni ccrit -> calibration Monte Carlo
    if (is.null(B) || B < 1) {
      stop("Si c_vec et ccrit sont NULL, il faut un B >= 1 pour la calibration Monte Carlo.")
    }
    if (verbose) {
      message("Calibration de c_vec(m) par simulation (global null)...")
    }
    # On calibre sur X non standardisé pour rester cohérent avec la dépendance brute.
    # Tu peux remplacer par la version standardisée si tu préfères.
    c_vec <- calibrate_cvec_mc(X, alpha = alpha, B = B)
    B_used <- B
    if (verbose) {
      message("c_vec calibré :")
      message(paste0("  c_vec = ", paste(round(c_vec, 3), collapse = ", ")))
    }
  }

  # S'assurer que c_vec est non décroissant (par sécurité)
  c_vec <- cummax(c_vec)

  ## ------------------------------------------------------------
  ## 3. Cascade : suppression séquentielle des prédicteurs
  ## ------------------------------------------------------------

  act <- seq_len(p)   # indices des colonnes de X encore actives
  steps <- list()
  step  <- 0

  repeat {
    m <- length(act)
    if (m == 0) {
      if (verbose) {
        message("Tous les prédicteurs ont été supprimés (modèle à intercept seul).")
      }
      break
    }

    step <- step + 1

    cols_used <- c(1, 1 + act) # 1 = intercept ; 1+act = colonnes prédicteurs actifs
    X_used    <- Xd[, cols_used, drop = FALSE]

    fit <- lm.fit(X_used, y)

    beta    <- fit$coefficients
    res     <- fit$residuals
    df_res  <- fit$df.residual
    sigma2  <- sum(res^2) / df_res

    XtX_inv <- solve(crossprod(X_used))
    se      <- sqrt(diag(XtX_inv) * sigma2)

    tvals   <- beta / se

    # Parties "prédicteurs" (sans intercept)
    beta_pred <- beta[-1]
    t_pred    <- tvals[-1]

    # Choix du candidat : plus petit |beta| (sur échelle standardisée)
    # k_local <- which.min(abs(beta_pred))

    # Choix du candidat : plus petit |T|
    k_local <- which.min(abs(t_pred))
    j_act   <- act[k_local]                  # index global dans X
    name_j  <- predictor_names[j_act]

    beta_j  <- beta_pred[k_local]
    t_j     <- t_pred[k_local]

    # Seuil local : procédure STEP-DOWN (max-T, Westfall-Young).
    # À la taille active m, le candidat testé (plus petit |t| parmi les m
    # survivants) est le m-ième plus grand |t| global ; sous le null global il
    # doit franchir le quantile du max sur (p - m + 1) hypothèses, soit
    # c_vec[p - m + 1] (et NON c_vec[m]). c_vec est croissant [c(1)..c(p)] ;
    # ainsi le dernier survivant (max|t|, m = 1) affronte c(p), le max complet.
    c_local <- c_vec[p - m + 1]

    # Décision : |t_j| > c_local -> rejet de H_j (on s'arrête)
    reject_Hj <- (abs(t_j) > c_local)

    if (verbose) {
      message(sprintf("Étape %d (m = %d) : candidat = %s | beta=%.4f (std), t=%.3f, c(m)=%.3f -> %s",
                      step, m, name_j, beta_j, t_j, c_local,
                      if (reject_Hj) "REJET (arrêt de la cascade)" else "ACCEPTÉ (supprimé)"))
    }

    steps[[step]] <- data.frame(
      step       = step,
      m_family   = m,
      predictor  = name_j,
      beta_std   = beta_j,
      t_value    = t_j,
      c_local    = c_local,
      removed    = !reject_Hj,
      stringsAsFactors = FALSE
    )

    if (reject_Hj) {
      # On s'arrête : on garde le modèle de l'étape actuelle sans supprimer ce prédicteur
      break
    } else {
      # On accepte H_j: beta_j = 0 -> on supprime j_act de la liste active
      act <- setdiff(act, j_act)
    }
  }

  steps_df <- do.call(rbind, steps)

  ## Prédicteurs restants à la fin de la cascade
  active_predictors <- predictor_names[act]

  ## Ajustement final avec les prédicteurs restants (sur X standardisé)
  if (length(act) == 0) {
    # modèle à intercept seul
    X_final <- matrix(1, nrow = n, ncol = 1)
    colnames(X_final) <- "(Intercept)"
  } else {
    cols_final <- c(1, 1 + act)
    X_final    <- Xd[, cols_final, drop = FALSE]
  }

  fit_final <- lm.fit(X_final, y)

  beta_final   <- fit_final$coefficients
  res_final    <- fit_final$residuals
  df_res_final <- fit_final$df.residual
  sigma2_final <- sum(res_final^2) / df_res_final
  XtX_inv_fin  <- solve(crossprod(X_final))
  se_final     <- sqrt(diag(XtX_inv_fin) * sigma2_final)
  t_final      <- beta_final / se_final

  result <- list(
    alpha             = alpha,
    B                 = B_used,
    c_vec             = c_vec,
    standardized      = standardize,
    predictor_names   = predictor_names,
    col_means         = col_means,
    col_sds           = col_sds,
    steps             = steps_df,
    active_indices    = act,
    active_predictors = active_predictors,
    final = list(
      X            = X_final,
      coefficients = beta_final,
      se           = se_final,
      t_values     = t_final,
      df_residual  = df_res_final,
      sigma2       = sigma2_final
    )
  )

  class(result) <- "fusionCascadeReg"
  return(result)
}

## ------------------------------------------------------------
## Méthode d'impression
## ------------------------------------------------------------

print.fusionCascadeReg <- function(x, ...) {
  cat("Fusion Cascade pour la régression linéaire (version c(m))\n")
  cat(sprintf("  alpha = %.3f, B = %s\n",
              x$alpha,
              ifelse(is.na(x$B), "fourni / non utilisé", as.character(x$B))))
  cat("  Standardisation des prédicteurs :", if (x$standardized) "oui" else "non", "\n")
  cat("  c_vec(m) (1..p) :\n")
  print(round(x$c_vec, 3))
  cat("\nÉtapes de la cascade :\n")
  print(x$steps, row.names = FALSE)
  cat("\nPrédicteurs encore actifs à la fin :\n")
  if (length(x$active_predictors) == 0) {
    cat("  (aucun : modèle à intercept seul)\n")
  } else {
    cat("  ", paste(x$active_predictors, collapse = ", "), "\n")
  }
  invisible(x)
}


## ------------------------------------------------------------
## Fusion Cascade pour la régression linéaire avec d(m)
## VERSION "HOLM-LIKE" : utilise des seuils de p-values
## Règle de décision : p-value ≤ alpha * d(m) / m
##
## Cette formulation est IDENTIQUE à l'ANOVA Fusion Cascade
## et permet une unification théorique complète des deux méthodes.
##
## - X : matrice (ou data.frame) de prédicteurs, sans intercept (n x p)
## - y : vecteur réponse (longueur n)
## - alpha : niveau global (utile si d_vec et dcrit sont NULL)
## - d_vec : vecteur de multiplicateurs d(m), m = 1..p (prioritaire si non NULL)
## - dcrit : multiplicateur global scalaire (optionnel).
##           S'il est fourni et d_vec est NULL, on met d_vec = rep(dcrit, p).
## - B : nb de simulations pour calibrer d_vec si d_vec et dcrit sont NULL
## - standardize : TRUE = standardisation des prédicteurs
## - verbose : messages d'avancement
## ------------------------------------------------------------

fusion_cascade_reg_alpha <- function(X, y,
                                      alpha = 0.05,
                                      d_vec = NULL,
                                      dcrit = NULL,
                                      B = 2000,
                                      standardize = TRUE,
                                      verbose = TRUE) {
  X <- as.matrix(X)
  y <- as.numeric(y)

  n <- nrow(X)
  p <- ncol(X)

  if (length(y) != n) {
    stop("X et y n'ont pas le même nombre de lignes.")
  }
  if (p < 1) {
    stop("Il faut au moins un prédicteur (p >= 1).")
  }

  predictor_names <- colnames(X)
  if (is.null(predictor_names)) {
    predictor_names <- paste0("X", seq_len(p))
    colnames(X) <- predictor_names
  }

  ## 1. Standardisation des prédicteurs (colonnes de X)
  col_means <- rep(0, p)
  col_sds   <- rep(1, p)

  if (standardize) {
    col_means <- colMeans(X)
    col_sds   <- apply(X, 2, sd)

    if (any(col_sds == 0)) {
      idx0 <- which(col_sds == 0)
      stop("Prédicteur(s) à variance nulle : ",
           paste(predictor_names[idx0], collapse = ", "))
    }

    X <- sweep(X, 2, col_means, FUN = "-")
    X <- sweep(X, 2, col_sds,   FUN = "/")
  }

  ## Matrice de régression avec intercept
  Xd <- cbind("(Intercept)" = 1, X)

  ## ------------------------------------------------------------
  ## 2. Détermination de d_vec(m)
  ## ------------------------------------------------------------

  B_used <- NA_integer_

  if (!is.null(d_vec)) {
    # d_vec fourni : vérifier la taille et la monotonie
    if (length(d_vec) < p) {
      stop("d_vec doit être de longueur au moins p.")
    }
    d_vec <- d_vec[1:p]
    if (verbose) {
      message("Vecteur d_vec fourni par l'utilisateur.")
    }
  } else if (!is.null(dcrit)) {
    # dcrit scalaire fourni : on le recycle en vecteur constant
    d_vec <- rep(dcrit, p)
    if (verbose) {
      message(sprintf("Multiplicateur scalaire fourni (dcrit = %.3f). Utilisation de d_vec constant.", dcrit))
    }
  } else {
    # Ni d_vec ni dcrit -> calibration Monte Carlo
    if (is.null(B) || B < 1) {
      stop("Si d_vec et dcrit sont NULL, il faut un B >= 1 pour la calibration Monte Carlo.")
    }
    if (verbose) {
      message("Calibration de d_vec(m) par simulation (global null)...")
    }
    # Appel à calibrate_dvec_mc (qui doit être sourcé depuis calibrate_c_reg.R)
    d_vec <- calibrate_dvec_mc(X, alpha = alpha, B = B)
    B_used <- B
    if (verbose) {
      message("d_vec calibré :")
      message(paste0("  d_vec = ", paste(round(d_vec, 3), collapse = ", ")))
    }
  }

  # S'assurer que d_vec est non décroissant (par sécurité)
  d_vec <- cummax(d_vec)

  ## ------------------------------------------------------------
  ## 3. Cascade : suppression séquentielle des prédicteurs
  ## ------------------------------------------------------------

  act <- seq_len(p)   # indices des colonnes de X encore actives
  steps <- list()
  step  <- 0

  repeat {
    m <- length(act)
    if (m == 0) {
      if (verbose) {
        message("Tous les prédicteurs ont été supprimés (modèle à intercept seul).")
      }
      break
    }

    step <- step + 1

    cols_used <- c(1, 1 + act) # 1 = intercept ; 1+act = colonnes prédicteurs actifs
    X_used    <- Xd[, cols_used, drop = FALSE]

    fit <- lm.fit(X_used, y)

    beta    <- fit$coefficients
    res     <- fit$residuals
    df_res  <- fit$df.residual
    sigma2  <- sum(res^2) / df_res

    XtX_inv <- solve(crossprod(X_used))
    se      <- sqrt(diag(XtX_inv) * sigma2)

    tvals   <- beta / se

    # Parties "prédicteurs" (sans intercept)
    beta_pred <- beta[-1]
    t_pred    <- tvals[-1]

    # Calcul des p-values bilatérales
    pval_pred <- 2 * pt(-abs(t_pred), df = df_res)

    # Choix du candidat : plus petit |T| (équivalent : plus petite p-value)
    k_local <- which.min(abs(t_pred))
    j_act   <- act[k_local]                  # index global dans X
    name_j  <- predictor_names[j_act]

    beta_j  <- beta_pred[k_local]
    t_j     <- t_pred[k_local]
    pval_j  <- pval_pred[k_local]

    # Seuil local en fonction de la taille de la famille actuelle m
    # Seuil de type Holm modifié : alpha * d(m) / m
    threshold_alpha <- alpha * d_vec[m] / m

    # Décision : p-value > threshold -> rejet de H_j (on garde le prédicteur, arrêt)
    reject_Hj <- (pval_j > threshold_alpha)

    if (verbose) {
      message(sprintf("Étape %d (m = %d) : candidat = %s | beta=%.4f (std), t=%.3f, p-val=%.4f, seuil=%.4f (alpha*d(m)/m) -> %s",
                      step, m, name_j, beta_j, t_j, pval_j, threshold_alpha,
                      if (reject_Hj) "REJET (arrêt de la cascade)" else "ACCEPTÉ (supprimé)"))
    }

    steps[[step]] <- data.frame(
      step       = step,
      m_family   = m,
      predictor  = name_j,
      beta_std   = beta_j,
      t_value    = t_j,
      p_value    = pval_j,
      threshold  = threshold_alpha,
      d_local    = d_vec[m],
      removed    = !reject_Hj,
      stringsAsFactors = FALSE
    )

    if (reject_Hj) {
      # On s'arrête : on garde le modèle de l'étape actuelle sans supprimer ce prédicteur
      break
    } else {
      # On accepte H_j: beta_j = 0 -> on supprime j_act de la liste active
      act <- setdiff(act, j_act)
    }
  }

  steps_df <- do.call(rbind, steps)

  ## Prédicteurs restants à la fin de la cascade
  active_predictors <- predictor_names[act]

  ## Ajustement final avec les prédicteurs restants (sur X standardisé)
  if (length(act) == 0) {
    # modèle à intercept seul
    X_final <- matrix(1, nrow = n, ncol = 1)
    colnames(X_final) <- "(Intercept)"
  } else {
    cols_final <- c(1, 1 + act)
    X_final    <- Xd[, cols_final, drop = FALSE]
  }

  fit_final <- lm.fit(X_final, y)

  beta_final   <- fit_final$coefficients
  res_final    <- fit_final$residuals
  df_res_final <- fit_final$df.residual
  sigma2_final <- sum(res_final^2) / df_res_final
  XtX_inv_fin  <- solve(crossprod(X_final))
  se_final     <- sqrt(diag(XtX_inv_fin) * sigma2_final)
  t_final      <- beta_final / se_final
  pval_final   <- 2 * pt(-abs(t_final), df = df_res_final)

  result <- list(
    alpha             = alpha,
    B                 = B_used,
    d_vec             = d_vec,
    standardized      = standardize,
    predictor_names   = predictor_names,
    col_means         = col_means,
    col_sds           = col_sds,
    steps             = steps_df,
    active_indices    = act,
    active_predictors = active_predictors,
    final = list(
      X            = X_final,
      coefficients = beta_final,
      se           = se_final,
      t_values     = t_final,
      p_values     = pval_final,
      df_residual  = df_res_final,
      sigma2       = sigma2_final
    )
  )

  class(result) <- "fusionCascadeRegAlpha"
  return(result)
}

## ------------------------------------------------------------
## Méthode d'impression pour fusion_cascade_reg_alpha
## ------------------------------------------------------------

print.fusionCascadeRegAlpha <- function(x, ...) {
  cat("Fusion Cascade pour la régression linéaire (version Holm-like avec d(m))\n")
  cat(sprintf("  alpha = %.3f, B = %s\n",
              x$alpha,
              ifelse(is.na(x$B), "fourni / non utilisé", as.character(x$B))))
  cat("  Standardisation des prédicteurs :", if (x$standardized) "oui" else "non", "\n")
  cat("  d_vec(m) (1..p) : multiplicateurs de Holm\n")
  print(round(x$d_vec, 3))
  cat("\nÉtapes de la cascade :\n")
  print(x$steps, row.names = FALSE)
  cat("\nPrédicteurs encore actifs à la fin :\n")
  if (length(x$active_predictors) == 0) {
    cat("  (aucun : modèle à intercept seul)\n")
  } else {
    cat("  ", paste(x$active_predictors, collapse = ", "), "\n")
  }
  invisible(x)
}
