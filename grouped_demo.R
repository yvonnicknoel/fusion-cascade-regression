## ============================================================
## Grouped Fusion Cascade -- demonstration with categorical predictors
## ------------------------------------------------------------
## Each predictor is a multi-level factor (a *group* of dummy columns).
## The grouped cascade removes, at each step, the group with the weakest
## joint evidence (largest drop F-test p-value) and halts when even the
## weakest remaining group is significant. Its stopping threshold p* is
## Monte Carlo calibrated under the global null so that, when all groups
## are null, the probability of retaining any group is alpha.
## Scored at the GROUP (factor) level against grouped competitors.
## ============================================================
suppressMessages({ library(parallel); library(dplyr); library(grpreg); library(MASS) })

NREP  <- as.integer(Sys.getenv("NREP",  "200"))
CORES <- as.integer(Sys.getenv("CORES", "60"))
RHO   <- as.numeric(Sys.getenv("RHO",   "0"))    # latent AR(1) correlation among factors
B_CAL <- 2000L                      # MC reps for null calibration of p*
alpha <- 0.05
LOG <- "grouped_demo.log"; cat("", file = LOG)
logm <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }

## ---- design: K factors, L levels, s_act active ----
make_groups <- function(K, L) lapply(seq_len(K), function(k) ((k-1)*(L-1) + 1):(k*(L-1)))

gen_X <- function(n, K, L, rho = 0) {
  if (rho == 0) {
    levs <- lapply(seq_len(K), function(k) sample.int(L, n, replace = TRUE))
  } else {                                           # correlated factors via latent AR(1)
    Sigma <- rho ^ abs(outer(seq_len(K), seq_len(K), "-"))
    Z <- MASS::mvrnorm(n, rep(0, K), Sigma)
    brk <- qnorm(seq(0, 1, length.out = L + 1))      # L equiprobable levels
    levs <- lapply(seq_len(K), function(k) as.integer(cut(Z[, k], brk, include.lowest = TRUE)))
  }
  cols <- lapply(levs, function(lev)
    model.matrix(~ factor(lev, levels = 1:L))[, -1, drop = FALSE])
  X <- do.call(cbind, cols); colnames(X) <- paste0("x", seq_len(ncol(X))); X
}
gen_data <- function(n, K, L, s_act, R2, groups, rho = 0) {
  X <- gen_X(n, K, L, rho)
  beta <- rep(0, ncol(X))
  for (g in groups[seq_len(s_act)]) beta[g] <- sample(c(-1, 1), length(g), TRUE)
  eta <- as.numeric(X %*% beta); v <- var(eta)
  sigma2 <- if (v == 0) 1 else v * (1 - R2) / R2
  list(X = X, y = eta + rnorm(n, 0, sqrt(sigma2)))
}

## ---- drop-test p-values for the active groups in the current model ----
group_pvals <- function(y, X, active_groups) {
  cols <- sort(unlist(active_groups)); Xa <- X[, cols, drop = FALSE]
  n <- length(y); P <- ncol(Xa); dfres <- n - P - 1
  fit <- lm.fit(cbind(1, Xa), y); rss <- sum(fit$residuals^2)
  sapply(active_groups, function(g) {
    keep <- setdiff(cols, g)
    rss0 <- if (length(keep) == 0) sum((y - mean(y))^2)
            else sum(lm.fit(cbind(1, X[, keep, drop = FALSE]), y)$residuals^2)
    Fg <- ((rss0 - rss) / length(g)) / (rss / dfres)
    pf(Fg, length(g), dfres, lower.tail = FALSE)
  })
}

## ---- calibrate p* under the global null for a fixed design X ----
calibrate_pstar <- function(X, groups, alpha, B) {
  n <- nrow(X)
  Tvals <- vapply(seq_len(B), function(b) {
    y <- rnorm(n); A <- groups; pmins <- Inf
    while (length(A) >= 1) {
      pv <- group_pvals(y, X, A); pw <- max(pv)     # weakest group's p-value
      pmins <- min(pmins, pw)
      A <- A[-which.max(pv)]                          # force-drop weakest
    }
    pmins
  }, numeric(1))
  as.numeric(quantile(Tvals, alpha, type = 8))        # P(retain any | H0) = alpha
}

## ---- the grouped Fusion Cascade ----
grouped_cascade <- function(y, X, groups, pstar) {
  gid <- seq_along(groups); A <- groups     # gid: original indices of active groups
  repeat {
    if (length(A) == 0) return(integer(0))
    pv <- group_pvals(y, X, A); j <- which.max(pv)   # weakest group
    if (pv[j] < pstar) return(gid)          # even weakest significant -> keep all remaining
    A <- A[-j]; gid <- gid[-j]              # drop weakest, continue
  }
}

## ---- competitors (group level) ----
sel_group_holm <- function(y, X, groups) {           # Holm on full-model group F-tests
  pv <- group_pvals(y, X, groups)
  which(p.adjust(pv, "holm") < alpha)
}
sel_group_bic <- function(y, X, groups) {            # greedy backward BIC on groups
  A <- seq_along(groups); n <- length(y)
  bic <- function(gs) { cols <- sort(unlist(groups[gs]))
    Xa <- if (length(cols)) cbind(1, X[, cols, drop = FALSE]) else matrix(1, n, 1)
    rss <- sum(lm.fit(Xa, y)$residuals^2); n*log(rss/n) + (ncol(Xa))*log(n) }
  cur <- bic(A)
  repeat {
    if (length(A) == 0) break
    cand <- sapply(seq_along(A), function(i) bic(A[-i]))
    if (min(cand) < cur) { cur <- min(cand); A <- A[-which.min(cand)] } else break
  }
  A
}
sel_group_lasso <- function(y, X, groups) {          # group LASSO (grpreg, CV)
  gv <- integer(ncol(X)); for (k in seq_along(groups)) gv[groups[[k]]] <- k
  cv <- cv.grpreg(X, y, group = gv, penalty = "grLasso", nfolds = 10)
  b <- as.numeric(coef(cv, lambda = cv$lambda.min))[-1]
  which(sapply(groups, function(g) any(b[g] != 0)))
}

gmetrics <- function(sel, active, K) {
  TP <- length(intersect(sel, active)); FP <- length(setdiff(sel, active))
  FN <- length(setdiff(active, sel))
  data.frame(gTPR = TP/length(active), gFPR = FP/(K - length(active)),
             n_groups = length(sel))
}

## ============================================================
## Simulation
## ============================================================
K <- 8; L <- 4; s_act <- 2; groups <- make_groups(K, L); active <- seq_len(s_act)
grid <- expand.grid(n = c(100, 200, 500), R2 = c(0.2, 0.5))
logm("Grouped demo: K=%d factors (%d levels), %d active; NREP=%d", K, L, s_act, NREP)

run_scen <- function(s) {
  n <- grid$n[s]; R2 <- grid$R2[s]
  set.seed(7000 + s)
  pstar <- calibrate_pstar(gen_X(n, K, L, RHO), groups, alpha, B_CAL)   # per-design calibration
  reps <- mclapply(seq_len(NREP), function(r) {
    set.seed(7000 + 1000*s + r)
    d <- gen_data(n, K, L, s_act, R2, groups, RHO)
    meth <- list(Fusion = grouped_cascade(d$y, d$X, groups, pstar),
                 Holm   = sel_group_holm(d$y, d$X, groups),
                 BIC    = sel_group_bic(d$y, d$X, groups),
                 gLASSO = sel_group_lasso(d$y, d$X, groups))
    do.call(rbind, lapply(names(meth), function(m)
      data.frame(method = m, gmetrics(meth[[m]], active, K))))
  }, mc.cores = CORES)
  out <- bind_rows(reps[vapply(reps, is.data.frame, logical(1))])
  out$n <- n; out$R2 <- R2
  logm("  n=%3d R2=%.1f  pstar=%.4f  done", n, R2, pstar)
  out
}
res <- bind_rows(lapply(seq_len(nrow(grid)), run_scen))
saveRDS(res, sprintf("grouped_demo_results_rho%.1f.rds", RHO))

## overall and by-n summaries at the group level
ov <- res %>% group_by(method) %>%
  summarise(gTPR = mean(gTPR), gFPR = mean(gFPR), n_groups = mean(n_groups), .groups="drop") %>%
  arrange(match(method, c("Holm","Fusion","gLASSO","BIC")))
write.csv2(ov, sprintf("grouped_demo_overall_rho%.1f.csv", RHO), row.names = FALSE)
logm("\n=== Group-level performance (rho=%.1f, mean over all scenarios) ===", RHO)
print(as.data.frame(ov), digits = 3)
