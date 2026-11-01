## ============================================================
## Robustness of the global-null guarantee BEYOND AR(1).
## Reports P(select any predictor | all beta = 0) under block,
## dense factor-model, and equicorrelated designs, using the
## transportable rho=0 surface. Also compares all competitors
## under equicorrelation, to show that the near-collinear regime
## is inhospitable to selection generally (only distribution-free
## Holm stays calibrated; AIC/BIC wildly anticonservative;
## penalized methods collapse to the empty model).
## Run: Rscript check_correlation_robustness.R
## ============================================================
suppressMessages({ library(MASS); library(parallel) })
source("cascade_reg.R"); source("get_critical_value.R"); source("method_wrappers.R")

n <- 200; p <- 20; B <- 1500; CORES <- min(40, detectCores() - 2)
alpha <- 0.05; SE <- sqrt(alpha * (1 - alpha) / B); tol <- alpha + 3 * SE

equi   <- function(p, r){ S <- matrix(r, p, p); diag(S) <- 1; S }
block  <- function(p, bs, r){ S <- diag(p); for (s in seq(1, p, by = bs)){ i <- s:min(s+bs-1, p); S[i,i] <- r }; diag(S) <- 1; S }
factorS<- function(p, q, load = 0.5, seed = 9){ set.seed(seed); L <- matrix(rnorm(p*q, 0, load), p, q)
           h <- rowSums(L^2); L <- L * sqrt(pmin(1, 0.85/h)); L %*% t(L) + diag(1 - rowSums(L^2)) }  # full rank
mkX <- function(S){ X <- MASS::mvrnorm(n, rep(0, nrow(S)), S); colnames(X) <- paste0("X", 1:ncol(X)); X }
sel <- function(expr) length(tryCatch(expr, error = function(e) integer(0))) > 0

## ---- (1) Fusion (rho=0 surface) across structures ----
cat("=== Fusion global-null FWER under non-AR(1) structures (rho=0 surface) ===\n")
structs <- list("Equicorr 0.3"=equi(p,0.3), "Equicorr 0.6"=equi(p,0.6), "Equicorr 0.9"=equi(p,0.9),
                "Block(4x5) 0.6"=block(p,5,0.6), "1-factor"=factorS(p,1), "3-factor"=factorS(p,3))
for (nm in names(structs)){
  set.seed(42); X <- mkX(structs[[nm]]); cc <- get_critical_value(1:p, n, p)
  hit <- mean(unlist(mclapply(1:B, function(b){ set.seed(5000+b); sel(select_fusion(X, rnorm(n), alpha, ccrit=cc)) }, mc.cores=CORES)))
  cat(sprintf("  %-16s FWER=%.3f  %s\n", nm, hit, if (hit <= tol) "OK" else "*** exceeds alpha ***"))
}

## ---- (2) all methods under equicorrelation (contextualising the breakdown) ----
cat("\n=== All methods, global-null FWER under equicorrelation ===\n")
for (r in c(0.6, 0.7, 0.9)){
  set.seed(42); X <- mkX(equi(p, r)); cc <- get_critical_value(1:p, n, p)
  rows <- mclapply(1:B, function(b){ set.seed(5000+b); y <- rnorm(n)
    c(Fusion=sel(select_fusion(X,y,alpha,ccrit=cc)), Holm=sel(select_holm(X,y,alpha)),
      AIC=sel(select_aic(X,y)), BIC=sel(select_bic(X,y)), Lasso=sel(select_lasso(X,y)),
      ElasticNet=sel(select_enet(X,y)), SCAD=sel(select_scad(X,y)), RelaxedLasso=sel(select_relaxed_lasso(X,y)))
  }, mc.cores=CORES)
  fw <- colMeans(do.call(rbind, rows))
  cat(sprintf("-- rho=%.1f --  ", r)); cat(paste(sprintf("%s=%.3f", names(fw), fw), collapse="  "), "\n")
}
cat(sprintf("\n(n=%d p=%d B=%d; tol=alpha+3SE=%.3f)\n", n, p, B, tol))
