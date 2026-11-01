## ============================================================
## SAPA spi worked example -- stability via repeated splits
## R = 50 random 70/30 train-test splits. For each split, fit all
## seven methods and record model size + out-of-sample R^2; also
## record each method's selected item set to measure (i) how often
## Fusion coincides with BIC, and (ii) how stable Fusion's item
## selection is. Repeats run in parallel (step() dominates).
## ============================================================
suppressMessages({ library(psychTools); library(parallel); library(dplyr) })
source("get_critical_value.R"); source("cascade_reg.R")
source("generate_data.R"); source("method_wrappers.R")

R     <- as.integer(Sys.getenv("R",     "50"))
CORES <- as.integer(Sys.getenv("CORES", "50"))
LOG <- "spi_cv.log"; cat("", file = LOG)
logm <- function(...) { m <- sprintf(...); cat(m,"\n"); cat(m,"\n",file=LOG,append=TRUE) }

data(spi)
demo  <- c("age","sex","health","education","p1edu","p2edu","wellness","exer","smoke","ER")
items <- setdiff(colnames(spi), demo)
d <- spi[!is.na(spi$health), ]
X <- as.matrix(d[, items]); y <- as.numeric(d$health)
n <- nrow(X); p <- ncol(X)
methods <- c("Fusion","Holm","AIC","BIC","Lasso","ElasticNet","SCAD","RelaxedLasso")
logm("SAPA spi: predict health from %d items; n=%d; %d repeated 70/30 splits", p, n, R)

run_rep <- function(r) {
  set.seed(1000 + r)
  tr <- sample(n, floor(0.7 * n))
  Xtr<-X[tr,]; ytr<-y[tr]; Xte<-X[-tr,]; yte<-y[-tr]
  ntr<-nrow(Xtr); base<-mean((yte-mean(ytr))^2)
  cc <- get_critical_value(1:p, ntr, p)
  pick <- function(m) switch(m,
    Fusion=select_fusion(Xtr,ytr,alpha=0.05,ccrit=cc), Holm=select_holm(Xtr,ytr,0.05),
    AIC=select_aic(Xtr,ytr), BIC=select_bic(Xtr,ytr), Lasso=select_lasso(Xtr,ytr),
    ElasticNet=select_enet(Xtr,ytr), SCAD=select_scad(Xtr,ytr),
    RelaxedLasso=select_relaxed_lasso(Xtr,ytr))
  sets <- list(); perf <- list()
  for (m in methods) {
    s <- tryCatch(pick(m), error=function(e) integer(0)); s <- s[!is.na(s)]
    sets[[m]] <- s
    mse <- prediction_mse(Xtr,ytr,Xte,yte,s)
    perf[[m]] <- data.frame(rep=r, method=m, size=length(s), mse=mse, r2=1-mse/base)
  }
  cat(sprintf("rep %d done\n", r), file=LOG, append=TRUE)
  list(perf=do.call(rbind, perf), sets=sets)
}

res <- mclapply(seq_len(R), run_rep, mc.cores=CORES, mc.preschedule=FALSE)
ok  <- vapply(res, is.list, logical(1)); res <- res[ok]
perf <- bind_rows(lapply(res, `[[`, "perf"))

## ---- performance summary ----
summ <- perf %>% group_by(method) %>%
  summarise(size_mean=mean(size), size_sd=sd(size),
            r2_mean=mean(r2), r2_sd=sd(r2), mse_mean=mean(mse), .groups="drop") %>%
  arrange(size_mean)
logm("\n=== Performance over %d splits (mean +/- sd) ===", length(res))
logm("%-11s %12s %14s", "method", "size", "test R2")
for (i in seq_len(nrow(summ))) logm("%-11s  %5.1f +/- %3.1f   %.3f +/- %.3f",
  summ$method[i], summ$size_mean[i], summ$size_sd[i], summ$r2_mean[i], summ$r2_sd[i])

## ---- Fusion vs BIC agreement, and Fusion vs others overlap ----
jac <- function(a,b) if(length(union(a,b))==0) 1 else length(intersect(a,b))/length(union(a,b))
agree_bic <- sapply(res, function(z) setequal(z$sets$Fusion, z$sets$BIC))
jac_bic   <- sapply(res, function(z) jac(z$sets$Fusion, z$sets$BIC))
logm("\n=== Fusion vs BIC across splits ===")
logm("exact-equal selection: %d/%d splits (%.0f%%);  mean Jaccard = %.2f",
     sum(agree_bic), length(res), 100*mean(agree_bic), mean(jac_bic))
logm("mean Jaccard(Fusion, .):")
for (m in setdiff(methods,"Fusion"))
  logm("  %-11s %.2f", m, mean(sapply(res, function(z) jac(z$sets$Fusion, z$sets[[m]]))))

## ---- Fusion item-selection stability ----
freq <- table(unlist(lapply(res, function(z) items[z$sets$Fusion])))
freq <- sort(freq, decreasing=TRUE)
dict <- tryCatch({ data(spi.dictionary, package="psychTools"); get("spi.dictionary") }, error=function(e) NULL)
logm("\n=== Items selected by Fusion in >=50%% of splits ===")
for (it in names(freq)[freq >= 0.5*length(res)]) {
  txt <- if(!is.null(dict) && it %in% rownames(dict)) dict[it,"item"] else ""
  logm("  %4.0f%%  %-8s %s", 100*freq[it]/length(res), it, txt)
}
saveRDS(list(perf=perf, summ=summ, agree_bic=agree_bic, freq=freq), "spi_cv_results.rds")
logm("\nDone (%d/%d splits ok).", length(res), R)
