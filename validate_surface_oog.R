source("calibrate_c_fast.R"); source("generate_data.R"); source("get_critical_value.R")
## Held-out (n,p): none is a calibration-grid node (n in {50,100,150,200,500,1000,2000,5000,10000},
## p in {5,10,20,30,40,60,100,150,250,350,500}). Genuine interpolation test.
grid <- rbind(c(75,15), c(300,15),c(300,80),c(300,200),
              c(750,15),c(750,80),c(750,200),c(750,400),
              c(1500,15),c(1500,80),c(1500,200),c(1500,400))
LOG <- "validation_oog.log"; cat("Out-of-grid validation of the 32-parameter surface (B=8000, rho=0)\n", file=LOG)
res <- data.frame(); allrel <- numeric(0)
for(i in 1:nrow(grid)){
  n <- grid[i,1]; p <- grid[i,2]
  set.seed(100+i); X <- gen_X_ar1(n,p,0)
  mc   <- calibrate_cvec_mc_fast(X, alpha=0.05, B=8000, seed=200+i)
  pred <- get_critical_value(1:p, n, p, 0.05)
  rel  <- abs(mc - pred)/pred
  allrel <- c(allrel, rel)
  res <- rbind(res, data.frame(n=n, p=p, med=100*median(rel), p95=100*quantile(rel,.95), max=100*max(rel)))
  cat(sprintf("done n=%4d p=%3d : median=%.2f%%  p95=%.2f%%  max=%.2f%%\n",
              n,p,100*median(rel),100*quantile(rel,.95),100*max(rel)), file=LOG, append=TRUE)
}
cat("\n=== PER-POINT (rel. error %, held-out) ===\n", file=LOG, append=TRUE)
capture.output(print(res, row.names=FALSE, digits=3), file=LOG, append=TRUE) -> ig
cat(sprintf("\n=== POOLED over all %d held-out c(k) values ===\n", length(allrel)), file=LOG, append=TRUE)
cat(sprintf("median %.2f%%   95th pct %.2f%%   max %.2f%%\n",
            100*median(allrel), 100*quantile(allrel,.95), 100*max(allrel)), file=LOG, append=TRUE)
saveRDS(list(res=res, allrel=allrel), "validation_oog.rds")
cat("VALIDATION DONE\n", file=LOG, append=TRUE)
