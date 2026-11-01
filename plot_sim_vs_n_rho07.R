## ============================================================
## Figure: behaviour under collinearity (rho=0.7) vs sample size.
## Black-and-white safe: every method has a distinct colour,
## line type AND point shape, so it is distinguishable in grayscale.
## No on-plot title/subtitle (the LaTeX caption carries the text).
## Run: Rscript plot_sim_vs_n_rho07.R  ->  plot_sim_vs_n_rho07.pdf
## ============================================================
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

d <- read.csv2("simulation_by_n_rho_aic.csv", stringsAsFactors = FALSE)
for (k in c("n","rho","TPR","FPR")) d[[k]] <- as.numeric(d[[k]])
d <- subset(d, rho == 0.7)
d$method <- recode(d$method, ElasticNet = "Elastic Net", Lasso = "LASSO",
                   RelaxedLasso = "Relaxed LASSO")
meth <- c("Fusion","Holm","AIC","BIC","LASSO","Elastic Net","Relaxed LASSO","SCAD")
d$method <- factor(d$method, levels = meth)

long <- pivot_longer(d[c("method","n","TPR","FPR")], c(TPR,FPR),
                     names_to = "metric", values_to = "value")
long$metric <- factor(long$metric, levels = c("TPR","FPR"),
                      labels = c("True positive rate","False positive rate"))

## Okabe-Ito (colour-blind safe; also separates in grayscale), plus
## distinct line types and shapes so the figure works in pure B&W.
cols <- c("Fusion"="#000000","Holm"="#0072B2","AIC"="#D55E00","BIC"="#009E73",
          "LASSO"="#56B4E9","Elastic Net"="#CC79A7","Relaxed LASSO"="#666666","SCAD"="#E69F00")
ltys <- c("Fusion"="solid","Holm"="dashed","AIC"="dotted","BIC"="dotdash",
          "LASSO"="longdash","Elastic Net"="twodash","Relaxed LASSO"="2262","SCAD"="4212")
shps <- c("Fusion"=16,"Holm"=1,"AIC"=17,"BIC"=15,"LASSO"=4,"Elastic Net"=8,
          "Relaxed LASSO"=6,"SCAD"=18)

p <- ggplot(long, aes(n, value, colour = method, linetype = method, shape = method)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 2.1, stroke = 0.7) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_x_log10(breaks = c(50,100,200,500,1000,2000)) +
  scale_colour_manual(values = cols) +
  scale_linetype_manual(values = ltys) +
  scale_shape_manual(values = shps) +
  labs(x = "Sample size n (log scale)", y = NULL,
       colour = NULL, linetype = NULL, shape = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", legend.key.width = unit(1.4, "cm"))

ggsave("plot_sim_vs_n_rho07.pdf", p, width = 9, height = 4.3)
cat("wrote plot_sim_vs_n_rho07.pdf (B&W-safe: colour + linetype + shape)\n")
