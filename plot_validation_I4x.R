## ============================================================
## Validation figure for the extended Model I4x
## Points = Monte Carlo calibrated c(m;n,p) (rho=0)
## Lines  = Model I4x prediction
## Faceted by n, coloured by p.  -> plot_validation_I4x.pdf
## Run: Rscript plot_validation_I4x.R
## ============================================================
suppressMessages({ library(ggplot2); library(dplyr) })
source("get_critical_value.R")

d <- read.csv("c_table_ext_full.csv", stringsAsFactors = FALSE)
d <- d[d$rho == 0, ]

## model prediction at every calibrated (m,n,p)
d$pred <- unlist(Map(function(m, n, p) get_critical_value(m, n, p),
                     d$m, d$n, d$p))

nlab <- function(n) factor(n, levels = sort(unique(n)), labels = paste0("n = ", sort(unique(n))))
d$np <- nlab(d$n)

## Ordinal p encoded by a viridis gradient: reproduces in grayscale as an
## ordered dark-to-light ramp (Psychological Methods prints in black & white).
pl <- ggplot(d, aes(m, c, colour = factor(p))) +
  geom_point(size = 0.35, alpha = 0.6) +
  geom_line(aes(y = pred, group = p), linewidth = 0.4) +
  facet_wrap(~ np, scales = "free_x", ncol = 3) +
  scale_colour_viridis_d(option = "viridis", end = 0.92) +
  labs(x = "Active predictors (m)", y = "Critical value  c(m; n, p)",
       colour = "p") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom")

ggsave("plot_validation_I4x.pdf", pl, width = 10, height = 9)

## quick fit summary printed alongside
rel <- 100 * abs(d$c - d$pred) / d$c
cat(sprintf("Validation over %d calibrated points: median rel = %.2f%%, 90th = %.2f%%, max = %.2f%%\n",
            nrow(d), median(rel), quantile(rel, .9), max(rel)))
cat("Wrote plot_validation_I4x.pdf\n")
