## ============================================================
## Operating-point figure for the SAPA health application:
## model size vs out-of-sample R^2 (mean over 50 splits, +/- sd).
## Black & white safe: black points, distinct shapes, text labels
## (no reliance on colour). -> plot_spi_operating_point.pdf
## ============================================================
suppressMessages({ library(ggplot2); library(dplyr); library(ggrepel) })
res <- readRDS("spi_cv_results.rds")
s <- res$perf %>% group_by(method) %>%
  summarise(size = mean(size), size_sd = sd(size),
            r2 = mean(r2), r2_sd = sd(r2), .groups = "drop")
s$method <- recode(s$method, ElasticNet = "Elastic Net", Lasso = "LASSO",
                   RelaxedLasso = "Relaxed LASSO")

## Labels placed by ggrepel: kept close to their dots, nudged apart only as
## much as needed to avoid overlap, with a short connector when moved.
p <- ggplot(s, aes(size, r2)) +
  geom_errorbar(aes(ymin = r2 - r2_sd, ymax = r2 + r2_sd), width = 0, colour = "grey55") +
  geom_errorbarh(aes(xmin = size - size_sd, xmax = size + size_sd), height = 0, colour = "grey55") +
  geom_point(size = 1.6, colour = "black") +
  geom_text_repel(aes(label = method), size = 3.1,
                  box.padding = 0.35, point.padding = 0.15,
                  min.segment.length = 0.3, segment.size = 0.3,
                  segment.colour = "grey60", max.overlaps = Inf, seed = 1) +
  labs(x = "Model size (number of selected items)",
       y = expression("Out-of-sample "*R^2)) +
  coord_cartesian(xlim = c(0, 62), ylim = c(0, 0.30)) +
  theme_bw(base_size = 11)

ggsave("plot_spi_operating_point.pdf", p, width = 6.2, height = 4.4)
cat("wrote plot_spi_operating_point.pdf\n")
print(as.data.frame(s[c("method","size","r2")]), digits = 3)
