## ============================================================
## Decision-mechanics figure: the effective |t| keep-threshold each
## selection rule imposes, shown analytically.
##   Panel A -- threshold vs sample size n (fixed p): scaling in n.
##   Panel B -- threshold vs family size m (fixed n,p): multiplicity
##              adaptivity along the cascade.
## Message: Fusion applies a multiplicity-aware, joint-null-calibrated
## threshold (adaptive in m, ~flat in n), unlike the penalty-anchored
## AIC/BIC (blind to m) or the worst-case union bound (Bonferroni/Holm).
## Black & white safe: line type + point shape + greyscale.
## Run: Rscript plot_decision_mechanics.R -> plot_decision_mechanics.pdf
## ============================================================
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork) })
source("get_critical_value.R")

alpha <- 0.05
z2 <- function(k) qnorm(1 - alpha / (2 * k))   # two-sided Bonferroni bar for k tests

meth_lv <- c("Fusion", "Bonferroni-Holm", "BIC", "AIC")
leg_order <- c("Bonferroni-Holm", "Fusion", "BIC", "AIC")   # longest label first
cols <- c("Fusion"="#000000", "Bonferroni-Holm"="#0072B2", "BIC"="#009E73", "AIC"="#D55E00")
ltys <- c("Fusion"="solid", "Bonferroni-Holm"="dashed", "BIC"="dotdash", "AIC"="dotted")
shps <- c("Fusion"=16, "Bonferroni-Holm"=6, "BIC"=15, "AIC"=17)

## ---------- Panel A: threshold vs n (fixed p) ----------
pA <- 20
nA <- c(50, 75, 100, 150, 200, 300, 500, 750, 1000, 1500, 2000)
fus <- lapply(nA, function(n) {
  cc <- get_critical_value(1:pA, n = n, p = pA, alpha = alpha)
  data.frame(n = n, lo = min(cc), hi = max(cc))
}) %>% bind_rows()

linesA <- bind_rows(
  data.frame(n = nA, value = sqrt(log(nA)),        method = "BIC"),
  data.frame(n = nA, value = sqrt(2),              method = "AIC"),
  data.frame(n = nA, value = z2(pA),               method = "Bonferroni-Holm")
)
linesA$method <- factor(linesA$method, levels = meth_lv)

gA <- ggplot() +
  geom_ribbon(data = fus, aes(n, ymin = lo, ymax = hi),
              fill = "grey80", colour = NA, alpha = 0.9) +
  geom_line(data = fus, aes(n, hi, colour = "Fusion", linetype = "Fusion"), linewidth = 0.6) +
  geom_line(data = fus, aes(n, lo, colour = "Fusion", linetype = "Fusion"), linewidth = 0.6) +
  geom_point(data = fus, aes(n, hi, colour = "Fusion", shape = "Fusion"), size = 1.9) +
  geom_line(data = linesA, aes(n, value, colour = method, linetype = method), linewidth = 0.6) +
  geom_point(data = linesA, aes(n, value, colour = method, shape = method), size = 1.9) +
  annotate("text", x = 95, y = 2.55, label = "Fusion c(k)\n(band over k)", size = 2.7, hjust = 0, lineheight = 0.9) +
  scale_x_log10(breaks = c(50,100,200,500,1000,2000)) +
  scale_colour_manual(values = cols, breaks = leg_order) +
  scale_linetype_manual(values = ltys, breaks = leg_order) +
  scale_shape_manual(values = shps, breaks = leg_order) +
  labs(x = "Sample size n (log scale)", y = "Effective |t| keep-threshold",
       colour = NULL, linetype = NULL, shape = NULL,
       subtitle = "(A) Scaling in n  (p = 20)") +
  coord_cartesian(ylim = c(1.2, 3.4)) +
  theme_bw(base_size = 10) + theme(legend.position = "none")

## ---------- Panel B: threshold vs m (fixed n,p) ----------
nB <- 200; pB <- 20; mB <- 1:pB
fusB <- data.frame(m = mB, value = get_critical_value(mB, n = nB, p = pB, alpha = alpha),
                   method = "Fusion")
linesB <- bind_rows(
  fusB,
  data.frame(m = mB, value = z2(mB),           method = "Bonferroni-Holm"),
  data.frame(m = mB, value = sqrt(log(nB)),     method = "BIC"),
  data.frame(m = mB, value = sqrt(2),           method = "AIC")
)
linesB$method <- factor(linesB$method, levels = meth_lv)

gB <- ggplot(linesB, aes(m, value, colour = method, linetype = method, shape = method)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.9) +
  scale_colour_manual(values = cols, breaks = leg_order) +
  scale_linetype_manual(values = ltys, breaks = leg_order) +
  scale_shape_manual(values = shps, breaks = leg_order) +
  labs(x = "Family size k (jointly tested coefficients)", y = "Effective |t| keep-threshold",
       colour = NULL, linetype = NULL, shape = NULL,
       subtitle = "(B) Multiplicity adaptivity  (n = 200, p = 20)") +
  coord_cartesian(ylim = c(1.2, 3.4)) +
  theme_bw(base_size = 10) +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "white", colour = "grey80"),
        legend.key = element_blank(),
        legend.key.size = unit(0.9, "lines"),
        legend.text = element_text(size = 8),
        legend.margin = margin(2, 4, 2, 2))

g <- gA | gB

ggsave("plot_decision_mechanics.pdf", g, width = 9, height = 4.4)
cat("wrote plot_decision_mechanics.pdf\n")
cat("\nPanel A Fusion band at selected n:\n"); print(fus, row.names = FALSE)
cat("\nPanel B thresholds vs m (n=200,p=20):\n")
print(pivot_wider(linesB, names_from = method, values_from = value), n = 20)
