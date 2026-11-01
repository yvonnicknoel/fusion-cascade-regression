## ============================================================
## The Fusion Cascade traces its own frontier: TPR vs FPR operating
## space, with the competitor methods as fixed points and the Fusion
## (and Holm) operating points as an alpha-indexed PATH through them.
##
## This materializes the Pareto-frontier claim that the manuscript
## previously made only in prose, and answers Reviewer 2's points 3-4:
## the procedure is a one-parameter family, not a single conservative
## point, and its dial carries a frequentist meaning.
##
## Inputs : simulation_overall_with_aic.csv (main study, fixed methods)
##          alpha_sweep_overall.csv         (Analysis B)
##          alpha_sweep_global_null.csv     (realized vs nominal alpha)
## Output : plot_alpha_frontier.pdf
##
## Black & white safe: shape and line type carry the information,
## no reliance on colour.
## ============================================================
suppressMessages({ library(ggplot2); library(dplyr); library(ggrepel); library(patchwork) })

fixed <- read.csv2("simulation_overall_with_aic.csv")
sweep <- read.csv2("alpha_sweep_overall.csv")
gn    <- read.csv2("alpha_sweep_global_null.csv")

nice <- function(x) dplyr::recode(x, ElasticNet = "Elastic Net", Lasso = "LASSO",
                                  RelaxedLasso = "Relaxed LASSO")
fixed$method <- nice(fixed$method)
sweep$method <- nice(sweep$method)

## Competitor points: everything except the two sweepable procedures,
## which are shown as paths instead.
comp <- fixed %>% filter(!method %in% c("Fusion", "Holm"))

fus  <- sweep %>% filter(method == "Fusion") %>% arrange(alpha)
holm <- sweep %>% filter(method == "Holm")   %>% arrange(alpha)

## ------------------------------------------------------------
## Panel A -- operating space
## ------------------------------------------------------------
pA <- ggplot() +
  ## competitor operating points (fixed tuning)
  geom_point(data = comp, aes(FPR, TPR), shape = 1, size = 2.2, colour = "grey30") +
  geom_text_repel(data = comp, aes(FPR, TPR, label = method), size = 2.9,
                  colour = "grey25", box.padding = 0.35, point.padding = 0.15,
                  min.segment.length = 0.3, segment.size = 0.25,
                  segment.colour = "grey65", max.overlaps = Inf, seed = 1) +
  ## Holm's alpha-path (the only competitor with a comparable dial)
  geom_path(data = holm, aes(FPR, TPR), linetype = "22", linewidth = 0.4,
            colour = "grey45") +
  geom_point(data = holm, aes(FPR, TPR), shape = 2, size = 1.7, colour = "grey45") +
  ## Fusion's alpha-path
  geom_path(data = fus, aes(FPR, TPR), linewidth = 0.7, colour = "black") +
  geom_point(data = fus, aes(FPR, TPR), shape = 16, size = 2.4, colour = "black") +
  geom_text_repel(data = fus, aes(FPR, TPR, label = sprintf("%.2f", alpha)),
                  size = 2.7, fontface = "bold", box.padding = 0.42,
                  nudge_x = -0.012, direction = "y",
                  min.segment.length = 0.25, segment.size = 0.25,
                  segment.colour = "grey60", max.overlaps = Inf, seed = 2) +
  labs(x = "False positive rate", y = "True positive rate",
       title = "A  Operating space") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 10, face = "bold"))

## ------------------------------------------------------------
## Panel B -- the dial means what it says: realized global-null
## selection rate vs nominal alpha (s = 0 scenarios)
## ------------------------------------------------------------
gp <- gn %>% filter(scope == "pooled")
gp$method <- nice(gp$method)

pB <- ggplot(gp, aes(alpha, realized, shape = method, linetype = method)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2.2) +
  scale_shape_manual(values = c(Fusion = 16, Holm = 2)) +
  scale_linetype_manual(values = c(Fusion = "solid", Holm = "22")) +
  annotate("text", x = 0.145, y = 0.185, label = "nominal", size = 2.9,
           colour = "grey45", angle = 38) +
  coord_equal(xlim = c(0, 0.22), ylim = c(0, 0.22)) +
  labs(x = expression("Nominal "*alpha),
       y = "Realized global-null selection rate",
       title = "B  Calibration of the dial") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),
        legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", colour = "grey80"),
        legend.key.size = unit(0.8, "lines"))

p <- pA + pB + plot_layout(widths = c(1.25, 1))
ggsave("plot_alpha_frontier.pdf", p, width = 9.0, height = 4.2)
cat("wrote plot_alpha_frontier.pdf\n\n")

cat("=== Fusion alpha-path (signal scenarios) ===\n")
print(as.data.frame(fus[c("alpha","TPR","FPR","FDR","F1","model_size","FWER_FP")]), digits = 3)
cat("\n=== Holm alpha-path ===\n")
print(as.data.frame(holm[c("alpha","TPR","FPR","FDR","F1","model_size","FWER_FP")]), digits = 3)
cat("\n=== Realized vs nominal at the global null ===\n")
print(as.data.frame(gp[c("method","alpha","realized")]), digits = 3)
