# ============================================================================
# 34 : supplementary figures requested during pre-submission review
#
#   Supplementary Figure S2 : how the analysis sets were assembled
#                             (reviewer point 13 - a flow diagram)
#   Supplementary Figure S3 : posterior dependence of the two association
#                             parameters and the posterior of the rate odds
#                             ratio (evidence for reviewer point 2f, which
#                             asked for support for the "banana-shaped"
#                             description and for the sampling claim)
#
#   Output : fig_flow_design.png , fig_posterior_dependence.png
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
source("D:/BaiduSyncdisk/OAI/Scripts/poma/jm_core.R")
suppressPackageStartupMessages({library(ggplot2); library(grid); library(patchwork)})

## ================================================== Supplementary Fig S2 ==
## every number below is taken from Sections 2.1, 2.4, 2.6 and 3.1 of the
## manuscript; nothing is recomputed here.
bx <- 1.05                       # x of the main column
ex <- 4.15                       # x of the exclusion box
box <- data.frame(
  y = c(9.0, 7.7, 6.2, 6.2, 4.7, 3.4, 2.1),
  x = c(bx, bx, bx, ex, bx, bx, bx),
  col = c("#2E6DA4", "#2E6DA4", "#1B7A3E", "#B03A2E", "#1B7A3E", "#1B7A3E", "#7B4EA3"),
  txt = c(
    "OAI POMA knee-replacement nested\ncase-control release\n225 pairs / 450 knee units, matched\non the index visit",
    "Chondrometrics cartilage thickness,\n6 compartment metrics\npre-index visits only (0, 12, 24,\n36, 48 months)",
    "Primary analysis set (Set A)\nboth knees with >= 2 pre-operative\nvisits and complete covariates\n191 pairs / 382 knees / 1248 observations",
    "Excluded: 34 pairs\nat least one knee with < 2\npre-operative visits, or\nincomplete covariates",
    "Analyses on Set A\nconditional logistic ladder M0-M7\nwithin-pair conditional joint model\nthree formal compartment contrasts",
    "Trajectory set\n389 knees with >= 2 visits\nlatent-class mixed model, K = 1-6",
    "Independent automatic severity score\n450 knee units -> 391 linked (87%)\n-> 367 with >= 2 time points\n-> 148 pairs / 296 knees"),
  stringsAsFactors = FALSE)

seg <- function(x1, y1, x2, y2, col = "grey30", lty = 1)
  annotate("segment", x = x1, xend = x2, y = y1, yend = y2, colour = col,
           linetype = lty, linewidth = .7,
           arrow = arrow(length = unit(0.20, "cm")))

gp <- ggplot() +
  seg(bx, box$y[1] - 0.40, bx, box$y[2] + 0.40) +      # 1 -> 2
  seg(bx, box$y[2] - 0.40, bx, box$y[3] + 0.40) +      # 2 -> 3 (Set A)
  seg(bx + 1.35, box$y[3], ex - 1.15, box$y[4], col = "#B03A2E") +   # 3 -> excluded
  seg(bx, box$y[3] - 0.40, bx, box$y[5] + 0.40) +      # 3 -> analyses
  seg(bx, box$y[5] - 0.40, bx, box$y[6] + 0.40) +      # 5 -> 6
  seg(bx, box$y[6] - 0.40, bx, box$y[7] + 0.40) +      # 6 -> 7
  geom_label(data = box, aes(x, y, label = txt, fill = col),
             colour = "white", size = 2.75, lineheight = 1.15, label.size = 0,
             label.r = unit(0.10, "lines"), alpha = .93) +
  scale_fill_identity() +
  coord_cartesian(xlim = c(-0.7, 5.3), ylim = c(0, 10.2), clip = "off") +
  theme_void(base_size = 10) +
  labs(title = "Assembly of the analysis sets") +
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold",
                                  margin = margin(b = 8)))
ggsave(file.path(BASE, "fig_flow_design.png"), gp, width = 8.2, height = 8.8, dpi = 300)

## ================================================== Supplementary Fig S3 ==
L <- read.csv(file.path(BASE, "poma_long.csv"))
U <- read.csv(file.path(BASE, "poma_analysis_long.csv"))
COVARS <- c("V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")
PRIM <- "cMFTC_ThCtAB_aMe"
sub <- L[!is.na(L[[PRIM]]), ]
lg  <- data.frame(unit = sub$unit, year = sub$years, y = sub[[PRIM]],
                  newstrata = sub$newstrata, case = sub$case)
sv  <- U[, c("unit", "newstrata", "case", "t_index_months", COVARS)]
sv  <- sv[sv$unit %in% lg$unit, ]
sv  <- sv[complete.cases(sv[, COVARS]), ]
lg  <- lg[lg$unit %in% sv$unit, ]
sv$t_idx <- sv$t_index_months / 12
PP  <- prune_pairs(lg, sv, 2)
dt  <- build(PP$lg, PP$sv)
f   <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 2026, verbose = FALSE)
d   <- f$draws
rho <- cor(d$aV, d$aS)
cat(sprintf("posterior correlation aV vs aS (this chain) = %.3f\n", rho))
## prefer the pooled value over the three dispersed chains when 32 has been run
corfile <- file.path(BASE, "mcmc_posterior_correlation.csv")
rho_lab <- if (file.exists(corfile)) {
  pp <- read.csv(corfile)
  sprintf("posterior correlation %.2f\n(pooled over three chains)",
          pp$cor_aV_aS[pp$chain == "pooled"])
} else sprintf("posterior correlation %.2f\n(this chain)", rho)

set.seed(9)
idx <- sample(seq_along(d$aS), 1500)
sc <- ggplot(data.frame(aV = d$aV[idx], aS = d$aS[idx]), aes(aV, aS)) +
  geom_point(colour = "#2E6DA4", size = .7, alpha = .35) +
  geom_density_2d(colour = "#B03A2E", linewidth = .35, bins = 6) +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey50") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey50") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 3.0,
           label = rho_lab) +
  theme_bw(base_size = 10) +
  labs(x = expression(a[V] ~ "(per mm latent level)"),
       y = expression(a[S] ~ "(per mm/yr latent thinning)"),
       title = "Joint posterior of the two association parameters")

dfd <- data.frame(OR = exp(-d$aS * d$SDslo))
dens <- ggplot(dfd, aes(OR)) +
  geom_density(fill = "#2E6DA4", alpha = .35, colour = "#1B3F66") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = median(dfd$OR), colour = "#B03A2E") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.4, size = 3.1,
           label = sprintf("median %.2f\n95%% CrI %.2f-%.2f", median(dfd$OR),
                           quantile(dfd$OR, .025), quantile(dfd$OR, .975))) +
  scale_x_log10(breaks = c(1, 2, 5, 10), limits = c(1, 12)) +
  theme_bw(base_size = 10) +
  labs(x = "OR per 1 SD faster thinning", y = "posterior density",
       title = "Sampled posterior of the rate odds ratio")

ggsave(file.path(BASE, "fig_posterior_dependence.png"), sc + dens,
       width = 9, height = 3.9, dpi = 300)
cat("-> fig_flow_design.png / fig_posterior_dependence.png\n")
