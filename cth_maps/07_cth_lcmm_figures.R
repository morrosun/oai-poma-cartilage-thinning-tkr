# ============================================================================
# 07 : P2 figures and the corrected class tables for the LCMM route
#
#   Replaces the stale artefacts left by the first (degenerate) LCMM run:
#     cth_lcmm_class_profile.csv  -> now the converged K = 3 classes
#     fig_cth_lcmm_trajectories.png
#     fig_cth_lcmm_classrate.png
#     cth_lcmm_class_vs_KL.csv
#   and adds
#     fig_cth_lcmm_model_selection.png   (BIC / AIC against K)
#
#   Inputs   cth_long.csv, zib_klinfo.csv, cth_lcmm_cache.rds, cth_lcmm_chain.rds
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/oai_imaging"
suppressPackageStartupMessages({library(lcmm); library(ggplot2)})

raw <- read.csv(file.path(BASE, "cth_long.csv"))
raw <- raw[!is.na(raw$t) & !is.na(raw$score), ]
nv <- table(raw$uid)
keep <- names(nv)[nv >= 2]
d0 <- raw[raw$uid %in% keep, ]
d0$unit <- as.integer(factor(d0$uid))
d0 <- d0[order(d0$unit, d0$t), ]
N <- max(d0$unit)
uid_of_unit <- d0$uid[match(seq_len(N), d0$unit)]

d3 <- readRDS(file.path(BASE, "cth_lcmm_k3_diag.rds"))
chain <- readRDS(file.path(BASE, "cth_lcmm_chain.rds"))
k3 <- d3$k3
sel <- chain$selection

## ------------------------------------------------ corrected class profile ---
mslope <- tapply(k3$ols, k3$cl, mean, na.rm = TRUE)
ord <- as.integer(names(sort(mslope)))                  # slowest -> fastest
prof <- do.call(rbind, lapply(ord, function(g) {
  x <- k3[k3$cl == g, ]
  data.frame(class = g, n = nrow(x), pct = 100 * nrow(x) / nrow(k3),
             baseline_mean = mean(x$b0), baseline_sd = sd(x$b0),
             slope_mean = mean(x$ols, na.rm = TRUE), slope_sd = sd(x$ols, na.rm = TRUE),
             visits = mean(x$nvis))
}))
write.csv(prof, file.path(BASE, "cth_lcmm_class_profile.csv"), row.names = FALSE)
cat("---- K = 3 class profile (corrected; replaces the degenerate single-class table) ----\n")
print(prof, row.names = FALSE, digits = 4)

## --------------------------------------------- external criterion: KL ------
klf <- read.csv(file.path(BASE, "zib_klinfo.csv"))
kk <- data.frame(unit = k3$unit, cl = k3$cl,
                 ID = as.integer(sub("_.*$", "", uid_of_unit[k3$unit])),
                 SIDE = sub("^.*_", "", uid_of_unit[k3$unit]))
mm <- merge(kk[kk$SIDE == "RIGHT", ], klf[, c("SubjectID", "KLGrade")],
            by.x = "ID", by.y = "SubjectID")
mm <- mm[!is.na(mm$KLGrade), ]
byc <- data.frame(cl = sort(unique(mm$cl)),
                  n = as.integer(table(mm$cl)),
                  mean_KL = as.numeric(tapply(mm$KLGrade, mm$cl, mean)),
                  sd_KL = as.numeric(tapply(mm$KLGrade, mm$cl, sd)))
byc <- merge(byc, prof[, c("class", "baseline_mean", "slope_mean")],
             by.x = "cl", by.y = "class")
write.csv(byc, file.path(BASE, "cth_lcmm_class_vs_KL.csv"), row.names = FALSE)
cat("\n---- baseline KLGrade by K = 3 class (OAIZIB right knees) ----\n")
print(byc, row.names = FALSE, digits = 3)
cat(sprintf("  Spearman (class mean KL vs class mean baseline) = %+.2f\n",
            cor(byc$mean_KL, byc$baseline_mean, method = "spearman")))
cat(sprintf("  Spearman (class mean KL vs class mean rate)     = %+.2f\n",
            cor(byc$mean_KL, byc$slope_mean, method = "spearman")))

## ------------------------------------------------------- figures ----------
p1 <- ggplot(sel, aes(K, BIC)) +
  geom_line(colour = "#1f4e79", linewidth = .9) +
  geom_point(colour = "#1f4e79", size = 2.4) +
  geom_text(aes(label = sprintf("%.1f%%", min_pct)), vjust = -1.1, size = 3) +
  scale_x_continuous(breaks = sel$K) +
  scale_y_continuous(expand = expansion(mult = c(.05, .12))) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank()) +
  labs(x = "number of latent classes K", y = "BIC",
       title = "No identifiable number of trajectory classes",
       subtitle = paste0("latent-class mixed model on CTh-Score, 8,278 knees; ",
                         "labels show the smallest class. BIC still falls at K = 6"))
ggsave(file.path(BASE, "fig_cth_lcmm_model_selection.png"), p1,
       width = 7.2, height = 4.0, dpi = 200)

dd <- merge(d0[, c("unit", "t", "score")],
            data.frame(unit = k3$unit, cl = k3$cl), by = "unit")
ms <- tapply(k3$ols, k3$cl, mean, na.rm = TRUE); o2 <- as.integer(names(sort(ms)))
lab <- sprintf("Class %d\n%+.2f pts/yr | base %.0f | n=%d (%.0f%%)",
               o2, ms[o2], tapply(k3$b0, k3$cl, mean)[o2],
               table(k3$cl)[o2], 100 * table(k3$cl)[o2] / nrow(k3))
dd$clf <- factor(dd$cl, levels = o2, labels = lab)
pm <- do.call(rbind, lapply(split(dd, dd$clf), function(x) {
  f <- lm(score ~ t, data = x); tt <- seq(0, 8, by = .25)
  data.frame(clf = x$clf[1], t = tt, score = predict(f, data.frame(t = tt)))
}))
p2 <- ggplot(dd, aes(t, score, group = unit)) +
  geom_line(colour = "grey82", linewidth = .18, alpha = .55) +
  geom_line(data = pm, aes(t, score, group = clf), colour = "#B03A2E",
            linewidth = 1.5, inherit.aes = FALSE) +
  facet_wrap(~ clf, nrow = 1) + theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 8), panel.grid.minor = element_blank()) +
  labs(x = "years from baseline", y = "CTh-Score (0-100, higher = worse)",
       title = "K = 3 latent classes: severity strata with unequal rates, not rate classes",
       subtitle = sprintf(paste0("the classes account for %.0f%% of the variance in the ",
                                 "severity LEVEL but only %.0f%% of the variance in the RATE"),
                          100 * 0.6366, 100 * 0.3693))
ggsave(file.path(BASE, "fig_cth_lcmm_trajectories.png"), p2,
       width = 11.5, height = 3.9, dpi = 200)

p3 <- ggplot(prof, aes(factor(class, levels = ord), slope_mean)) +
  geom_col(fill = "#2E6DA4", width = .6) +
  geom_errorbar(aes(ymin = slope_mean - slope_sd, ymax = slope_mean + slope_sd),
                width = .15, colour = "grey30") +
  geom_text(aes(label = sprintf("%+.2f (%.0f%%; n=%d)", slope_mean, pct, n)),
            vjust = -2.6, size = 3) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(.05, .22))) +
  labs(x = "K = 3 class (slowest to fastest)", y = "class mean rate (points / year)",
       title = "Class rates differ, but within-class spread is as large as the gap",
       subtitle = "bars = class mean, error bars = +/- 1 SD of the individual per-knee rates")
ggsave(file.path(BASE, "fig_cth_lcmm_classrate.png"), p3,
       width = 7.6, height = 3.9, dpi = 200)

cat("\n-> cth_lcmm_class_profile.csv / cth_lcmm_class_vs_KL.csv")
cat(" / fig_cth_lcmm_model_selection.png / fig_cth_lcmm_trajectories.png")
cat(" / fig_cth_lcmm_classrate.png\n")
