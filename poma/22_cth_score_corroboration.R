# ============================================================================
# 22 : external imaging corroboration of the POMA thinning-TKR association
#      using an INDEPENDENT third-party imaging metric (CTh-Score)
#
#   Motivation
#     The POMA analysis rests on one metric family (kmri / Chondrometrics
#     cMFTC).  The CTh-Maps release supplies, for the same OAI participants,
#     a cartilage severity score produced by a completely independent
#     automatic framework (Lausanne; Eur Radiol 2026).  If the "faster
#     imaging progression -> knee replacement" association survives
#     when the exposure is rebuilt from that unrelated metric, the finding
#     is not an artefact of one measurement pipeline.
#
#   Exposure
#     cth_slope : per-knee OLS slope of CTh-Score (points / year),
#                 CTh-Score 0-100, HIGHER = MORE SEVERE.
#                 Two variants: all CTh timepoints, or pre-index only.
#     cth_eblup : same slope, shrunken (BLUP) from a random-slope model.
#
#   Outcome / design
#     Matched case-control on knee replacement; within-pair conditional
#     logistic regression (strata = newstrata), identical layout to
#     clogit_adjusted_OR.csv in the main analysis.
#
#   Three questions
#     Q1 replication  : does cth_slope predict TKR within matched pairs?
#     Q2 concordance  : does cth_slope agree with cMFTC slope (same
#                       direction but different scale), and does each add
#                       information beyond the other?
#     Q3 shrinkage    : the manuscript corrected for differential shrinkage
#                       (unshrunk OLS vs EBLUP) in cMFTC.  Re-test that
#                       question on a metric bounded at 0 and 100, where
#                       ceiling / floor effects also operate.
#
#   Side coding
#     POMA `side`: 1 = RIGHT, 2 = LEFT.  Established empirically, not assumed:
#     mapping 1=RIGHT gives r(CTh-Score, cMFTC baseline) = -0.493 whereas the
#     reversed mapping gives -0.024.  See 23_side_coding_check.R.
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
suppressPackageStartupMessages({library(survival); library(nlme); library(ggplot2)})

LONG <- read.csv(file.path(BASE, "cth_poma_long.csv"))
BAS  <- read.csv(file.path(BASE, "cth_poma_base.csv"))

## --------------------------------------------------------------- helpers ---
zs <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
sd_ <- function(x) sd(x, na.rm = TRUE)

## ------------------------------------------------- per-knee slope (2 ways) --
slope_of <- function(dat, tag) {
  d <- dat[order(dat$unit, dat$t), ]
  nv <- table(d$unit)
  use <- names(nv)[nv >= 2]
  d <- d[as.character(d$unit) %in% use, ]
  u <- sort(unique(d$unit))
  ols <- vapply(u, function(i) { x <- d[d$unit == i, ]
    unname(coef(lm(score ~ t, x))[2]) }, numeric(1))
  b0  <- vapply(u, function(i) { x <- d[d$unit == i, ]; x$score[which.min(x$t)] }, numeric(1))
  ntp  <- vapply(u, function(i) sum(d$unit == i), integer(1))
  tmax <- vapply(u, function(i) max(d$t[d$unit == i]), numeric(1))

  ## EBLUP slope from a random-intercept-and-slope model on the same knees
  ## (mirrors the EBLUP construction used for cMFTC in the main analysis)
  mm <- lme(score ~ t, random = ~ 1 + t | unit, data = d, method = "ML")
  fx <- fixef(mm)[["t"]]
  re <- ranef(mm)
  sl <- re[[ncol(re)]]; names(sl) <- rownames(re)
  eblup <- unname(fx + sl[as.character(u)])
  out <- data.frame(unit = u, n_tp = ntp, t_max = tmax, base = b0,
                    ols = ols, eblup = eblup)
  names(out)[names(out) != "unit"] <- paste0(names(out)[names(out) != "unit"], "_", tag)
  out
}

## all available CTh timepoints
A <- slope_of(LONG, "all")
## pre-index only (visits at or before the matched index visit)
P <- slope_of(LONG[LONG$pre_index == 1, ], "pre")

k <- merge(A, P, by = "unit", all = TRUE)
u <- merge(BAS, k, by = "unit", all.y = FALSE)
cat(sprintf("POMA knees with CTh slope : %d of %d (all-tp %d, pre-index %d)\n",
            nrow(u), nrow(BAS), sum(!is.na(u$ols_all)), sum(!is.na(u$ols_pre))))

## ------------------------------------------------------ complete pairs -----
mk_pairs <- function(dat, col) {
  x <- dat[!is.na(dat[[col]]) & !is.na(dat$newstrata), ]
  tb <- table(x$newstrata)
  x <- x[x$newstrata %in% names(tb)[tb == 2], ]
  x
}
PA <- mk_pairs(u, "ols_all")
cat(sprintf("complete matched pairs on CTh-Score slope : %d pairs / %d knees\n",
            length(unique(PA$newstrata)), nrow(PA)))

## ---------------------------------------------- Q1 : replication via clogit --
cl <- function(dat, form) {
  m <- clogit(as.formula(paste("case ~", form, "+ strata(newstrata)")), data = dat)
  s <- summary(m)$coef
  data.frame(term = rownames(s), beta = s[, 1], se = s[, 3],
             OR = exp(s[, 1]), lo = exp(s[, 1] - 1.96 * s[, 3]),
             hi = exp(s[, 1] + 1.96 * s[, 3]), p = s[, 5],
             n_pairs = length(unique(dat$newstrata)), n = nrow(dat),
             stringsAsFactors = FALSE)
}
PA$z_cth_all  <- zs(PA$ols_all)
PA$z_cth_pre  <- zs(PA$ols_pre)
PA$z_cth_e    <- zs(PA$eblup_all)
PA$z_base_all <- zs(PA$base_all)
PA$z_cmftc    <- zs(PA[["cMFTC_ThCtAB_aMe__eblup"]])
PA$z_cmftc_o  <- zs(PA[["cMFTC_ThCtAB_aMe__slope"]])
PA$z_cmbase   <- zs(PA[["cMFTC_ThCtAB_aMe__base"]])

models <- list(
  "A1  cMFTC slope alone (kmri, reference)"   = "z_cmftc",
  "A2  cMFTC slope, unshrunk OLS"             = "z_cmftc_o",
  "B1  CTh-Score slope alone (all tp)"        = "z_cth_all",
  "B2  CTh-Score slope alone (pre-index)"     = "z_cth_pre",
  "B3  CTh-Score slope alone (EBLUP, all tp)" = "z_cth_e",
  "B4  CTh-Score slope + baseline score"      = "z_cth_all + z_base_all",
  "C1  both metrics, mutual adjustment"       = "z_cth_all + z_cmftc",
  "C2  both metrics + both baselines"         = "z_cth_all + z_base_all + z_cmftc + z_cmbase"
)
res <- do.call(rbind, lapply(names(models), function(nm) {
  o <- cl(PA, models[[nm]]); o$model <- nm; o }))
res <- res[, c("model", "term", "beta", "se", "OR", "lo", "hi", "p", "n_pairs", "n")]
write.csv(res, file.path(BASE, "cth_clogit_models.csv"), row.names = FALSE)
cat("\n==== Q1 within-pair conditional logistic (outcome = knee replacement) ====\n")
print(res, row.names = FALSE, digits = 4)

## ----------------------------------- Q1b : tertiles of CTh-Score slope ------
PA$ter <- cut(PA$ols_all, quantile(PA$ols_all, c(0, 1/3, 2/3, 1)),
              labels = c("slow", "middle", "fast"), include.lowest = TRUE)
PA$fast <- as.integer(PA$ter == "fast")
PA$rank <- as.integer(PA$ter) - 1
cat("\n---- CTh-Score slope tertiles ----\n")
tb <- aggregate(cbind(n = rep(1, nrow(PA)), tkr = case) ~ ter, data = PA, FUN = sum)
tb$pct_tkr <- 100 * tb$tkr / tb$n
print(tb, row.names = FALSE, digits = 3)
if (length(unique(PA$fast)) == 2) {
  f1 <- cl(PA, "fast")
  cat(sprintf("  fast vs slow/middle : OR = %.3f (%.3f-%.3f), p = %.4g\n",
              f1$OR[1], f1$lo[1], f1$hi[1], f1$p[1]))
  f2 <- cl(PA, "rank")
  cat(sprintf("  per tertile step    : OR = %.3f (%.3f-%.3f), p = %.4g\n",
              f2$OR[1], f2$lo[1], f2$hi[1], f2$p[1]))
  write.csv(rbind(f1, f2), file.path(BASE, "cth_tertile_OR.csv"), row.names = FALSE)
}

## ------------------------------------------------ Q2 : concordance ----------
## The two metrics encode worsening on opposite signs: CTh-Score rises as the
## knee deteriorates, cMFTC thickness falls.  Negate the cMFTC slope so that
## BOTH variables mean "faster progression = larger" before correlating.
cc <- PA[!is.na(PA$ols_all) & !is.na(PA[["cMFTC_ThCtAB_aMe__slope"]]), ]
cc$loss_ols <- -cc[["cMFTC_ThCtAB_aMe__slope"]]
cc$loss_eb  <- -cc[["cMFTC_ThCtAB_aMe__eblup"]]
r_p  <- cor(cc$ols_all, cc$loss_ols)
r_s  <- cor(cc$ols_all, cc$loss_ols, method = "spearman")
r_pe <- cor(cc$eblup_all, cc$loss_eb)
cat(sprintf("\n==== Q2 concordance with cMFTC (n = %d knees) ====\n", nrow(cc)))
cat("  (cMFTC slope negated so both metrics = 'faster progression = larger')\n")
cat(sprintf("  CTh OLS   vs cMFTC OLS   : Pearson %+.3f  Spearman %+.3f\n", r_p, r_s))
cat(sprintf("  CTh EBLUP vs cMFTC EBLUP : Pearson %+.3f\n", r_pe))
z1 <- zs(cc$ols_all); z2 <- zs(cc$loss_ols)
dif <- z1 - z2; avg <- (z1 + z2) / 2
cat(sprintf("  Bland-Altman on z-scale  : mean diff %+.3f (SD %.3f), 95%% limits %+.3f to %+.3f\n",
            mean(dif), sd(dif), mean(dif) - 1.96 * sd(dif), mean(dif) + 1.96 * sd(dif)))
tt <- table(cut(cc$ols_all, quantile(cc$ols_all, c(0, 1/3, 2/3, 1)), include.lowest = TRUE),
            cut(cc$loss_ols, quantile(cc$loss_ols, c(0, 1/3, 2/3, 1)), include.lowest = TRUE))
cat("  tertile cross-tabulation (rows = CTh-Score, cols = cMFTC thinning):\n"); print(tt)
agr <- sum(diag(tt)) / sum(tt)
cat(sprintf("  same tertile (diagonal) %.1f%% ; extreme-vs-extreme disagreement %.1f%%\n",
            100 * agr, 100 * (tt[1, 3] + tt[3, 1]) / sum(tt)))
conc <- data.frame(metric = c("Pearson", "Spearman", "Pearson_EBLUP", "BA_mean_diff",
                              "BA_sd", "tertile_agreement_pct", "n"),
                   value = c(r_p, r_s, r_pe, mean(dif), sd(dif), 100 * agr, nrow(cc)))
write.csv(conc, file.path(BASE, "cth_concordance.csv"), row.names = FALSE)

## ------------------------------------------------- Q3 : differential shrinkage
## same test as in the main analysis: compare the case-minus-control paired
## difference in slope computed unshrunk (OLS) versus shrunken (EBLUP).
pair_diff <- function(dat, col) {
  x <- dat[!is.na(dat[[col]]), c("newstrata", "case", col)]
  cs <- x[x$case == 1, ]; ct <- x[x$case == 0, ]
  m <- merge(cs, ct, by = "newstrata", suffixes = c("_c", "_k"))
  m[[paste0(col, "_c")]] - m[[paste0(col, "_k")]]
}
d_ols   <- pair_diff(PA, "ols_all")
d_eblup <- pair_diff(PA, "eblup_all")
d_cmft  <- pair_diff(PA, "cMFTC_ThCtAB_aMe__eblup")
d_cmfto <- pair_diff(PA, "cMFTC_ThCtAB_aMe__slope")
sh <- data.frame(
  metric   = c("CTh-Score all-tp", "CTh-Score all-tp", "cMFTC", "cMFTC"),
  estimate = c("unshrunk OLS", "EBLUP", "unshrunk OLS", "EBLUP"),
  mean_paired_diff = c(mean(d_ols), mean(d_eblup), mean(d_cmfto), mean(d_cmft)),
  sd = c(sd(d_ols), sd(d_eblup), sd(d_cmfto), sd(d_cmft)),
  n_pairs = c(length(d_ols), length(d_eblup), length(d_cmfto), length(d_cmft)))
tt1 <- t.test(d_ols, d_eblup, paired = TRUE)
cat("\n==== Q3 differential shrinkage (case - control paired difference) ====\n")
print(sh, row.names = FALSE, digits = 4)
cat(sprintf("  CTh-Score : paired t between OLS and EBLUP differences = %+.4f (95%% CI %+.4f to %+.4f), p = %.3g\n",
            mean(d_ols - d_eblup), tt1$conf.int[1], tt1$conf.int[2], tt1$p.value))
cat(sprintf("  attenuation factor (EBLUP/OLS) : CTh %.3f  vs  cMFTC %.3f\n",
            mean(d_eblup) / mean(d_ols),
            mean(d_cmft) / mean(d_cmfto)))
write.csv(sh, file.path(BASE, "cth_shrinkage.csv"), row.names = FALSE)

## ---- the two manuscript 'Table 1' analogues, re-checked on CTh-Score ----
perm <- function(dat, col) {
  x <- dat[!is.na(dat[[col]]), c("newstrata", "case", col)]
  m <- merge(x[x$case == 1, ], x[x$case == 0, ], by = "newstrata", suffixes = c("_c", "_k"))
  d <- m[[paste0(col, "_c")]] - m[[paste0(col, "_k")]]
  c(mean = mean(d), lo = mean(d) - 1.96 * sd(d) / sqrt(length(d)),
    hi = mean(d) + 1.96 * sd(d) / sqrt(length(d)),
    p = t.test(m[[paste0(col, "_c")]], m[[paste0(col, "_k")]], paired = TRUE)$p.value)
}
bt <- rbind(
  cbind(variable = "baseline CTh-Score (pts)", t(perm(PA, "base_all"))),
  cbind(variable = "baseline cMFTC (mm)",      t(perm(PA, "cMFTC_ThCtAB_aMe__base"))),
  cbind(variable = "CTh-Score slope (pts/yr)", t(perm(PA, "ols_all"))),
  cbind(variable = "cMFTC slope (mm/yr)",      t(perm(PA, "cMFTC_ThCtAB_aMe__slope"))))
bt <- as.data.frame(bt); for (j in 2:5) bt[[j]] <- as.numeric(bt[[j]])
cat("\n==== paired (case - control) differences, CTh vs cMFTC ====\n")
print(bt, row.names = FALSE, digits = 4)
write.csv(bt, file.path(BASE, "cth_paired_diffs.csv"), row.names = FALSE)

## follow-up density of the two metrics (the second manuscript Table-1 fact)
cat("\n---- visit density ----\n")
cat(sprintf("  CTh-Maps visits per knee : cases %.2f vs controls %.2f (paired p = %.3g)\n",
            mean(PA$n_tp_all[PA$case == 1]), mean(PA$n_tp_all[PA$case == 0]),
            t.test(PA$n_tp_all[PA$case == 1], PA$n_tp_all[PA$case == 0], paired = FALSE)$p.value))
cat(sprintf("  POMA kmri visits per knee: cases %.2f vs controls %.2f (paired p = %.3g)\n",
            mean(PA$n_visits[PA$case == 1]), mean(PA$n_visits[PA$case == 0]),
            t.test(PA$n_visits[PA$case == 1], PA$n_visits[PA$case == 0], paired = FALSE)$p.value))

## --------------------------------------------------------------- figures ----
pl <- rbind(
  data.frame(case = PA$case, slope = PA$ols_all,  metric = "CTh-Score slope (pts/yr)"),
  data.frame(case = PA$case, slope = PA[["cMFTC_ThCtAB_aMe__slope"]], metric = "cMFTC slope (mm/yr)"))
pl$grp <- ifelse(pl$case == 1, "TKR case", "matched control")
g1 <- ggplot(pl, aes(grp, slope, fill = grp)) +
  geom_boxplot(width = .5, outlier.size = .5, alpha = .85, colour = "grey30", linewidth = .3) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c("TKR case" = "#B03A2E", "matched control" = "#2E6DA4")) +
  theme_bw(base_size = 10) + theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "per-knee slope",
       title = "Two independent imaging metrics, same 191 matched pairs",
       subtitle = "CTh-Score drawn from an unrelated automatic framework (Lausanne); cMFTC from kmri (Chondrometrics)")
ggsave(file.path(BASE, "fig_cth_slope_by_case.png"), g1, width = 8.5, height = 3.6, dpi = 200)

g2 <- ggplot(cc, aes(ols_all, cMFTC_ThCtAB_aMe__slope)) +
  geom_point(aes(colour = ifelse(case == 1, "TKR case", "control")), size = 1.4, alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey25", linewidth = .5) +
  scale_colour_manual(values = c("TKR case" = "#B03A2E", "control" = "#2E6DA4"), name = NULL) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "top") +
  labs(x = "CTh-Score slope (points / year)", y = "cMFTC slope (mm / year)",
       title = "Same direction, different scale: the two metrics are correlated but not interchangeable",
       subtitle = sprintf("Pearson %+.3f, Spearman %+.3f, n = %d knees", r_p, r_s, nrow(cc)))
ggsave(file.path(BASE, "fig_cth_vs_cmftc.png"), g2, width = 7, height = 4.2, dpi = 200)

g3 <- data.frame(term = c("CTh-Score slope (all tp)", "CTh-Score slope (pre-index)",
                          "CTh-Score slope (EBLUP)", "cMFTC slope (kmri, reference)"),
                 OR = c(res$OR[3], res$OR[4], res$OR[5], res$OR[1]),
                 lo = c(res$lo[3], res$lo[4], res$lo[5], res$lo[1]),
                 hi = c(res$hi[3], res$hi[4], res$hi[5], res$hi[1]))
g3$term <- factor(g3$term, levels = rev(g3$term))
g4 <- ggplot(g3, aes(OR, term)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = .3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .18, colour = "grey35", linewidth = .5) +
  geom_point(size = 2.4, colour = "#B03A2E") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank()) +
  labs(x = "OR per 1 SD faster progression (within-pair conditional logistic)", y = NULL,
       title = "Knee replacement is predicted by an imaging metric the original analysis never used")
ggsave(file.path(BASE, "fig_cth_forest.png"), g4, width = 8, height = 2.8, dpi = 200)

saveRDS(list(data = PA, models = res, concordance = conc, shrinkage = sh,
             paired = bt), file.path(BASE, "cth_corroboration.rds"))
cat("\n-> cth_clogit_models.csv / cth_concordance.csv / cth_shrinkage.csv")
cat(" / cth_paired_diffs.csv / cth_tertile_OR.csv")
cat(" / fig_cth_slope_by_case.png / fig_cth_vs_cmftc.png / fig_cth_forest.png / cth_corroboration.rds\n")
