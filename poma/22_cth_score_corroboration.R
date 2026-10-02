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
#     cth_eblup : same slope, shrunken (BLUP) from a random-slope model.
#
#   PRIMARY vs SENSITIVITY TRUNCATION  (revised 2026-10-02)
#     The PRIMARY exposure is the PRE-INDEX slope: visits at or before the
#     matched index visit, i.e. before the knee replacement for every case and
#     before the index date for every matched control.  The same rule is applied
#     to cases and controls; no case contributes a post-operative measurement.
#     An ALL-TIMEPOINT slope uses measurements taken after the outcome and
#     therefore cannot be described as a pre-operative or prognostic exposure.
#     It is retained for reference only and is labelled SENSITIVITY throughout.
#
#   Outcome / design
#     Matched case-control on knee replacement; within-pair conditional
#     logistic regression (strata = newstrata), identical layout to
#     clogit_adjusted_OR.csv in the main analysis.
#
#   Three questions
#     Q1 reproduction : does cth_slope predict TKR within matched pairs?
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
#
#   Wording: this is CROSS-PIPELINE CORROBORATION, not independent external
#     validation.  The automatic CTh-Score pipeline is independent of kmri, but
#     the knees are the SAME OAI participants, so nothing here validates the
#     model in a new population, and nothing here is a prospective prediction.
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

## pre-index only = PRIMARY (visits at or before the matched index visit)
P <- slope_of(LONG[LONG$pre_index == 1, ], "pre")
## all available CTh timepoints = SENSITIVITY ONLY (may include post-outcome visits)
A <- slope_of(LONG, "all")

k <- merge(A, P, by = "unit", all = TRUE)
u <- merge(BAS, k, by = "unit", all.y = FALSE)
cat(sprintf("POMA knees with CTh slope : %d of %d (pre-index PRIMARY %d, all-tp %d)\n",
            nrow(u), nrow(BAS), sum(!is.na(u$ols_pre)), sum(!is.na(u$ols_all))))

## ------------------------------------------------------ complete pairs -----
mk_pairs <- function(dat, col) {
  x <- dat[!is.na(dat[[col]]) & !is.na(dat$newstrata), ]
  tb <- table(x$newstrata)
  x <- x[x$newstrata %in% names(tb)[tb == 2], ]
  x
}
PA <- mk_pairs(u, "ols_pre")
cat(sprintf("complete matched pairs on the PRE-INDEX CTh-Score slope : %d pairs / %d knees\n",
            length(unique(PA$newstrata)), nrow(PA)))
cat(sprintf("  (all-timepoint set for comparison : %d pairs)\n",
            length(unique(mk_pairs(u, "ols_all")$newstrata))))

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
PA$z_cth_pre   <- zs(PA$ols_pre)        # PRIMARY
PA$z_cth_pre_e <- zs(PA$eblup_pre)      # PRIMARY, shrunken
PA$z_base_pre  <- zs(PA$base_pre)
PA$z_cth_all   <- zs(PA$ols_all)        # SENSITIVITY
PA$z_cth_e     <- zs(PA$eblup_all)      # SENSITIVITY
PA$z_base_all  <- zs(PA$base_all)       # SENSITIVITY
PA$z_cmftc    <- zs(PA[["cMFTC_ThCtAB_aMe__eblup"]])
PA$z_cmftc_o  <- zs(PA[["cMFTC_ThCtAB_aMe__slope"]])
PA$z_cmbase   <- zs(PA[["cMFTC_ThCtAB_aMe__base"]])

## NB "z_cmftc" is the NEGATIVE of the cMFTC thinning slope (a fall in thickness
## is a rise in the variable), so its OR is below 1; every CTh-Score slope OR is
## above 1.  Reporting both directions in one table is deliberate and is stated
## in the legend.
models <- list(
  "A1  cMFTC slope alone (kmri, reference)"          = "z_cmftc",
  "A2  cMFTC slope, unshrunk OLS"                    = "z_cmftc_o",
  "B1  CTh-Score slope, PRE-INDEX, OLS   [primary]"  = "z_cth_pre",
  "B2  CTh-Score slope, PRE-INDEX, EBLUP [primary]"  = "z_cth_pre_e",
  "B3  CTh-Score slope + baseline (pre-index)"       = "z_cth_pre + z_base_pre",
  "C1  both metrics, mutual adjustment (pre-index)"  = "z_cth_pre + z_cmftc",
  "C2  both metrics + both baselines (pre-index)"    = "z_cth_pre + z_base_pre + z_cmftc + z_cmbase",
  "S1  SENS cMFTC slope + baseline"                  = "z_cmftc + z_cmbase",
  "S2  SENS CTh-Score slope, ALL TIMEPOINTS"         = "z_cth_all",
  "S3  SENS CTh-Score slope, all tp, EBLUP"          = "z_cth_e",
  "S4  SENS CTh-Score + baseline, all tp"            = "z_cth_all + z_base_all",
  "S5  SENS both metrics, all tp"                    = "z_cth_all + z_cmftc",
  "S6  SENS both metrics + baselines, all tp"        = "z_cth_all + z_base_all + z_cmftc + z_cmbase"
)
res <- do.call(rbind, lapply(names(models), function(nm) {
  o <- cl(PA, models[[nm]]); o$model <- nm; o }))
res <- res[, c("model", "term", "beta", "se", "OR", "lo", "hi", "p", "n_pairs", "n")]
write.csv(res, file.path(BASE, "cth_clogit_models.csv"), row.names = FALSE)
cat("\n==== Q1 within-pair conditional logistic (outcome = knee replacement) ====\n")
cat("     OR > 1 = faster CTh-Score progression (or thinner cMFTC) -> TKR\n")
cat("     A1/A2 cMFTC rows are OR < 1 because z_cmftc is minus the thinning slope\n")
print(res, row.names = FALSE, digits = 4)
cat("\n---- headline (PRE-INDEX, primary) ----\n")
for (lab in c("B1  CTh-Score slope, PRE-INDEX, OLS   [primary]",
              "B2  CTh-Score slope, PRE-INDEX, EBLUP [primary]")) {
  rr <- res[res$model == lab, ]
  cat(sprintf("  %-48s OR %5.3f (%5.3f-%5.3f)  p = %.3g\n",
              lab, rr$OR[1], rr$lo[1], rr$hi[1], rr$p[1]))
}
cat("\n---- C1 mutual adjustment, the number quoted in the manuscript ----\n")
r1 <- res[res$model == "C1  both metrics, mutual adjustment (pre-index)", ]
for (k2 in seq_len(nrow(r1)))
  cat(sprintf("  %-10s OR %5.3f (%5.3f-%5.3f)  p = %.3g\n",
              r1$term[k2], r1$OR[k2], r1$lo[k2], r1$hi[k2], r1$p[k2]))
cat(sprintf("  cMFTC term inverted for the manuscript : %.3f (%.3f-%.3f)\n",
            1 / r1$OR[r1$term == "z_cmftc"], 1 / r1$hi[r1$term == "z_cmftc"],
            1 / r1$lo[r1$term == "z_cmftc"]))

## ----------------------------------- Q1b : tertiles of CTh-Score slope ------
PA$ter <- cut(PA$ols_pre, quantile(PA$ols_pre, c(0, 1/3, 2/3, 1)),
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
cc <- PA[!is.na(PA$ols_pre) & !is.na(PA[["cMFTC_ThCtAB_aMe__slope"]]), ]
cc$loss_ols <- -cc[["cMFTC_ThCtAB_aMe__slope"]]
cc$loss_eb  <- -cc[["cMFTC_ThCtAB_aMe__eblup"]]
r_p  <- cor(cc$ols_pre, cc$loss_ols)
r_s  <- cor(cc$ols_pre, cc$loss_ols, method = "spearman")
r_pe <- cor(cc$eblup_pre, cc$loss_eb)
## the all-timepoint version is kept as a sensitivity figure
r_p_all <- cor(cc$ols_all, cc$loss_ols)
cat(sprintf("\n==== Q2 concordance with cMFTC (n = %d knees) ====\n", nrow(cc)))
cat("  (cMFTC slope negated so both metrics = 'faster progression = larger')\n")
cat(sprintf("  CTh PRE-INDEX OLS vs cMFTC OLS : Pearson %+.3f  Spearman %+.3f\n", r_p, r_s))
cat(sprintf("  CTh PRE-INDEX EBLUP vs cMFTC   : Pearson %+.3f\n", r_pe))
cat(sprintf("  [SENS] CTh all-tp OLS vs cMFTC : Pearson %+.3f\n", r_p_all))
z1 <- zs(cc$ols_pre); z2 <- zs(cc$loss_ols)
dif <- z1 - z2; avg <- (z1 + z2) / 2
cat(sprintf("  Bland-Altman on z-scale  : mean diff %+.3f (SD %.3f), 95%% limits %+.3f to %+.3f\n",
            mean(dif), sd(dif), mean(dif) - 1.96 * sd(dif), mean(dif) + 1.96 * sd(dif)))
tt <- table(cut(cc$ols_pre, quantile(cc$ols_pre, c(0, 1/3, 2/3, 1)), include.lowest = TRUE),
            cut(cc$loss_ols, quantile(cc$loss_ols, c(0, 1/3, 2/3, 1)), include.lowest = TRUE))
cat("  tertile cross-tabulation (rows = CTh-Score, cols = cMFTC thinning):\n"); print(tt)
agr <- sum(diag(tt)) / sum(tt)
cat(sprintf("  same tertile (diagonal) %.1f%% ; extreme-vs-extreme disagreement %.1f%%\n",
            100 * agr, 100 * (tt[1, 3] + tt[3, 1]) / sum(tt)))
conc <- data.frame(metric = c("Pearson_pre", "Spearman_pre", "Pearson_EBLUP_pre",
                              "Pearson_alltp_SENS",
                              "BA_mean_diff", "BA_sd", "tertile_agreement_pct", "n"),
                   value = c(r_p, r_s, r_pe, r_p_all,
                             mean(dif), sd(dif), 100 * agr, nrow(cc)))
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
d_ols   <- pair_diff(PA, "ols_pre")
d_eblup <- pair_diff(PA, "eblup_pre")
d_cmft  <- pair_diff(PA, "cMFTC_ThCtAB_aMe__eblup")
d_cmfto <- pair_diff(PA, "cMFTC_ThCtAB_aMe__slope")
sh <- data.frame(
  metric   = c("CTh-Score pre-index", "CTh-Score pre-index", "cMFTC", "cMFTC"),
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
  cbind(variable = "baseline CTh-Score, pre-index (pts)", t(perm(PA, "base_pre"))),
  cbind(variable = "baseline cMFTC (mm)",      t(perm(PA, "cMFTC_ThCtAB_aMe__base"))),
  cbind(variable = "CTh-Score slope, pre-index (pts/yr)", t(perm(PA, "ols_pre"))),
  cbind(variable = "cMFTC slope (mm/yr)",      t(perm(PA, "cMFTC_ThCtAB_aMe__slope"))))
bt <- as.data.frame(bt); for (j in 2:5) bt[[j]] <- as.numeric(bt[[j]])
cat("\n==== paired (case - control) differences, CTh vs cMFTC ====\n")
print(bt, row.names = FALSE, digits = 4)
write.csv(bt, file.path(BASE, "cth_paired_diffs.csv"), row.names = FALSE)

## follow-up density of the two metrics (the second manuscript Table-1 fact)
cat("\n---- visit density ----\n")
cat(sprintf("  CTh-Maps PRE-INDEX visits/knee : cases %.2f vs controls %.2f (p = %.3g)\n",
            mean(PA$n_tp_pre[PA$case == 1]), mean(PA$n_tp_pre[PA$case == 0]),
            t.test(PA$n_tp_pre[PA$case == 1], PA$n_tp_pre[PA$case == 0], paired = FALSE)$p.value))
cat(sprintf("  CTh-Maps visits/knee, all tp   : cases %.2f vs controls %.2f (p = %.3g)  [SENS]\n",
            mean(PA$n_tp_all[PA$case == 1]), mean(PA$n_tp_all[PA$case == 0]),
            t.test(PA$n_tp_all[PA$case == 1], PA$n_tp_all[PA$case == 0], paired = FALSE)$p.value))
cat(sprintf("  POMA kmri visits per knee      : cases %.2f vs controls %.2f (paired p = %.3g)\n",
            mean(PA$n_visits[PA$case == 1]), mean(PA$n_visits[PA$case == 0]),
            t.test(PA$n_visits[PA$case == 1], PA$n_visits[PA$case == 0], paired = FALSE)$p.value))
cat("  NOTE: the all-timepoint count is asymmetric because controls keep being")
cat("  imaged after the index date while cases leave the imaging stream at TKR.")
cat("  The PRE-INDEX count is the one that matters for the primary exposure and")
cat("  is reported above; a residual difference in pre-index dose is a limitation.")
write.csv(data.frame(
  metric = c("CTh pre-index", "CTh all-tp", "kmri POMA"),
  visits_case = c(mean(PA$n_tp_pre[PA$case == 1]), mean(PA$n_tp_all[PA$case == 1]),
                  mean(PA$n_visits[PA$case == 1])),
  visits_control = c(mean(PA$n_tp_pre[PA$case == 0]), mean(PA$n_tp_all[PA$case == 0]),
                     mean(PA$n_visits[PA$case == 0]))),
  file.path(BASE, "cth_visit_density.csv"), row.names = FALSE)

## --------------------------------------------------------------- figures ----
pl <- rbind(
  data.frame(case = PA$case, slope = PA$ols_pre,  metric = "CTh-Score slope, pre-index (pts/yr)"),
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

g2 <- ggplot(cc, aes(ols_pre, cMFTC_ThCtAB_aMe__slope)) +
  geom_point(aes(colour = ifelse(case == 1, "TKR case", "control")), size = 1.4, alpha = .8) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey25", linewidth = .5) +
  scale_colour_manual(values = c("TKR case" = "#B03A2E", "control" = "#2E6DA4"), name = NULL) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "top") +
  labs(x = "CTh-Score slope, pre-index (points / year)", y = "cMFTC slope (mm / year)",
       title = "Same direction, different scale: the two metrics are correlated but not interchangeable",
       subtitle = sprintf("Pearson %+.3f, Spearman %+.3f, n = %d knees", r_p, r_s, nrow(cc)))
ggsave(file.path(BASE, "fig_cth_vs_cmftc.png"), g2, width = 7, height = 4.2, dpi = 200)

pick <- function(lab) {
  i <- which(res$model == lab & res$term %in% c("z_cth_pre", "z_cth_pre_e", "z_cth_all", "z_cmftc"))
  if (!length(i)) i <- which(res$model == lab)[1]
  i[1]
}
iB1 <- pick("B1  CTh-Score slope, PRE-INDEX, OLS   [primary]")
iB2 <- pick("B2  CTh-Score slope, PRE-INDEX, EBLUP [primary]")
iS2 <- pick("S2  SENS CTh-Score slope, ALL TIMEPOINTS")
iA1 <- pick("A1  cMFTC slope alone (kmri, reference)")
g3 <- data.frame(term = c("CTh-Score slope, PRE-INDEX (primary)",
                          "CTh-Score slope, pre-index, EBLUP",
                          "CTh-Score slope, all timepoints (sensitivity)",
                          "cMFTC slope (kmri, reference)"),
                 OR = res$OR[c(iB1, iB2, iS2, iA1)],
                 lo = res$lo[c(iB1, iB2, iS2, iA1)],
                 hi = res$hi[c(iB1, iB2, iS2, iA1)])
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
