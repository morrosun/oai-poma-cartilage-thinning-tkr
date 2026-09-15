# ============================================================================
# POMA nested case-control : conditional logistic, crude and covariate-adjusted
#   exposure = EBLUP thinning rate of MFTC mean cartilage thickness (mm/yr)
#   z-standardised within the analysis sample; positive = faster thinning
#   strata  = newstrata (225 matched knee pairs)
# ============================================================================
suppressWarnings(suppressMessages({
  library(survival); library(dplyr); library(ggplot2); library(tidyr)
}))

OUT <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
d <- read.csv(file.path(OUT, "poma_analysis_long.csv"),
              stringsAsFactors = FALSE)

sink(file.path(OUT, "clogit_adjusted_report.txt"), split = TRUE)
cat("================================================================================\n")
cat("POMA : conditional logistic regression, crude + covariate adjusted\n")
cat("================================================================================\n\n")

PRIM   <- "MFTC_ThCtAB_aMe"
d$slope <- -d[[paste0(PRIM, "__eblup")]]      # thinning rate, + = faster loss
d$base  <- d[[paste0(PRIM, "__base")]]

# ---- analysis set: COMPLETE matched pairs with slope + all covariates -------
need <- c("slope", "base", "V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")
complete_pair <- function(dat, vars) {
  ok <- dat %>%
    group_by(newstrata) %>%
    mutate(n_ok = sum(complete.cases(across(all_of(vars)))),
           n_tot = n()) %>%
    ungroup() %>%
    filter(n_ok == n_tot, n_tot == 2,
           complete.cases(across(all_of(vars))))
  ok
}

a <- complete_pair(d, need)
cat(sprintf("complete matched pairs with all covariates : %d / %d\n",
            n_distinct(a$newstrata), n_distinct(d$newstrata)))
cat(sprintf("knees in analysis                         : %d\n\n", nrow(a)))

# ---- standardise within analysis set ---------------------------------------
z <- function(x) as.numeric(scale(x))
a <- a %>% mutate(
  z_loss   = z(slope),
  loss01   = slope / 0.1,                       # per 0.1 mm/yr
  z_base   = z(base),
  age5     = V00AGE / 5,
  bmi5     = P01BMI / 5,
  kl       = V00XRKL,
  womac5   = WOMAC_pain / 5,
  pase10   = V00PASE / 10
)

cat("--- baseline characteristics of the analysis set (matched pairs) ---\n")
tab1 <- a %>% group_by(case) %>% summarise(
  n = n(),
  age = sprintf("%.1f (%.1f)", mean(V00AGE, na.rm = TRUE), sd(V00AGE, na.rm = TRUE)),
  female_pct = sprintf("%.1f", 100 * mean(female)),
  bmi = sprintf("%.1f (%.1f)", mean(P01BMI, na.rm = TRUE), sd(P01BMI, na.rm = TRUE)),
  kl  = sprintf("%.2f (%.2f)", mean(kl, na.rm = TRUE), sd(kl, na.rm = TRUE)),
  womac = sprintf("%.1f (%.1f)", mean(WOMAC_pain, na.rm = TRUE), sd(WOMAC_pain, na.rm = TRUE)),
  thin_rate = sprintf("%.4f (%.4f)", mean(slope), sd(slope)),
  .groups = "drop")
print(as.data.frame(t(tab1)), quote = FALSE)

# ---- Table 1 style paired comparison ---------------------------------------
cat("\n--- paired within-stratum comparisons (case - control) ---\n")
pairs_tab <- a %>%
  select(newstrata, case, V00AGE, P01BMI, kl, WOMAC_pain, V00PASE, slope) %>%
  pivot_wider(names_from = case, values_from = c(V00AGE, P01BMI, kl, WOMAC_pain, V00PASE, slope)) %>%
  mutate(d_age = V00AGE_1 - V00AGE_0,
         d_bmi = P01BMI_1 - P01BMI_0,
         d_kl  = kl_1 - kl_0,
         d_wom = WOMAC_pain_1 - WOMAC_pain_0,
         d_pase= V00PASE_1 - V00PASE_0,
         d_slope = slope_1 - slope_0)
for (v in c("d_age", "d_bmi", "d_kl", "d_wom", "d_pase", "d_slope")) {
  x <- pairs_tab[[v]]; x <- x[is.finite(x)]
  tt <- t.test(x)
  cat(sprintf("  %-8s mean %+8.3f  95%%CI [%+8.3f, %+8.3f]  paired t P = %s\n",
              v, mean(x), tt$conf.int[1], tt$conf.int[2],
              format.pval(tt$p.value, digits = 3)))
}

# ---- conditional logistic models -------------------------------------------
fit_clogit <- function(form, dat, lab) {
  m <- clogit(as.formula(form), data = dat, method = "exact")
  s <- summary(m)
  co <- s$coefficients
  # main exposure = first term
  i <- 1
  data.frame(model = lab,
             term = rownames(co)[i],
             beta = co[i, 1], se = co[i, 3],
             OR = exp(co[i, 1]),
             lo = exp(co[i, 1] - 1.96 * co[i, 3]),
             hi = exp(co[i, 1] + 1.96 * co[i, 3]),
             p = co[i, 5],
             n_pairs = m$nevent,
             n = length(m$residuals) + m$nevent,
             loglik = as.numeric(logLik(m)),
             stringsAsFactors = FALSE)
}

mods <- list(
  list("M0 crude",              "case ~ z_loss + strata(newstrata)"),
  list("M1 + age(/5y)",         "case ~ z_loss + age5 + strata(newstrata)"),
  list("M2 + sex",              "case ~ z_loss + female + strata(newstrata)"),
  list("M3 + BMI(/5)",          "case ~ z_loss + bmi5 + strata(newstrata)"),
  list("M4 + KL grade",         "case ~ z_loss + kl + strata(newstrata)"),
  list("M5 + WOMAC pain(/5)",   "case ~ z_loss + womac5 + strata(newstrata)"),
  list("M6 full adjusted",      "case ~ z_loss + age5 + female + bmi5 + kl + womac5 + strata(newstrata)"),
  list("M7 full + PASE",        "case ~ z_loss + age5 + female + bmi5 + kl + womac5 + pase10 + strata(newstrata)")
)

cat("\n\n================================================================================\n")
cat("A. effect of the THINNING RATE, sequentially adjusted (OR per 1 SD)\n")
cat("================================================================================\n")
resA <- bind_rows(lapply(mods, function(m) fit_clogit(m[[2]], a, m[[1]])))
print(resA %>% mutate(across(c(beta, se, OR, lo, hi), ~round(.x, 4)),
                      p = signif(p, 3)) %>%
        select(model, OR, lo, hi, p, n_pairs), row.names = FALSE)

cat("\n--- full coefficient table for M6 ---\n")
print(summary(clogit(case ~ z_loss + age5 + female + bmi5 + kl + womac5 +
                       strata(newstrata), data = a, method = "exact"))$coefficients)

cat("\n--- and for M0 (crude) ---\n")
print(summary(clogit(case ~ z_loss + strata(newstrata), data = a,
                     method = "exact"))$coefficients)

# ---- level vs rate ----------------------------------------------------------
cat("\n\n================================================================================\n")
cat("B. BASELINE LEVEL vs RATE OF CHANGE  (which carries the signal?)\n")
cat("================================================================================\n")
mods2 <- list(
  list("level only  (z_base)",  "case ~ z_base + strata(newstrata)"),
  list("rate only   (z_loss)",  "case ~ z_loss + strata(newstrata)"),
  list("level + rate",          "case ~ z_base + z_loss + strata(newstrata)"),
  list("level + rate + covars", "case ~ z_base + z_loss + age5 + female + bmi5 + kl + strata(newstrata)")
)
resB <- bind_rows(lapply(mods2, function(m) fit_clogit(m[[2]], a, m[[1]])))
print(resB %>% mutate(across(c(beta, se, OR, lo, hi), ~round(.x, 4)),
                      p = signif(p, 3)) %>%
        select(model, term, OR, lo, hi, p, n_pairs, loglik), row.names = FALSE)

m_lev <- clogit(case ~ z_base + strata(newstrata), data = a, method = "exact")
m_rat <- clogit(case ~ z_loss + strata(newstrata), data = a, method = "exact")
m_bot <- clogit(case ~ z_base + z_loss + strata(newstrata), data = a, method = "exact")
lr1 <- 2 * (as.numeric(logLik(m_rat)) - as.numeric(logLik(m_lev)))
cat(sprintf("\nLR test, rate-only vs level-only (non-nested, descriptive): dLogLik = %.2f\n", lr1))
cat(sprintf("LR test, adding rate to level  : chi2 = %.2f, df=1, P = %s\n",
            2 * (as.numeric(logLik(m_bot)) - as.numeric(logLik(m_lev))),
            format.pval(pchisq(2 * (as.numeric(logLik(m_bot)) - as.numeric(logLik(m_lev))),
                               1, lower.tail = FALSE), digits = 3)))
cat(sprintf("LR test, adding level to rate  : chi2 = %.2f, df=1, P = %s\n",
            2 * (as.numeric(logLik(m_bot)) - as.numeric(logLik(m_rat))),
            format.pval(pchisq(2 * (as.numeric(logLik(m_bot)) - as.numeric(logLik(m_rat))),
                               1, lower.tail = FALSE), digits = 3)))

# ---- equivalence of clogit and stratified Cox (design point) ----------------
cat("\n\n================================================================================\n")
cat("C. DESIGN CHECK: stratified Cox == conditional logistic ?\n")
cat("================================================================================\n")
cx <- coxph(Surv(t_index_months, case) ~ z_loss + strata(newstrata), data = a)
cg <- clogit(case ~ z_loss + strata(newstrata), data = a, method = "exact")
cat(sprintf("  coxph  beta = %+.6f (SE %.6f)\n", coef(cx), sqrt(vcov(cx)[1, 1])))
cat(sprintf("  clogit beta = %+.6f (SE %.6f)\n", coef(cg), sqrt(vcov(cg)[1, 1])))
cat("  -> identical because the index time is shared within every matched set,\n")
cat("     so survival time carries no within-stratum information.\n")

# ---- forest plot ------------------------------------------------------------
p <- ggplot(resA, aes(x = model, y = OR, ymin = lo, ymax = hi)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(size = .6, colour = "#C44E52") +
  coord_flip() + theme_bw(base_size = 11) +
  labs(x = "", y = "OR per 1 SD faster MFTC thinning",
       title = "POMA: thinning rate and odds of knee replacement",
       subtitle = "conditional logistic, 225 matched pairs (complete-case n shown)") +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(OUT, "fig_forest_adjusted.png"), p, width = 8, height = 4.2,
       dpi = 200)

write.csv(resA, file.path(OUT, "clogit_adjusted_OR.csv"), row.names = FALSE)
write.csv(resB, file.path(OUT, "clogit_level_vs_rate.csv"), row.names = FALSE)
cat("\n-> clogit_adjusted_OR.csv / clogit_level_vs_rate.csv / fig_forest_adjusted.png\n")
sink()
