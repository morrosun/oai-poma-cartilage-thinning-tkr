# ============================================================================
# POMA : conditional logistic, crude + covariate adjusted  --  PRIMARY METRIC
#   exposure = EBLUP thinning rate of cMFTC central medial weight-bearing
#              cartilage thickness (mm/yr)  <<< the pre-specified primary metric
#   z-standardised within the analysis sample; positive = faster thinning
#   strata  = newstrata ; analysis set = complete matched pairs (== Set A)
#   Purpose: give the manuscript a covariate-adjustment table that uses the
#            SAME metric and SAME pair set as the joint model / gradient tests.
# ============================================================================
suppressWarnings(suppressMessages({
  library(survival); library(dplyr); library(ggplot2); library(tidyr)
}))

OUT <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
d <- read.csv(file.path(OUT, "poma_analysis_long.csv"), stringsAsFactors = FALSE)

sink(file.path(OUT, "clogit_cMFTC_seq_report.txt"), split = TRUE)
cat("================================================================================\n")
cat("POMA : conditional logistic (PRIMARY METRIC = cMFTC), crude + adjusted\n")
cat("================================================================================\n\n")

PRIM   <- "cMFTC_ThCtAB_aMe"
d$slope <- -d[[paste0(PRIM, "__eblup")]]
d$base  <- d[[paste0(PRIM, "__base")]]

need <- c("slope", "base", "V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")
complete_pair <- function(dat, vars) {
  dat %>%
    group_by(newstrata) %>%
    mutate(n_ok  = sum(complete.cases(across(all_of(vars)))),
           n_tot = n()) %>%
    ungroup() %>%
    filter(n_ok == n_tot, n_tot == 2,
           complete.cases(across(all_of(vars))))
}
a <- complete_pair(d, need)
cat(sprintf("complete matched pairs with all covariates : %d / %d\n",
            n_distinct(a$newstrata), n_distinct(d$newstrata)))
cat(sprintf("knees in analysis                         : %d\n\n", nrow(a)))

z <- function(x) as.numeric(scale(x))
a <- a %>% mutate(
  z_loss = z(slope), loss01 = slope / 0.1, z_base = z(base),
  age5 = V00AGE / 5, bmi5 = P01BMI / 5, kl = V00XRKL,
  womac5 = WOMAC_pain / 5, pase10 = V00PASE / 10)

SDB <- sd(a$slope); cat(sprintf("SD of EBLUP thinning rate (cMFTC) : %.4f mm/yr\n\n", SDB))

cat("--- baseline characteristics of the analysis set (matched pairs) ---\n")
tab1 <- a %>% group_by(case) %>% summarise(
  n = n(),
  age   = sprintf("%.1f (%.1f)", mean(V00AGE, na.rm = TRUE), sd(V00AGE, na.rm = TRUE)),
  female_pct = sprintf("%.1f", 100 * mean(female)),
  bmi   = sprintf("%.1f (%.1f)", mean(P01BMI, na.rm = TRUE), sd(P01BMI, na.rm = TRUE)),
  kl    = sprintf("%.2f (%.2f)", mean(kl, na.rm = TRUE), sd(kl, na.rm = TRUE)),
  womac = sprintf("%.1f (%.1f)", mean(WOMAC_pain, na.rm = TRUE), sd(WOMAC_pain, na.rm = TRUE)),
  cmftc_thin_rate = sprintf("%.4f (%.4f)", mean(slope), sd(slope)),
  .groups = "drop")
print(as.data.frame(t(tab1)), quote = FALSE)

cat("\n--- paired within-stratum comparisons (case - control) ---\n")
pt <- a %>% select(newstrata, case, V00AGE, P01BMI, kl, WOMAC_pain, V00PASE, slope) %>%
  pivot_wider(names_from = case,
              values_from = c(V00AGE, P01BMI, kl, WOMAC_pain, V00PASE, slope)) %>%
  mutate(d_age = V00AGE_1 - V00AGE_0, d_bmi = P01BMI_1 - P01BMI_0,
         d_kl = kl_1 - kl_0, d_wom = WOMAC_pain_1 - WOMAC_pain_0,
         d_pase = V00PASE_1 - V00PASE_0, d_slope = slope_1 - slope_0)
for (v in c("d_age","d_bmi","d_kl","d_wom","d_pase","d_slope")) {
  x <- pt[[v]]; x <- x[is.finite(x)]; tt <- t.test(x)
  cat(sprintf("  %-8s mean %+8.3f  95%%CI [%+8.3f, %+8.3f]  paired t P = %s\n",
              v, mean(x), tt$conf.int[1], tt$conf.int[2],
              format.pval(tt$p.value, digits = 3)))
}

fit_clogit <- function(form, dat, lab) {
  m <- clogit(as.formula(form), data = dat, method = "exact"); s <- summary(m); co <- s$coefficients
  data.frame(model = lab, term = rownames(co)[1], beta = co[1,1], se = co[1,3],
             OR = exp(co[1,1]), lo = exp(co[1,1] - 1.96*co[1,3]),
             hi = exp(co[1,1] + 1.96*co[1,3]), p = co[1,5],
             n_pairs = m$nevent, loglik = as.numeric(logLik(m)), stringsAsFactors = FALSE)
}

cat("\n\n================================================================================\n")
cat("A. effect of the cMFTC thinning RATE, sequentially adjusted (OR per 1 SD)\n")
cat("   and re-expressed per 0.1 mm/yr (OR/0.1 = exp(beta*0.1/SD))\n")
cat("================================================================================\n")
mods <- list(
  list("M0 crude",            "case ~ z_loss + strata(newstrata)"),
  list("M1 + age(/5y)",       "case ~ z_loss + age5 + strata(newstrata)"),
  list("M2 + sex",            "case ~ z_loss + female + strata(newstrata)"),
  list("M3 + BMI(/5)",        "case ~ z_loss + bmi5 + strata(newstrata)"),
  list("M4 + KL grade",       "case ~ z_loss + kl + strata(newstrata)"),
  list("M5 + WOMAC pain(/5)", "case ~ z_loss + womac5 + strata(newstrata)"),
  list("M6 full adjusted",    "case ~ z_loss + age5 + female + bmi5 + kl + womac5 + strata(newstrata)"),
  list("M7 full + PASE",      "case ~ z_loss + age5 + female + bmi5 + kl + womac5 + pase10 + strata(newstrata)"))
resA <- bind_rows(lapply(mods, function(m) fit_clogit(m[[2]], a, m[[1]])))
resA$OR_01 <- exp(resA$beta * 0.1 / SDB)
resA$OR_01_lo <- exp((resA$beta - 1.96*resA$se) * 0.1 / SDB)
resA$OR_01_hi <- exp((resA$beta + 1.96*resA$se) * 0.1 / SDB)
print(resA %>% mutate(across(c(beta,se,OR,lo,hi,OR_01,OR_01_lo,OR_01_hi), ~round(.x,4)),
                      p = signif(p,3)) %>%
        select(model, OR, lo, hi, OR_01, OR_01_lo, OR_01_hi, p, n_pairs), row.names = FALSE)

cat("\n--- full coefficient table for M6 ---\n")
print(summary(clogit(case ~ z_loss + age5 + female + bmi5 + kl + womac5 + strata(newstrata),
                     data = a, method = "exact"))$coefficients)
cat("\n--- and for M0 (crude) ---\n")
print(summary(clogit(case ~ z_loss + strata(newstrata), data = a, method = "exact"))$coefficients)

cat("\n\n================================================================================\n")
cat("B. BASELINE LEVEL vs RATE OF CHANGE  (cMFTC)\n")
cat("================================================================================\n")
mods2 <- list(
  list("level only  (z_base)",  "case ~ z_base + strata(newstrata)"),
  list("rate only   (z_loss)",  "case ~ z_loss + strata(newstrata)"),
  list("level + rate",          "case ~ z_base + z_loss + strata(newstrata)"),
  list("level + rate + covars", "case ~ z_base + z_loss + age5 + female + bmi5 + kl + strata(newstrata)"))
resB <- bind_rows(lapply(mods2, function(m) fit_clogit(m[[2]], a, m[[1]])))
print(resB %>% mutate(across(c(beta,se,OR,lo,hi), ~round(.x,4)), p = signif(p,3)) %>%
        select(model, term, OR, lo, hi, p, n_pairs, loglik), row.names = FALSE)

m_lev <- clogit(case ~ z_base + strata(newstrata), data = a, method = "exact")
m_rat <- clogit(case ~ z_loss + strata(newstrata), data = a, method = "exact")
m_bot <- clogit(case ~ z_base + z_loss + strata(newstrata), data = a, method = "exact")
cat(sprintf("\nLR test, rate-only vs level-only (non-nested, descriptive): dLogLik = %.2f\n",
            2*(as.numeric(logLik(m_rat)) - as.numeric(logLik(m_lev)))))
cat(sprintf("LR test, adding rate to level  : chi2 = %.2f, df=1, P = %s\n",
            2*(as.numeric(logLik(m_bot)) - as.numeric(logLik(m_lev))),
            format.pval(pchisq(2*(as.numeric(logLik(m_bot)) - as.numeric(logLik(m_lev))),
                               1, lower.tail = FALSE), digits = 3)))
cat(sprintf("LR test, adding level to rate  : chi2 = %.2f, df=1, P = %s\n",
            2*(as.numeric(logLik(m_bot)) - as.numeric(logLik(m_rat))),
            format.pval(pchisq(2*(as.numeric(logLik(m_bot)) - as.numeric(logLik(m_rat))),
                               1, lower.tail = FALSE), digits = 3)))

cat("\n\n================================================================================\n")
cat("C. DESIGN CHECK: stratified Cox == conditional logistic ?\n")
cat("================================================================================\n")
cx <- coxph(Surv(t_index_months, case) ~ z_loss + strata(newstrata), data = a)
cg <- clogit(case ~ z_loss + strata(newstrata), data = a, method = "exact")
cat(sprintf("  coxph  beta = %+.6f (SE %.6f)\n", coef(cx), sqrt(vcov(cx)[1,1])))
cat(sprintf("  clogit beta = %+.6f (SE %.6f)\n", coef(cg), sqrt(vcov(cg)[1,1])))

write.csv(resA, file.path(OUT, "clogit_cMFTC_seq_OR.csv"), row.names = FALSE)
write.csv(resB, file.path(OUT, "clogit_cMFTC_level_vs_rate.csv"), row.names = FALSE)
cat("\n-> clogit_cMFTC_seq_OR.csv / clogit_cMFTC_level_vs_rate.csv\n")
sink()
