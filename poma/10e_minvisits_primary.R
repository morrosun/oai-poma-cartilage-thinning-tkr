# ============================================================================
# 10e : SENSITIVITY -- restrict to matched pairs whose BOTH knees have
#       >= J MRI visits (J = 2 / 3 / 4), PRIMARY METRIC = cMFTC.
#
#   Why: the manuscript's primary analysis set requires only >= 2 visits per
#   knee, so ~1/3 of knees contribute exactly 2 points and their EBLUP /
#   latent slope is the least stable.  Requiring >= J visits for BOTH knees
#   of a pair removes the sparsest trajectories and is the natural robustness
#   check for the (joint model AND conditional logistic) rate association.
#
#   All estimators (joint model + conditional logistic) run on the IDENTICAL
#   knees at each J, exactly as in 10_jm_multimetric.R, so the two paths stay
#   comparable and the J = 2 row reproduces the primary Set A numbers.
#
#   NOTE the analysis SET changes with J (pairs are dropped), therefore the SD
#   of the EBLUP slope changes too: OR/0.1 is re-derived with the SD of the
#   CURRENT set, never with the Set-A SD.  (See skill gotcha: dose conversion
#   must use sd(eblup) of the analysis set actually used.)
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
source("D:/BaiduSyncdisk/OAI/Scripts/poma/jm_core.R")
suppressPackageStartupMessages(library(survival))

L <- read.csv(file.path(BASE, "poma_long.csv"))
U <- read.csv(file.path(BASE, "poma_analysis_long.csv"))
COVARS <- c("V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")
PRIM <- "cMFTC_ThCtAB_aMe"

make_input <- function(metric, min_visits, require_covars) {
  sub <- L[!is.na(L[[metric]]), ]
  lg  <- data.frame(unit = sub$unit, year = sub$years, y = sub[[metric]],
                    newstrata = sub$newstrata, case = sub$case)
  sv  <- U[, c("unit", "newstrata", "case", "t_index_months", COVARS, "V00PASE")]
  sv  <- sv[sv$unit %in% lg$unit, ]
  if (require_covars) sv <- sv[complete.cases(sv[, COVARS]), ]
  lg  <- lg[lg$unit %in% sv$unit, ]
  sv$t_idx <- sv$t_index_months / 12
  prune_pairs(lg, sv, min_visits)
}

sink(file.path(BASE, "sens_minvisits_primary_report.txt"), split = TRUE)
cat("================================================================================\n")
cat("POMA : SENSITIVITY -- >= J visits on BOTH knees; PRIMARY METRIC = cMFTC\n")
cat("================================================================================\n")

rows <- list(); keep <- list()

for (J in c(2, 3, 4)) {
  cat(sprintf("\n\n################  min_visits = %d  ################\n", J))
  A <- make_input(PRIM, J, TRUE)
  dt <- build(A$lg, A$sv)
  uA <- U[match(A$sv$unit, U$unit), ]
  nv <- table(A$lg$unit)
  cat(sprintf("pairs %d | knees %d | rows %d | visits per knee: min %d / median %.1f / max %d\n",
              length(dt$cas), dt$n, nrow(A$lg), min(nv), median(as.numeric(nv)), max(nv)))
  if (length(dt$cas) < 50) { cat("  -> fewer than 50 pairs, SKIP fitting\n"); next }

  ## ---- joint model (unadjusted + within-pair KL/pain adjusted) ----
  f  <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 2026, verbose = FALSE)
  r  <- report(f, sprintf("JM  min_visits=%d : cMFTC", J))
  Zd <- cbind(kl_diff    = uA$V00XRKL[dt$cas]    - uA$V00XRKL[dt$ctl],
              womac_diff = uA$WOMAC_pain[dt$cas] - uA$WOMAC_pain[dt$ctl])
  fa <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 2028,
               Zd = Zd, verbose = FALSE)
  ra <- report(fa, sprintf("JM  min_visits=%d : cMFTC + within-pair KL/pain", J))

  ## ---- conditional logistic on the IDENTICAL knees ----
  sl  <- uA[[paste0(PRIM, "__eblup")]]
  sdz <- sd(sl, na.rm = TRUE)
  dfc <- data.frame(case = uA$case, strata = uA$newstrata,
                    zz = (sl - mean(sl, na.rm = TRUE)) / sdz,
                    kl = uA$V00XRKL, wom5 = uA$WOMAC_pain / 5,
                    age5 = uA$V00AGE / 5, bmi5 = uA$P01BMI / 5,
                    female = uA$female, pase10 = uA$V00PASE / 10)
  cg0 <- clogit(case ~ zz + strata(strata), data = dfc, method = "exact")
  cg1 <- clogit(case ~ zz + age5 + female + bmi5 + kl + wom5 + strata(strata),
                data = dfc, method = "exact")
  s0 <- summary(cg0)$coef; s1 <- summary(cg1)$coef
  b0 <- s0["zz", 1]; e0 <- s0["zz", 3]
  b1 <- s1["zz", 1]; e1 <- s1["zz", 3]
  ## NOTE b0/b1 are NEGATIVE (raw EBLUP slope: negative = thinning), the OR is
  ## exp(-b); the interval is therefore exp(-b -/+ 1.96*se) -- do NOT write
  ## exp(-(b - 1.96*se)) which silently swaps the bounds.
  cat(sprintf("  [clogit crude] OR/SD %.3f (%.3f-%.3f)  OR/0.1 %.3f (%.3f-%.3f)  P=%.3g\n",
              exp(-b0), exp(-b0 - 1.96*e0), exp(-b0 + 1.96*e0),
              exp(-0.1*b0/sdz), exp(-0.1*b0/sdz - 0.196*e0/sdz), exp(-0.1*b0/sdz + 0.196*e0/sdz),
              s0["zz", 5]))
  cat(sprintf("  [clogit adj  ] OR/SD %.3f (%.3f-%.3f)  OR/0.1 %.3f (%.3f-%.3f)  P=%.3g\n",
              exp(-b1), exp(-b1 - 1.96*e1), exp(-b1 + 1.96*e1),
              exp(-0.1*b1/sdz), exp(-0.1*b1/sdz - 0.196*e1/sdz), exp(-0.1*b1/sdz + 0.196*e1/sdz),
              s1["zz", 5]))

  ## ---- level vs rate inside the restricted set (does the rate still absorb?) ----
  base_z <- scale(uA[[paste0(PRIM, "__base")]])
  dfc$zb <- as.numeric(base_z)
  mL <- clogit(case ~ zb + strata(strata), data = dfc, method = "exact")
  mR <- clogit(case ~ zz + strata(strata), data = dfc, method = "exact")
  mB <- clogit(case ~ zb + zz + strata(strata), data = dfc, method = "exact")
  chi_rate <- 2 * (as.numeric(logLik(mB)) - as.numeric(logLik(mL)))
  chi_lev  <- 2 * (as.numeric(logLik(mB)) - as.numeric(logLik(mR)))
  cat(sprintf("  [level vs rate] add-rate chi2 = %.2f (P=%.3g) ; add-level chi2 = %.2f (P=%.3g)\n",
              chi_rate, pchisq(chi_rate, 1, lower.tail = FALSE),
              chi_lev,  pchisq(chi_lev, 1, lower.tail = FALSE)))

  d <- f$draws; da <- fa$draws
  rows[[as.character(J)]] <- data.frame(
    min_visits = J, n_pairs = length(dt$cas), n_knees = dt$n, n_rows = nrow(A$lg),
    visits_per_knee = nrow(A$lg) / dt$n, sd_eblup = sdz,
    JM_OR01 = median(exp(-0.1 * d$aS)),
    JM_OR01_lo = quantile(exp(-0.1 * d$aS), .025),
    JM_OR01_hi = quantile(exp(-0.1 * d$aS), .975),
    JM_ORSD = median(exp(-d$aS * d$SDslo)),
    JM_ORSD_lo = quantile(exp(-d$aS * d$SDslo), .025),
    JM_ORSD_hi = quantile(exp(-d$aS * d$SDslo), .975),
    JM_aS = median(d$aS), JM_aS_lo = quantile(d$aS, .025),
    JM_aS_hi = quantile(d$aS, .975),
    JM_OR01_adj = median(exp(-0.1 * da$aS)),
    JM_OR01_adj_lo = quantile(exp(-0.1 * da$aS), .025),
    JM_OR01_adj_hi = quantile(exp(-0.1 * da$aS), .975),
    clogit_ORSD = exp(-b0), clogit_ORSD_lo = exp(-b0 - 1.96*e0), clogit_ORSD_hi = exp(-b0 + 1.96*e0),
    clogit_p = s0["zz", 5],
    clogit_OR01 = exp(-0.1*b0/sdz), clogit_OR01_lo = exp(-0.1*b0/sdz - 0.196*e0/sdz),
    clogit_OR01_hi = exp(-0.1*b0/sdz + 0.196*e0/sdz),
    clogitAdj_ORSD = exp(-b1), clogitAdj_ORSD_lo = exp(-b1 - 1.96*e1),
    clogitAdj_ORSD_hi = exp(-b1 + 1.96*e1), clogitAdj_p = s1["zz", 5],
    clogitAdj_OR01 = exp(-0.1*b1/sdz), clogitAdj_OR01_lo = exp(-0.1*b1/sdz - 0.196*e1/sdz),
    clogitAdj_OR01_hi = exp(-0.1*b1/sdz + 0.196*e1/sdz),
    LR_add_rate = chi_rate, LR_add_rate_p = pchisq(chi_rate, 1, lower.tail = FALSE),
    LR_add_level = chi_lev, LR_add_level_p = pchisq(chi_lev, 1, lower.tail = FALSE),
    stringsAsFactors = FALSE)
  keep[[as.character(J)]] <- list(JM = f, JMadj = fa)
}

tab <- do.call(rbind, rows)
write.csv(tab, file.path(BASE, "sens_minvisits_primary.csv"), row.names = FALSE)
saveRDS(keep, file.path(BASE, "sens_minvisits_primary_fit.rds"))

cat("\n\n================ SENSITIVITY SUMMARY : cMFTC rate association ================\n")
print(tab[, c("min_visits", "n_pairs", "n_knees", "sd_eblup",
              "JM_OR01", "JM_OR01_lo", "JM_OR01_hi",
              "JM_ORSD", "JM_ORSD_lo", "JM_ORSD_hi",
              "clogit_ORSD", "clogit_ORSD_lo", "clogit_ORSD_hi", "clogit_p",
              "clogitAdj_ORSD", "clogitAdj_p")], row.names = FALSE)
cat("\n-> sens_minvisits_primary.csv / sens_minvisits_primary_report.txt\n")
sink()
