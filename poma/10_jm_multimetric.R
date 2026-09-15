# ============================================================================
# 10 : multi-metric joint model -> medial vs lateral gradient
#
#   Primary metric = cMFTC (central medial femur+tibia, weight-bearing).
#   A-priori LATERAL negative controls: cLF, cLFTC.
#   A monotone medial-to-lateral gradient in the RATE association is the
#   strongest single piece of internal-validity evidence in this design.
#
#   ANALYSIS SETS
#     A (primary) : complete matched pairs in which BOTH knees have >= 2 MRI
#                   visits and complete baseline covariates -> the thinning
#                   slope is estimable for every knee and BOTH estimators
#                   (joint model and conditional logistic) run on the SAME
#                   knees, so they are directly comparable.
#     B (gain)    : all complete matched pairs (>= 1 visit).  Only the joint
#                   model can use these, because for a single-visit knee the
#                   EBLUP slope collapses onto the population mean and an
#                   EBLUP-based conditional logistic is attenuated.  Set B is
#                   therefore reported as a data-gain sensitivity analysis.
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
source("D:/BaiduSyncdisk/OAI/Scripts/poma/jm_core.R")
suppressPackageStartupMessages(library(survival))

L <- read.csv(file.path(BASE, "poma_long.csv"))
U <- read.csv(file.path(BASE, "poma_analysis_long.csv"))

COVARS <- c("V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")

METRICS <- list(
  list(m = "cMFTC_ThCtAB_aMe", lab = "cMFTC central medial (weight-bearing)", reg = "Medial"),
  list(m = "MFTC_ThCtAB_aMe",  lab = "MFTC medial (whole compartment)",       reg = "Medial"),
  list(m = "MT_ThCtAB_aMe",    lab = "MT medial tibia",                       reg = "Medial"),
  list(m = "cMF_ThCtAB_aMe",   lab = "cMF central medial femur",              reg = "Medial"),
  list(m = "cLFTC_ThCtAB_aMe", lab = "cLFTC central lateral (weight-bearing)", reg = "Lateral"),
  list(m = "cLF_ThCtAB_aMe",   lab = "cLF central lateral femur",             reg = "Lateral")
)

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

## OR + CI for a log-hazard coefficient, converted to a target dose
orci <- function(b, se, dose) {
  z <- -dose                    # faster thinning = lower slope
  c(OR = exp(z * b), lo = exp(z * (b - 1.96 * se)), hi = exp(z * (b + 1.96 * se)))
}

res <- list(); keep_fits <- list()

for (k in seq_along(METRICS)) {
  mm <- METRICS[[k]]
  cat(sprintf("\n\n########## [%d/%d] %s ##########\n", k, length(METRICS), mm$lab))

  ## ---------------- SET A : primary, both estimators identical set ---------
  A <- make_input(mm$m, 2, TRUE)
  dtA <- build(A$lg, A$sv)
  uA <- U[match(A$sv$unit, U$unit), ]
  cat(sprintf("SET A : %d pairs / %d knees / %d rows (%.2f visits per knee)\n",
              length(dtA$cas), dtA$n, nrow(A$lg), nrow(A$lg) / dtA$n))

  fA <- fit_jm(dtA, NITER = 30000, BURN = 5000, THIN = 5, seed = 2026, verbose = FALSE)
  rA <- report(fA, paste0("JM SET A : ", mm$lab))

  ## within-pair covariate differences (only these are identifiable)
  Zd <- cbind(kl_diff    = uA$V00XRKL[dtA$cas]    - uA$V00XRKL[dtA$ctl],
              womac_diff = uA$WOMAC_pain[dtA$cas] - uA$WOMAC_pain[dtA$ctl])
  fAadj <- fit_jm(dtA, NITER = 30000, BURN = 5000, THIN = 5, seed = 2028,
                  Zd = Zd, verbose = FALSE)
  report(fAadj, paste0("JM SET A + within-pair KL/pain : ", mm$lab))

  ## conditional logistic on the SAME knees
  sl <- uA[[paste0(mm$m, "__eblup")]]
  sdz <- sd(sl, na.rm = TRUE)
  dfc <- data.frame(case = uA$case, strata = uA$newstrata,
                    zz = (sl - mean(sl, na.rm = TRUE)) / sdz,
                    loss = -sl,
                    kl = uA$V00XRKL, wom5 = uA$WOMAC_pain / 5,
                    age5 = uA$V00AGE / 5, bmi5 = uA$P01BMI / 5,
                    female = uA$female)
  cg0 <- clogit(case ~ zz + strata(strata), data = dfc, method = "exact")
  cg1 <- clogit(case ~ zz + age5 + female + bmi5 + kl + wom5 + strata(strata),
                data = dfc, method = "exact")
  s0 <- summary(cg0)$coef; s1 <- summary(cg1)$coef
  b0 <- s0["zz", 1]; e0 <- s0["zz", 3]
  b1 <- s1["zz", 1]; e1 <- s1["zz", 3]
  ## NOTE b0/b1 are NEGATIVE (raw EBLUP slope, negative = thinning) so the OR is
  ## exp(-b); the interval is exp(-b -/+ 1.96*se).  Writing exp(-(b - 1.96*se))
  ## silently SWAPS the bounds -- fixed 2026-09-14 (same defect found in 10e).
  cat(sprintf("  [clogit crude]  OR/SD %.3f (%.3f-%.3f)   OR/0.1mm/yr %.3f\n",
              exp(-b0), exp(-b0 - 1.96*e0), exp(-b0 + 1.96*e0),
              exp(-0.1*b0/sdz)))
  cat(sprintf("  [clogit adj  ]  OR/SD %.3f (%.3f-%.3f)   OR/0.1mm/yr %.3f\n",
              exp(-b1), exp(-b1 - 1.96*e1), exp(-b1 + 1.96*e1),
              exp(-0.1*b1/sdz)))

  ## ---------------- SET B : all complete pairs, joint model only ----------
  Bst <- make_input(mm$m, 1, FALSE)
  dtB <- build(Bst$lg, Bst$sv)
  fB <- fit_jm(dtB, NITER = 30000, BURN = 5000, THIN = 5, seed = 2027, verbose = FALSE)
  cat(sprintf("SET B : %d pairs / %d knees\n", length(dtB$cas), dtB$n))

  d <- fA$draws
  row <- data.frame(
    metric = mm$m, label = mm$lab, region = mm$reg,
    n_pairs_A = length(dtA$cas), n_pairs_B = length(dtB$cas),
    sd_eblup = sdz,
    sd_latent_slope = median(d$SDslo),
    beta1_JMA = median(d$beta1),
    aS_A = median(d$aS), aS_lo = quantile(d$aS, .025), aS_hi = quantile(d$aS, .975),
    JM_A_OR01 = median(exp(-0.1 * d$aS)),
    JM_A_OR01_lo = quantile(exp(-0.1 * d$aS), .025),
    JM_A_OR01_hi = quantile(exp(-0.1 * d$aS), .975),
    JM_A_ORSD = median(exp(-d$aS * d$SDslo)),
    JM_A_ORSD_lo = quantile(exp(-d$aS * d$SDslo), .025),
    JM_A_ORSD_hi = quantile(exp(-d$aS * d$SDslo), .975),
    JM_Aadj_OR01 = median(exp(-0.1 * fAadj$draws$aS)),
    JM_Aadj_OR01_lo = quantile(exp(-0.1 * fAadj$draws$aS), .025),
    JM_Aadj_OR01_hi = quantile(exp(-0.1 * fAadj$draws$aS), .975),
    JM_Aadj_ORSD = median(exp(-fAadj$draws$aS * fAadj$draws$SDslo)),    clogit_A_OR01 = exp(-0.1 * b0 / sdz),
    clogit_A_OR01_lo = exp(-0.1 * b0 / sdz - 0.196 * e0 / sdz),
    clogit_A_OR01_hi = exp(-0.1 * b0 / sdz + 0.196 * e0 / sdz),
    clogit_A_ORSD = exp(-b0), clogit_A_ORSD_lo = exp(-b0 - 1.96 * e0),
    clogit_A_ORSD_hi = exp(-b0 + 1.96 * e0),
    clogit_A_p = s0["zz", 5],
    clogitAdj_A_OR01 = exp(-0.1 * b1 / sdz),
    clogitAdj_A_OR01_lo = exp(-0.1 * b1 / sdz - 0.196 * e1 / sdz),
    clogitAdj_A_OR01_hi = exp(-0.1 * b1 / sdz + 0.196 * e1 / sdz),
    clogitAdj_A_ORSD = exp(-b1), clogitAdj_A_p = s1["zz", 5],
    JM_B_OR01 = median(exp(-0.1 * fB$draws$aS)),
    JM_B_OR01_lo = quantile(exp(-0.1 * fB$draws$aS), .025),
    JM_B_OR01_hi = quantile(exp(-0.1 * fB$draws$aS), .975),
    JM_B_ORSD = median(exp(-fB$draws$aS * fB$draws$SDslo)),
    stringsAsFactors = FALSE)
  ## level (per 1 SD of latent level) for the level-vs-rate contrast
  row$JM_A_ORSDlevel <- median(exp(-d$aV * d$SDval))
  row$JM_A_ORSDlevel_lo <- quantile(exp(-d$aV * d$SDval), .025)
  row$JM_A_ORSDlevel_hi <- quantile(exp(-d$aV * d$SDval), .975)

  res[[k]] <- row
  keep_fits[[mm$m]] <- list(A = fA, Aadj = fAadj, B = fB)
  saveRDS(keep_fits[[mm$m]], file.path(BASE, paste0("jmfit_", mm$m, ".rds")))
}

tab <- do.call(rbind, res)
write.csv(tab, file.path(BASE, "jm_multimetric.csv"), row.names = FALSE)

## -------- level-vs-rate and Bayes factors for the primary metric ----------
prim <- keep_fits[["cMFTC_ThCtAB_aMe"]]
cat("\n\n########## PRIMARY cMFTC : level-only model (aS = 0) ##########\n")
A <- make_input("cMFTC_ThCtAB_aMe", 2, TRUE); dtA <- build(A$lg, A$sv)
f0 <- fit_jm(dtA, NITER = 30000, BURN = 5000, THIN = 5, seed = 2029,
             aS_fixed = 0, verbose = FALSE)
report(f0, "JM SET A : cMFTC, latent LEVEL only (aS = 0)")
bfS <- sd_bf01(prim$A$draws$aS); bfV <- sd_bf01(prim$A$draws$aV)
cat(sprintf("\n  Savage-Dickey  H0: aS = 0 -> BF01 = %.4g  (BF10 = %.3g)\n", bfS, 1 / bfS))
cat(sprintf("  Savage-Dickey  H0: aV = 0 -> BF01 = %.4g  (BF10 = %.3g)\n", bfV, 1 / bfV))

## ------------------------------------------------------------- figure -----
suppressPackageStartupMessages({library(ggplot2)})
mk <- function(est, lo, hi, est2, lo2, hi2) {
  rbind(
    data.frame(label = tab$label, region = tab$region, est = tab[[est]],
               lo = tab[[lo]], hi = tab[[hi]], type = "Joint model"),
    data.frame(label = tab$label, region = tab$region, est = tab[[est2]],
               lo = tab[[lo2]], hi = tab[[hi2]], type = "Conditional logistic"))
}
pd1 <- mk("JM_A_OR01", "JM_A_OR01_lo", "JM_A_OR01_hi",
          "clogit_A_OR01", "clogit_A_OR01_lo", "clogit_A_OR01_hi")
pd1$scale <- "per 0.1 mm/yr faster thinning"
pd2 <- mk("JM_A_ORSD", "JM_A_ORSD_lo", "JM_A_ORSD_hi",
          "clogit_A_ORSD", "clogit_A_ORSD_lo", "clogit_A_ORSD_hi")
pd2$scale <- "per 1 SD faster thinning"
pd <- rbind(pd1, pd2)

ordl <- c(as.character(tab$label[tab$region == "Medial"]),
          as.character(tab$label[tab$region == "Lateral"]))
pd$label <- factor(pd$label, levels = rev(ordl))
pd$type <- factor(pd$type, levels = c("Joint model", "Conditional logistic"))
pd$region <- factor(pd$region, levels = c("Medial", "Lateral"))

g <- ggplot(pd, aes(x = est, y = label, colour = region, shape = type)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .2, linewidth = .5,
                 position = position_dodge(width = .55)) +
  geom_point(size = 2.4, position = position_dodge(width = .55)) +
  facet_wrap(~ scale, scales = "free_x") +
  scale_x_log10() +
  scale_colour_manual(values = c(Medial = "#C0392B", Lateral = "#2E6DA4")) +
  scale_shape_manual(values = c("Joint model" = 16, "Conditional logistic" = 1)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA)) +
  labs(x = "OR for knee replacement (same 191 matched pairs)",
       y = NULL,
       title = "Medial-to-lateral gradient in the pre-operative thinning-rate association",
       subtitle = "Within-pair conditional joint model (filled) vs conditional logistic on the identical knees (open)")
ggsave(file.path(BASE, "fig_jm_medial_lateral.png"), g,
       width = 10, height = 4.6, dpi = 200)

saveRDS(list(fits = keep_fits, level_only = f0, table = tab,
             bf = c(aS = bfS, aV = bfV)),
        file.path(BASE, "jm_multimetric_fit.rds"))
cat("\n-> jm_multimetric.csv / fig_jm_medial_lateral.png / jm_multimetric_fit.rds\n")
cat("\n================ GRADIENT SUMMARY (Set A, 191 pairs) ================\n")
print(tab[, c("label", "region", "JM_A_ORSD", "JM_A_ORSD_lo", "JM_A_ORSD_hi",
              "clogit_A_ORSD", "clogit_A_ORSD_lo", "clogit_A_ORSD_hi",
              "JM_A_OR01", "clogit_A_OR01", "JM_B_ORSD")], row.names = FALSE)
