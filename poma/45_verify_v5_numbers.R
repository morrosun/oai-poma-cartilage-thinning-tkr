# ============================================================================
# 45 : audit every number that v5 of the manuscript introduces or changes
#
#   Nothing is recomputed from models here: the point is to confirm that what
#   the manuscript now says is what the stored result objects actually contain.
#   v5 replaced the sampler and re-ran the joint model, the gradient test and
#   the cross-pipeline corroboration, so the numbers below are the *audited*
#   values and every one of them must be traceable to a file on disk.
#
#   Run with:  Rscript 45_verify_v5_numbers.R > _audit_v5.log 2>&1
# ============================================================================
BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
q3 <- function(x) sprintf("%.3f (%.3f, %.3f)", median(x), quantile(x, .025), quantile(x, .975))

cat("========================================================================\n")
cat("audit of the numbers in Manuscript_EN_OAC_v5 (post-sampler-audit)\n")
cat("========================================================================\n")

## ---- 0. provenance: which sampler produced each stored fit -----------------
cat("\n[0] provenance\n")
cat("  shared sampler      : Scripts/poma/jm_core.R      (corrected)\n")
cat("  legacy sampler kept : Scripts/poma/jm_core_legacy_sampler.R\n")
cat("  validators          : Scripts/poma/40 / 41 / 42 / 43 / 44\n")
cat("  this audit          : Scripts/poma/45\n")

## ---- 1. Abstract + Section 3.2/3.3 : the primary conditional-logistic OR ---
cat("\n[1] Abstract & Sections 3.2-3.3 : conditional logistic on the 191 pairs\n")
mm <- read.csv(file.path(BASE, "jm_multimetric.csv"))
cm <- mm[mm$metric == "cMFTC_ThCtAB_aMe", ]
cat(sprintf("  clogit OR/SD = %.2f (%.2f-%.2f), clogit OR/0.1 = %.2f (%.2f-%.2f), P = %.3g\n",
            cm$clogit_A_ORSD, cm$clogit_A_ORSD_lo, cm$clogit_A_ORSD_hi,
            cm$clogit_A_OR01, cm$clogit_A_OR01_lo, cm$clogit_A_OR01_hi, cm$clogit_A_p))
cat("  manuscript (abstract, 3.2, Table 2 M0) : 2.08 (1.51 to 2.86), P = 7.2e-06\n")
cat(sprintf("  clogitAdj OR/SD = %.2f (%.2f-%.2f), P = %.3g\n",
            cm$clogitAdj_A_ORSD, cm$clogitAdj_A_ORSD_lo, cm$clogitAdj_A_ORSD_hi, cm$clogitAdj_A_p))
cat("  manuscript (3.2, Table 2 M6)          : 2.16 (1.55 to 3.03)\n")

## ---- 2. Table 2 panel B : the joint model on the same 191 pairs -----------
cat("\n[2] Table 2 panel B : joint model, unadjusted and + within-pair covariates\n")
f <- readRDS(file.path(BASE, "jmfit_cMFTC_ThCtAB_aMe.rds"))
for (nm in c("A", "Aadj")) {
  d <- f[[nm]]$draws
  cat(sprintf("  %-5s : OR/SD %s | OR/0.1 %s\n", nm,
              q3(exp(-d$aS * d$SDslo)), q3(exp(-0.1 * d$aS))))
}
cat("  manuscript unadjusted : OR/SD 3.04 (1.63-9.47), OR/0.1 1.88 (1.32-3.68)\n")
cat("  manuscript adjusted   : OR/SD 4.14 (1.91-31.75), OR/0.1 2.25 (1.44-7.47)\n")

## ---- 3. Table 3 : the six compartments, both scales -----------------------
cat("\n[3] Table 3 : compartment-specific joint-model associations\n")
for (i in seq_len(nrow(mm))) {
  cat(sprintf("  %-34s n=%3d  OR/SD %.2f (%.2f-%.2f) | OR/0.1 %.2f (%.2f-%.2f)\n",
              substr(mm$label[i], 1, 34), mm$n_pairs_A[i],
              mm$JM_A_ORSD[i], mm$JM_A_ORSD_lo[i], mm$JM_A_ORSD_hi[i],
              mm$JM_A_OR01[i], mm$JM_A_OR01_lo[i], mm$JM_A_OR01_hi[i]))
}

## ---- 4. Table 3 lower panels : the three formal contrasts -----------------
cat("\n[4] Table 3 : formal contrast tests (bivariate / composite / within-knee)\n")
g <- read.csv(file.path(BASE, "gradient_contrast_test.csv"))
for (i in seq_len(nrow(g)))
  cat(sprintf("  %-28s medial %.3f (%.3f-%.3f) | lateral %.3f (%.3f-%.3f) | ratio %.3f (%.3f-%.3f) p=%.4f\n",
              g$contrast[i], g$OR_medial[i], g$OR_medial_lo[i], g$OR_medial_hi[i],
              g$OR_lateral[i], g$OR_lateral_lo[i], g$OR_lateral_hi[i],
              g$OR_ratio[i], g$OR_ratio_lo[i], g$OR_ratio_hi[i], g$p[i]))
lg <- readLines(file.path(BASE, "10b_gradient_test_rerun.log"))
cat("  -- within-knee difference lines from 10b_gradient_test_rerun.log --\n")
for (l in lg) if (grepl("OR / 1 SD faster|P\\(aS<0\\)|DIFF_", l)) cat("   ", trimws(l), "\n")
cat("  manuscript cMFTC-cLF row : 2.15 (1.33-4.11) / 1.51 (1.16-2.15), P=0.999, clogit 1.79 (1.36-2.36)\n")
cat("  manuscript MFTC-cLF  row : 0.89 (0.60-1.28) / 0.94 (0.77-1.13), P=0.271, clogit 1.06 (0.85-1.31)\n")

## ---- 5. Table 6 : cross-pipeline corroboration ----------------------------
cat("\n[5] Table 6 : cross-pipeline corroboration (148 pairs / 296 knees)\n")
ct <- read.csv(file.path(BASE, "cth_clogit_models.csv"))
for (i in seq_len(nrow(ct))) {
  r <- ct[i, ]
  cat(sprintf("  %-46s %-11s OR %.2f (%.2f-%.2f)  inverted %.2f (%.2f-%.2f)  P %.3g\n",
              substr(r$model, 1, 46), r$term, r$OR, r$lo, r$hi,
              1 / r$OR, 1 / r$hi, 1 / r$lo, r$p))
}
co <- read.csv(file.path(BASE, "cth_concordance.csv"))
cat(sprintf("  concordance : Pearson_pre %.3f | Spearman %.3f | diagonal %.1f%% | n=%d\n",
            co$value[co$metric == "Pearson_pre"], co$value[co$metric == "Spearman_pre"],
            co$value[co$metric == "tertile_agreement_pct"], co$value[co$metric == "n"]))
sh <- read.csv(file.path(BASE, "cth_shrinkage.csv"))
for (i in seq_len(nrow(sh)))
  cat(sprintf("  shrinkage   : %-20s %-12s mean paired diff %+.3f\n",
              sh$metric[i], sh$estimate[i], sh$mean_paired_diff[i]))
cat("  manuscript : Pearson 0.511, Spearman 0.356, diagonal 46.6%, extreme discordance 12.8%\n")
cat("  manuscript : attenuation CTh-Score 0.585 (+3.567 -> +2.087), cMFTC 0.619 (-0.150 -> -0.093)\n")

## ---- 6. Supplementary Results : sampler adequacy of the reported chain -----
cat("\n[6] Supplementary Results : sampler adequacy (primary cMFTC chain)\n")
cd <- read.csv(file.path(BASE, "chain_diagnostics_primary.csv"))
print(cd, row.names = FALSE, digits = 5)
r <- read.csv(file.path(BASE, "mcmc_multichain_rhat.csv"))
print(r, row.names = FALSE, digits = 4)
cat(sprintf("  worst R-hat over all parameters = %.4f (manuscript says 1.0002)\n", max(r$Rhat)))
cat(sprintf("  cross-chain ESS for aS = %.0f (manuscript says 4 451)\n",
            r$ESS_cross_chain[r$parameter == "aS"]))
cc <- read.csv(file.path(BASE, "mcmc_posterior_correlation.csv"))
cat(sprintf("  pooled cor(aV,aS) = %.3f ; per chain %s (manuscript says -0.58, -0.56 to -0.60)\n",
            cc$cor_aV_aS[cc$chain == "pooled"],
            paste(sprintf("%.2f", cc$cor_aV_aS[cc$chain != "pooled"]), collapse = ", ")))

## ---- 7. Supplementary Methods : the sampler audit evidence -----------------
cat("\n[7] Supplementary Methods + Table S4 : sampler audit evidence\n")
di <- read.csv(file.path(BASE, "defect1_isolation.csv"))
print(di, row.names = FALSE, digits = 6)
cat("  manuscript : corrected kernel matches the exact reference to within 1.6 MC-SE\n")
cat("               the v4 kernel deviates by 9 (E[eta]), 13 (E[plogis]) and 25 (E[log plogis]) MC-SE\n")
ab <- read.csv(file.path(BASE, "sampler_ablation.csv"))
agg <- aggregate(cbind(ORSD, aS) ~ label, data = ab, FUN = mean)
print(agg, row.names = FALSE, digits = 4)
cat("  manuscript : switching the random-effect correction alone moves OR per SD 2.60 -> 3.04\n")

## ---- 8. Supplementary Results : prior sensitivity of the joint model ------
cat("\n[8] Supplementary Results : prior sensitivity of the joint model\n")
pj <- read.csv(file.path(BASE, "prior_sensitivity_joint.csv"))
print(pj[pj$term == "OR per 1 SD faster thinning", ], row.names = FALSE)
cat("  manuscript : 2.56 / 2.92 / 3.04 / 3.10 for N(0,5^2), N(0,10^2), N(0,20^2), N(0,50^2)\n")
cat("               inverse-Wishart scale 0.5 / 2 against 1 -> 3.038 / 3.036 against 3.037\n")
pb <- read.csv(file.path(BASE, "prior_sensitivity_bclogit.csv"))
print(pb[pb$prior == "N(0, 2^2)", c("term", "OR", "lo", "hi", "P_OR_gt1")],
      row.names = FALSE, digits = 4)

## ---- 9. Table S1 : the set-sensitivity rows -------------------------------
cat("\n[9] Table S1 : set sensitivity (J = number of pre-index visits required)\n")
sv <- read.csv(file.path(BASE, "sens_minvisits_primary.csv"))
keep <- c("min_visits", "n_pairs", "n_knees", "visits_per_knee",
          "JM_OR01", "JM_OR01_lo", "JM_OR01_hi", "JM_ORSD", "JM_ORSD_lo", "JM_ORSD_hi",
          "clogit_ORSD", "clogit_ORSD_lo", "clogit_ORSD_hi",
          "clogitAdj_ORSD", "clogitAdj_ORSD_lo", "clogitAdj_ORSD_hi",
          "LR_add_rate", "LR_add_level")
print(sv[, keep], row.names = FALSE, digits = 3)
cat("  manuscript J=2 : 1.88 (1.32-3.68) / 3.04 (1.63-9.47) / clogit 2.08 (1.51-2.86)\n")
cat("  manuscript J=3 : 1.79 (1.26-3.45) / 2.76 (1.51-8.39) / clogit 2.13 (1.47-3.07), adj 2.30 (1.53-3.45)\n")

## ---- 10. cross-pipeline visit density --------------------------------------
cat("\n[10] Supplementary Results : pre-index visit density in the score pipeline\n")
vd <- read.csv(file.path(BASE, "cth_visit_density.csv"))
print(vd, row.names = FALSE, digits = 4)
cat("  manuscript : 3.57 vs 4.22 pre-index measurements per knee\n")

cat("\n========================================================================\n")
cat("reading of the audit: every value quoted in v5 must appear above\n")
cat("========================================================================\n")
