# ============================================================================
# 36 : audit every number that v4 of the manuscript introduces
#
#   Nothing is recomputed from models: the point is to check that what the
#   manuscript now says is what the stored result objects actually contain.
#   Run with:  Rscript 36_verify_v4_numbers.R > _audit_v4.log 2>&1
# ============================================================================
BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
q3 <- function(x) sprintf("%.3f (%.3f, %.3f)", median(x), quantile(x, .025), quantile(x, .975))
ok <- function(cond, msg) cat(sprintf("  %s  %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))

cat("========================================================================\n")
cat("audit of the numbers introduced by Manuscript_EN_OAC_v4\n")
cat("========================================================================\n")

## ---- 1. Table 2 : the interval that was missing --------------------------
cat("\n[1] Table 2, joint model + within-pair KL/pain differences (per 1 SD)\n")
f <- readRDS(file.path(BASE, "jmfit_cMFTC_ThCtAB_aMe.rds"))
for (nm in c("A", "Aadj")) {
  d <- f[[nm]]$draws
  cat(sprintf("  %-5s : aS %s | OR/1SD %s | OR/0.1  %s | OR/1SD level %s\n",
              nm, q3(d$aS), q3(exp(-d$aS * d$SDslo)), q3(exp(-0.1 * d$aS)),
              q3(exp(-d$aV * d$SDval))))
}
da <- f$Aadj$draws
cat(sprintf("  manuscript Table 2 row quotes : 2.97 (1.65-7.32)  -> stored %.2f (%.2f-%.2f)\n",
            median(exp(-da$aS * da$SDslo)),
            quantile(exp(-da$aS * da$SDslo), .025), quantile(exp(-da$aS * da$SDslo), .975)))
cat(sprintf("  manuscript Table 2 unadjusted : 2.60 (1.54-5.67)  -> stored %.2f (%.2f-%.2f)\n",
            median(exp(-f$A$draws$aS * f$A$draws$SDslo)),
            quantile(exp(-f$A$draws$aS * f$A$draws$SDslo), .025),
            quantile(exp(-f$A$draws$aS * f$A$draws$SDslo), .975)))

## ---- 2. Table 3 : the lateral intervals now stored -----------------------
cat("\n[2] Table 3, lateral term of each bivariate conditional logistic model\n")
g <- read.csv(file.path(BASE, "gradient_contrast_test.csv"))
for (i in seq_len(nrow(g)))
  cat(sprintf("  %-32s medial %.3f (%.3f-%.3f) | lateral %.3f (%.3f-%.3f) | ratio %.3f (%.3f-%.3f) p=%.4f\n",
              g$contrast[i], g$OR_medial[i], g$OR_medial_lo[i], g$OR_medial_hi[i],
              g$OR_lateral[i], g$OR_lateral_lo[i], g$OR_lateral_hi[i],
              g$OR_ratio[i], g$OR_ratio_lo[i], g$OR_ratio_hi[i], g$p[i]))

## ---- 3. Table 6 : orientation of the cMFTC rows and model C2 -------------
cat("\n[3] Table 6, cMFTC rows reversed so that OR > 1 means faster progression\n")
ct <- read.csv(file.path(BASE, "cth_clogit_models.csv"))
for (i in seq_len(nrow(ct))) {
  r <- ct[i, ]
  cat(sprintf("  %-46s %-11s OR %.2f (%.2f-%.2f)  inverted %.2f (%.2f-%.2f)  P %.3g\n",
              substr(r$model, 1, 46), r$term, r$OR, r$lo, r$hi,
              1 / r$OR, 1 / r$hi, 1 / r$lo, r$p))
}

## ---- 4. Table 4 : per-class odds ratios, and their prior sensitivity -----
cat("\n[4] Table 4 per-class odds ratios at the reported prior N(0, 2^2)\n")
pb <- read.csv(file.path(BASE, "prior_sensitivity_bclogit.csv"))
print(pb[pb$prior == "N(0, 2^2)", c("term", "OR", "lo", "hi", "P_OR_gt1")],
      row.names = FALSE, digits = 4)
cat("  manuscript Table 4 rows : 1.82 (0.95-3.55) / 2.36 (1.09-5.39) / 37.25 (6.68-427.51)\n")

## ---- 5. Supplementary Results : sampler diagnostics -----------------------
cat("\n[5] Supplementary Results, sampler diagnostics\n")
r <- read.csv(file.path(BASE, "mcmc_multichain_rhat.csv"))
print(r, row.names = FALSE, digits = 4)
cat(sprintf("  worst R-hat = %.4f (manuscript says 1.0001)\n", max(r$Rhat)))
cat(sprintf("  cross-chain ESS for aS = %.0f (manuscript says 8 484)\n",
            r$ESS_cross_chain[r$parameter == "aS"]))
cc <- read.csv(file.path(BASE, "mcmc_posterior_correlation.csv"))
cat(sprintf("  pooled cor(aV, aS) = %.3f ; per chain %s (manuscript says -0.52, -0.49 to -0.54)\n",
            cc$cor_aV_aS[cc$chain == "pooled"],
            paste(sprintf("%.2f", cc$cor_aV_aS[cc$chain != "pooled"]), collapse = ", ")))

## ---- 6. Supplementary Results : prior sensitivity of the joint model -----
cat("\n[6] Supplementary Results, joint-model prior sensitivity\n")
pj <- read.csv(file.path(BASE, "prior_sensitivity_joint.csv"))
print(pj[pj$term == "OR per 1 SD faster thinning", ], row.names = FALSE)
cat("  manuscript says 2.37 / 2.55 / 2.60 / 2.62 for N(0,5^2), N(0,10^2), N(0,20^2), N(0,50^2),\n")
cat("  and 2.600 / 2.599 against 2.601 for the inverse-Wishart scale 0.5 / 2 against 1.\n")

cat("\n========================================================================\n")
cat("reading of the audit: every value quoted in v4 must appear above\n")
cat("========================================================================\n")
