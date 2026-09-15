# ============================================================================
# 10b : formal test of the medial-to-lateral gradient
#
#   Comparing two separate ORs with non-overlapping CIs is NOT a test of a
#   difference, and the posteriors of two separate joint models are correlated
#   (same knees).  Three proper contrasts are computed here, all on the same
#   191 complete matched pairs used in 06 / 10:
#
#   (1) bivariate conditional logistic : case ~ z_medial + z_lateral + strata.
#       The contrast exp(b_med - b_lat) with SE from the full covariance is a
#       direct test of "medial change adds beyond lateral change".
#   (2) composite medial vs composite lateral score (mean of standardised
#       thinning rates within each region).
#   (3) within-knee medial-minus-lateral difference run through the SAME
#       within-pair conditional joint model.  A significant rate association
#       for the difference is compartment-specific evidence: if a global
#       (both-compartment) thinning process drove the outcome the difference
#       would carry no signal.
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
source("D:/BaiduSyncdisk/OAI/Scripts/poma/jm_core.R")
suppressPackageStartupMessages(library(survival))

L <- read.csv(file.path(BASE, "poma_long.csv"))
U <- read.csv(file.path(BASE, "poma_analysis_long.csv"))
COVARS <- c("V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")

## within-knee medial-minus-lateral differences, evaluated at the same visit
L$DIFF_MFTC  <- L$MFTC_ThCtAB_aMe  - L$cLFTC_ThCtAB_aMe
L$DIFF_cMFTC <- L$cMFTC_ThCtAB_aMe - L$cLF_ThCtAB_aMe

MED <- c("cMFTC_ThCtAB_aMe", "MFTC_ThCtAB_aMe", "MT_ThCtAB_aMe", "cMF_ThCtAB_aMe")
LAT <- c("cLFTC_ThCtAB_aMe", "cLF_ThCtAB_aMe")

set_A <- function(ycol) {
  sub <- L[!is.na(L[[ycol]]), ]
  lg  <- data.frame(unit = sub$unit, year = sub$years, y = sub[[ycol]],
                    newstrata = sub$newstrata, case = sub$case)
  sv  <- U[, c("unit", "newstrata", "case", "t_index_months", COVARS)]
  sv  <- sv[sv$unit %in% lg$unit, ]
  sv  <- sv[complete.cases(sv[, COVARS]), ]
  lg  <- lg[lg$unit %in% sv$unit, ]
  sv$t_idx <- sv$t_index_months / 12
  prune_pairs(lg, sv, 2)
}
zloss <- function(x) (mean(x, na.rm = TRUE) - x) / sd(x, na.rm = TRUE)   # + = faster thinning

cat("========================================================================\n")
cat("formal medial-vs-lateral contrasts, Set A (191 complete matched pairs)\n")
cat("========================================================================\n")

## ---------------------------------------------------------------- (1) -----
A <- set_A(MED[1]); uA <- U[match(A$sv$unit, U$unit), ]
cat(sprintf("Set A : %d pairs / %d knees\n\n", length(unique(A$sv$newstrata)), nrow(uA)))

z <- data.frame(case = uA$case, strata = uA$newstrata)
for (m in c(MED, LAT)) z[[m]] <- zloss(uA[[paste0(m, "__eblup")]])

cat("---- (1) bivariate conditional logistic : medial vs lateral ----\n")
cat(sprintf("%-42s %9s %9s %9s %9s %10s\n", "contrast",
            "OR medial", "lo", "hi", "OR ratio", "p (diff)"))
pairs_to_test <- list(c("cMFTC_ThCtAB_aMe", "cLF_ThCtAB_aMe"),
                      c("cMFTC_ThCtAB_aMe", "cLFTC_ThCtAB_aMe"),
                      c("MFTC_ThCtAB_aMe",  "cLF_ThCtAB_aMe"),
                      c("MFTC_ThCtAB_aMe",  "cLFTC_ThCtAB_aMe"))
res1 <- list()
for (pr in pairs_to_test) {
  f <- clogit(as.formula(paste("case ~", pr[1], "+", pr[2], "+ strata(strata)")),
              data = z, method = "exact")
  s <- summary(f)$coef; V <- vcov(f)
  b <- s[1, 1]; e <- s[1, 3]
  d <- s[1, 1] - s[2, 1]; sd_d <- sqrt(V[1, 1] + V[2, 2] - 2 * V[1, 2])
  lab <- sprintf("%s vs %s", sub("_ThCtAB_aMe", "", pr[1]), sub("_ThCtAB_aMe", "", pr[2]))
  cat(sprintf("%-42s %9.3f %9.3f %9.3f %9.3f %10.4f\n", lab,
              exp(b), exp(b - 1.96 * e), exp(b + 1.96 * e),
              exp(d), 2 * pnorm(-abs(d / sd_d))))
  res1[[lab]] <- data.frame(contrast = lab, OR_medial = exp(b),
                            OR_medial_lo = exp(b - 1.96 * e), OR_medial_hi = exp(b + 1.96 * e),
                            OR_lateral = exp(s[2, 1]),
                            OR_ratio = exp(d), OR_ratio_lo = exp(d - 1.96 * sd_d),
                            OR_ratio_hi = exp(d + 1.96 * sd_d), p = 2 * pnorm(-abs(d / sd_d)))
}

## ---------------------------------------------------------------- (2) -----
cat("\n---- (2) composite medial vs composite lateral score ----\n")
z$zM <- rowMeans(z[, MED]); z$zL <- rowMeans(z[, LAT])
f <- clogit(case ~ zM + zL + strata(strata), data = z, method = "exact")
s <- summary(f)$coef; V <- vcov(f)
d <- s[1, 1] - s[2, 1]; sd_d <- sqrt(V[1, 1] + V[2, 2] - 2 * V[1, 2])
cat(sprintf("  composite MEDIAL  OR/SD = %.3f (%.3f-%.3f)\n", exp(s[1, 1]),
            exp(s[1, 1] - 1.96 * s[1, 3]), exp(s[1, 1] + 1.96 * s[1, 3])))
cat(sprintf("  composite LATERAL OR/SD = %.3f (%.3f-%.3f)\n", exp(s[2, 1]),
            exp(s[2, 1] - 1.96 * s[2, 3]), exp(s[2, 1] + 1.96 * s[2, 3])))
cat(sprintf("  medial/lateral OR ratio  = %.3f (%.3f-%.3f), p = %.4f\n",
            exp(d), exp(d - 1.96 * sd_d), exp(d + 1.96 * sd_d), 2 * pnorm(-abs(d / sd_d))))
res1[["composite"]] <- data.frame(contrast = "composite medial vs lateral",
                                  OR_medial = exp(s[1, 1]), OR_medial_lo = exp(s[1, 1] - 1.96*s[1,3]),
                                  OR_medial_hi = exp(s[1, 1] + 1.96*s[1,3]),
                                  OR_lateral = exp(s[2, 1]),
                                  OR_ratio = exp(d), OR_ratio_lo = exp(d - 1.96*sd_d),
                                  OR_ratio_hi = exp(d + 1.96*sd_d),
                                  p = 2 * pnorm(-abs(d / sd_d)))

## ---------------------------------------------------------------- (3) -----
cat("\n---- (3) within-knee medial-minus-lateral difference, joint model ----\n")
for (dc in c("DIFF_cMFTC", "DIFF_MFTC")) {
  B <- set_A(dc); dt <- build(B$lg, B$sv)
  fit <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 2030, verbose = FALSE)
  cat(sprintf("\n  %s : %d pairs / %d knees\n", dc, length(dt$cas), dt$n))
  report(fit, paste0("JM on ", dc, " (medial minus lateral thinning)"))
  eb <- uA[[paste0(sub("DIFF_", "", dc) , "_ThCtAB_aMe__eblup")]]
  ## clogit on the EBLUP difference, same knees
  Bl <- U[match(B$sv$unit, U$unit), ]
  if (dc == "DIFF_cMFTC") ed <- Bl$cMFTC_ThCtAB_aMe__eblup - Bl$cLF_ThCtAB_aMe__eblup
  else                    ed <- Bl$MFTC_ThCtAB_aMe__eblup  - Bl$cLFTC_ThCtAB_aMe__eblup
  f2 <- clogit(case ~ zl + strata(newstrata),
               data = data.frame(case = Bl$case, newstrata = Bl$newstrata,
                                 zl = zloss(ed)), method = "exact")
  s2 <- summary(f2)$coef
  cat(sprintf("  [clogit on the same difference] OR/SD = %.3f (%.3f-%.3f), p = %.4f\n",
              exp(s2[1, 1]), exp(s2[1, 1] - 1.96 * s2[1, 3]),
              exp(s2[1, 1] + 1.96 * s2[1, 3]), s2[1, 5]))
}

out <- do.call(rbind, res1)
write.csv(out, file.path(BASE, "gradient_contrast_test.csv"), row.names = FALSE)
cat("\n-> gradient_contrast_test.csv\n")
