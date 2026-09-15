# ============================================================================
# 10d : verify the joint model's LONGITUDINAL submodel against nlme::lme
#       on the PRIMARY metric (cMFTC) and on the SAME Set A used everywhere else.
#       This is the mutual-validation check quoted in the manuscript Results.
# ============================================================================
BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
source("D:/BaiduSyncdisk/OAI/Scripts/poma/jm_core.R")
suppressPackageStartupMessages({ library(nlme); library(survival) })

L <- read.csv(file.path(BASE, "poma_long.csv"))
U <- read.csv(file.path(BASE, "poma_analysis_long.csv"))
COVARS <- c("V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")

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

for (m in c("cMFTC_ThCtAB_aMe", "MFTC_ThCtAB_aMe")) {
  A <- make_input(m, 2, TRUE)
  lg <- A$lg
  cat(sprintf("\n================ %s : Set A = %d knees / %d rows ================\n",
              m, length(unique(lg$unit)), nrow(lg)))
  f <- lme(y ~ year, random = ~ year | unit, data = lg, method = "ML")
  b  <- fixef(f); sdv <- as.numeric(VarCorr(f)[, "StdDev"])
  cat(sprintf("  lme   beta0 %.6f  beta1 %+.6f | resid SD %.6f | SD(u0) %.6f SD(u1) %.6f\n",
              b[1], b[2], f$sigma, sdv[1], sdv[2]))
  cat(sprintf("  lme   beta1 as printed to 3 dp : %+.3f ; resid SD to 3 dp : %.3f\n",
              round(b[2], 3), round(f$sigma, 3)))
}
cat("\n-> compare these with the joint-model longitudinal estimates in jm_multimetric.log\n")
