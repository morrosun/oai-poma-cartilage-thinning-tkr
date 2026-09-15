# ============================================================================
#  31_cfg_poma_lcmm.R
#  Config for Scripts/common/30_class_continuum_check.R
#
#  The POMA cMFTC LCMM four-class solution -- the original hard-coded case.
#  This config exists to prove that the generic engine reproduces
#  28_lcmm_continuum_check.R exactly.  Numbers are compared in
#  Scripts/common/33_compare_poma.R
#
#  Run:
#    Rscript Scripts/common/30_class_continuum_check.R Scripts/poma/31_cfg_poma_lcmm.R
# ============================================================================

FIT <- readRDS("D:/BaiduSyncdisk/OAI/Analysis/poma_pilot/lcmm_fit.rds")

CFG <- list(
  label  = "POMA cMFTC LCMM 4-class solution (cMFTC, matched TKR pairs)",
  tag    = "poma_generic",
  outdir = "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot",

  ## ---- (a) the per-knee table: class + the continuous axis ----
  unit      = FIT$knee,          # 389 knees with a class
  unit_col  = "unit",
  class_col = "cl",

  ## the axis the class might merely be thresholding.  `eblup` is the fitted
  ## per-knee thickness slope in mm/yr (negative = thinning); `ols_slope` is
  ## the unshrunk per-knee OLS slope carried along as a robustness axis.
  exposures = c("eblup", "ols_slope"),
  exposure_labels = c(
    eblup     = "cMFTC thinning rate, EBLUP (mm/yr, + = faster)",
    ols_slope = "cMFTC thinning rate, OLS unshrunk (mm/yr, + = faster)"),
  signs = c(eblup = -1, ols_slope = -1),

  ## severity level carried along for the eta-squared contrast
  levels       = "b0",
  level_labels = c(b0 = "baseline cartilage thickness"),

  rapid_k = 2,                   # the two fastest classes = "the rapid tail"

  ## ---- (b) the outcome table: matched case-control pairs ----
  out       = FIT$pair_data,     # 382 knees = 191 complete pairs
  out_model = "clogit",
  out_cols  = list(unit = "unit", case = "case", strata = "strata")
)
