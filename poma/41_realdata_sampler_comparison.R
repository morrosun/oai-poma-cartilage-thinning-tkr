# ============================================================================
# 41 : does the corrected sampler change the REAL POMA results?
#
#   40_validate_sampler.R shows, on simulated data, that the corrected kernel and
#   the legacy kernel target the same posterior.  That is necessary but not
#   sufficient: the manuscript numbers were produced by the legacy kernel on the
#   real 191 matched pairs, and defect (3) only bites when Zd is supplied.
#
#   Here both kernels are run on the REAL Set A of the primary metric, with and
#   without the within-pair covariate differences, four seeds each.  Four seeds
#   give an honest Monte-Carlo error, so "no material difference" can be stated
#   as  |median_corrected - median_legacy|  vs  the seed-to-seed spread, rather
#   than as an eyeball comparison of two chains.
#
#   -> sampler_comparison_realdata.csv
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
SDIR <- "D:/BaiduSyncdisk/OAI/Scripts/poma"
source(file.path(SDIR, "jm_core.R"))                                   # corrected
LEG <- new.env()
sys.source(file.path(SDIR, "jm_core_legacy_sampler.R"), envir = LEG)   # legacy

L <- read.csv(file.path(BASE, "poma_long.csv"))
U <- read.csv(file.path(BASE, "poma_analysis_long.csv"))
COVARS <- c("V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")

## --- exactly the Set A construction used by 10_jm_multimetric.R --------------
metric <- "cMFTC_ThCtAB_aMe"
sub <- L[!is.na(L[[metric]]), ]
lg  <- data.frame(unit = sub$unit, year = sub$years, y = sub[[metric]],
                  newstrata = sub$newstrata, case = sub$case)
sv  <- U[, c("unit", "newstrata", "case", "t_index_months", COVARS, "V00PASE")]
sv  <- sv[sv$unit %in% lg$unit, ]
sv  <- sv[complete.cases(sv[, COVARS]), ]
lg  <- lg[lg$unit %in% sv$unit, ]
sv$t_idx <- sv$t_index_months / 12
A <- prune_pairs(lg, sv, 2)
dtA <- build(A$lg, A$sv)
uA  <- U[match(A$sv$unit, U$unit), ]
cat(sprintf("Set A : %d pairs / %d knees / %d rows\n",
            length(dtA$cas), dtA$n, nrow(A$lg)))

Zd <- cbind(kl_diff    = uA$V00XRKL[dtA$cas]    - uA$V00XRKL[dtA$ctl],
            womac_diff = uA$WOMAC_pain[dtA$cas] - uA$WOMAC_pain[dtA$ctl])

## --- one row of summaries per fit ------------------------------------------
summ_fit <- function(f, sampler, model, seed) {
  d <- f$draws
  q <- function(x) c(med = median(x), lo = quantile(x, .025), hi = quantile(x, .975))
  o1  <- q(exp(-0.1 * d$aS))
  oSD <- q(exp(-d$aS * d$SDslo))
  oLv <- q(exp(-d$aV * d$SDval))
  data.frame(sampler = sampler, model = model, seed = seed,
             aS_med = median(d$aS), aS_lo = quantile(d$aS, .025),
             aS_hi = quantile(d$aS, .975),
             aV_med = median(d$aV),
             beta1 = median(d$beta1),
             sd_latent_slope = median(d$SDslo),
             OR01 = o1["med"], OR01_lo = o1["lo"], OR01_hi = o1["hi"],
             ORSD = oSD["med"], ORSD_lo = oSD["lo"], ORSD_hi = oSD["hi"],
             ORSDlev = oLv["med"], ORSDlev_lo = oLv["lo"], ORSDlev_hi = oLv["hi"],
             gamma_kl = if (is.null(d$gamma)) NA_real_ else median(d$gamma[, 1]),
             gamma_wom = if (is.null(d$gamma)) NA_real_ else median(d$gamma[, 2]),
             acc_re = f$acc_b, stringsAsFactors = FALSE)
}

SEEDS <- c(2026, 2027, 2028, 2029)
res <- list(); k <- 0
for (sd_i in SEEDS) {
  for (md in c("crude", "adjusted")) {
    Z <- if (md == "adjusted") Zd else NULL
    cat(sprintf("  corrected / %-8s / seed %d ...\n", md, sd_i)); flush.console()
    fC <- fit_jm(dtA, NITER = 30000, BURN = 5000, THIN = 5, seed = sd_i,
                 Zd = Z, verbose = FALSE)
    k <- k + 1; res[[k]] <- summ_fit(fC, "corrected", md, sd_i)
    cat(sprintf("  legacy    / %-8s / seed %d ...\n", md, sd_i)); flush.console()
    fL <- LEG$fit_jm(dtA, NITER = 30000, BURN = 5000, THIN = 5, seed = sd_i,
                     Zd = Z, verbose = FALSE)
    k <- k + 1; res[[k]] <- summ_fit(fL, "legacy", md, sd_i)
  }
}
tab <- do.call(rbind, res)
write.csv(tab, file.path(BASE, "sampler_comparison_realdata.csv"), row.names = FALSE)

## --- report: seed spread vs sampler difference ------------------------------
cat("\n================ REAL DATA: corrected vs legacy ================\n")
for (md in c("crude", "adjusted")) {
  cat(sprintf("\n--- %s ---\n", md))
  for (v in c("aS_med", "aV_med", "OR01", "ORSD", "ORSDlev", "gamma_kl", "gamma_wom")) {
    sc <- tab[tab$model == md, ]
    cc <- sc[[v]][sc$sampler == "corrected"]; ll <- sc[[v]][sc$sampler == "legacy"]
    if (all(is.na(cc))) next
    seed_sd <- mean(c(sd(cc), sd(ll)), na.rm = TRUE)
    bias    <- abs(mean(cc) - mean(ll))
    cat(sprintf("  %-14s corrected %8.4f (seed sd %.4f) | legacy %8.4f (seed sd %.4f) | diff %7.4f  %s\n",
                v, mean(cc), sd(cc), mean(ll), sd(ll), bias,
                ifelse(bias < 2 * seed_sd, "within MC error", "EXCEEDS MC error")))
  }
}
cat("\n-> sampler_comparison_realdata.csv\n")
