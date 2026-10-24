# ============================================================================
#  Put the joint model and the conditional logistic regression on ONE scale.
#
#  Problem identified in review: the two paths report "per 1 SD" effects that
#  are NOT the same SD.
#    - conditional logistic : the SD of the observed EBLUP rate in Set A
#                             (cMFTC 0.1411 mm/yr)
#    - joint model          : the SD of the LATENT random slope, which differs
#                             from the EBLUP SD because the EBLUP is shrunk,
#                             and which differs by a factor of 4.5 across the
#                             six regions (0.039 - 0.176 mm/yr)
#
#  Consequence: "per 1 SD" odds ratios cannot be ranked ACROSS regions, and
#  the joint-model 3.04 cannot be set beside the conditional-logistic 2.08 as
#  if they measured the same contrast.
#
#  Fix: recompute the joint-model effect at the FIXED reference SD used by the
#  conditional logistic regression (0.141132 mm/yr for cMFTC), directly from
#  the posterior draws, so the two paths become comparable, and report the
#  per-region latent-slope SD alongside every per-SD effect.
#
#  exp(-aS * s) is monotone in aS for fixed s, so credible-interval endpoints
#  transform exactly; the point estimate is the posterior MEDIAN of the
#  transformed draws (not exp(-median(aS) * s)).
# ============================================================================
options(warn = 1)
B <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot/"
p <- function(...) cat(..., "\n", sep = "")

met <- read.csv(paste0(B, "jm_multimetric.csv"), stringsAsFactors = FALSE)

## the reference SD: within-set SD of the EBLUP rate actually used as the
## exposure in the conditional logistic regression, per metric
ref_sd <- setNames(met$sd_eblup, met$metric)

summ <- function(v) c(est = median(v),
                      lo  = unname(quantile(v, 0.025)),
                      hi  = unname(quantile(v, 0.975)))

rows <- list()
for (m in met$metric) {
  fl <- paste0(B, "jmfit_", m, ".rds")
  if (!file.exists(fl)) { p("missing fit for ", m); next }
  x  <- readRDS(fl)

  for (part in c("A", "Aadj")) {
    if (is.null(x[[part]])) next
    dr <- x[[part]]$draws
    if (is.null(dr$aS)) { p("no aS draws for ", m, " / ", part); next }
    aS <- as.numeric(dr$aS)
    s  <- ref_sd[[m]]

    or_fix <- exp(-aS * s)                        # fixed reference SD
    or_01  <- exp(-aS * 0.1)                      # per 0.1 mm/yr
    sd_lat <- if (!is.null(dr$SDslo)) mean(as.numeric(dr$SDslo)) else NA
    or_own <- if (!is.null(dr$SDslo)) exp(-aS * as.numeric(dr$SDslo)) else NA

    rows[[length(rows) + 1]] <- data.frame(
      metric = m, model = part,
      ref_sd = round(s, 6),
      sd_latent_slope = round(sd_lat, 6),
      OR_fixedSD = round(summ(or_fix)["est"], 3),
      OR_fixedSD_lo = round(summ(or_fix)["lo"], 3),
      OR_fixedSD_hi = round(summ(or_fix)["hi"], 3),
      OR_per0p1 = round(summ(or_01)["est"], 3),
      OR_per0p1_lo = round(summ(or_01)["lo"], 3),
      OR_per0p1_hi = round(summ(or_01)["hi"], 3),
      OR_ownSD = if (is.na(sd_lat)) NA else round(summ(or_own)["est"], 3),
      P_OR_gt_1 = round(mean(or_fix > 1), 4),
      stringsAsFactors = FALSE)
  }
}

out <- do.call(rbind, rows)
print(out, row.names = FALSE)
write.csv(out, paste0(B, "jm_fixed_sd_effect.csv"), row.names = FALSE)

p("\n=============================== SUMMARY ===============================")
pr <- out[out$metric == met$metric[1], ]
for (i in seq_len(nrow(pr)))
  p(sprintf("%-5s  OR per fixed SD %.4f = %.2f (%.2f-%.2f)   | OR per 0.1 = %.2f (%.2f-%.2f)   | OR per OWN latent SD %.4f = %s",
            pr$model[i], pr$ref_sd[i], pr$OR_fixedSD[i], pr$OR_fixedSD_lo[i],
            pr$OR_fixedSD_hi[i], pr$OR_per0p1[i], pr$OR_per0p1_lo[i],
            pr$OR_per0p1_hi[i], pr$sd_latent_slope[i], pr$OR_ownSD[i]))
