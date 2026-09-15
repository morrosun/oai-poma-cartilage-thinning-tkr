# ============================================================================
# 13 : methodological result 2 -- an EBLUP slope used as an exposure must be
#      accompanied by a differential-shrinkage control
#
#   Shrinkage maps a noisy per-knee slope to  EBLUP_i = mu + lam_i (OLS_i - mu),
#   so a within-pair contrast of EBLUP slopes is NOT a within-pair contrast of
#   slopes: it equals lam_c*OLS_c - lam_j*OLS_j - mu*(lam_c - lam_j).  The last
#   term vanishes only if the two members of a pair shrink equally, i.e. only if
#   they carry the same amount of information.
#
#   In POMA they do not: cases have 3.50 pre-index visits and controls 3.03
#   (paired mean difference +0.47), so lam_c > lam_j systematically.  A
#   simulation with a KNOWN truth separates the consequences:
#
#     S0  null (aS = 0), observed visit asymmetry preserved
#           -> does differential shrinkage manufacture an association?
#     S1  POMA design (aS = calibrated)
#           -> attenuation of the EBLUP coefficient relative to the truth, and
#              the same quantity for the unshrunk OLS slope
#     S2  asymmetry exaggerated (controls truncated to 2 visits)
#           -> how fast the EBLUP coefficient drifts once the two arms shrink
#              by different amounts
# ============================================================================
source("D:/BaiduSyncdisk/OAI/Scripts/sim/sim_core.R")

BASE   <- "D:/BaiduSyncdisk/OAI"
OUT    <- file.path(BASE, "Analysis/sim")
design <- read_design(file.path(OUT, "sim_visit_schedule.csv"))
CAL    <- readRDS(file.path(OUT, "sim_as.rds"))
aS <- CAL$aS; CONST <- CAL$const

LOG <- list(); P <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG[[length(LOG)+1]] <<- s }

REPS <- as.integer(Sys.getenv("SIM_REPS", "400"))
P("==========================================================================")
P("SIMULATION 2 : differential shrinkage  (true aS = %.3f, %d replicates)", aS, REPS)
P("==========================================================================")

## -------------------------------------------------------- analytic shrinkage --
shrinkage_factors <- function(dat, t2, s2) {
  sx <- tapply(dat$t, dat$unit, function(t) sum((t - mean(t))^2))
  k1 <- dat[!duplicated(dat$unit), c("unit", "pair", "case")]
  k1$sxx <- as.numeric(sx[as.character(k1$unit)])
  k1$lam <- t2 / (t2 + s2 / k1$sxx)
  k1
}

one_rep <- function(seed, aS_use, visits_ctl = NULL) {
  d <- gen_data(design, NULL, CONST$tau, CONST$mu1, CONST$b0mu, CONST$b0sd,
                CONST$sigma_e, aS = aS_use, seed = seed, visits_ctl = visits_ctl)
  kn <- add_truth(knee_slopes(d), d)
  m  <- lme(y ~ t, random = ~ 1 + t | unit, data = d,
            control = lmeControl(opt = "optim"))
  vc <- as.numeric(VarCorr(m)[, 2])^2       # 1 = intercept, 2 = slope, 3 = residual
  sf <- shrinkage_factors(d, vc[2], vc[3])
  sf <- sf[match(kn$unit, sf$unit), ]

  ft <- clogit_fit(kn, "slope_true")
  fe <- clogit_fit(kn, "slope_eblup")
  fo <- clogit_fit(kn, "slope_ols")

  ## paired differences in LOSS units (positive = cases thinner); rows are
  ## ordered case-first within pair by knee_slopes()
  ag <- aggregate(cbind(slope_true, slope_eblup, slope_ols) ~ pair, kn,
                  function(x) -(x[1] - x[2]))

  c(b_true = ft[["beta"]], s_true = ft[["se"]],
    b_eblup = fe[["beta"]], s_eblup = fe[["se"]],
    b_ols = fo[["beta"]], s_ols = fo[["se"]],
    pd_true = mean(ag$slope_true), pd_eblup = mean(ag$slope_eblup),
    pd_ols = mean(ag$slope_ols),
    lam_case = mean(sf$lam[sf$case == 1]), lam_ctl = mean(sf$lam[sf$case == 0]),
    v_ols = var(kn$slope_ols), v_eblup = var(kn$slope_eblup),
    sd_true = sd(kn$slope_true))
}

run <- function(aS_use, visits_ctl, seed0) {
  L <- lapply(seq_len(REPS), function(i)
    try(one_rep(seed0 + i, aS_use, visits_ctl), silent = TRUE))
  L <- L[!vapply(L, inherits, TRUE, "try-error")]
  if (!length(L)) return(NULL)
  nm <- names(L[[1]])
  M <- t(vapply(L, function(z) as.numeric(z[nm]), numeric(length(nm))))
  colnames(M) <- nm
  M
}

report <- function(res, tag, truth) {
  P("")
  P("  %s", tag)
  P("  %-30s %9s %9s %9s %9s", "quantity", "mean", "sd",
    if (truth == 0) "type-I" else "bias", if (truth == 0) "" else "ratio")
  for (v in c("b_true", "b_eblup", "b_ols")) {
    b <- res[, v]; s <- res[, sub("^b_", "s_", v)]
    ok <- is.finite(b) & is.finite(s)
    if (truth == 0) {
      rej <- if (any(ok)) mean(abs(b[ok] / s[ok]) > 1.96) else NA
      P("  %-30s %9.3f %9.3f %9.3f %9s", v, mean(b, na.rm = TRUE),
        sd(b, na.rm = TRUE), rej, "")
    } else {
      P("  %-30s %9.3f %9.3f %9.3f %9.3f", v, mean(b, na.rm = TRUE),
        sd(b, na.rm = TRUE), mean(b, na.rm = TRUE) - truth,
        mean(b, na.rm = TRUE) / truth)
    }
  }
  P("  %-30s %9.4f %9.4f", "paired diff, TRUE slope (loss)", mean(res[, "pd_true"]), sd(res[, "pd_true"]))
  P("  %-30s %9.4f %9.4f", "paired diff, EBLUP slope",      mean(res[, "pd_eblup"]), sd(res[, "pd_eblup"]))
  P("  %-30s %9.4f %9.4f", "paired diff, OLS slope",        mean(res[, "pd_ols"]), sd(res[, "pd_ols"]))
  P("  %-30s %9.4f %9.4f", "shrinkage factor lam, cases",   mean(res[, "lam_case"]), sd(res[, "lam_case"]))
  P("  %-30s %9.4f %9.4f", "shrinkage factor lam, controls",mean(res[, "lam_ctl"]),  sd(res[, "lam_ctl"]))
  P("  %-30s %9.5f %9.5f", "Var(OLS slope)",                mean(res[, "v_ols"]), sd(res[, "v_ols"]))
  P("  %-30s %9.5f %9.5f", "Var(EBLUP slope)",              mean(res[, "v_eblup"]), sd(res[, "v_eblup"]))
  P("  %-30s %9.5f %9.5f", "Var(true slope)",               mean(res[, "sd_true"]^2), sd(res[, "sd_true"]^2))
  P("  %-30s %9.5f", "Var(OLS) - Var(true)  [= E Var(est)]",
    mean(res[, "v_ols"]) - mean(res[, "sd_true"]^2))
}

P("")
P("REAL POMA reference values (primary set, 191 pairs):")
P("  conditional logistic, EBLUP slope  beta = +5.180   (sd of EBLUP slope 0.1411)")
P("  conditional logistic, OLS   slope  beta = +2.578   (sd of OLS   slope 0.2523)")
P("  paired difference in LOSS units    EBLUP 0.0743 | OLS 0.1189")
P("  marginal (unstratified)            EBLUP beta = +1.893")

P("\n--- S0 : NULL (aS = 0), observed visit asymmetry preserved ---")
S0 <- run(0, NULL, 10000);  report(S0, "S0  aS = 0", 0)

P("\n--- S1 : POMA design (aS = %.3f) ---", aS)
S1 <- run(aS, NULL, 20000); report(S1, "S1  POMA design", aS)

P("\n--- S2 : controls truncated to 2 visits (asymmetry exaggerated) ---")
S2 <- run(aS, 2, 30000);    report(S2, "S2  differential visit density", aS)

saveRDS(list(S0 = S0, S1 = S1, S2 = S2, aS = aS, log = LOG),
        file.path(OUT, "sim2_shrinkage.rds"))
writeLines(unlist(LOG), file.path(OUT, "sim2_shrinkage_log.txt"))
P("\nwritten -> Analysis/sim/sim2_shrinkage.rds")
