# ============================================================================
# 55b : RELIABILITY OF THE PER-KNEE OLS SLOPE, BY DESIGN
#
#   The OLS slope of a short trajectory is attenuated by within-knee
#   measurement error (regression dilution), and the attenuation is stronger
#   the fewer visits a knee contributes:
#
#       lambda = Var(true slope) / Var(observed OLS slope)
#              = D[2,2] / sd(OLS slope)^2
#
#   with D[2,2] from the random-slope LME fitted to the SAME rows as the
#   design.  lambda is reported so that the reader can see whether a
#   difference between two designs could be dilution rather than a real
#   change in the association; the estimates themselves are NOT corrected
#   (lambda is too unstable when every trajectory has two points).
#
#   Depends only on the LME, so it runs in seconds and does not re-fit the
#   joint model.
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
suppressPackageStartupMessages(library(nlme))

L <- read.csv(file.path(BASE, "poma_long.csv"))
U <- read.csv(file.path(BASE, "poma_analysis_long.csv"))
PRIM <- "cMFTC_ThCtAB_aMe"
COVARS <- c("V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")

sv0 <- U[, c("unit", "newstrata", "case", "t_index_months", COVARS)]
sv0 <- sv0[complete.cases(sv0[, COVARS]), ]
lg0 <- L[!is.na(L[[PRIM]]), c("unit", "months", PRIM)]
names(lg0)[3] <- "y"
lg0 <- lg0[lg0$unit %in% sv0$unit, ]
sv0 <- sv0[sv0$unit %in% lg0$unit, ]

mm <- merge(lg0, sv0[, c("unit", "t_index_months")], by = "unit")
mm$rel <- mm$months - mm$t_index_months
mm <- mm[mm$rel < 0, ]

## the pair sets actually used by 55_window_lag_sensitivity.R
pre_n <- table(mm$unit)
dens  <- data.frame(unit = names(pre_n), n = as.integer(pre_n), stringsAsFactors = FALSE)
dens  <- merge(dens, sv0[, c("unit", "newstrata", "case")], by = "unit")
wd    <- reshape(dens[, c("newstrata", "case", "n")], idvar = "newstrata",
                 timevar = "case", direction = "wide")
wd    <- wd[complete.cases(wd), ]
EQ    <- wd$newstrata[wd$n.1 == wd$n.0]

uok <- function(sub, need) { tb <- table(sub$unit); names(tb)[tb >= need] }

DESIGN <- list(
  primary   = list(units = sv0$unit,                       rows = quote(rep(TRUE, nrow(mm))), need = 2),
  cw24      = list(units = uok(mm[mm$rel >= -24, ], 2),     rows = quote(mm$rel >= -24),        need = 2),
  cw3       = list(units = uok(mm[mm$rel >= -36, ], 3),     rows = quote(mm$rel >= -36),        need = 3),
  cw3_ref   = list(units = uok(mm[mm$rel >= -36, ], 3),     rows = quote(rep(TRUE, nrow(mm))),  need = 3),
  lag24     = list(units = uok(mm[mm$rel <= -24, ], 2),     rows = quote(mm$rel <= -24),        need = 2),
  lag24_ref = list(units = uok(mm[mm$rel <= -24, ], 2),     rows = quote(rep(TRUE, nrow(mm))),  need = 2),
  lag36     = list(units = uok(mm[mm$rel <= -36, ], 2),     rows = quote(mm$rel <= -36),        need = 2),
  dens      = list(units = dens$unit[dens$newstrata %in% EQ], rows = quote(rep(TRUE, nrow(mm))), need = 2)
)
need <- sapply(DESIGN, `[[`, "need")
tab <- read.csv(file.path(BASE, "sens_window_lag.csv"), stringsAsFactors = FALSE)

out <- list()
for (nm in names(DESIGN)) {
  D <- DESIGN[[nm]]
  keep <- eval(D$rows) & (mm$unit %in% D$units)
  sub  <- mm[keep, c("unit", "months", "y")]
  sub$year <- sub$months / 12
  tb <- table(sub$unit)
  sub <- sub[sub$unit %in% names(tb)[tb >= need[nm]], ]

  ## restrict to complete pairs of the corresponding sensitivity design
  svS <- sv0[sv0$unit %in% sub$unit, ]
  pn  <- table(svS$newstrata)
  ## a design specific pair set is only available through its CSV row count;
  ## for the referent designs the units are inherited, so use the parent's
  ## pair set when the row is marked as a referent
  parent <- if (nm == "cw3_ref") "cw3" else if (nm == "lag24_ref") "lag24" else nm
  np <- tab$n_pairs[match(nm, tab$design)]

  m  <- lme(y ~ year, random = ~ year | unit, data = sub,
            control = lmeControl(opt = "optim"))
  D22 <- as.numeric(VarCorr(m)[2, 1])          # Var of the random slope
  sl  <- tapply(seq_len(nrow(sub)), sub$unit, function(ii)
    coef(lm(y ~ year, data = sub[ii, ]))[2])
  sdo <- sd(-as.numeric(sl))
  ## sampling variance of an OLS slope averaged over the observed visit counts
  nv  <- as.integer(table(sub$unit))
  Tt  <- nv * (nv^2 - 1) / 12                   # sum of squared deviances of t
  s2  <- as.numeric(VarCorr(m)[3, 1])
  samp <- mean(s2 / Tt)
  out[[nm]] <- data.frame(
    design = nm, n_pairs = np, n_knees = length(unique(sub$unit)),
    mean_visits = mean(nv), D22 = D22, sd_ols = sdo,
    lambda_lme = D22 / (D22 + samp),
    lambda_obs = D22 / sdo^2,
    stringsAsFactors = FALSE)
  ## dilution-corrected odds ratio: log(OR_true) = log(OR_obs) / lambda
  orr <- as.numeric(tab$clogit_ols_OR01[match(nm, tab$design)])
  lam <- out[[nm]]$lambda_obs
  out[[nm]]$or01_ols <- orr
  out[[nm]]$or01_corrected <- if (is.finite(lam) && lam > 0.05) exp(log(orr) / lam) else NA_real_
  cat(sprintf("%-10s pairs %3d knees %3d visits %.2f | SD(true) %.3f SD(ols) %.3f | lambda(lme) %.3f lambda(obs) %.3f\n",
              nm, np, length(unique(sub$unit)), mean(nv),
              sqrt(D22), sdo, out[[nm]]$lambda_lme, out[[nm]]$lambda_obs))
}

res <- do.call(rbind, out)
write.csv(res, file.path(BASE, "sens_wl_reliability.csv"), row.names = FALSE)
cat("\n-> sens_wl_reliability.csv\n")
