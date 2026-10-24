# ============================================================================
# 55 : WINDOW-ALIGNMENT, VISIT-DENSITY AND LAG SENSITIVITY
#      (review batch 3, item 9)
#
#   Three criticisms that the paired design alone does not answer.
#
#   (A) OBSERVATION WINDOW.  Set A uses every pre-index visit, so the width of
#       the exposure window (12-48 months) and the number of visits per knee
#       (2-5) both vary, and cases and controls differ in visit density
#       (3.50 vs 3.03 pre-index visits).  If the rate association were an
#       artefact of differential precision it should move when both are held
#       equal.
#         cw3  : the three most recent pre-index visits of EVERY knee, i.e.
#                three points over a 24-month span for all of them
#         dens : the 139 pairs whose case and control contribute the same
#                number of pre-index visits
#         cw24 : the two most recent visits (a 24-month lookback) -- perfectly
#                balanced, but every trajectory then has two points and the
#                joint model is not identified (no within-knee residual
#                degrees of freedom), so this design is conditional-logistic
#                only and is reported as such.
#
#   (B) LAG.  The index visit is the visit at which the case knee was
#       ascertained, so the knee replacement happened between the LAST analysed
#       visit (index-12) and the index visit: the final pre-index year still
#       overlaps the run-up to surgery, and terminal accelerated thinning could
#       drive the association by reverse causation.
#         lag24 : keep only visits at index-24 or earlier -> the last retained
#                 MRI sits >= 24 months before the index visit
#         lag36 : keep only visits at index-36 or earlier -> the last retained
#                 MRI sits >= 36 months before the index visit
#
#   (C) SAME-PAIR REFERENTS.  A restriction also changes WHICH pairs are
#       analysable, so comparing a restricted design with the primary estimate
#       confounds the window effect with a sample effect.  For the designs that
#       carry the argument (cw3, lag24) the SAME pairs are refitted with ALL
#       their pre-index visits (rows "*_ref"):
#           window effect = restricted vs same-pair referent
#           sample effect = referent vs primary
#
#   Estimators, so that every row is on one scale:
#     * conditional logistic (survival::clogit, exact) on the LME slope
#       (fixef + EBLUP), the same exposure used in Table 2 -> OR per 0.1 mm/yr
#       obtained as exp(b * 0.1 / sd), not from a standardised coefficient;
#     * the within-pair conditional joint model (jm_core.R), crude and adjusted
#       for the within-pair KL and WOMAC-pain differences.
#
#   Usage:  Rscript 55_window_lag_sensitivity.R [set1,set2,...]
#           no argument = every design.  Each design writes its own CSV, so a
#           long run may be split across invocations without losing work.
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
source("D:/BaiduSyncdisk/OAI/Scripts/poma/jm_core.R")
suppressPackageStartupMessages({ library(survival); library(nlme) })

L  <- read.csv(file.path(BASE, "poma_long.csv"))
U  <- read.csv(file.path(BASE, "poma_analysis_long.csv"))
PRIM <- "cMFTC_ThCtAB_aMe"
COVARS <- c("V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")

## ---------------------------------------------------------------- base data
sv0 <- U[, c("unit", "newstrata", "case", "t_index_months", COVARS,
             "cMFTC_ThCtAB_aMe__base")]
sv0 <- sv0[complete.cases(sv0[, COVARS]), ]
lg0 <- L[!is.na(L[[PRIM]]), c("unit", "months", PRIM)]
names(lg0)[3] <- "y"
lg0 <- lg0[lg0$unit %in% sv0$unit, ]
sv0 <- sv0[sv0$unit %in% lg0$unit, ]
sv0$t_idx <- sv0$t_index_months / 12

mm <- merge(lg0, sv0[, c("unit", "t_index_months")], by = "unit")
mm$rel <- mm$months - mm$t_index_months
mm <- mm[mm$rel < 0, ]

pre_n <- table(mm$unit)
dens  <- data.frame(unit = names(pre_n), n = as.integer(pre_n), stringsAsFactors = FALSE)
dens  <- merge(dens, sv0[, c("unit", "newstrata", "case")], by = "unit")
wd    <- reshape(dens[, c("newstrata", "case", "n")], idvar = "newstrata",
                 timevar = "case", direction = "wide")
wd    <- wd[complete.cases(wd), ]
EQ_PAIRS <- wd$newstrata[wd$n.1 == wd$n.0]

units_ok <- function(sub, need = 2) { tb <- table(sub$unit); names(tb)[tb >= need] }

## ------------------------------------------------------------------ fitting
fit_one <- function(name, family, keep_units, rows_ok, label, jm = TRUE, minv = 2) {
  lgS <- mm[mm$unit %in% keep_units & rows_ok(mm$rel), c("unit", "months", "y")]
  lgS$year <- lgS$months / 12
  ok  <- units_ok(lgS, minv)
  lgS <- lgS[lgS$unit %in% ok, ]
  svS <- sv0[sv0$unit %in% ok, ]
  svS$pair <- NULL

  A  <- prune_pairs(lgS[, c("unit", "year", "y")], svS, min_visits = minv)
  dt <- build(A$lg, A$sv)
  np <- length(dt$cas)
  if (np < 50) { cat(sprintf("\n[%s] only %d pairs -> SKIP\n", name, np)); return(NULL) }
  uA <- svS[match(A$sv$unit, svS$unit), ]

  ## ---- descriptive ----
  nv <- table(A$lg$unit)
  nv_case <- mean(nv[as.character(uA$unit[uA$case == 1])])
  nv_ctl  <- mean(nv[as.character(uA$unit[uA$case == 0])])
  span <- tapply(A$lg$year, A$lg$unit, function(x) diff(range(x)) * 12)
  gap  <- uA$t_index_months - tapply(A$lg$year, A$lg$unit, max)[as.character(uA$unit)] * 12

  ## ---- exposure: LME slope (fixef + EBLUP), the estimator behind Table 2 ----
  m  <- lme(y ~ year, random = ~ year | unit, data = A$lg,
            control = lmeControl(opt = "optim"))
  re <- ranef(m)
  eb <- fixef(m)[2] + re[, 2]
  names(eb) <- rownames(re)          # re[,2] drops the names of the group levels
  loss_eb <- -as.numeric(eb[as.character(uA$unit)])
  ## naive per-knee OLS slope, for the dilution check only
  loss_ol <- -as.numeric(tapply(seq_len(nrow(A$lg)), A$lg$unit, function(ii)
    coef(lm(y ~ year, data = A$lg[ii, ]))[2])[as.character(uA$unit)])

  cl <- function(loss, adj) {
    sdz <- sd(loss)
    d <- data.frame(case = uA$case, strata = uA$newstrata,
                    zz = (loss - mean(loss)) / sdz,
                    kl = uA$V00XRKL, wom5 = uA$WOMAC_pain / 5,
                    age5 = uA$V00AGE / 5, bmi5 = uA$P01BMI / 5, female = uA$female)
    f <- if (adj) clogit(case ~ zz + age5 + female + bmi5 + kl + wom5 + strata(strata),
                         data = d, method = "exact")
         else     clogit(case ~ zz + strata(strata), data = d, method = "exact")
    s <- summary(f)$coef
    k <- 0.1 / sdz
    b <- unname(s["zz", 1]); se <- unname(s["zz", 3]); pv <- unname(s["zz", 5])
    c(OR = exp(b * k), lo = exp(b * k - 1.96 * se * k), hi = exp(b * k + 1.96 * se * k),
      p = pv, sd = unname(sdz), b_per_sd = exp(b))
  }
  c0 <- cl(loss_eb, FALSE); c1 <- cl(loss_eb, TRUE)
  o0 <- cl(loss_ol, FALSE)

  ## ---- joint model ----
  q <- function(x, p) as.numeric(quantile(x, p))
  if (jm) {
    f  <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 2026, verbose = FALSE)
    report(f, sprintf("JM  %s : cMFTC", label))
    Zd <- cbind(kl_diff    = uA$V00XRKL[dt$cas]    - uA$V00XRKL[dt$ctl],
                womac_diff = uA$WOMAC_pain[dt$cas] - uA$WOMAC_pain[dt$ctl])
    fa <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 2028,
                 Zd = Zd, verbose = FALSE)
    report(fa, sprintf("JM  %s : cMFTC + within-pair KL/pain", label))
    jm_or <- q(exp(-0.1 * f$draws$aS), .5);   jm_lo <- q(exp(-0.1 * f$draws$aS), .025)
    jm_hi <- q(exp(-0.1 * f$draws$aS), .975); jm_p  <- mean(f$draws$aS < 0)
    ja_or <- q(exp(-0.1 * fa$draws$aS), .5);  ja_lo <- q(exp(-0.1 * fa$draws$aS), .025)
    ja_hi <- q(exp(-0.1 * fa$draws$aS), .975); ja_p <- mean(fa$draws$aS < 0)
  } else {
    jm_or <- jm_lo <- jm_hi <- jm_p <- ja_or <- ja_lo <- ja_hi <- ja_p <- NA_real_
  }

  g <- function(x, k) as.numeric(x[[k]])
  row <- data.frame(
    design = name, family = family, label = label,
    n_pairs = np, n_knees = dt$n,
    visits_case = nv_case, visits_control = nv_ctl,
    span_median = median(span), span_min = min(span), span_max = max(span),
    gap_median = median(gap), gap_min = min(gap), gap_max = max(gap),
    sd_loss_eblup = g(c0, "sd"), sd_loss_ols = g(o0, "sd"),
    JM_OR01 = jm_or, JM_OR01_lo = jm_lo, JM_OR01_hi = jm_hi, JM_P = jm_p,
    JM_OR01_adj = ja_or, JM_OR01_adj_lo = ja_lo, JM_OR01_adj_hi = ja_hi,
    JM_P_adj = ja_p,
    clogit_OR01 = g(c0, "OR"), clogit_OR01_lo = g(c0, "lo"),
    clogit_OR01_hi = g(c0, "hi"), clogit_p = g(c0, "p"),
    clogit_ORSD = g(c0, "b_per_sd"),
    clogit_OR01_adj = g(c1, "OR"), clogit_OR01_adj_lo = g(c1, "lo"),
    clogit_OR01_adj_hi = g(c1, "hi"), clogit_adj_p = g(c1, "p"),
    clogit_ols_OR01 = g(o0, "OR"), clogit_ols_OR01_lo = g(o0, "lo"),
    clogit_ols_OR01_hi = g(o0, "hi"),
    row.names = NULL, stringsAsFactors = FALSE)

  out <- file.path(BASE, sprintf("sens_wl_%s.csv", name))
  write.csv(row, out, row.names = FALSE)

  cat(sprintf("\n---- %s ----\n", name))
  cat(sprintf("pairs %d | knees %d | visits/knee case %.2f control %.2f | span median %.0f mo (%0.f-%0.f) | last visit to index median %.0f mo (min %.0f)\n",
              np, dt$n, nv_case, nv_ctl, median(span), min(span), max(span),
              median(gap), min(gap)))
  cat(sprintf("EBLUP clogit OR/0.1 %.3f (%.3f-%.3f) P=%.3g | adj %.3f (%.3f-%.3f) P=%.3g\n",
              c0["OR"], c0["lo"], c0["hi"], c0["p"],
              c1["OR"], c1["lo"], c1["hi"], c1["p"]))
  cat(sprintf("OLS   clogit OR/0.1 %.3f (%.3f-%.3f)  [regression-dilution reference]\n",
              o0["OR"], o0["lo"], o0["hi"]))
  if (jm) cat(sprintf("JM      OR/0.1 %.3f (%.3f-%.3f) P(aS<0)=%.3f | adj %.3f (%.3f-%.3f) P(aS<0)=%.3f\n",
                      jm_or, jm_lo, jm_hi, jm_p, ja_or, ja_lo, ja_hi, ja_p))
  cat(sprintf("-> %s\n", basename(out)))
  row
}

## ------------------------------------------------------------------ designs
all_pre <- function(rel) rep(TRUE, length(rel))
win24   <- function(rel) rel >= -24
win36   <- function(rel) rel >= -36
lag24f  <- function(rel) rel <= -24
lag36f  <- function(rel) rel <= -36

DES <- list(
  primary   = list(fam = "primary",      u = NULL,        rows = all_pre, jm = TRUE,
                   lab = "primary: all pre-index visits"),
  cw24      = list(fam = "common window", u = "win24e2",  rows = win24,   jm = FALSE,
                   lab = "two most recent visits, 24-month lookback (clogit only)"),
  cw3       = list(fam = "common window", u = "win36e3",  rows = win36,   jm = TRUE,
                   lab = "three most recent visits, 36-month lookback"),
  cw3_ref   = list(fam = "referent",      u = "win36e3",  rows = all_pre, jm = TRUE,
                   lab = "cw3 pairs, all pre-index visits"),
  cw3_far   = list(fam = "common window", u = "win36e3",  rows = lag24f,  jm = TRUE,
                   minv = 3,
                   lab = "cw3 pairs, same 3-visit window shifted back 12 months"),
  lag24     = list(fam = "lag",           u = "lag24",    rows = lag24f,  jm = TRUE,
                   lab = "visits at index-24 or earlier (last MRI >= 24 months before index)"),
  lag24_ref = list(fam = "referent",      u = "lag24",    rows = all_pre, jm = TRUE,
                   lab = "lag24 pairs, all pre-index visits"),
  lag36     = list(fam = "lag",           u = "lag36",    rows = lag36f,  jm = TRUE,
                   lab = "visits at index-36 or earlier (last MRI >= 36 months before index)"),
  dens      = list(fam = "density",       u = "eqdensity", rows = all_pre, jm = TRUE,
                   lab = "pairs with equal pre-index visit counts")
)

units_for <- function(key) {
  if (is.null(key))        return(sv0$unit)
  if (key == "win24e2")    return(units_ok(mm[mm$rel >= -24, ], 2))
  if (key == "win36e3")    return(units_ok(mm[mm$rel >= -36, ], 3))
  if (key == "lag24")      return(units_ok(mm[mm$rel <= -24, ], 2))
  if (key == "lag36")      return(units_ok(mm[mm$rel <= -36, ], 2))
  if (key == "eqdensity")  return(dens$unit[dens$newstrata %in% EQ_PAIRS])
  stop("unknown unit key: ", key)
}

ARG  <- commandArgs(trailingOnly = TRUE)
want <- if (length(ARG)) strsplit(ARG[1], ",")[[1]] else names(DES)
want <- intersect(want, names(DES))

sink(file.path(BASE, "sens_window_lag_report.txt"), split = TRUE, append = length(ARG) > 0)
cat("================================================================================\n")
cat("POMA : WINDOW-ALIGNMENT, VISIT-DENSITY AND LAG SENSITIVITY -- cMFTC\n")
cat("================================================================================\n")
cat(sprintf("pre-index rows %d | knee units %d | complete pairs %d\n",
            nrow(mm), length(unique(mm$unit)), nrow(wd)))
cat(sprintf("pre-index visits per knee: case %.2f vs control %.2f (paired Wilcoxon P = %.3g)\n",
            mean(wd$n.1), mean(wd$n.0), wilcox.test(wd$n.1, wd$n.0, paired = TRUE)$p.value))
cat(sprintf("pairs with equal case/control visit counts: %d of %d\n", length(EQ_PAIRS), nrow(wd)))

res <- list()
for (nm in want) {
  D <- DES[[nm]]
  cat(sprintf("\n\n################  %s  ################\n", nm))
  r <- try(fit_one(nm, D$fam, units_for(D$u), D$rows, D$lab, D$jm,
                   if (is.null(D$minv)) 2 else D$minv), silent = TRUE)
  if (inherits(r, "try-error")) cat("  ERROR:", conditionMessage(attr(r, "condition")), "\n") else res[[nm]] <- r
}

if (length(res)) {
  tab <- do.call(rbind, res)
  prev <- file.path(BASE, "sens_window_lag.csv")
  if (file.exists(prev)) {
    old <- read.csv(prev, stringsAsFactors = FALSE)
    old <- old[!old$design %in% tab$design, , drop = FALSE]
    if (nrow(old)) {
      for (m in setdiff(names(tab), names(old))) old[[m]] <- NA
      tab <- rbind(tab, old[, names(tab), drop = FALSE])
    }
  }
  tab <- tab[match(intersect(names(DES), tab$design), tab$design), , drop = FALSE]
  rownames(tab) <- NULL
  write.csv(tab, file.path(BASE, "sens_window_lag.csv"), row.names = FALSE)

  cat("\n\n================ SUMMARY ================\n")
  print(tab[, c("design", "n_pairs", "visits_case", "visits_control", "span_median",
                "gap_median", "clogit_OR01", "clogit_OR01_lo", "clogit_OR01_hi",
                "JM_OR01", "JM_OR01_lo", "JM_OR01_hi")], row.names = FALSE)
}
cat("\n-> sens_window_lag.csv / sens_window_lag_report.txt\n")
sink()
