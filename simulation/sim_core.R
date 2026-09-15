# ============================================================================
# sim_core.R : shared engine for the POMA-calibrated simulation studies
#
#   Design fidelity
#     The visit schedule of every simulated pair is resampled (without
#     replacement where possible) from the 191 observed POMA pairs, so the
#     number and timing of pre-index visits, the index visit (t*) and the
#     case-versus-control visit-density asymmetry are exactly those of the
#     real study (cases 3.50 visits, controls 3.03).
#
#   Generating model (identical in form to the one the manuscript fits)
#     y_ij      = b0_i + b1_i * t_ij + e_ij ,          e ~ N(0, sigma_e^2)
#     logit P(i is the case in pair s) =
#                 aV * [ (b0 + b1 t*)_case - (b0 + b1 t*)_ctl ]
#               + aS * [ (-b1)_case      - (-b1)_ctl      ]
#     i.e. the hazard depends on the current thickness LEVEL (aV) and on the
#     thinning RATE (aS, positive = faster thinning -> higher risk).  Both
#     members of a pair share t*, so within-pair constant terms cancel and
#     this is exactly the within-pair conditional likelihood of jm_core.R.
#
#   Calibrated constants come from 10_calibrate_poma.py / pilot_report.txt:
#     tau = 0.1779 mm/yr (SD of the true slope), sigma_e = 0.2317 mm,
#     mean slope -0.1455, mean baseline 3.38 mm (SD 1.55).
# ============================================================================

suppressPackageStartupMessages({
  library(nlme); library(survival); library(MASS)
})

## ------------------------------------------------------------------ design --
read_design <- function(path) {
  S <- read.csv(path, stringsAsFactors = FALSE)
  S$mlist <- lapply(strsplit(S$months, ";"), function(x) as.numeric(x) / 12)
  split(S, S$pair)
}

## --------------------------------------------------------------- generator --
## npair     : number of pairs to generate (<= length(design) -> no resampling)
## visits_ctl: optional integer; force controls to have at most this many visits
##             (used to create a differential visit-density scenario)
## zsd, aZ   : optional PAIR-LEVEL confounder.  A pair-level latent Z_p ~ N(0,1)
##             shifts both the true slope (SD zsd) and the hazard (log-HR aZ).
##             It is removed by within-pair differencing but NOT by an
##             unstratified analysis -- this is the quantity that makes matched
##             designs require a within-stratum analysis.
gen_data <- function(design, npair = NULL, tau, mu1, b0mu, b0sd, sigma_e,
                     aS, aV = 0, seed = 1,
                     visits_case = NULL, visits_ctl = NULL,
                     zsd = 0, aZ = 0) {
  set.seed(seed)
  D <- design
  if (is.null(npair) || npair >= length(D)) {
    idx <- seq_along(D); npair <- length(D)
  } else {
    idx <- sample(length(D), npair)
  }
  tau_w <- sqrt(max(tau^2 - zsd^2, 1e-8))
  rows <- vector("list", npair); k <- 0L
  for (s in seq_len(npair)) {
    d <- D[[idx[s]]]
    tstar <- unique(d$t_index)[1] / 12
    sc <- d$mlist[d$role == 1][[1]]           # case visit times (years)
    sj <- d$mlist[d$role == 0][[1]]           # control visit times
    if (!is.null(visits_case)) sc <- head(sc, visits_case)
    if (!is.null(visits_ctl))  sj <- head(sj, visits_ctl)
    if (length(sc) < 2 || length(sj) < 2) { next }

    zp <- rnorm(1)
    b0 <- rnorm(2, b0mu, b0sd)
    b1 <- mu1 + zsd * zp + rnorm(2, 0, tau_w)
    eta <- aV * (b0 + b1 * tstar) + aS * (-b1) + aZ * zp
    is_case1 <- runif(1) < plogis(eta[1] - eta[2])
    role <- if (is_case1) c(1L, 0L) else c(0L, 1L)

    for (i in 1:2) {
      tt <- if (role[i] == 1L) sc else sj
      k <- k + 1L
      rows[[k]] <- data.frame(
        pair = s, unit = k, case = role[i], t = tt, tstar = tstar,
        b0 = b0[i], b1 = b1[i], zp = zp,
        y = b0[i] + b1[i] * tt + rnorm(length(tt), 0, sigma_e)
      )
    }
  }
  do.call(rbind, rows[seq_len(k)])
}

## -------------------------------------------------------------- estimators --
## returns one row per knee: OLS slope, EBLUP slope (random intercept+slope),
## and the number of visits
knee_slopes <- function(dat) {
  sp <- split(dat, dat$unit)
  ols <- vapply(sp, function(g) {
    if (nrow(g) < 2 || diff(range(g$t)) == 0) return(NA_real_)
    unname(coef(lm(y ~ t, g))[2])
  }, numeric(1))
  m <- try(lme(y ~ t, random = ~ 1 + t | unit, data = dat,
               control = lmeControl(opt = "optim")), silent = TRUE)
  if (inherits(m, "try-error")) {
    eb <- rep(NA_real_, length(sp)); names(eb) <- names(sp)
  } else {
    re <- ranef(m)
    eb <- as.numeric(fixef(m)[["t"]]) + re[, "t"]     # keep the rownames
    names(eb) <- rownames(re)
  }
  k1 <- dat[!duplicated(dat$unit), c("unit", "pair", "case", "tstar")]
  k1$n_visits <- as.integer(table(dat$unit)[as.character(k1$unit)])
  k1$slope_ols <- as.numeric(ols[as.character(k1$unit)])
  k1$slope_eblup <- as.numeric(eb[as.character(k1$unit)])
  k1[order(k1$pair, -k1$case), ]
}

## within-pair conditional logistic regression (paired likelihood)
## sign = -1 puts the exposure into LOSS units (positive = thinner)
clogit_fit <- function(kn, var, sign = -1) {
  kn$xb <- sign * kn[[var]]
  kn$str <- factor(kn$pair)
  fit <- try(clogit(case ~ xb + strata(str), data = kn), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(beta = NA, se = NA))
  s <- summary(fit)
  c(beta = unname(s$coef[1, 1]), se = unname(s$coef[1, 3]))
}

## same but with a level term added (for the level-versus-rate argument)
clogit_fit2 <- function(kn, varrate, varlevel) {
  kn$zr <- -kn[[varrate]]; kn$zl <- kn[[varlevel]]
  kn$str <- factor(kn$pair)
  fit <- try(clogit(case ~ zr + zl + strata(str), data = kn), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(b1 = NA, s1 = NA, b2 = NA, s2 = NA))
  s <- summary(fit)
  c(b1 = unname(s$coef[1, 1]), s1 = unname(s$coef[1, 3]),
    b2 = unname(s$coef[2, 1]), s2 = unname(s$coef[2, 3]))
}

## oracle: conditional logistic on the TRUE slope
add_truth <- function(kn, dat) {
  tr <- dat[!duplicated(dat$unit), c("unit", "b1")]
  kn$slope_true <- tr$b1[match(kn$unit, tr$unit)]
  kn
}

## plain unstratified (population) partial likelihood, used to show what is
## lost when the matching is ignored
cox_fit <- function(kn, var, sign = -1) {
  kn$xb <- sign * kn[[var]]
  fit <- try(coxph(Surv(tstar, case) ~ xb, data = kn), silent = TRUE)
  if (inherits(fit, "try-error")) return(c(beta = NA, se = NA))
  s <- summary(fit)
  c(beta = unname(s$coef[1, 1]), se = unname(s$coef[1, 3]))
}

## ---------------------------------------------------------------- reporting --
## x : numeric vector, possibly containing NA for failed / non-estimable fits
summ <- function(x, truth = NA_real_) {
  n_all <- length(x); ok <- is.finite(x)
  fail <- sum(!ok)
  if (!any(ok)) {
    return(data.frame(mean = NA, sd = NA, lo = NA, hi = NA,
                      bias = NA, rmse = NA, cover = NA, n = 0L, fail = fail))
  }
  z <- x[ok]
  data.frame(mean = mean(z), sd = sd(z),
             lo = mean(z) - 1.96 * sd(z) / sqrt(length(z)),
             hi = mean(z) + 1.96 * sd(z) / sqrt(length(z)),
             bias = if (is.na(truth)) NA_real_ else mean(z) - truth,
             rmse = if (is.na(truth)) NA_real_ else sqrt(mean((z - truth)^2)),
             cover = NA_real_, n = length(z), fail = fail)
}

## coverage needs the per-replicate SE, so it gets its own summary
summ_cover <- function(beta, se, truth) {
  ok <- is.finite(beta) & is.finite(se)
  if (!any(ok)) return(c(cover = NA, n = 0, fail = length(beta)))
  lo <- beta[ok] - 1.96 * se[ok]; hi <- beta[ok] + 1.96 * se[ok]
  c(cover = mean(lo <= truth & truth <= hi), n = sum(ok),
    fail = sum(!ok))
}
