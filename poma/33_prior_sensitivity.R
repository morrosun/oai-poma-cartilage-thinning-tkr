# ============================================================================
# 33 : prior sensitivity (reviewer points 2d and 7c)
#
#   Two priors carry informative weight in this paper and neither had been
#   checked:
#
#   (A) the N(0, sd^2) prior on the log-OR of the BAYESIAN WITHIN-PAIR
#       CONDITIONAL LOGISTIC used for the per-class odds ratios of Table 4,
#       reported with sd = 2.  Class 4 is almost completely separated
#       (24 of 25 knees replaced), so this prior is doing real work.
#   (B) the N(0, prior_sd^2) prior on the joint-model association parameters
#       aV and aS (prior_sd = 20), together with the inverse-Wishart prior on
#       the random-effects covariance D (scale diag(1.2, 0.01), 4 df).
#
#   Everything else in the pipeline is either a conditional-logistic exact
#   likelihood or an empirical-Bayes plug-in and has no prior at all.
#
#   Output : prior_sensitivity_bclogit.csv , prior_sensitivity_joint.csv
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
source("D:/BaiduSyncdisk/OAI/Scripts/poma/jm_core.R")
suppressPackageStartupMessages(library(survival))

## =========================================================== (A) ========
## identical helpers to 11_trajectory_lcmm.R (copied, not sourced, so that
## re-running this script cannot refit the LCMM)
log1pexp  <- function(x) ifelse(x > 30, x + log1p(exp(-x)), log1p(exp(pmin(x, 30))))
logplogis <- function(x) -log1pexp(-x)
bclogit <- function(Xp, niter = 40000, burn = 8000, thin = 8, prior = 2, seed = 1) {
  set.seed(seed); p <- ncol(Xp); a <- rep(0, p); s <- rep(0.4, p)
  keep <- seq(burn + 1, niter, by = thin); out <- matrix(0, length(keep), p)
  k <- 0; acc <- 0; win <- 0
  lp <- function(a) sum(logplogis(drop(Xp %*% a))) - 0.5 * sum(a^2) / prior^2
  cur <- lp(a)
  for (it in 1:niter) {
    prop <- a + rnorm(p) * s; lpn <- lp(prop)
    if (log(runif(1)) < lpn - cur) { a <- prop; cur <- lpn; acc <- acc + 1 }
    win <- win + 1
    if (win == 500) { ar <- acc / 500; s <- s * ifelse(ar < 0.25, 0.8, ifelse(ar > 0.45, 1.25, 1))
                      acc <- 0; win <- 0 }
    if (k < length(keep) && it == keep[k + 1]) { k <- k + 1; out[k, ] <- a }
  }
  out
}
pair_contrast <- function(dat, cols) {
  dat <- dat[order(dat$strata, -dat$case), ]
  ii <- split(seq_len(nrow(dat)), dat$strata)
  ii <- ii[vapply(ii, length, 1L) == 2]
  as.matrix(dat[unlist(lapply(ii, `[`, 1)), cols, drop = FALSE]) -
    as.matrix(dat[unlist(lapply(ii, `[`, 2)), cols, drop = FALSE])
}

rd <- readRDS(file.path(BASE, "lcmm_fit.rds"))
an <- rd$pair_data
dmy <- grep("^cv", names(an), value = TRUE)
cat(sprintf("pair data : %d knees, %d class contrasts (%s)\n",
            nrow(an), length(dmy), paste(dmy, collapse = ", ")))
Xp <- pair_contrast(an, dmy)
cat(sprintf("contrast matrix : %d pairs x %d columns\n\n", nrow(Xp), ncol(Xp)))

## reference (slowest class) -> class labels follow the ordinal used in Table 4
ordcls <- levels(an$clfac)

outA <- data.frame()
PS <- c(1, 2, 4, 10)          # prior SDs on the log-OR scale
for (ps in PS) {
  post <- bclogit(Xp, prior = ps, seed = 11)
  for (j in seq_along(dmy)) {
    outA <- rbind(outA, data.frame(
      prior = sprintf("N(0, %g^2)", ps), term = sprintf("class %s vs %s",
                                                        ordcls[j + 1], ordcls[1]),
      OR = exp(median(post[, j])),
      lo = exp(quantile(post[, j], .025)),
      hi = exp(quantile(post[, j], .975)),
      P_OR_gt1 = mean(post[, j] > 0)))
  }
  cat(sprintf("---- bclogit prior N(0,%g^2) ----\n", ps))
  print(outA[outA$prior == sprintf("N(0, %g^2)", ps), ], row.names = FALSE, digits = 4)
  cat("\n")
}
write.csv(outA, file.path(BASE, "prior_sensitivity_bclogit.csv"), row.names = FALSE)

## =========================================================== (B) ========
L <- read.csv(file.path(BASE, "poma_long.csv"))
U <- read.csv(file.path(BASE, "poma_analysis_long.csv"))
COVARS <- c("V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")
PRIM <- "cMFTC_ThCtAB_aMe"
sub <- L[!is.na(L[[PRIM]]), ]
lg  <- data.frame(unit = sub$unit, year = sub$years, y = sub[[PRIM]],
                  newstrata = sub$newstrata, case = sub$case)
sv  <- U[, c("unit", "newstrata", "case", "t_index_months", COVARS)]
sv  <- sv[sv$unit %in% lg$unit, ]
sv  <- sv[complete.cases(sv[, COVARS]), ]
lg  <- lg[lg$unit %in% sv$unit, ]
sv$t_idx <- sv$t_index_months / 12
PP  <- prune_pairs(lg, sv, 2)
dt  <- build(PP$lg, PP$sv)
cat(sprintf("\nSet A : %d pairs / %d knees\n", dt$npair, dt$n))

q3 <- function(x) sprintf("%.3f (%.3f, %.3f)", median(x), quantile(x, .025), quantile(x, .975))
rowjm <- function(lab, spec, f) {
  d <- f$draws
  data.frame(specification = spec, term = lab,
             estimate = q3(switch(lab,
               "OR per 1 SD faster thinning" = exp(-d$aS * d$SDslo),
               "OR per 0.1 mm/yr faster"     = exp(-0.1 * d$aS),
               "OR per 1 SD thinner (level)" = exp(-d$aV * d$SDval))),
             check.names = FALSE)
}

outB <- data.frame()
## ---- (B1) SD of the normal prior on (aV, aS)
for (ps in c(5, 10, 20, 50)) {
  f <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 2026,
              prior_sd = ps, verbose = FALSE)
  sp <- sprintf("log-OR prior N(0, %g^2)", ps)
  for (lab in c("OR per 1 SD faster thinning", "OR per 0.1 mm/yr faster",
                "OR per 1 SD thinner (level)")) outB <- rbind(outB, rowjm(lab, sp, f))
  cat(sprintf("%-28s OR/SD rate %s | OR/0.1 %s | OR/SD level %s\n", sp,
              q3(exp(-f$draws$aS * f$draws$SDslo)), q3(exp(-0.1 * f$draws$aS)),
              q3(exp(-f$draws$aV * f$draws$SDval))))
}
## ---- (B2) inverse-Wishart prior on the random-effects covariance D
for (sc in c(0.5, 1, 2)) {
  f <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 2026,
              iw_scale = sc, verbose = FALSE)
  sp <- sprintf("D ~ IW(scale x %g, 4 df)", sc)
  for (lab in c("OR per 1 SD faster thinning", "OR per 0.1 mm/yr faster",
                "OR per 1 SD thinner (level)")) outB <- rbind(outB, rowjm(lab, sp, f))
  cat(sprintf("%-28s OR/SD rate %s | OR/0.1 %s | OR/SD level %s\n", sp,
              q3(exp(-f$draws$aS * f$draws$SDslo)), q3(exp(-0.1 * f$draws$aS)),
              q3(exp(-f$draws$aV * f$draws$SDval))))
}
## ---- (B3) joint change: weaker prior on the rate association plus a
##           different D prior, i.e. the least favourable corner
f <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 2026,
            prior_sd = 5, iw_scale = 2, verbose = FALSE)
sp <- "combination: N(0, 5^2) + IW(scale x 2)"
for (lab in c("OR per 1 SD faster thinning", "OR per 0.1 mm/yr faster",
              "OR per 1 SD thinner (level)")) outB <- rbind(outB, rowjm(lab, sp, f))
cat(sprintf("%-28s OR/SD rate %s | OR/0.1 %s | OR/SD level %s\n", sp,
            q3(exp(-f$draws$aS * f$draws$SDslo)), q3(exp(-0.1 * f$draws$aS)),
            q3(exp(-f$draws$aV * f$draws$SDval))))

write.csv(outB, file.path(BASE, "prior_sensitivity_joint.csv"), row.names = FALSE)
cat("\n-> prior_sensitivity_bclogit.csv / prior_sensitivity_joint.csv\n")
