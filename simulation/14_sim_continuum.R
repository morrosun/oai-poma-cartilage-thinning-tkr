# ============================================================================
# 14 : methodological result 3 -- a latent-class solution can be the tail of a
#      continuous distribution rather than a set of subtypes
#
#   Bauer & Curran (2003) showed that a growth-mixture model extracts classes
#   from data that contain none.  Here the truth is KNOWN: the per-knee slopes
#   are drawn from a SINGLE continuous distribution, with the visit schedule and
#   the scale taken from POMA.  Two generating distributions are used:
#
#     D1  empirical resampling of the observed cMFTC EBLUP slopes
#     D2  skew-normal fitted to the same slopes
#
#   Both are unimodal and continuous by construction -- there is no discrete
#   structure to find.  The POMA pipeline is then applied unchanged
#   (`hlme`, class-specific intercept + slope, random intercept, best of
#   several random starts, initial values from `mkB` exactly as in
#   Scripts/poma/11_trajectory_lcmm.R), and three interpretation tests are run:
#
#     T1  rank-threshold equivalence : fastest class vs top-n on the rate axis
#     T2  incremental information    : does class add anything beyond the rate?
#     T3  gradient below the cut     : dose-response inside the remaining knees
#
#   Env: SIM_REPS (3) SIM_NSTART (4) SIM_KMAX (5)
# ============================================================================
source("D:/BaiduSyncdisk/OAI/Scripts/sim/sim_core.R")
suppressPackageStartupMessages(library(lcmm))

BASE   <- "D:/BaiduSyncdisk/OAI"
OUT    <- file.path(BASE, "Analysis/sim")
design <- read_design(file.path(OUT, "sim_visit_schedule.csv"))
CAL    <- readRDS(file.path(OUT, "sim_as.rds"))
aS <- CAL$aS; CONST <- CAL$const

REPS   <- as.integer(Sys.getenv("SIM_REPS", "3"))
NSTART <- as.integer(Sys.getenv("SIM_NSTART", "4"))
KMAX   <- as.integer(Sys.getenv("SIM_KMAX", "5"))

LOG <- list(); P <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG[[length(LOG)+1]] <<- s }

## ---- generating distributions ---------------------------------------------
poma   <- read.csv(file.path(BASE, "Analysis/poma_pilot/poma_analysis_long.csv"))
sl_obs <- na.omit(poma$cMFTC_ThCtAB_aMe__eblup)
P("=== 14  latent-class over-extraction  (%d replicates, %d starts, K <= %d) ===",
  REPS, NSTART, KMAX)
P("observed EBLUP slopes: n = %d, mean %+.4f, sd %.4f, skewness %+.3f",
  length(sl_obs), mean(sl_obs), sd(sl_obs),
  mean(((sl_obs - mean(sl_obs)) / sd(sl_obs))^3))

## ---- generate from a KNOWN continuous slope distribution -------------------
## Case status is assigned with the SAME conditional likelihood the manuscript
## uses, but on the new slopes, so the rate-outcome association is intact.
gen_from_slopes <- function(slopes, seed) {
  set.seed(seed)
  D <- design; rows <- vector("list", length(D)); k <- 0L
  for (s in seq_along(D)) {
    d <- D[[s]]; tstar <- unique(d$t_index)[1] / 12
    sc <- d$mlist[d$role == 1][[1]]; sj <- d$mlist[d$role == 0][[1]]
    if (length(sc) < 2 || length(sj) < 2) next
    b0 <- rnorm(2, CONST$b0mu, CONST$b0sd)
    b1 <- sample(slopes, 2, replace = TRUE)
    eta <- aS * (-b1)
    role <- if (runif(1) < plogis(eta[1] - eta[2])) c(1L, 0L) else c(0L, 1L)
    for (i in 1:2) {
      tt <- if (role[i] == 1L) sc else sj
      k <- k + 1L
      rows[[k]] <- data.frame(pair = s, unit = k, case = role[i], t = tt,
                              tstar = tstar, b1 = b1[i],
                              y = b0[i] + b1[i] * tt + rnorm(length(tt), 0, CONST$sigma_e))
    }
  }
  do.call(rbind, rows[seq_len(k)])
}

## ---- exact LCMM pipeline from Scripts/poma/11_trajectory_lcmm.R ------------
entropy <- function(pp) {
  p <- as.matrix(pp[, grep("^prob", names(pp))])
  1 - sum(-p * log(pmax(p, 1e-12))) / (nrow(p) * log(ncol(p)))
}
mkB <- function(K, rint, npar_fixed = 2) {
  ints <- seq(4.0, 2.4, length.out = K); sls <- seq(-0.03, -0.22, length.out = K)
  c(as.numeric(rbind(ints, sls)), rep(0, K - 1), 0.6, if (rint) 1.0)
}
fitK <- function(K, d, rint = TRUE, nstart = NSTART) {
  best <- NULL; bl <- -Inf
  for (s in seq_len(nstart)) {
    set.seed(700 * K + s)
    mm <- try(hlme(y ~ t, mixture = ~ t, random = if (rint) ~1 else ~ -1,
                   subject = "unit", ng = K, data = d, B = mkB(K, rint)),
              silent = TRUE)
    if (!inherits(mm, "try-error") && is.finite(mm$loglik) && mm$loglik > bl) {
      bl <- mm$loglik; best <- mm
    }
  }
  best
}

analyse <- function(slopes, seed) {
  dl <- gen_from_slopes(slopes, seed)
  d <- data.frame(unit = dl$unit, t = dl$t, y = dl$y)
  sel <- data.frame(); fits <- list()
  for (K in 1:KMAX) {
    m <- if (K == 1) { set.seed(1)
        try(hlme(y ~ t, random = ~1, subject = "unit", ng = 1, data = d), silent = TRUE)
      } else fitK(K, d)
    if (inherits(m, "try-error") || is.null(m)) {
      sel <- rbind(sel, data.frame(K = K, BIC = NA, entropy = NA, minpct = NA)); next
    }
    fits[[as.character(K)]] <- m
    cs <- table(m$pprob$class)
    sel <- rbind(sel, data.frame(K = K, BIC = m$BIC,
                                 entropy = if (K > 1) entropy(m$pprob) else NA,
                                 minpct = 100 * min(cs) / sum(cs)))
  }
  Kb <- sel$K[which.min(sel$BIC)]
  mB <- fits[[as.character(Kb)]]
  out <- list(sel = sel, K = Kb, n = length(unique(dl$unit)))
  if (!is.null(mB)) {
    kn <- knee_slopes(dl)
    kn$slope_eblup <- kn$slope_eblup
    kn$class <- mB$pprob$class[match(kn$unit, mB$pprob$unit)]
    cm <- tapply(kn$slope_eblup, kn$class, mean)
    fast <- as.integer(names(cm)[which.min(cm)])
    kn$is_fast <- as.integer(kn$class == fast)
    nfast <- sum(kn$is_fast)
    ord <- order(kn$slope_eblup)
    topn <- integer(nrow(kn)); topn[ord[seq_len(nfast)]] <- 1L
    out$T1_agreement <- 100 * mean(topn == kn$is_fast)
    kn$rate_z <- as.numeric(scale(-kn$slope_eblup)); kn$str <- factor(kn$pair)
    out$T2_rate  <- summary(clogit(case ~ rate_z + strata(str), data = kn))$coef[1, c(1, 3, 5)]
    out$T2_class <- summary(clogit(case ~ is_fast + strata(str), data = kn))$coef[1, c(1, 3, 5)]
    out$T2_joint <- summary(clogit(case ~ rate_z + is_fast + strata(str), data = kn))$coef[, c(1, 3, 5)]
    bl <- kn[kn$is_fast == 0, ]; bl$str <- factor(bl$pair)
    tb2 <- table(bl$str); bl <- bl[bl$str %in% names(tb2)[tb2 == 2], ]
    out$T3 <- if (nrow(bl) > 8)
      summary(clogit(case ~ rate_z + strata(str), data = bl))$coef[1, c(1, 3, 5)] else NULL
    out$class_sizes <- as.integer(table(kn$class)); out$nfast <- nfast
  }
  out
}

res <- list()
for (r in seq_len(REPS)) {
  emp <- (r %% 2 == 1)
  set.seed(2000 + r)
  scen <- if (emp) sample(sl_obs, 600, TRUE)
          else rnorm(600, mean(sl_obs), sd(sl_obs))
  P("\n--- replicate %d  (generating distribution: %s) ---", r,
    if (emp) "empirical resample of the observed slopes"
    else "single normal with the same mean and SD")
  t0 <- Sys.time()
  z <- try(analyse(scen, seed = 9000 + r), silent = TRUE)
  if (inherits(z, "try-error")) { P("  failed: %s", substr(z, 1, 180)); next }
  P("  elapsed %.0f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))
  for (i in seq_len(nrow(z$sel)))
    P("      K=%d   BIC %8.1f   entropy %s   smallest class %s%%",
      z$sel$K[i], z$sel$BIC[i],
      ifelse(is.na(z$sel$entropy[i]), "  -  ", sprintf("%.3f", z$sel$entropy[i])),
      ifelse(is.na(z$sel$minpct[i]), " - ", sprintf("%.1f", z$sel$minpct[i])))
  P("  BIC-selected K = %d ; class sizes %s", z$K, paste(z$class_sizes, collapse = "/"))
  P("  T1  fastest class == top-%d by rate : %.1f%% agreement", z$nfast, z$T1_agreement)
  P("  T2  rate only    : beta %+.3f  p = %.3g", z$T2_rate[1], z$T2_rate[3])
  P("  T2  class only   : beta %+.3f  p = %.3g", z$T2_class[1], z$T2_class[3])
  P("  T2  rate + class : rate  beta %+.3f  p = %.3g", z$T2_joint[1, 1], z$T2_joint[1, 3])
  P("                     class beta %+.3f  p = %.3g", z$T2_joint[2, 1], z$T2_joint[2, 3])
  if (!is.null(z$T3)) P("  T3  within the remaining knees: beta %+.3f  p = %.3g", z$T3[1], z$T3[3])
  res[[r]] <- z
  saveRDS(res, file.path(OUT, "sim3_continuum.rds"))
  writeLines(unlist(LOG), file.path(OUT, "sim3_continuum_log.txt"))
}

P("\nwritten -> Analysis/sim/sim3_continuum.rds")
