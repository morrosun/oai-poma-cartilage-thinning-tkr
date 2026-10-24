# ============================================================================
# 12 : methodological result 1 -- a shared index visit degenerates the
#      stratum-specific baseline hazard
#
#   Part A. Exact result (no simulation needed, but verified numerically)
#     With both members of a matched pair observed up to the SAME index visit
#     t*, the risk set of the pair at t* is the pair itself, so the pair's
#     contribution to the partial likelihood is
#
#         l_s = (g_s + eta_case) - log( exp(g_s + eta_case) + exp(g_s + eta_ctl) )
#             = log plogis(eta_case - eta_ctl),     g_s = h_{0s}(t*)
#
#     with eta_i = aV * (thickness level at t*) + aS * (thinning rate) the
#     within-pair linear predictor of the model that is actually fitted.
#
#     The stratum-specific baseline hazard h_{0s}(t) enters numerator and
#     denominator alike and cancels.  Hence the log-likelihood is exactly
#     invariant to the value of the baseline at the index visit and the Fisher
#     information with respect to it is the zero matrix.  The stratum-specific
#     baseline cancels from this CONDITIONAL EVENT likelihood.  That does not
#     establish its identifiability in a full joint-model likelihood, and it does
#     not by itself explain the implementation-specific failures recorded in the
#     software probe (see 07, 12b and Supplementary Table S3E).  Nothing here is
#     a statement about the posterior of any particular implementation.
#
#     v21 correction: this check previously multiplied aS (a coefficient per
#     mm/yr of thinning) by the marker LEVEL at the index visit, giving |eta| up
#     to 40 mm and a total of -1223.08 over 191 pairs.  The invariance itself was
#     unaffected, but that total was the log-likelihood of the mis-scaled
#     predictor rather than of the analysis model.  Eta is now the model's
#     within-pair linear predictor, and the diagnostic lines report the number of
#     pairs, the per-pair mean and the range of plogis(eta_case - eta_ctl) so
#     that the magnitude of the total can be checked at a glance.
#
#     v24 correction: the A2 information matrix took dltS = eta_case - eta_ctl,
#     which already contains aS, and then multiplied it by a again with the
#     starting value a = aS, so the Hessian was evaluated at aS^2 * rate
#     difference instead of aS * rate difference.  dltS is now the raw
#     case-minus-control thinning-rate difference and a is the coefficient per
#     mm/yr, so a = aS reproduces the generating predictor exactly.
#
#   Part B. Consequence for the paired contrast
#     The same algebra shows what IS identified: only within-pair contrasts
#     survive.  Anything that is constant within a pair -- the baseline hazard,
#     the pair-level component of a covariate, the matching variables --
#     cancels.  An unstratified analysis therefore pools risk sets across
#     pairs and re-introduces exactly the confounding that matching removed.
#
#   Part C. Estimator comparison over replicates (calibrated to POMA)
#     conditional likelihood on the TRUE slope (oracle)
#     conditional likelihood on the EBLUP slope (what the paper uses)
#     conditional likelihood on the unshrunk OLS slope
#     unstratified (marginal) partial likelihood on the EBLUP slope
# ============================================================================
source("D:/BaiduSyncdisk/OAI/Scripts/sim/sim_core.R")

BASE   <- "D:/BaiduSyncdisk/OAI"
OUT    <- file.path(BASE, "Analysis/sim")
design <- read_design(file.path(OUT, "sim_visit_schedule.csv"))
CAL    <- readRDS(file.path(OUT, "sim_as.rds"))
aS     <- CAL$aS; CONST <- CAL$const

cat("==========================================================================\n")
cat("SIMULATION 1 : shared index visit and the baseline hazard\n")
cat("  aS (true rate effect on the log-odds) =", round(aS, 4), "\n")
cat("  target observed EBLUP coefficient    =", CAL$target, "\n")
cat("==========================================================================\n\n")

LOG <- list()
P <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG[[length(LOG) + 1]] <<- s }

## =========================================================== Part A =========
P("--- Part A : the stratum baseline hazard is EXACTLY unidentifiable ---")

## A1: numeric invariance check.  Draw a simulated data set, then evaluate the
## pair-level log-likelihood at arbitrary stratum-specific baselines.  Within a
## stratum both members are evaluated at the SAME index visit, so the baseline
## value h_{0s}(t*) enters numerator and denominator alike and cancels.
##
## The linear predictor is the one the fitted model actually has,
##   eta_i = aV * (level at t*) + aS * (thinning rate),
## evaluated here with aV = 0 because gen_data() below is called with aV = 0.
## (v21: the earlier version of this check multiplied aS -- a coefficient per
## mm/yr of thinning -- by the marker LEVEL at the index visit.  Invariance was
## unaffected, but |eta| reached 40 mm and the reported total (-1223.08) was the
## log-likelihood of that mis-scaled predictor, not of the analysis model.)
d0 <- gen_data(design, NULL, CONST$tau, CONST$mu1, CONST$b0mu, CONST$b0sd,
               CONST$sigma_e, aS = aS, seed = 7)
d1 <- d0[!duplicated(d0$unit), ]               # one row per knee
cs <- which(d1$case == 1); js <- which(d1$case == 0)
stopifnot(length(cs) == length(js), length(cs) == length(unique(d0$pair)))
eta_c <- aS * (-d1$b1[cs])                     # aV = 0 -> rate term only
eta_j <- aS * (-d1$b1[js])
pairll <- function(ec, ej, g)                  # g = h_{0s}(t*), one value/pair
  (g + ec) - log(exp(g + ec) + exp(g + ej))
set.seed(99)
g_s <- rnorm(length(cs), 0, 3)                 # arbitrary per-stratum baseline value
ll0 <- sum(pairll(eta_c, eta_j, rep(0, length(cs))))
ll1 <- sum(pairll(eta_c, eta_j, g_s))
## Third check: two structurally different baseline FUNCTIONS of time that agree
## at the index visit.  Each is evaluated at t*, which is the only time the
## pair likelihood ever visits; f2 is deliberately not constant away from t*.
ts <- d1$tstar[cs]
f1 <- function(tt, ts) rep(0, length(tt))          # constant baseline
f2 <- function(tt, ts) 2.5 * (tt - ts) / 12        # linear, zero exactly at t*
ll2 <- sum(pairll(eta_c, eta_j, f2(ts, ts)))       # = f1(ts) = 0 at the visit
P("  A1  log-likelihood, baseline h_{0s}(t*) = 0        : %.10f", ll0)
P("  A1  log-likelihood, arbitrary h_{0s}(t*)        : %.10f", ll1)
P("  A1  log-likelihood, functions agreeing at t*    : %.10f  (f1 %.10f / f2 %.10f)",
  ll2, sum(pairll(eta_c, eta_j, f1(ts, ts))), ll2)
P("  A1  f1 and f2 differ away from t*: f2(t*+6) = %.3f (first pair)",
  f2(ts[1] + 6, ts[1]))
P("  A1  absolute change (both comparisons)          : %.3e  %.3e  <- machine zero",
  abs(ll1 - ll0), abs(ll2 - ll0))
P("  A1  pairs summed / mean contribution per pair   : %d / %.4f",
  length(cs), ll0 / length(cs))
P("  A1  plogis(eta_case - eta_control) : min %.3f  median %.3f  max %.3f",
  min(plogis(eta_c - eta_j)), median(plogis(eta_c - eta_j)),
  max(plogis(eta_c - eta_j)))
stopifnot(abs(ll1 - ll0) < 1e-8, abs(ll2 - ll0) < 1e-8,
          abs(ll0 - sum(plogis(eta_c - eta_j, log.p = TRUE))) < 1e-8,
          all(plogis(eta_c - eta_j) > 0.001),
          abs(f2(ts[1] + 6, ts[1])) > 0.1)   # the two functions really differ

## A2: Fisher information.  Each stratum baseline is a three-coefficient
## function of time, evaluated at t*; ALL THREE coefficients enter the
## likelihood, so a null direction is produced by cancellation and not by a
## coefficient never being used.  Stack (a, gamma_1..gamma_S) and require every
## gamma direction to be null.  S is kept small only to keep the numeric
## Hessian cheap -- the result is algebraic and does not depend on S.
## `a` is the association coefficient per mm/yr of thinning, so the regressor is
## the RAW case-minus-control thinning-rate difference (not eta, which already
## carries aS).  At a = aS this reproduces the generating predictor exactly.
## Tolerance: 1e-5 x largest absolute eigenvalue, five orders of magnitude above
## the rounding level of the finite-difference Hessian; the largest direction it
## excludes is printed so the gap is visible.
S <- 12; K <- 3
rate  <- -d1$b1                                     # thinning rate, mm/yr
dltS  <- (rate[cs] - rate[js])[seq_len(S)]           # raw within-pair difference
stopifnot(max(abs(aS * dltS - (eta_c - eta_j)[seq_len(S)])) < 1e-10)
tsS  <- d1$tstar[cs][seq_len(S)]                    # index visit of each stratum
BAS  <- cbind(1, tsS / 12, (tsS / 12)^2)            # baseline basis in time
mk_ll <- function(par) {
  a <- par[1]; g <- matrix(par[-1], nrow = S, ncol = K)
  off <- as.vector(rowSums(g * BAS))                # h_{0s}(t*) for the stratum
  sum((off + a * dltS) - log(exp(off + a * dltS) + exp(off)))
}
par0 <- c(aS, rep(0, S * K))
H <- matrix(0, length(par0), length(par0))
h <- 1e-4
for (i in seq_along(par0)) for (j in seq_along(par0)) {
  pp <- par0; pp[i] <- pp[i] + h; pp[j] <- pp[j] + h; fpp <- mk_ll(pp)
  pm <- par0; pm[i] <- pm[i] - h; pm[j] <- pm[j] - h; fmm <- mk_ll(pm)
  pq <- par0; pq[i] <- pq[i] + h; pq[j] <- pq[j] - h; fpq <- mk_ll(pq)
  qp <- par0; qp[i] <- qp[i] - h; qp[j] <- qp[j] + h; fqp <- mk_ll(qp)
  H[i, j] <- -(fpp - fpq - fqp + fmm) / (4 * h * h)
}
ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
tol <- max(abs(ev)) * 1e-5
n_null <- sum(abs(ev) <= tol)
rnd <- max(abs(ev)[abs(ev) <= tol])            # largest null-direction eigenvalue
P("  A2  information parameters : %d  (1 association + %d strata x %d basis)",
  length(par0), S, K)
P("  A2  null directions (|eigen| <= 1e-5 x max) : %d of %d", n_null, length(ev))
P("  A2  largest |eigenvalue|            : %.4g", max(abs(ev)))
P("  A2  largest null-direction |eigen|  : %.3e  (rounding level)", rnd)

## machine-readable copy of panel (A) for the number audit
write.csv(data.frame(
  n_pairs            = length(cs),
  ll_zero_baseline   = ll0,
  ll_arbitrary       = ll1,
  ll_other_function  = ll2,
  abs_change         = max(abs(ll1 - ll0), abs(ll2 - ll0)),
  mean_per_pair      = ll0 / length(cs),
  plogis_min         = min(plogis(eta_c - eta_j)),
  plogis_median      = median(plogis(eta_c - eta_j)),
  plogis_max         = max(plogis(eta_c - eta_j)),
  n_params           = length(par0),
  null_directions    = n_null,
  max_abs_eigen      = max(abs(ev)),
  next_abs_eigen     = rnd,
  tol_rel            = 1e-5
), file.path(OUT, "sim1a_invariance.csv"), row.names = FALSE)

## A3: what the unstratified analysis actually pools.  Diagnostic on the REAL
## 191-pair design: how many knees are at risk at an index visit, and how many
## of them belong to the matched pair that fails there.
Sd <- read.csv(file.path(OUT, "sim_visit_schedule.csv"), stringsAsFactors = FALSE)
Sd$t_idx <- Sd$t_index / 12
P("  A3  risk-set composition in an UNSRATIFIED partial likelihood:")
for (tt in sort(unique(Sd$t_idx))) {
  n <- sum(Sd$t_idx >= tt)
  P("        index visit %g yr : risk set %3d knees, 2 of them from the index pair (%.2f%%)",
    tt, n, 200 / n)
}

## =========================================================== Part C =========
run_rep <- function(r, nrep, zsd, aZ, tag) {
  cols <- c("b_true", "b_eblup", "b_ols", "b_marg_eblup", "b_marg_true")
  out <- matrix(NA_real_, nrep, length(cols) * 2,
                dimnames = list(NULL, c(cols, sub("^b_", "se_", cols))))
  for (i in seq_len(nrep)) {
    d <- gen_data(design, NULL, CONST$tau, CONST$mu1, CONST$b0mu, CONST$b0sd,
                  CONST$sigma_e, aS = aS, aV = 0, seed = r + i,
                  zsd = zsd, aZ = aZ)
    kn <- add_truth(knee_slopes(d), d)
    out[i, c("b_true", "se_true")]          <- clogit_fit(kn, "slope_true")
    out[i, c("b_eblup", "se_eblup")]        <- clogit_fit(kn, "slope_eblup")
    out[i, c("b_ols", "se_ols")]            <- clogit_fit(kn, "slope_ols")
    out[i, c("b_marg_eblup", "se_marg_eblup")] <- cox_fit(kn, "slope_eblup")
    out[i, c("b_marg_true", "se_marg_true")]   <- cox_fit(kn, "slope_true")
  }
  attr(out, "tag") <- tag
  out
}

REPS <- as.integer(Sys.getenv("SIM_REPS", "250"))
P("\n--- Part C : estimator comparison (%d replicates per scenario) ---", REPS)

scen <- list(
  list(tag = "C1  no pair-level confounder",       zsd = 0.00, aZ = 0),
  list(tag = "C2  pair-level confounder present",  zsd = 0.10, aZ = 1.0)
)
LBL <- c(b_true = "conditional, TRUE slope (oracle)",
         b_eblup = "conditional, EBLUP slope",
         b_ols = "conditional, OLS slope",
         b_marg_eblup = "marginal, EBLUP slope",
         b_marg_true = "marginal, TRUE slope")
store <- list()
for (sc in scen) {
  M <- run_rep(0, REPS, sc$zsd, sc$aZ, sc$tag)
  store[[sc$tag]] <- M
  P("")
  P("  %s", sc$tag)
  P("  %-34s %8s %8s %8s %8s", "estimator", "mean", "bias", "sd", "coverage")
  for (v in names(LBL)) {
    b <- M[, v]; s <- M[, sub("^b_", "se_", v)]
    ok <- is.finite(b) & is.finite(s)
    cv <- if (any(ok)) mean(b[ok] - 1.96 * s[ok] <= aS & aS <= b[ok] + 1.96 * s[ok]) else NA
    P("  %-34s %8.3f %8.3f %8.3f %8.3f", LBL[[v]], mean(b, na.rm = TRUE),
      mean(b, na.rm = TRUE) - aS, sd(b, na.rm = TRUE), cv)
  }
}

saveRDS(list(log = LOG, scen = store, aS = aS),
        file.path(OUT, "sim1_degeneracy.rds"))
writeLines(unlist(LOG), file.path(OUT, "sim1_degeneracy_log.txt"))
P("\nwritten -> Analysis/sim/sim1_degeneracy.rds")
