# ============================================================================
# 12 : methodological result 1 -- a shared index visit degenerates the
#      stratum-specific baseline hazard
#
#   Part A. Exact result (no simulation needed, but verified numerically)
#     With both members of a matched pair observed up to the SAME index visit
#     t*, the risk set of the pair at t* is the pair itself, so the pair's
#     contribution to the partial likelihood is
#
#         l_s = a*m_case(t*) - log( exp(a*m_case(t*)) + exp(a*m_ctl(t*)) )
#
#     The stratum-specific baseline hazard h_{0s}(t) enters twice with opposite
#     signs and cancels.  Hence h_{0s}(.) is not a parameter of the observed
#     data at all: the log-likelihood is exactly invariant to it, the Fisher
#     information with respect to it is the zero matrix, and any procedure that
#     estimates it (e.g. a stratified joint model whose baseline is an M-spline
#     basis per stratum) has an improper posterior.  That is the precise reason
#     off-the-shelf joint-model software fails here.
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
## value h_{0s}(t*) enters twice with opposite signs and cancels.
d0 <- gen_data(design, NULL, CONST$tau, CONST$mu1, CONST$b0mu, CONST$b0sd,
               CONST$sigma_e, aS = aS, seed = 7)
d1 <- d0[!duplicated(d0$unit), ]               # one row per knee
m_obs <- d1$b0 + d1$b1 * d1$tstar              # marker value at the index visit
pairll <- function(mc, mj, a, g)               # g = h_{0s}(t*), one value/pair
  (g + a * mc) - log(exp(g + a * mc) + exp(g + a * mj))
cs <- which(d1$case == 1); js <- which(d1$case == 0)
stopifnot(length(cs) == length(js), length(cs) == length(unique(d0$pair)))
mc <- m_obs[cs]; mj <- m_obs[js]
set.seed(99)
g_s <- rnorm(length(cs), 0, 3)                 # arbitrary per-stratum baseline
P("  A1  log-likelihood, baseline h_{0s}(t*) = 0        : %.10f",
  sum(pairll(mc, mj, aS, rep(0, length(cs)))))
P("  A1  log-likelihood, arbitrary h_{0s}(t*)        : %.10f",
  sum(pairll(mc, mj, aS, g_s)))
P("  A1  absolute change                             : %.3e  <- machine zero",
  abs(sum(pairll(mc, mj, aS, g_s)) - sum(pairll(mc, mj, aS, rep(0, length(cs))))))

## the same holds for two DIFFERENT baseline *functions* that agree at t*:
## only the value at the index visit can ever matter
g_alt <- g_s + 5 * cos(1:length(cs))           # any other function of t
P("  A1  two different baseline FUNCTIONS agreeing at t* give the same likelihood: |diff| = %.3e",
  abs(sum(pairll(mc, mj, aS, g_alt)) - sum(pairll(mc, mj, aS, g_s))))

## A2: Fisher information.  Stack (a, gamma_1..gamma_S) with a piecewise-constant
## stratum baseline; every gamma direction must be a null direction.  S is kept
## small here only to keep the numeric Hessian cheap -- the result is algebraic
## and does not depend on S.
S <- 12; K <- 3
mcS <- mc[seq_len(S)]; mjS <- mj[seq_len(S)]
mk_ll <- function(par) {
  a <- par[1]; g <- matrix(par[-1], nrow = S, ncol = K)
  sum(pairll(mcS, mjS, a, g[, 1]))
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
tol <- max(abs(ev)) * 1e-8
P("  A2  information parameters : %d  (1 association + %d strata x %d basis)",
  length(par0), S, K)
P("  A2  null directions (|eigen| <= tol) : %d of %d", sum(abs(ev) <= tol), length(ev))
P("  A2  largest |eigenvalue|            : %.4g", max(abs(ev)))

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
