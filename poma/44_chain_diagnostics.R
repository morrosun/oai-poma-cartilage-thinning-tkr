# ============================================================================
# 44 : single-chain sampling adequacy for the PRIMARY chain (reported in 2.5.2)
#
#   Replaces the v4 figures, which were computed from the legacy kernel.  All
#   quantities are recomputed here from the corrected sampler on the primary
#   cMFTC Set A fit (seed 2026), so that the numbers quoted in the manuscript
#   come from the same chain as the reported odds ratio.
#
#   -> chain_diagnostics_primary.csv
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
source("D:/BaiduSyncdisk/OAI/Scripts/poma/jm_core.R")

L <- read.csv(file.path(BASE, "poma_long.csv"))
U <- read.csv(file.path(BASE, "poma_analysis_long.csv"))
COVARS <- c("V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")
metric <- "cMFTC_ThCtAB_aMe"
sub <- L[!is.na(L[[metric]]), ]
lg  <- data.frame(unit = sub$unit, year = sub$years, y = sub[[metric]],
                  newstrata = sub$newstrata, case = sub$case)
sv  <- U[, c("unit", "newstrata", "case", "t_index_months", COVARS, "V00PASE")]
sv  <- sv[sv$unit %in% lg$unit, ]
sv  <- sv[complete.cases(sv[, COVARS]), ]
lg  <- lg[lg$unit %in% sv$unit, ]
sv$t_idx <- sv$t_index_months / 12
A <- prune_pairs(lg, sv, 2); dt <- build(A$lg, A$sv)

f <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 2026, verbose = FALSE)
d <- f$draws
N <- length(d$aV)

## autocorrelation time (initial positive sequence estimator)
act <- function(x) {
  x <- x - mean(x); n <- length(x); v <- sum(x^2) / n
  s <- 0
  for (lag in 1:min(500, n - 1)) {
    r <- sum(x[1:(n - lag)] * x[(lag + 1):n]) / (n * v)
    if (r <= 0) break
    s <- s + r
  }
  1 + 2 * s
}
lag1 <- function(x) {
  x <- x - mean(x)
  sum(x[1:(N - 1)] * x[2:N]) / sum(x^2)
}

orsd <- exp(-d$aS * d$SDslo)
ess_aS <- N / act(d$aS)
half <- split(seq_len(N), rep(1:2, each = N / 2))
or_h1 <- median(exp(-d$aS[half[[1]]] * d$SDslo[half[[1]]]))
or_h2 <- median(exp(-d$aS[half[[2]]] * d$SDslo[half[[2]]]))

res <- data.frame(
  quantity = c("draws retained", "OR per SD (median)",
               "OR per SD 95% CrI lower", "OR per SD 95% CrI upper",
               "ESS aS", "autocorrelation time aS", "lag-1 autocorrelation aS",
               "MC-SE of aS (in units of posterior SD)", "MC-SE on log-OR scale",
               "first-half OR per SD", "second-half OR per SD",
               "MH acceptance, random effects", "MH acceptance, association",
               "fraction of iterations with proposal scale capped"),
  value = c(N, median(orsd), quantile(orsd, .025), quantile(orsd, .975),
            ess_aS, act(d$aS), lag1(d$aS),
            sd(d$aS) / sqrt(ess_aS) / sd(d$aS),
            sd(-d$aS * d$SDslo) / sqrt(N / act(d$aS)),
            or_h1, or_h2, f$acc_b, f$acc_a, f$a_scale_capped))
write.csv(res, file.path(BASE, "chain_diagnostics_primary.csv"), row.names = FALSE)

cat("=========== PRIMARY CHAIN SAMPLING ADEQUACY (corrected sampler) ===========\n")
cat(sprintf("OR per SD              : %.4f (%.3f, %.3f)\n",
            median(orsd), quantile(orsd, .025), quantile(orsd, .975)))
cat(sprintf("ESS(aS)                : %.0f  of %d draws\n", ess_aS, N))
cat(sprintf("autocorrelation time   : %.2f   (lag-1 %.3f)\n", act(d$aS), lag1(d$aS)))
cat(sprintf("MC-SE / posterior SD   : %.4f\n", sd(d$aS) / sqrt(ess_aS) / sd(d$aS)))
cat(sprintf("MC-SE on log-OR scale  : %.4f\n", sd(-d$aS * d$SDslo) / sqrt(ess_aS)))
cat(sprintf("halves                 : OR/SD %.3f vs %.3f\n", or_h1, or_h2))
cat(sprintf("acceptance             : random effects %.3f | association %.3f\n",
            f$acc_b, f$acc_a))
cat(sprintf("proposal scale capped  : %.4f of iterations\n", f$a_scale_capped))
cat("\n-> chain_diagnostics_primary.csv\n")
