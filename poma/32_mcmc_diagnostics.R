# ============================================================================
# 32 : multi-chain convergence diagnostics for the joint model
#
#   Reviewer point 2g/2h asked whether 30 000 / 5 000 burn-in / thin-5 is
#   enough and asked for a convergence statistic.  Section 3.1 of the review
#   response could only report single-chain ESS / MC-SE / half-chain stability,
#   because the retained objects were single chains.  This script supplies a
#   real Gelman-Rubin R-hat from three DISPERSED-start chains on the primary
#   analysis (cMFTC, Set A, identical construction to 10_jm_multimetric.R),
#   together with a cross-chain effective sample size.
#
#   Dispersed starts are obtained through fit_jm(init_aVS = ), an optional
#   argument that leaves the default path bit-for-bit unchanged.
#
#   Output : mcmc_multichain_rhat.csv , mcmc_multichain_rhat.log
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
source("D:/BaiduSyncdisk/OAI/Scripts/poma/jm_core.R")

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
cat(sprintf("Set A : %d pairs / %d knees\n", dt$npair, dt$n))

## ---------------------------------------------------- helpers --------------
## standard Gelman-Rubin potential scale reduction (draws x chains matrix)
rhat <- function(M) {
  n <- nrow(M); m <- ncol(M)
  W  <- mean(apply(M, 2, var))
  B  <- n * var(colMeans(M))
  vh <- ((n - 1) / n) * W + B / n
  sqrt(vh / W)
}
## per-chain ESS by initial-positive-sequence autocorrelation, summed
ess_chain <- function(x) {
  n <- length(x)
  s <- 0
  for (lag in 1:(n - 1)) {
    r <- cor(x[1:(n - lag)], x[(lag + 1):n])
    if (!is.finite(r) || r < 0.05) break
    s <- s + r
  }
  n / (1 + 2 * s)
}
ess_total <- function(M) sum(apply(M, 2, ess_chain))

summ3 <- function(x) sprintf("%.3f (%.3f, %.3f)", median(x),
                             quantile(x, .025), quantile(x, .975))

## ---------------------------------------------------- run three chains -----
## dispersed associative-parameter starts: at zero, and two over-dispersed
## points of opposite sign (|a| far outside the posterior support)
STARTS <- list(c(0, 0), c(2, -2), c(-2, 2))
SEEDS  <- c(101, 202, 303)

chains <- vector("list", 3)
for (j in 1:3) {
  t0 <- Sys.time()
  cat(sprintf("\n---- chain %d : seed %d, aV0 = %+.2f, aS0 = %+.2f ----\n",
              j, SEEDS[j], STARTS[[j]][1], STARTS[[j]][2]))
  f <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = SEEDS[j],
              init_aVS = STARTS[[j]], verbose = FALSE)
  d <- f$draws
  cat(sprintf("  elapsed %.1f s | MH acc (random effects) %.2f\n",
              as.numeric(difftime(Sys.time(), t0, units = "secs")), f$acc_b))
  cat("  aS (per mm/yr latent)    :", summ3(d$aS), "\n")
  cat("  OR / 1 SD faster         :", summ3(exp(-d$aS * d$SDslo)), "\n")
  cat("  OR / 0.1 mm/yr faster    :", summ3(exp(-0.1 * d$aS)), "\n")
  chains[[j]] <- d
}

## ---------------------------------------------------- diagnostics ----------
pars <- c("aS", "aV", "beta1", "s2", "SDslo", "SDval")
res <- data.frame()
for (p in pars) {
  M <- sapply(chains, function(d) d[[p]])
  colnames(M) <- paste0("chain", 1:3)
  res <- rbind(res, data.frame(parameter = p,
                               Rhat = rhat(M), ESS_cross_chain = ess_total(M),
                               q025 = quantile(M, .025), median = median(M),
                               q975 = quantile(M, .975),
                               check.names = FALSE))
}
## derived quantities, computed draw-wise inside each chain
ORSD <- sapply(chains, function(d) exp(-d$aS * d$SDslo))
OR01 <- sapply(chains, function(d) exp(-0.1 * d$aS))
ORSDlev <- sapply(chains, function(d) exp(-d$aV * d$SDval))
for (nm in c("OR_per_SD_rate", "OR_per_0.1mm_yr", "OR_per_SD_level")) {
  M <- switch(nm, OR_per_SD_rate = ORSD, OR_per_0.1mm_yr = OR01, OR_per_SD_level = ORSDlev)
  res <- rbind(res, data.frame(parameter = nm,
                               Rhat = rhat(M), ESS_cross_chain = ess_total(M),
                               q025 = quantile(M, .025), median = median(M),
                               q975 = quantile(M, .975),
                               check.names = FALSE))
}

cat("\n================ multi-chain diagnostics (3 dispersed chains) ============\n")
print(res, row.names = FALSE, digits = 4)
cat(sprintf("\nworst R-hat = %.4f  (threshold 1.01)\n", max(res$Rhat)))
cat(sprintf("OR per SD of thinning rate, pooled 3 chains : %s\n",
            summ3(as.vector(ORSD))))

## posterior dependence between the level and rate association parameters,
## the quantity behind the "banana-shaped" description
CV <- sapply(chains, function(d) cor(d$aV, d$aS))
pooled <- cor(unlist(lapply(chains, function(d) d$aV)),
              unlist(lapply(chains, function(d) d$aS)))
cat(sprintf("posterior cor(aV, aS) : per chain %s ; pooled %.3f\n",
            paste(sprintf("%.3f", CV), collapse = " / "), pooled))
write.csv(data.frame(chain = c(1:3, "pooled"), cor_aV_aS = c(CV, pooled)),
          file.path(BASE, "mcmc_posterior_correlation.csv"), row.names = FALSE)

write.csv(res, file.path(BASE, "mcmc_multichain_rhat.csv"), row.names = FALSE)
cat("\n-> mcmc_multichain_rhat.csv / mcmc_posterior_correlation.csv\n")
