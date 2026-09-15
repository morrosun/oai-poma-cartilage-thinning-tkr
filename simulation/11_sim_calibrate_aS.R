# ============================================================================
# 11 : calibrate the true rate effect aS so that the simulation reproduces the
#      observed POMA estimate when the exposure is the EBLUP slope
#
#   Target (real data, primary set 191 pairs):
#       within-pair conditional logistic, exposure = cMFTC EBLUP slope
#       beta = +5.180 per mm/yr of loss   (se 1.154)  ->  OR/SD = 2.08
#
#   aS is the coefficient of (-b1) in the generating hazard.  Because the
#   EBLUP slope is a shrunken version of the true slope, E[beta_hat_EBLUP]
#   is attenuated relative to aS, so aS must be solved numerically.
# ============================================================================
source("D:/BaiduSyncdisk/OAI/Scripts/sim/sim_core.R")

BASE <- "D:/BaiduSyncdisk/OAI"
design <- read_design(file.path(BASE, "Analysis/sim/sim_visit_schedule.csv"))

CONST <- list(tau = 0.1779, mu1 = -0.1455, b0mu = 3.3822, b0sd = 1.5476,
              sigma_e = 0.2317)
TARGET <- 5.180          # observed EBLUP-based conditional-logistic coefficient

est_beta <- function(aS, R = 60, seed0 = 1000) {
  b <- numeric(R)
  for (r in seq_len(R)) {
    d <- gen_data(design, NULL, CONST$tau, CONST$mu1, CONST$b0mu, CONST$b0sd,
                  CONST$sigma_e, aS = aS, seed = seed0 + r)
    kn <- knee_slopes(d)
    b[r] <- clogit_fit(kn, "slope_eblup")[["beta"]]
  }
  c(mean = mean(b, na.rm = TRUE), sd = sd(b, na.rm = TRUE),
    fail = sum(!is.finite(b)))
}

cat("=== calibrating aS (target mean beta_EBLUP =", TARGET, ") ===\n")
grid <- c(5, 8, 11, 14, 18)
res <- data.frame()
for (a in grid) {
  e <- est_beta(a)
  cat(sprintf("  aS = %5.2f -> beta_EBLUP %+7.3f (sd %.3f, fail %d)\n",
              a, e[["mean"]], e[["sd"]], e[["fail"]]))
  res <- rbind(res, data.frame(aS = a, beta = e[["mean"]]))
}
# linear interpolation of aS on the response (near-linear over this range)
fit <- lm(aS ~ beta, data = res)
aS_hat <- unname(predict(fit, data.frame(beta = TARGET)))
cat(sprintf("\nfirst-pass aS = %.3f\n", aS_hat))

for (it in 1:3) {
  e <- est_beta(aS_hat, R = 150, seed0 = 5000 + 1000 * it)
  cat(sprintf("iter %d: aS = %.3f -> beta_EBLUP %+7.3f (sd %.3f, fail %d)\n",
              it, aS_hat, e[["mean"]], e[["sd"]], e[["fail"]]))
  if (abs(e[["mean"]] - TARGET) < 0.12) break
  aS_hat <- aS_hat * TARGET / e[["mean"]]
}

saveRDS(list(aS = aS_hat, target = TARGET, const = CONST),
        file.path(BASE, "Analysis/sim/sim_as.rds"))
cat(sprintf("\nFROZEN aS = %.4f\n", aS_hat))
