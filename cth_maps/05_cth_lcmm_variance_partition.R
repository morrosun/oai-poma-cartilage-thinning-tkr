# ============================================================================
# 05 : how much of the LEVEL and how much of the RATE do the LCMM classes explain?
#
#   A latent-class model of a trajectory can, in principle, separate knees on two
#   axes: where they start (level) and how fast they change (rate).  This script
#   quantifies which axis the converged solutions actually use, with eta-squared
#   (the share of between-knee variance in each variable that the class labels
#   account for).
#
#   It also prints the coefficient names of the K = 3 model so that the parameter
#   layout of lcmm::hlme can be read off rather than assumed.
#
#   Inputs   cth_lcmm_k3_diag.rds   (04)
#            cth_lcmm_k2_fit.rds    (04)
#   Outputs  cth_lcmm_variance_partition.csv
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/oai_imaging"
suppressPackageStartupMessages(library(lcmm))

eta2 <- function(y, g) {
  ok <- is.finite(y)
  y <- y[ok]; g <- g[ok]
  ybar <- mean(y)
  ss_tot <- sum((y - ybar)^2)
  ss_within <- sum(vapply(split(y, g), function(z) sum((z - mean(z))^2), numeric(1)))
  1 - ss_within / ss_tot
}

d3 <- readRDS(file.path(BASE, "cth_lcmm_k3_diag.rds"))
d2 <- readRDS(file.path(BASE, "cth_lcmm_k2_fit.rds"))
k3 <- d3$k3; m3 <- d3$m3
k2 <- d2$knee

cat("================ K = 3 model, as reported by lcmm ================\n")
print(summary(m3))
cat("\nraw parameter vector (`best`), in this order:\n")
print(m3$best)
cat("class-membership intercepts (K-1)\n")
cat("-> longitudinal intercepts, one per class\n")
cat("-> longitudinal slopes, one per class\n")
cat("-> random-intercept variance\n")
cat("-> residual standard error\n")
cat("note: the class-membership model is parameterised relative to the LAST class\n")

K3 <- length(unique(k3$cl))
npar3 <- length(m3$best)
off_int <- (K3 - 1) + 1
imp <- data.frame(class = seq_len(K3),
                  model_intercept = m3$best[off_int + seq_len(K3) - 1],
                  model_slope = m3$best[off_int + K3 + seq_len(K3) - 1],
                  pct = 100 * as.numeric(table(k3$cl)) / nrow(k3),
                  obs_baseline_mean = as.numeric(tapply(k3$b0, k3$cl, mean)),
                  obs_slope_mean = as.numeric(tapply(k3$ols, k3$cl, mean, na.rm = TRUE)))
cat("\n---- model-implied vs observed class profiles (K = 3) ----\n")
print(imp, row.names = FALSE, digits = 3)

cat("\n---- variance partition by the class labels ----\n")
res <- rbind(
  data.frame(model = "LCMM K=2", variable = "baseline level",
             eta2 = eta2(k2$b0, k2$cl), n_classes = length(unique(k2$cl))),
  data.frame(model = "LCMM K=2", variable = "progression rate (OLS slope)",
             eta2 = eta2(k2$ols, k2$cl), n_classes = length(unique(k2$cl))),
  data.frame(model = "LCMM K=3", variable = "baseline level",
             eta2 = eta2(k3$b0, k3$cl), n_classes = length(unique(k3$cl))),
  data.frame(model = "LCMM K=3", variable = "progression rate (OLS slope)",
             eta2 = eta2(k3$ols, k3$cl), n_classes = length(unique(k3$cl)))
)
print(res, row.names = FALSE, digits = 4)
cat(sprintf("\n  K=3 : the classes account for %.0f%% of the variance in the LEVEL but only %.0f%% of the variance in the RATE\n",
            100 * res$eta2[3], 100 * res$eta2[4]))
cat(sprintf("  K=2 : %.0f%% of the LEVEL vs %.0f%% of the RATE\n",
            100 * res$eta2[1], 100 * res$eta2[2]))
cat("  a trajectory classification only earns its name if the second number is the larger one\n")
write.csv(res, file.path(BASE, "cth_lcmm_variance_partition.csv"), row.names = FALSE)
write.csv(imp, file.path(BASE, "cth_lcmm_k3_model_profile.csv"), row.names = FALSE)

## ---- is the "fast" class anything more than an upper-tail cutoff on the rate? ----
## If a nominal latent class reproduces, knee for knee, the set you would pick by
## simply taking the top n knees by observed rate, then the class carries no
## information beyond a rank threshold -- it is a dichotomisation, not a discovery.
cat("\n---- test: is the fastest class just a rank cutoff on the rate? ----\n")
cut_tab <- data.frame()
for (nm in c("K=2", "K=3")) {
  kk <- if (nm == "K=2") k2 else k3
  ms <- tapply(kk$ols, kk$cl, mean, na.rm = TRUE)
  fast <- as.integer(names(which.max(ms)))
  sel <- kk$cl == fast
  ns <- sum(sel)
  thr <- sort(kk$ols, decreasing = TRUE)[ns]      # top-n by observed rate
  topn <- kk$ols >= thr
  agree <- mean(sel == topn)                      # simple agreement rate
  jac <- sum(sel & topn) / sum(sel | topn)        # Jaccard
  cat(sprintf("  %s : fastest class n=%d ; taking the top %d knees by observed rate gives\n",
              nm, ns, ns))
  cat(sprintf("        agreement %.1f%%, Jaccard %.2f ; cutoff %.2f pts/yr vs class minimum %.2f\n",
              100 * agree, jac, thr, min(kk$ols[sel])))
  cat(sprintf("        class mean rate %+.2f vs top-%d mean rate %+.2f\n",
              mean(kk$ols[sel]), ns, mean(kk$ols[topn])))
  cut_tab <- rbind(cut_tab, data.frame(model = nm, class_n = ns, cutoff = thr,
                                       class_min = min(kk$ols[sel]),
                                       agreement = agree, jaccard = jac,
                                       class_mean_rate = mean(kk$ols[sel]),
                                       topn_mean_rate = mean(kk$ols[topn])))
}
write.csv(cut_tab, file.path(BASE, "cth_lcmm_class_vs_cutoff.csv"), row.names = FALSE)

cat("\n---- composition of the fastest class in each solution ----\n")
for (nm in c("K=2", "K=3")) {
  kk <- if (nm == "K=2") k2 else k3
  ms <- tapply(kk$ols, kk$cl, mean, na.rm = TRUE)
  fast <- as.integer(names(which.max(ms)))
  x <- kk[kk$cl == fast, ]
  cat(sprintf("  %s : fastest class %d holds %d knees (%.1f%%), mean rate %+.2f pts/yr (sd %.2f);\n",
              nm, fast, nrow(x), 100 * nrow(x) / nrow(kk), ms[[which.max(ms)]], sd(x$ols)))
  cat(sprintf("        rate range within it %.2f to %.2f pts/yr -- the class is not rate-homogeneous\n",
              min(x$ols), max(x$ols)))
}
cat("\n-> cth_lcmm_variance_partition.csv\n")
