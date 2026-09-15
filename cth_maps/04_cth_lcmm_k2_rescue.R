# ============================================================================
# 04 : LCMM diagnostics for P2
#
#   (i)  What separates the three classes of the converged K = 3 solution --
#        baseline severity level, or rate of progression?  The random intercept
#        already absorbs between-knee level differences, so this decides whether
#        the "classes" carry any slope information at all.
#
#   (ii) K = 2 could not be reached from any of six k-means-seeded starts -- every
#        one collapsed to a solution with an empty class.  Here K = 2 is re-started
#        from the CONVERGED K = 3 solution (dropping one mixture component), which
#        guarantees a populated two-class start.
#
#   Reads the cache written by 01_cth_trajectory_lcmm.R; nothing written here is
#   used by 01, so this script is safe to run at any time.
#   Outputs  cth_lcmm_k3_diagnostic.csv
#            cth_lcmm_k2_rescue.csv
#            cth_lcmm_k3_diag.rds / cth_lcmm_k2_fit.rds
# ============================================================================

BASE  <- "D:/BaiduSyncdisk/OAI/Analysis/oai_imaging"
CTH   <- file.path(BASE, "cth_long.csv")
CACHE <- file.path(BASE, "cth_lcmm_cache.rds")

suppressPackageStartupMessages(library(lcmm))

## ------------------------------------------------------- rebuild the data ---
## identical construction to 01, so `unit` numbering matches the cached fits
raw <- read.csv(CTH)
raw <- raw[!is.na(raw$t) & !is.na(raw$score), ]
nv <- table(raw$uid)
keep <- names(nv)[nv >= 2]
d0 <- raw[raw$uid %in% keep, ]
d0$unit <- as.integer(factor(d0$uid))
d0 <- d0[order(d0$unit, d0$t), ]
N <- max(d0$unit)

b0   <- vapply(seq_len(N), function(i) { x <- d0[d0$unit == i, ]; x$score[which.min(x$t)] }, numeric(1))
ols  <- vapply(seq_len(N), function(i) { x <- d0[d0$unit == i, ]
  if (length(unique(x$t)) < 2) NA_real_ else unname(coef(lm(score ~ t, x))[2]) }, numeric(1))
nvis <- as.integer(table(d0$unit))
sd_ols <- sd(ols, na.rm = TRUE); sd_b0 <- sd(b0)

cc <- readRDS(CACHE)
fits <- cc$fits
cat(sprintf("cache holds K = %s ; cache N = %d, data N = %d\n",
            paste(sub("^K", "", names(fits)), collapse = ", "), cc$N, N))

m3 <- fits[["K3"]]
stopifnot(!is.null(m3))

## ------------------------------------------------------ K = 3 diagnostic ----
pp3 <- m3$pprob
pp3$unit <- as.integer(pp3$unit)
pcols <- grep("^prob", names(pp3), value = TRUE)
ord <- order(pp3$unit)                      # pprob is already one row per knee
stopifnot(identical(pp3$unit[ord], seq_len(N)))
k3 <- data.frame(unit = seq_len(N),
                 cl = pp3$class[ord],
                 b0 = b0, ols = ols, nvis = nvis)
k3 <- cbind(k3, pp3[ord, pcols, drop = FALSE])

prof <- do.call(rbind, lapply(sort(unique(k3$cl)), function(g) {
  x <- k3[k3$cl == g, ]
  data.frame(class = g, n = nrow(x), pct = 100 * nrow(x) / nrow(k3),
             baseline_mean = mean(x$b0), baseline_sd = sd(x$b0),
             slope_mean = mean(x$ols, na.rm = TRUE), slope_sd = sd(x$ols, na.rm = TRUE),
             visits = mean(x$nvis))
}))
cat("\n---- K = 3 class profile (unshrunk OLS slope, baseline CTh-Score) ----\n")
print(prof, row.names = FALSE, digits = 4)

lev_rng <- diff(range(prof$baseline_mean)) / sd_b0     # in baseline SD units
rate_rng <- diff(range(prof$slope_mean)) / sd_ols      # in slope SD units
cat("\n")
cat(sprintf("between-class LEVEL range : %.1f points = %.2f baseline SD\n",
            diff(range(prof$baseline_mean)), lev_rng))
cat(sprintf("between-class RATE  range : %.2f pts/yr = %.2f slope SD\n",
            diff(range(prof$slope_mean)), rate_rng))
cat(sprintf("-> the classes separate %.1fx more strongly on level than on rate\n",
            lev_rng / rate_rng))
cat(sprintf("Spearman (class baseline vs class slope) = %+.2f\n",
            cor(prof$baseline_mean, prof$slope_mean, method = "spearman")))
cat(sprintf("mean within-class slope sd / overall slope sd = %.2f  (1.0 = classes explain nothing)\n",
            mean(prof$slope_sd) / sd_ols))
cat(sprintf("mean within-class baseline sd / overall baseline sd = %.2f\n",
            mean(prof$baseline_sd) / sd_b0))

fx <- m3$best
stopifnot(length(fx) == 10)
cat("\nK = 3 model coefficients\n")
cat(sprintf("  fixed intercept %+.3f, fixed slope %+.4f\n", fx[1], fx[2]))
cat(sprintf("  mixture intercepts (class 1 / 2 / 3) : %s\n",
            paste(sprintf("%+.2f", c(0, fx[3], fx[5])), collapse = "  ")))
cat(sprintf("  mixture slopes     (class 1 / 2 / 3) : %s\n",
            paste(sprintf("%+.4f", c(0, fx[4], fx[6])), collapse = "  ")))
cat(sprintf("  proportions %s ; RE variance %.1f ; residual SD %.3f\n",
            paste(sprintf("%.3f", plogis(c(fx[7], fx[8]))), collapse = " / "),
            fx[9], fx[10]))
write.csv(prof, file.path(BASE, "cth_lcmm_k3_diagnostic.csv"), row.names = FALSE)

## ---------------------------------------------- K = 2 rescue from K = 3 ------
## B layout for ng classes, mixture = ~ t, random = ~ 1 :
##   fixed intercept, fixed t,
##   (mixture intercept, mixture t) for classes 2..K,
##   K - 1 class proportions (logit scale),
##   random-intercept variance, residual SD
##   -> length 3K + 1
b2 <- function(mi, mt) c(fx[1], fx[2], mi, mt, 0, fx[9], fx[10])

starts <- list("keep class 2"           = b2(fx[3], fx[4]),
               "keep class 3"           = b2(fx[5], fx[6]),
               "average classes 2 and 3" = b2((fx[3] + fx[5]) / 2, (fx[4] + fx[6]) / 2))

cat("\n---- K = 2 re-started from the converged K = 3 solution ----\n")
res <- data.frame(start = character(), loglik = numeric(), BIC = numeric(),
                  entropy = numeric(), sizes = character(), ok = logical())
best2 <- NULL
for (nm in names(starts)) {
  mm <- try(hlme(score ~ t, mixture = ~ t, random = ~ 1, subject = "unit",
                 ng = 2, data = d0, B = starts[[nm]], maxiter = 200), silent = TRUE)
  if (inherits(mm, "try-error")) {
    cat(sprintf("  %-24s : error  %s\n", nm, substr(attr(mm, "condition")$message, 1, 70)))
    res <- rbind(res, data.frame(start = nm, loglik = NA_real_, BIC = NA_real_,
                                 entropy = NA_real_, sizes = NA_character_, ok = FALSE))
    next
  }
  cs <- table(mm$pprob$class)
  okk <- length(cs) == 2 && min(cs) >= 20 && is.finite(mm$loglik) && mm$loglik < 1e8
  P <- as.matrix(mm$pprob[, grep("^prob", names(mm$pprob)), drop = FALSE])
  en <- 1 - sum(-P * log(pmax(P, 1e-12))) / (nrow(P) * log(2))
  cat(sprintf("  %-24s : loglik %.1f  BIC %.1f  entropy %.3f  sizes %s  %s\n",
              nm, mm$loglik, mm$BIC, en, paste(cs, collapse = "/"),
              if (okk) "ACCEPTED" else "rejected (degenerate)"))
  res <- rbind(res, data.frame(start = nm, loglik = mm$loglik, BIC = mm$BIC,
                               entropy = en, sizes = paste(cs, collapse = "/"), ok = okk))
  if (okk && (is.null(best2) || mm$loglik > best2$loglik)) best2 <- mm
}

if (!is.null(best2)) {
  pp <- best2$pprob; pp$unit <- as.integer(pp$unit); ord2 <- order(pp$unit)
  kk <- data.frame(unit = seq_len(N), cl = pp$class[ord2], b0 = b0, ols = ols)
  pr <- do.call(rbind, lapply(sort(unique(kk$cl)), function(g) {
    x <- kk[kk$cl == g, ]
    data.frame(class = g, n = nrow(x), pct = 100 * nrow(x) / nrow(kk),
               baseline_mean = mean(x$b0), slope_mean = mean(x$ols, na.rm = TRUE),
               slope_sd = sd(x$ols, na.rm = TRUE))
  }))
  cat("\n  accepted K = 2 profile:\n"); print(pr, row.names = FALSE, digits = 4)
  cat(sprintf("  level gap %.1f points (%.2f SD) vs rate gap %.2f pts/yr (%.2f SD)\n",
              abs(diff(pr$baseline_mean)), abs(diff(pr$baseline_mean)) / sd_b0,
              abs(diff(pr$slope_mean)), abs(diff(pr$slope_mean)) / sd_ols))
  res$best <- FALSE
  okix <- which(res$ok)
  res$best[okix[which.max(res$loglik[okix])]] <- TRUE
  saveRDS(list(fit = best2, profile = pr, knee = kk), file.path(BASE, "cth_lcmm_k2_fit.rds"))
} else {
  cat("\n  no admissible K = 2 solution from any of the three starts\n")
}
write.csv(res, file.path(BASE, "cth_lcmm_k2_rescue.csv"), row.names = FALSE)
saveRDS(list(k3 = k3, m3 = m3, profile3 = prof), file.path(BASE, "cth_lcmm_k3_diag.rds"))
cat("\n-> cth_lcmm_k3_diagnostic.csv / cth_lcmm_k2_rescue.csv / cth_lcmm_k3_diag.rds\n")
