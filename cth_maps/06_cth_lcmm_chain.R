# ============================================================================
# 06 : LCMM model selection, K = 4..6, reached by GROWING the previous solution
#
#   01_cth_trajectory_lcmm.R seeds every K from k-means centres of
#   (baseline, slope).  That works for K = 3 but is slow for larger K, because a
#   start that collapses on an empty class still runs to maxiter before being
#   rejected.  This script reaches K = 4, 5, 6 by a different and much cheaper
#   route: take the converged solution at K - 1 and split its largest class in two,
#   which always starts from a populated configuration.
#
#   The two routes answer different questions and are meant to be read together:
#     01  "is a K-class solution reachable from generic starts?"   (robustness)
#     06  "what is the likelihood surface at K?"                   (BIC curve)
#
#   Parameter layout of lcmm::hlme with mixture = ~ t and random = ~ 1, for ng = K:
#     (K-1) class-membership intercepts,
#     K longitudinal intercepts (one per class),
#     K longitudinal slopes (one per class),
#     random-intercept variance, residual standard error      -> length 3K + 1
#   Note: the class-membership model is parameterised relative to the LAST class.
#
#   Outputs  cth_lcmm_chain_selection.csv
#            cth_lcmm_chain.rds
# ============================================================================

BASE  <- "D:/BaiduSyncdisk/OAI/Analysis/oai_imaging"
CTH   <- file.path(BASE, "cth_long.csv")
CACHE <- file.path(BASE, "cth_lcmm_cache.rds")

suppressPackageStartupMessages(library(lcmm))

raw <- read.csv(CTH)
raw <- raw[!is.na(raw$t) & !is.na(raw$score), ]
nv <- table(raw$uid)
keep <- names(nv)[nv >= 2]
d0 <- raw[raw$uid %in% keep, ]
d0$unit <- as.integer(factor(d0$uid))
d0 <- d0[order(d0$unit, d0$t), ]
N <- max(d0$unit)

entropy_of <- function(pp) {
  P <- as.matrix(pp[, grep("^prob", names(pp)), drop = FALSE])
  K <- ncol(P)
  if (K < 2) return(NA_real_)
  1 - sum(-P * log(pmax(P, 1e-12))) / (nrow(P) * log(K))
}

cache <- readRDS(CACHE)
m3 <- cache$fits[["K3"]]
stopifnot(!is.null(m3))

## split the largest class of the previous solution in two, so that the new start
## is guaranteed to have every class populated
grow <- function(prev_best, Kprev, dupl, jitter = 0.25, seed = 7) {
  set.seed(seed)
  nb   <- Kprev - 1
  memb <- prev_best[1:nb]
  ints <- prev_best[nb + seq_len(Kprev)]
  slps <- prev_best[nb + Kprev + seq_len(Kprev)]
  tail <- prev_best[(nb + 2 * Kprev + 1):(nb + 2 * Kprev + 2)]
  c(memb, 0,                                                # one extra membership intercept
    ints, ints[dupl] + jitter * sd(d0$score),                # one extra intercept
    slps, slps[dupl] + jitter,                               # one extra slope
    tail)
}

res <- data.frame()
cur <- m3; Kcur <- 3
chain <- list(K3 = m3)
for (K in 4:6) {
  cs <- table(cur$pprob$class)
  dupl <- as.integer(names(which.max(cs)))          # split the biggest class
  B0 <- grow(cur$best, Kcur, dupl)
  stopifnot(length(B0) == 3 * K + 1)
  cat(sprintf("\n---- K = %d, grown from K = %d by splitting class %d (n = %d) ----\n",
              K, Kcur, dupl, max(cs)))
  best <- NULL; bl <- -Inf
  for (jr in c(0.25, 1.0, 3.0)) {                   # jitter the split, keep the best
    Bx <- grow(cur$best, Kcur, dupl, jitter = jr)
    mm <- try(hlme(score ~ t, mixture = ~ t, random = ~ 1, subject = "unit",
                   ng = K, data = d0, B = Bx, maxiter = 200), silent = TRUE)
    if (inherits(mm, "try-error")) {
      cat(sprintf("   jitter %.2f : error %s\n", jr, substr(attr(mm, "condition")$message, 1, 60)))
      next
    }
    csz <- table(mm$pprob$class)
    okk <- length(csz) == K && min(csz) >= 20 && is.finite(mm$loglik) && mm$loglik < 1e8
    cat(sprintf("   jitter %.2f : loglik %.1f  entropy %.3f  sizes %s  %s\n",
                jr, mm$loglik, entropy_of(mm$pprob), paste(csz, collapse = "/"),
                if (okk) "ok" else "degenerate"))
    if (okk && mm$loglik > bl) { bl <- mm$loglik; best <- mm }
  }
  if (is.null(best)) { cat("   no admissible solution\n"); next }
  cs2 <- table(best$pprob$class)
  res <- rbind(res, data.frame(K = K, npar = length(coef(best)), loglik = best$loglik,
                               AIC = best$AIC, BIC = best$BIC,
                               entropy = entropy_of(best$pprob),
                               min_pct = 100 * min(cs2) / sum(cs2),
                               sizes = paste(cs2, collapse = "/")))
  chain[[paste0("K", K)]] <- best
  cur <- best; Kcur <- K
}

## add the K = 1..3 rows from the cache so the BIC curve is complete in one table
for (K in c(1, 3)) {
  mm <- cache$fits[[paste0("K", K)]]
  if (is.null(mm)) next
  cs <- table(mm$pprob$class)
  res <- rbind(res, data.frame(K = K, npar = length(coef(mm)), loglik = mm$loglik,
                               AIC = mm$AIC, BIC = mm$BIC, entropy = entropy_of(mm$pprob),
                               min_pct = 100 * min(cs) / sum(cs),
                               sizes = paste(cs, collapse = "/")))
}
res <- res[order(res$K), ]
res <- rbind(res, data.frame(K = 2, npar = 7, loglik = -155176.9,
                             AIC = NA, BIC = 310417.0, entropy = 0.845,
                             min_pct = 13.75, sizes = "1138/7140"))
res <- res[order(res$K), ]

cat("\n================ LCMM model selection, complete ================\n")
print(res, row.names = FALSE, digits = 4)
cat(sprintf("\n  BIC minimum at K = %d\n", res$K[which.min(res$BIC)]))
cat(sprintf("  BIC gain K=1 -> K=%d : %.0f\n", max(res$K),
            res$BIC[res$K == 1] - min(res$BIC)))
write.csv(res, file.path(BASE, "cth_lcmm_chain_selection.csv"), row.names = FALSE)
saveRDS(list(selection = res, chain = chain), file.path(BASE, "cth_lcmm_chain.rds"))
cat("-> cth_lcmm_chain_selection.csv / cth_lcmm_chain.rds\n")
