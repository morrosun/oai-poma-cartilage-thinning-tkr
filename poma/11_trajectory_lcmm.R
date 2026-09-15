# ============================================================================
# 11 : trajectory classification of pre-operative cartilage thinning (LCMM)
#
#   Question: does the POMA cohort contain a "rapid thinning" trajectory class,
#   and is membership of that class associated with subsequent knee replacement?
#
#   Model   : hlme() latent-class mixed model, class-specific intercept + slope,
#             random intercept per knee.  Pure group-based trajectory models
#             (random = ~ -1) were also fitted but their BIC never turns over
#             (K = 1..5: 4839 -> 2804), i.e. the no-random-effect assumption is
#             untenable here; the LCMM BIC has a clear minimum at K = 4.
#
#   Two-stage logic is unbiased by construction: the class model uses ONLY the
#   longitudinal trajectory and never the outcome.
#
#   Complete separation: the fastest-thinning class reaches ~96% KR, so the
#   exact conditional likelihood diverges.  Per-class effects are therefore
#   estimated with a Bayesian within-pair conditional likelihood under a
#   weakly informative N(0, 2^2) prior on the log-OR, which returns finite
#   intervals and is the standard remedy for separation.
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
source("D:/BaiduSyncdisk/OAI/Scripts/poma/jm_core.R")
suppressPackageStartupMessages({library(lcmm); library(survival); library(ggplot2)})

PRIM <- "cMFTC_ThCtAB_aMe"
L <- read.csv(file.path(BASE, "poma_long.csv"))
U <- read.csv(file.path(BASE, "poma_analysis_long.csv"))

## --------------------------------------------------------------- data -----
sub <- L[!is.na(L[[PRIM]]), ]
d0 <- data.frame(uid = sub$unit, t = sub$years, y = sub[[PRIM]],
                 case = sub$case, strata = sub$newstrata)
nv <- table(d0$uid); keep_uid <- names(nv)[nv >= 2]
d <- d0[d0$uid %in% keep_uid, ]
d$unit <- as.integer(factor(d$uid)); d <- d[order(d$unit, d$t), ]
N <- max(d$unit)
uid_of_unit <- d$uid[match(seq_len(N), d$unit)]
uU <- U[match(uid_of_unit, U$unit), , drop = FALSE]
stopifnot(nrow(uU) == N)
## first (baseline) observed value of the metric, used as a level covariate
b0 <- vapply(seq_len(N), function(i) { x <- d[d$unit == i, ]; x$y[which.min(x$t)] }, numeric(1))
d$b0 <- b0[d$unit]
ols_slope <- vapply(seq_len(N), function(i) {
  x <- d[d$unit == i, ]; unname(coef(lm(y ~ t, x))[2]) }, numeric(1))

cat(sprintf("metric %s : %d knee-visits, %d knees, %d with >= 2 visits\n",
            PRIM, nrow(d0), length(nv), N))

entropy <- function(pp) { p <- as.matrix(pp[, grep("^prob", names(pp))])
  1 - sum(-p * log(pmax(p, 1e-12))) / (nrow(p) * log(ncol(p))) }
mkB <- function(K, rint, npar_fixed = 2) {
  ints <- seq(4.0, 2.4, length.out = K); sls <- seq(-0.03, -0.22, length.out = K)
  c(as.numeric(rbind(ints, sls)), rep(0, K - 1), 0.6, if (rint) 1.0)
}
fitK <- function(K, rint, nstart = 8) {
  best <- NULL; bl <- -Inf
  for (s in seq_len(nstart)) {
    set.seed(700 * K + s)
    mm <- try(hlme(y ~ t, mixture = ~ t, random = if (rint) ~1 else ~ -1,
                   subject = "unit", ng = K, data = d, B = mkB(K, rint)), silent = TRUE)
    if (!inherits(mm, "try-error") && is.finite(mm$loglik) && mm$loglik > bl) { bl <- mm$loglik; best <- mm }
  }
  best
}

## ------------------------------------------------- model selection ---------
CACHE <- file.path(BASE, "lcmm_fits_cache.rds")
sel <- data.frame(); fits <- list()
if (file.exists(CACHE)) {
  cc <- readRDS(CACHE)
  if (identical(cc$N, N) && identical(cc$PRIM, PRIM)) {
    sel <- cc$sel; fits <- cc$fits
    cat("\n[model selection restored from cache]\n")
  }
}
if (length(fits) == 0) {
cat("\n=========================== model selection ===========================\n")
cat(sprintf("%-6s %-4s %5s %10s %10s %10s %8s %7s  %s\n",
            "model", "K", "npar", "loglik", "AIC", "BIC", "entropy", "min%", "class sizes"))
for (rint in c(TRUE, FALSE)) {
  tag <- if (rint) "LCMM" else "GBTM"
  for (K in 1:5) {
    m <- if (K == 1) { set.seed(1)
      hlme(y ~ t, random = if (rint) ~1 else ~ -1, subject = "unit", ng = 1, data = d)
    } else fitK(K, rint)
    if (is.null(m)) { cat(sprintf("%-6s %-4d  FAILED\n", tag, K)); next }
    fits[[paste(tag, K)]] <- m
    cs <- table(m$pprob$class); en <- if (K > 1) entropy(m$pprob) else NA
    sel <- rbind(sel, data.frame(model = tag, K = K, npar = length(coef(m)), loglik = m$loglik,
                                 AIC = m$AIC, BIC = m$BIC, entropy = en,
                                 min_pct = 100 * min(cs) / sum(cs), sizes = paste(cs, collapse = "/")))
    cat(sprintf("%-6s %-4d %5d %10.1f %10.1f %10.1f %8.3f %7.1f  %s\n", tag, K,
                length(coef(m)), m$loglik, m$AIC, m$BIC, ifelse(is.na(en), NaN, en),
                100 * min(cs) / sum(cs), paste(cs, collapse = "/")))
  }
}
saveRDS(list(sel = sel, fits = fits, N = N, PRIM = PRIM), CACHE)
}
write.csv(sel, file.path(BASE, "lcmm_model_selection.csv"), row.names = FALSE)
lcm <- sel[sel$model == "LCMM", ]; Kbest <- lcm$K[which.min(lcm$BIC)]
mB <- fits[[paste("LCMM", Kbest)]]
cat(sprintf("\nBIC-selected LCMM : K = %d (BIC %.1f, entropy %.3f, smallest class %.1f%%)\n",
            Kbest, mB$BIC, entropy(mB$pprob), 100 * min(table(mB$pprob$class)) / nrow(mB$pprob)))

## -------------------------------------------------------- class definition
pp <- mB$pprob; stopifnot(identical(as.integer(pp[[1]]), seq_len(N)))
P <- as.matrix(pp[, grep("^prob", names(pp))]); colnames(P) <- paste0("P", seq_len(ncol(P)))
hard <- max.col(P, ties.method = "first")
Kk <- ncol(P)

knee <- data.frame(unit = seq_len(N), uid = uid_of_unit, cl = hard,
                   pmax = P[cbind(seq_len(N), hard)], ols_slope = ols_slope, b0 = b0,
                   eblup = uU[[paste0(PRIM, "__eblup")]], case = uU$case,
                   strata = uU$newstrata, kl = uU$V00XRKL, index_yr = uU$t_index_months / 12,
                   nvis = as.integer(table(d$unit)))
mslope <- tapply(knee$ols_slope, knee$cl, mean)
ord <- as.integer(names(sort(mslope, decreasing = TRUE)))   # slowest -> fastest

cat("\n---- trajectory class profile (slowest to fastest thinning) ----\n")
prof <- do.call(rbind, lapply(ord, function(g) {
  x <- knee[knee$cl == g, ]
  data.frame(class = g, n = nrow(x), pct = 100 * nrow(x) / nrow(knee),
             slope_mm_yr = mean(x$ols_slope), slope_sd = sd(x$ols_slope),
             baseline_mm = mean(x$b0), eblup_slope = mean(x$eblup),
             pct_KR = 100 * mean(x$case), mean_KL = mean(x$kl),
             visits = mean(x$nvis), index_yr = mean(x$index_yr), mean_pmax = mean(x$pmax))
}))
print(prof, row.names = FALSE, digits = 4)
write.csv(prof, file.path(BASE, "lcmm_class_profile.csv"), row.names = FALSE)
cat(sprintf("\n  Spearman (class thinning rate vs class %%KR) = %+.2f\n",
            cor(prof$slope_mm_yr, prof$pct_KR, method = "spearman")))
cat(sprintf("  Spearman (class baseline thickness vs class %%KR) = %+.2f\n",
            cor(prof$baseline_mm, prof$pct_KR, method = "spearman")))

## --------------------------------------- Bayesian within-pair cond. logit --
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
  Z <- as.matrix(dat[unlist(lapply(ii, `[`, 1)), cols, drop = FALSE]) -
       as.matrix(dat[unlist(lapply(ii, `[`, 2)), cols, drop = FALSE])
  Z
}

an <- knee[!is.na(knee$strata), ]
gt <- table(an$strata); an <- an[an$strata %in% names(gt)[gt == 2], ]
an$clfac <- factor(an$cl, levels = ord)          # ord[1] = reference = slowest
REF <- ord[1]; an$rapid <- an$cl %in% ord[(Kk - 1):Kk]
cat(sprintf("\n---- outcome link : %d complete matched pairs / %d knees ----\n",
            length(unique(an$strata)), nrow(an)))
cat(sprintf("reference class = %d (slowest, slope %+.4f mm/yr, KR %.1f%%)\n",
            REF, mslope[as.character(REF)], 100 * mean(knee$case[knee$cl == REF])))
cat(sprintf("'rapid thinning' class(es) = %s (n = %d; KR %.1f%%)\n",
            paste(ord[(Kk - 1):Kk], collapse = " + "), sum(knee$cl %in% ord[(Kk - 1):Kk]),
            100 * mean(knee$case[knee$cl %in% ord[(Kk - 1):Kk]])))

## per-class contrasts, Bayesian (separation-safe)
dmy <- model.matrix(~ clfac - 1, data = an)[, -1, drop = FALSE]
colnames(dmy) <- paste0("cv", seq_len(ncol(dmy)))
an <- cbind(an, dmy)
Xp <- pair_contrast(an, colnames(dmy))
post <- bclogit(Xp, seed = 11)
pc <- data.frame(contrast = paste0("class ", ord[-1], " vs ", REF),
                 class = ord[-1],
                 slope_mm_yr = mslope[as.character(ord[-1])],
                 pct_KR_all = sapply(ord[-1], function(g) 100 * mean(knee$case[knee$cl == g])),
                 OR = exp(apply(post, 2, median)),
                 lo = exp(apply(post, 2, quantile, .025)),
                 hi = exp(apply(post, 2, quantile, .975)),
                 P_OR_gt1 = colMeans(post > 0))
cat("\n[per-class OR, Bayesian within-pair conditional logistic, prior N(0,2^2)]\n")
print(cbind(pc[, c("contrast", "class")], round(pc[, !names(pc) %in% c("contrast", "class")], 4)),
      row.names = FALSE)
write.csv(pc, file.path(BASE, "lcmm_class_OR.csv"), row.names = FALSE)

## binary 'rapid thinning' indicator and ordinal trend, hard classification
cgr <- clogit(case ~ rapid + strata(strata), data = an, method = "exact")
sr <- summary(cgr)$coef
cat(sprintf("\n[binary rapid-thinning indicator, exact conditional logistic]\n  OR = %.3f (%.3f-%.3f), p = %.4g\n",
            exp(sr[1, 1]), exp(sr[1, 1] - 1.96 * sr[1, 3]), exp(sr[1, 1] + 1.96 * sr[1, 3]), sr[1, 5]))
an$rank <- as.numeric(factor(an$cl, levels = ord))
cgt <- clogit(case ~ rank + strata(strata), data = an, method = "exact")
st <- summary(cgt)$coef
cat(sprintf("[ordinal trend across classes] OR per one class step = %.3f (%.3f-%.3f), p = %.4g\n",
            exp(st[1, 1]), exp(st[1, 1] - 1.96 * st[1, 3]), exp(st[1, 1] + 1.96 * st[1, 3]), st[1, 5]))
cat(sprintf("  implied OR, fastest vs slowest class = %.1f\n", exp(st[1, 1])^(Kk - 1)))
cat(sprintf("  Kendall tau (class rank vs thinning rate) = %+.3f\n",
            cor(an$rank, -an$eblup, method = "kendall")))

## -------------------- MI over the posterior class probabilities ------------
mid <- matrix(NA_real_, 20, 2); miv <- matrix(NA_real_, 20, 2)
for (r in 1:20) {
  set.seed(9000 + r)
  Pm <- P[match(an$unit, seq_len(N)), , drop = FALSE]; Pm <- Pm / rowSums(Pm)
  dr <- vapply(seq_len(nrow(Pm)), function(i) sample(seq_len(Kk), 1, prob = Pm[i, ]), numeric(1))
  a2 <- an
  a2$clfac <- factor(dr, levels = ord); a2$rapid <- dr %in% ord[(Kk - 1):Kk]
  a2$rank <- as.numeric(factor(dr, levels = ord))
  f1 <- try(clogit(case ~ rapid + strata(strata), data = a2, method = "exact"), silent = TRUE)
  f2 <- try(clogit(case ~ rank  + strata(strata), data = a2, method = "exact"), silent = TRUE)
  if (!inherits(f1, "try-error")) { c1 <- summary(f1)$coef; mid[r, 1] <- c1[1, 1]; miv[r, 1] <- c1[1, 3]^2 }
  if (!inherits(f2, "try-error")) { c2 <- summary(f2)$coef; mid[r, 2] <- c2[1, 1]; miv[r, 2] <- c2[1, 3]^2 }
}
miout <- do.call(rbind, lapply(1:2, function(j) {
  ok <- !is.na(mid[, j])
  Q <- mid[ok, j]; Uu <- miv[ok, j]; r <- sum(ok)
  Bv <- var(Q); Tv <- mean(Uu) + (1 + 1/r) * Bv
  df <- (r - 1) * (1 + mean(Uu) / ((1 + 1/r) * Bv))^2; tc <- qt(.975, df)
  qb <- mean(Q)
  data.frame(term = c("rapid class vs rest", "OR per class step")[j],
             OR = exp(qb), lo = exp(qb - tc * sqrt(Tv)), hi = exp(qb + tc * sqrt(Tv)),
             p = 2 * pt(-abs(qb / sqrt(Tv)), df), m_used = r)
}))
cat("\n[multiple imputation over posterior class probabilities, m = 20, Rubin]\n")
print(cbind(miout[, "term", drop = FALSE], round(miout[, !names(miout) %in% "term"], 4)),
      row.names = FALSE)
write.csv(miout, file.path(BASE, "lcmm_class_OR_mi.csv"), row.names = FALSE)

## checkpoint: keep the main results even if a sensitivity below fails
saveRDS(list(selection = sel, K = Kbest, model = mB, fits = fits, knee = knee,
             profile = prof, per_class = pc, mi = miout, rapid_hard = sr,
             trend_hard = st, bclogit_post = post, pair_data = an),
        file.path(BASE, "lcmm_fit.rds"))

## ---- sensitivity 1 : K = 3 ----------------------------------------------
cat("\n---- sensitivity : K = 3 ----\n")
try({
if (!is.null(fits[["LCMM 3"]])) {
  p3 <- fits[["LCMM 3"]]$pprob
  P3 <- as.matrix(p3[, grep("^prob", names(p3))])
  k3 <- data.frame(unit = seq_len(N), cl = max.col(P3, ties.method = "first"),
                   ols_slope = ols_slope, eblup = uU[[paste0(PRIM, "__eblup")]],
                   case = uU$case, strata = uU$newstrata)
  q3 <- tapply(k3$ols_slope, k3$cl, mean); o3 <- as.integer(names(sort(q3, decreasing = TRUE)))
  a3 <- k3[!is.na(k3$strata), ]; g3 <- table(a3$strata)
  a3 <- a3[a3$strata %in% names(g3)[g3 == 2], ]
  a3$clfac <- factor(a3$cl, levels = o3)
  print(round(summary(clogit(case ~ clfac + strata(strata), data = a3, method = "exact"))$conf.int, 4))
  for (g in o3) cat(sprintf("  class %d : n = %3d  slope %+.4f mm/yr  baseline %.2f mm  KR %.1f%%\n",
                            g, sum(k3$cl == g), q3[as.character(g)],
                            mean(b0[k3$cl == g]), 100 * mean(k3$case[k3$cl == g])))
  a3$rank <- as.numeric(factor(a3$cl, levels = o3))
  s3 <- summary(clogit(case ~ rank + strata(strata), data = a3, method = "exact"))$coef
  cat(sprintf("  ordinal trend OR per class step = %.3f (%.3f-%.3f), p = %.4g\n",
              exp(s3[1, 1]), exp(s3[1, 1] - 1.96 * s3[1, 3]), exp(s3[1, 1] + 1.96 * s3[1, 3]), s3[1, 5]))
}
})

## ---- sensitivity 2 : is the class trend just baseline thickness? ---------
## The classes differ in BASELINE thickness as well as in rate (POMA is enriched
## for knees with denuded cartilage).  Adjust the ordinal class trend for the
## within-pair difference in baseline thickness.  This is the direct test of
## "classes are level artefacts".
cat("\n---- sensitivity : class trend adjusted for baseline thickness ----\n")
try({
  an$zb0 <- (an$b0 - mean(an$b0)) / sd(an$b0)
  s2a <- summary(clogit(case ~ rank + zb0 + strata(strata), data = an,
                        method = "exact"))$coef
  s2b <- summary(clogit(case ~ zb0 + strata(strata), data = an,
                        method = "exact"))$coef
  cat(sprintf("  baseline thickness alone        OR/SD = %.3f (%.3f-%.3f), p = %.4g\n",
              exp(s2b[1, 1]), exp(s2b[1, 1] - 1.96 * s2b[1, 3]),
              exp(s2b[1, 1] + 1.96 * s2b[1, 3]), s2b[1, 5]))
  cat(sprintf("  class trend + baseline thickness OR/step = %.3f (%.3f-%.3f), p = %.4g\n",
              exp(s2a[1, 1]), exp(s2a[1, 1] - 1.96 * s2a[1, 3]),
              exp(s2a[1, 1] + 1.96 * s2a[1, 3]), s2a[1, 5]))
  cat(sprintf("  baseline thickness in that model        = %.3f (%.3f-%.3f), p = %.4g\n",
              exp(s2a[2, 1]), exp(s2a[2, 1] - 1.96 * s2a[2, 3]),
              exp(s2a[2, 1] + 1.96 * s2a[2, 3]), s2a[2, 5]))
  cat(sprintf("  Spearman (class slope vs class %%KR) = %+.2f ; (class baseline vs class %%KR) = %+.2f\n",
              cor(prof$slope_mm_yr, prof$pct_KR, method = "spearman"),
              cor(prof$baseline_mm, prof$pct_KR, method = "spearman")))
})

## ---- sensitivity 3 : time aligned on the index visit ---------------------
cat("\n---- sensitivity : time aligned on the index visit ----\n")
try({
d$trel <- d$t - (U$t_index_months[match(d$uid, U$unit)] / 12)
mrel <- local({
  best <- NULL; bl <- -Inf
  for (s in 1:8) { set.seed(4242 + s)
    mm <- try(hlme(y ~ trel, mixture = ~ trel, random = ~1, subject = "unit",
                   ng = Kbest, data = d, B = mkB(Kbest, TRUE)), silent = TRUE)
    if (!inherits(mm, "try-error") && is.finite(mm$loglik) && mm$loglik > bl) { bl <- mm$loglik; best <- mm } }
  best })
if (is.null(mrel)) cat("  failed\n") else {
  Pr <- as.matrix(mrel$pprob[, grep("^prob", names(mrel$pprob))])
  kr <- data.frame(unit = seq_len(N), cl = max.col(Pr, ties.method = "first"),
                   ols_slope = ols_slope, eblup = uU[[paste0(PRIM, "__eblup")]],
                   case = uU$case, strata = uU$newstrata)
  qr <- tapply(kr$ols_slope, kr$cl, mean); orr <- as.integer(names(sort(qr, decreasing = TRUE)))
  ar <- kr[!is.na(kr$strata), ]; gr <- table(ar$strata); ar <- ar[ar$strata %in% names(gr)[gr == 2], ]
  ar$rank <- as.numeric(factor(ar$cl, levels = orr))
  sr2 <- summary(clogit(case ~ rank + strata(strata), data = ar, method = "exact"))$coef
  cat(sprintf("  BIC %.1f (calendar axis %.1f); class sizes %s\n", mrel$BIC, mB$BIC,
              paste(table(kr$cl), collapse = "/")))
  cat(sprintf("  ordinal trend OR per class step = %.3f (%.3f-%.3f), p = %.4g\n",
              exp(sr2[1, 1]), exp(sr2[1, 1] - 1.96 * sr2[1, 3]), exp(sr2[1, 1] + 1.96 * sr2[1, 3]), sr2[1, 5]))
  for (g in orr) cat(sprintf("  class %d : n = %3d  slope %+.4f mm/yr  KR %.1f%%\n",
                             g, sum(kr$cl == g), qr[as.character(g)], 100 * mean(kr$case[kr$cl == g])))
}
})

## ------------------------------------------------------------- figure -----
labv <- sprintf("Class %d\n%+.3f mm/yr | %.1f mm | %.0f%% KR", ord,
                mslope[as.character(ord)], prof$baseline_mm, prof$pct_KR)
knee$clf <- factor(knee$cl, levels = ord, labels = labv)
dd <- merge(d[, c("unit", "t", "y")], knee[, c("unit", "clf")], by = "unit")
pm <- do.call(rbind, lapply(split(dd, dd$clf), function(x) {
  f <- lm(y ~ t, data = x); tt <- seq(0, 4, by = .5)
  data.frame(clf = x$clf[1], t = tt, y = predict(f, newdata = data.frame(t = tt))) }))
gp1 <- ggplot(dd, aes(t, y, group = unit)) +
  geom_line(colour = "grey78", linewidth = .22, alpha = .75) +
  geom_line(data = pm, aes(t, y, group = clf), colour = "#B03A2E", linewidth = 1.4, inherit.aes = FALSE) +
  facet_wrap(~ clf, nrow = 1) + theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 8), panel.grid.minor = element_blank()) +
  labs(x = "years from baseline", y = "cMFTC mean thickness (mm)",
       title = "Pre-operative cMFTC thinning trajectories: latent classes",
       subtitle = sprintf("Latent-class mixed model with random intercept, K = %d (BIC-selected); red line = within-class fitted mean", Kbest))
ggsave(file.path(BASE, "fig_lcmm_trajectories.png"), gp1, width = 11, height = 3.9, dpi = 200)

bt <- data.frame(clf = levels(knee$clf), n = as.integer(table(knee$clf)),
                 pct = 100 * as.numeric(tapply(knee$case, knee$clf, mean)),
                 slope = mslope[as.character(ord)])
gp2 <- ggplot(bt, aes(clf, pct)) +
  geom_col(fill = "#2E6DA4", width = .6) +
  geom_text(aes(label = sprintf("%.1f%% (n=%d)", pct, n)), vjust = -0.4, size = 3) +
  theme_bw(base_size = 10) + ylim(0, max(bt$pct) * 1.2) +
  theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "% undergoing knee replacement",
       title = "Knee replacement rate rises monotonically with the class thinning rate")
ggsave(file.path(BASE, "fig_lcmm_outcome.png"), gp2, width = 7, height = 3.4, dpi = 200)

saveRDS(list(selection = sel, K = Kbest, model = mB, fits = fits, knee = knee,
             profile = prof, per_class = pc, mi = miout,
             rapid_hard = sr, trend_hard = st, bclogit_post = post,
             pair_data = an), file.path(BASE, "lcmm_fit.rds"))
cat("\n-> lcmm_model_selection.csv / lcmm_class_profile.csv / lcmm_class_OR.csv")
cat(" / lcmm_class_OR_mi.csv / fig_lcmm_trajectories.png / fig_lcmm_outcome.png / lcmm_fit.rds\n")
