# ============================================================================
# 01 : trajectory classification of CTh-Score (LCMM)  —— P2
#
#   Data : OAI CTh-Maps / CTh-Score dataset (Zenodo 10.5281/zenodo.18745638)
#          4,340 subjects, 8,581 knees, 45,343 knee-visits, 7 timepoints
#          (0/12/24/36/48/72/96 months).  CTh-Score = 0-100 cartilage
#          osteoarthritis severity; HIGHER = MORE SEVERE.
#
#   Question: do CTh-Score trajectories over 8 years resolve into discrete
#             progression classes, and is the per-knee linear slope a faithful
#             summary of "how fast a knee gets worse"?
#
#   Model   : hlme() latent-class mixed model, class-specific intercept + slope,
#             random intercept per knee.  Same specification as
#             11_trajectory_lcmm.R used for the POMA manuscript, so the two
#             analyses are directly comparable.
#
#   Note    : this dataset carries NO clinical covariates.  Baseline KL is
#             available for a 431-knee subset through the OAIZIB-CM release
#             (right knees only) and is used here purely as an external
#             criterion, not as a covariate.
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/oai_imaging"
CTH  <- "D:/BaiduSyncdisk/OAI/Analysis/oai_imaging/cth_long.csv"
KLF  <- "D:/BaiduSyncdisk/OAI/Analysis/oai_imaging/zib_klinfo.csv"
dir.create(BASE, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({library(lcmm); library(ggplot2)})

## ------------------------------------------------------------------ data ---
## cth_long.csv : ID, SIDE, Timepoint, t (years), score (0-100), uid = ID_SIDE
raw <- read.csv(CTH)
raw <- raw[!is.na(raw$t) & !is.na(raw$score), ]

nv <- table(raw$uid)
keep <- names(nv)[nv >= 2]
d0 <- raw[raw$uid %in% keep, ]
d0$unit <- as.integer(factor(d0$uid))
d0 <- d0[order(d0$unit, d0$t), ]
N <- max(d0$unit)
uid_of_unit <- d0$uid[match(seq_len(N), d0$unit)]

cat(sprintf("CTh-Score : %d knee-visits, %d knees total, %d knees with >= 2 visits\n",
            nrow(raw), length(nv), N))
cat(sprintf("follow-up length (years): median %.1f, IQR %.1f-%.1f, max %.1f\n",
            median(d0$t), quantile(d0$t, .25), quantile(d0$t, .75), max(d0$t)))
cat(sprintf("visits per knee: mean %.2f, distribution %s\n",
            nrow(d0) / N, paste(names(table(table(d0$unit))),
                                table(table(d0$unit)), sep = ":", collapse = " / ")))

## baseline (first observed) score and unshrunk per-knee OLS slope
b0 <- vapply(seq_len(N), function(i) { x <- d0[d0$unit == i, ]; x$score[which.min(x$t)] }, numeric(1))
ols <- vapply(seq_len(N), function(i) { x <- d0[d0$unit == i, ]
  if (length(unique(x$t)) < 2) NA_real_ else unname(coef(lm(score ~ t, x))[2]) }, numeric(1))

## ------------------------------------------------------- model selection ---
entropy <- function(pp) { P <- as.matrix(pp[, grep("^prob", names(pp))])
  1 - sum(-P * log(pmax(P, 1e-12))) / (nrow(P) * log(ncol(P))) }

## ---- starting values -------------------------------------------------------
## hlme's B vector order for `random = ~ 1` is
##   fixed intercept, fixed t,
##   (mixture intercept, mixture t) x (K-1),
##   K-1 class proportions, residual variance, random-intercept variance
## i.e. length 3K + 1.
##
## Grid / heuristic starting values are NOT usable here.  CTh-Score has a floor
## at 0 and its between-knee level spread (random-intercept SD ~ 23) dwarfs the
## within-knee residual (SD ~ 4-6), so any hand-picked spread collapses the
## solution onto one class.  Seeding the class means with k-means centres of
## (baseline score, per-knee OLS slope) guarantees every starting class is
## populated, and turns the fit from divergent into a 3-5 s convergence.
m1 <- hlme(score ~ t, random = ~ 1, subject = "unit", ng = 1, data = d0, maxiter = 200)
cat(sprintf("\nK = 1 reference : intercept %.2f, slope %.3f, RE var %.1f, resid SD %.3f\n",
            m1$best[1], m1$best[2], m1$best[3], m1$best[4]))

## A single k-means seed is NOT reliable : for some K the k-means centres happen
## to sit in a configuration from which hlme stalls (non-finite loglik, or a
## solution with an empty class).  Each K is therefore retried over a small set
## of k-means seeds and the best converged, non-degenerate fit is kept.
mkB <- function(K, seed = 11) {
  set.seed(seed)
  km <- kmeans(cbind(b0, ols), centers = K, nstart = 50)
  ord <- order(km$centers[, 2])
  ints <- km$centers[ord, 1]; sls <- km$centers[ord, 2]
  c(as.numeric(rbind(ints, sls)), rep(0, K - 1), m1$best[3], m1$best[4])
}

LAST_ERR <- NULL
fitK <- function(K, seeds = c(11, 29, 53, 97, 151, 211)) {
  best <- NULL; bl <- -Inf; msgs <- character(0)
  for (sd_ in seeds) {
    mm <- try(hlme(score ~ t, mixture = ~ t, random = ~ 1, subject = "unit",
                   ng = K, data = d0, B = mkB(K, seed = sd_),
                   maxiter = 200), silent = TRUE)
    if (inherits(mm, "try-error")) {
      msgs <- c(msgs, sprintf("seed %d: %s", sd_, attr(mm, "condition")$message))
      next
    }
    if (!is.finite(mm$loglik) || mm$loglik > 1e8) {
      msgs <- c(msgs, sprintf("seed %d: non-finite loglik", sd_)); next
    }
    if (length(table(mm$pprob$class)) < K) {
      msgs <- c(msgs, sprintf("seed %d: empty class", sd_)); next
    }
    if (mm$loglik > bl) { bl <- mm$loglik; best <- mm }
  }
  if (is.null(best)) {
    LAST_ERR <<- paste(unique(msgs), collapse = " | ")
    return(NULL)
  }
  LAST_ERR <<- NULL
  cat(sprintf("      (K=%d converged; %d candidate starts tried, best loglik %.1f)\n",
              K, length(seeds), bl))
  best
}

CACHE <- file.path(BASE, "cth_lcmm_cache.rds")
sel <- data.frame(); fits <- list()
if (file.exists(CACHE)) {
  cc <- readRDS(CACHE)
  if (identical(cc$N, N)) { sel <- cc$sel; fits <- cc$fits
    cat(sprintf("\n[model selection restored from cache : K = %s]\n",
                paste(sub("^K", "", names(fits)), collapse = ", "))) }
}
cat("\n=========================== model selection ===========================\n")
cat(sprintf("%-6s %5s %10s %10s %10s %8s %7s  %s\n",
            "K", "npar", "loglik", "AIC", "BIC", "entropy", "min%", "class sizes"))
for (K in 1:6) {
  key <- paste0("K", K)
  if (!is.null(fits[[key]])) {
    cat(sprintf("%-6d  [cached]\n", K)); next
  }
  m <- if (K == 1) m1 else fitK(K)
  if (is.null(m)) {
    cat(sprintf("%-6d  FAILED  %s\n", K,
                substr(if (is.null(LAST_ERR)) "" else LAST_ERR, 1, 150))); next
  }
  fits[[key]] <- m
  cs <- table(m$pprob$class); en <- if (K > 1) entropy(m$pprob) else NA
  sel <- rbind(sel, data.frame(K = K, npar = length(coef(m)), loglik = m$loglik,
                               AIC = m$AIC, BIC = m$BIC, entropy = en,
                               min_pct = 100 * min(cs) / sum(cs),
                               sizes = paste(cs, collapse = "/")))
  cat(sprintf("%-6d %5d %10.1f %10.1f %10.1f %8.3f %7.1f  %s\n", K,
              length(coef(m)), m$loglik, m$AIC, m$BIC,
              ifelse(is.na(en), NaN, en), 100 * min(cs) / sum(cs),
              paste(cs, collapse = "/")))
  saveRDS(list(sel = sel, fits = fits, N = N), CACHE)
}
sel <- sel[order(sel$K), , drop = FALSE]
write.csv(sel, file.path(BASE, "cth_lcmm_model_selection.csv"), row.names = FALSE)
Kbest <- sel$K[which.min(sel$BIC)]
mB <- fits[[paste0("K", Kbest)]]
cat(sprintf("\nBIC-selected LCMM : K = %d (BIC %.1f, entropy %.3f, smallest class %.1f%%)\n",
            Kbest, mB$BIC, entropy(mB$pprob),
            100 * min(table(mB$pprob$class)) / nrow(mB$pprob)))

## ------------------------------------------------------ class definition ---
pp <- mB$pprob
P <- as.matrix(pp[, grep("^prob", names(pp))]); colnames(P) <- paste0("P", seq_len(ncol(P)))
hard <- max.col(P, ties.method = "first")
Kk <- ncol(P)
nvis <- as.integer(table(d0$unit))

knee <- data.frame(unit = seq_len(N), uid = uid_of_unit, cl = hard,
                   pmax = P[cbind(seq_len(N), hard)], ols = ols, b0 = b0, nvis = nvis)
mslope <- tapply(knee$ols, knee$cl, mean, na.rm = TRUE)
ord <- as.integer(names(sort(mslope, decreasing = FALSE)))   # slowest -> fastest

cat("\n---- trajectory class profile (slowest to fastest progression) ----\n")
prof <- do.call(rbind, lapply(ord, function(g) {
  x <- knee[knee$cl == g, ]
  data.frame(class = g, n = nrow(x), pct = 100 * nrow(x) / nrow(knee),
             slope_pts_yr = mean(x$ols, na.rm = TRUE), slope_sd = sd(x$ols, na.rm = TRUE),
             baseline = mean(x$b0), visits = mean(x$nvis), mean_pmax = mean(x$pmax))
}))
print(prof, row.names = FALSE, digits = 4)
write.csv(prof, file.path(BASE, "cth_lcmm_class_profile.csv"), row.names = FALSE)
cat(sprintf("  Spearman (class slope vs class baseline) = %+.2f\n",
            cor(prof$slope_pts_yr, prof$baseline, method = "spearman")))
cat(sprintf("  overall slope: mean %+.3f pts/yr, median %+.3f, sd %.3f, range %+.2f to %+.2f\n",
            mean(knee$ols, na.rm = TRUE), median(knee$ols, na.rm = TRUE),
            sd(knee$ols, na.rm = TRUE), min(knee$ols, na.rm = TRUE), max(knee$ols, na.rm = TRUE)))
cat(sprintf("  knees worsening (slope > 0) %.1f%%\n", 100 * mean(knee$ols > 0, na.rm = TRUE)))

## --------------------------------------- external criterion : baseline KL ---
## OAIZIB-CM carries KLGrade for 507 RIGHT knees only -> right-knee classes only
sub <- read.csv(KLF)                    # SubjectID, SIDE, KLGrade, Gender, Age, BMI
knee$ID <- as.integer(sub("_.*$", "", knee$uid)); knee$SIDE <- sub("^.*_", "", knee$uid)
kl <- merge(knee[knee$SIDE == "RIGHT", ], sub[, c("SubjectID", "KLGrade", "Age", "BMI")],
            by.x = "ID", by.y = "SubjectID")
if (nrow(kl) > 0) {
  kl <- kl[!is.na(kl$KLGrade), ]
  cat(sprintf("\n---- external criterion: baseline KLGrade (right knees, n = %d) ----\n", nrow(kl)))
  tab <- aggregate(cbind(mean_slope = ols, mean_base = b0) ~ cl + KLGrade, data = kl, FUN = mean)
  print(tab, row.names = FALSE, digits = 3)
  byc <- aggregate(cbind(KL = KLGrade) ~ cl, data = kl, FUN = mean)
  print(byc, row.names = FALSE, digits = 3)
  cat(sprintf("  Spearman (class mean KL vs class slope) = %+.2f\n",
              cor(byc$KL, mslope[as.character(byc$cl)], method = "spearman")))
  cat(sprintf("  Spearman (per-knee KL vs per-knee slope) = %+.3f\n",
              cor(kl$KLGrade, kl$ols, method = "spearman", use = "complete.obs")))
  write.csv(byc, file.path(BASE, "cth_lcmm_class_vs_KL.csv"), row.names = FALSE)
}

## --------------------------------------------------------------- figure ----
labv <- sprintf("Class %d\n%+.2f pts/yr | base %.0f | n=%d",
                ord, prof$slope_pts_yr, prof$baseline, prof$n)
knee$clf <- factor(knee$cl, levels = ord, labels = labv)
dd <- merge(d0[, c("unit", "t", "score")], knee[, c("unit", "clf")], by = "unit")
pm <- do.call(rbind, lapply(split(dd, dd$clf), function(x) {
  f <- lm(score ~ t, data = x); tt <- seq(0, 8, by = .25)
  data.frame(clf = x$clf[1], t = tt, score = predict(f, newdata = data.frame(t = tt))) }))
gp1 <- ggplot(dd, aes(t, score, group = unit)) +
  geom_line(colour = "grey80", linewidth = .18, alpha = .6) +
  geom_line(data = pm, aes(t, score, group = clf), colour = "#B03A2E",
            linewidth = 1.4, inherit.aes = FALSE) +
  facet_wrap(~ clf, nrow = 1) + theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 8), panel.grid.minor = element_blank()) +
  labs(x = "years from baseline", y = "CTh-Score (0-100, higher = worse)",
       title = "CTh-Score progression trajectories: latent classes",
       subtitle = sprintf("Latent-class mixed model, random intercept, K = %d (BIC-selected); red = within-class fitted mean", Kbest))
ggsave(file.path(BASE, "fig_cth_lcmm_trajectories.png"), gp1, width = 12, height = 3.9, dpi = 200)

gp2 <- ggplot(prof, aes(factor(class, levels = ord), slope_pts_yr)) +
  geom_col(fill = "#2E6DA4", width = .6) +
  geom_text(aes(label = sprintf("%+.2f (%.1f%%; n=%d)", slope_pts_yr, pct, n)),
            vjust = -0.4, size = 3) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = expansion(mult = c(0, .18))) +
  labs(x = "CTh-Score progression class (slowest to fastest)",
       y = "class mean slope (points / year)",
       title = "Class-specific rates of CTh-Score progression")
ggsave(file.path(BASE, "fig_cth_lcmm_classrate.png"), gp2, width = 8, height = 3.4, dpi = 200)

saveRDS(list(selection = sel, K = Kbest, model = mB, fits = fits, knee = knee,
             profile = prof, ols = ols, b0 = b0), file.path(BASE, "cth_lcmm_fit.rds"))
cat("\n-> cth_lcmm_model_selection.csv / cth_lcmm_class_profile.csv")
cat(" / fig_cth_lcmm_trajectories.png / fig_cth_lcmm_classrate.png / cth_lcmm_fit.rds\n")
