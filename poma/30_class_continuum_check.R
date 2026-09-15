#!/usr/bin/env Rscript
# ============================================================================
#  30_class_continuum_check.R
#  Generic "is the latent class anything more than a threshold on a continuous
#  axis?" checker for latent-class / cluster solutions.
#
#  Generalises Scripts/poma/28_lcmm_continuum_check.R, which was hard-coded to
#  the POMA cMFTC LCMM fit.  Everything dataset-specific now lives in a CONFIG
#  file; this engine only needs
#
#      (a) a per-unit table : unit id, class label, one or more continuous
#                             exposure axes (plus optional "level" axes), and
#      (b) optionally a per-unit outcome table : case / time+event / strata.
#
#  Usage
#      Rscript 30_class_continuum_check.R <config.R>
#
#  The config is plain R and must assign `CFG`, a list:
#    label          chr   human-readable label printed in the log
#    tag            chr   prefix for every output file
#    outdir         chr   output directory (created if absent)
#    unit           df    one row per unit; needs a unit id, a class label and
#                         the exposure columns
#    unit_col       chr   name of the unit-id column      (default "unit")
#    class_col      chr   name of the class column        (default "class")
#    exposures      chr   continuous exposure columns to test in turn
#    exposure_labels named chr  optional display names, names = column
#    signs           named num  +1 / -1 per exposure, applied so that a LARGER
#                               value means WORSE (faster / more severe).
#                               Missing entries default to +1.
#    levels         chr   optional "level" columns for the eta^2 contrast
#    level_labels   named chr  optional display names, names = column
#    rapid_k        int   how many fastest classes count as "the rapid tail"
#                         (default 2)
#    bins           int   how many equal-size bins for the dose-response step
#                         (default 5 = quintiles)
#    out            df    optional; tests 3-5 need it
#    out_cols       list  unit = <col>, case = <col>, strata = <col|NULL>,
#                         time = <col|NULL>, event = <col|NULL>
#    out_model      chr   "clogit" (needs strata) | "cox" (needs time+event)
#                         | "logit" (needs case only)
#    out_covars     chr   optional further columns of `out` added to every
#                         model; an "adjusted" model set is then reported too
#    out_covars_by  list  optional, names = exposure -> covariate vector, so an
#                         exposure that IS a covariate (e.g. age) can be left
#                         unadjusted while the others are adjusted
#
#  Tests
#   (1) rank-threshold equivalence : fastest class(es) vs top-n on the axis,
#                                    whole partition vs same-size axis bins
#   (2) variance partition         : eta-squared of class on the axis vs on the
#                                    "level" axes
#   (3) information increment      : axis-only / class-only / both, mutually
#                                    adjusted, LRT in both directions
#   (4) graded below the threshold : is the axis still predictive among the
#                                    units that are all NON-rapid?
#   (5) quintile dose-response     : event % per quintile, per-step and
#                                    factor contrasts
#
#  Convention : each exposure is multiplied by signs[[e]] before use, so the
#               engine always works with `rate` = "larger means worse".  All
#               reported numbers are on that signed scale.
# ============================================================================
options(warn = 1)
suppressMessages({ library(survival) })

ARGS <- commandArgs(trailingOnly = TRUE)
if (length(ARGS) < 1L)
  stop("usage: Rscript 30_class_continuum_check.R <config.R>")
source(ARGS[1])
if (!exists("CFG")) stop("config file must define `CFG`")

TAG     <- CFG$tag
OUTDIR  <- CFG$outdir
LABEL   <- if (!is.null(CFG$label)) CFG$label else TAG
RAPID_K <- if (!is.null(CFG$rapid_k)) as.integer(CFG$rapid_k) else 2L
UCOL    <- if (!is.null(CFG$unit_col))  CFG$unit_col  else "unit"
CCOL    <- if (!is.null(CFG$class_col)) CFG$class_col else "class"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

LOG <- character(0)
p <- function(...) {
  s <- paste(vapply(list(...), function(z) paste0(z, collapse = " "), ""),
             collapse = " ")
  LOG <<- c(LOG, s); cat(s, "\n", sep = "")
}
w <- function(df, f) {
  utils::write.csv(df, file.path(OUTDIR, f), row.names = FALSE)
  p("  [written] ", f)
}

## ====================================================== generic helpers =====
agree <- function(a, b)
  data.frame(n_a = sum(a), n_b = sum(b), both = sum(a & b),
             a_only = sum(a & !b), b_only = sum(b & !a),
             pct_agreement = 100 * mean(a == b),
             jaccard = sum(a & b) / sum(a | b))

adj_rand <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  sa <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
  si <- sum(choose(as.vector(tab), 2))
  (si - sa * sb / choose(n, 2)) / ((sa + sb) / 2 - sa * sb / choose(n, 2))
}

eta2 <- function(v, g) {
  a <- summary(aov(v ~ factor(g)))[[1]]
  a[["Sum Sq"]][1] / sum(a[["Sum Sq"]])
}

LL    <- function(m) as.numeric(logLik(m))
lrt_p <- function(m_big, m_small, df = 1)
  pchisq(2 * (LL(m_big) - LL(m_small)), df, lower.tail = FALSE)

tidy_fit <- function(m, nm) {
  s <- summary(m)$coefficients
  ci <- exp(confint(m))
  data.frame(model = nm, term = names(coef(m)),
             est = round(exp(coef(m)), 3),
             lo = round(ci[, 1], 3), hi = round(ci[, 2], 3),
             p = signif(s[, "Pr(>|z|)"], 3), row.names = NULL)
}

fit_one <- function(dat, terms, covars, fam) {
  rhs <- paste(c(terms, covars), collapse = " + ")
  if (fam == "clogit") {
    survival::clogit(
      as.formula(paste(CFG$out_cols$case, "~", rhs,
                       "+ strata(", CFG$out_cols$strata, ")")), data = dat)
  } else if (fam == "cox") {
    survival::coxph(
      as.formula(paste("Surv(", CFG$out_cols$time, ",", CFG$out_cols$event,
                       ") ~", rhs)), data = dat)
  } else {
    glm(as.formula(paste(CFG$out_cols$case, "~", rhs)),
        data = dat, family = binomial())
  }
}

## signed rate, class codes, and the fastest-to-slowest class order
signed <- function(v, e) {
  s <- if (!is.null(CFG$signs) && e %in% names(CFG$signs)) CFG$signs[[e]] else 1
  as.numeric(s) * as.numeric(v)
}
order_by_mean <- function(rate, cl)
  as.integer(names(sort(tapply(rate, cl, mean), decreasing = TRUE)))

## ============================================================ header ========
p("########################################################################")
p("#  ", LABEL)
p("#  tag = ", TAG, "  |  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
p("########################################################################")

U <- CFG$unit
## class codes are assigned once, from the unit table, and then MATCHED by
## label everywhere else -- so a class that happens to be absent from the
## outcome table can never shift the coding of the others.  Blank / whitespace
## labels (e.g. a class column that is only defined for a subgroup, read back
## from CSV as "") count as missing, never as a class.
clean_class <- function(x) {
  v <- as.character(x)
  v[is.na(v) | !nzchar(trimws(v))] <- NA_character_
  v
}
cls_lv <- sort(unique(clean_class(U[[CCOL]])))
cls_lv <- cls_lv[!is.na(cls_lv)]
cls_code <- function(x) match(clean_class(x), cls_lv)
U$.cl   <- cls_code(U[[CCOL]])
U$.unit <- as.character(U[[UCOL]])
if (anyNA(U$.cl)) {
  p("note: dropped ", sum(is.na(U$.cl)), " units with a missing class label")
  U <- U[!is.na(U$.cl), ]
}

EXPS   <- CFG$exposures
ELAB   <- if (!is.null(CFG$exposure_labels)) CFG$exposure_labels else
  setNames(EXPS, EXPS)
LEVS   <- if (!is.null(CFG$levels)) as.character(CFG$levels) else character(0)
LLAB   <- if (!is.null(CFG$level_labels)) CFG$level_labels else
  setNames(LEVS, LEVS)
lab_of <- function(nm, tab) if (!is.null(tab[[nm]])) as.character(tab[[nm]]) else nm

p("unit table : ", nrow(U), " units / ", length(cls_lv), " classes")
p("class sizes: ", paste(sprintf("%s=%d", cls_lv, as.integer(table(U$.cl))),
                         collapse = "  "))
p("exposures  : ", paste(sprintf("%s (sign %s)", EXPS,
                                 ifelse(EXPS %in% names(CFG$signs),
                                        CFG$signs[EXPS], 1)), collapse = "; "))
if (length(LEVS)) p("levels     : ", paste(sapply(LEVS, lab_of, LLAB), collapse = "; "))

SUMMARY <- list()

## ======================================================= per exposure =======
for (e in EXPS) {
  lab <- ELAB[[e]]
  p("\n")
  p("========================================================================")
  p("===  exposure : ", lab, "  (", e, ")", sep = "")
  p("========================================================================")

  U$rate <- signed(U[[e]], e)
  ok <- is.finite(U$rate)
  if (!all(ok)) p("  note: dropped ", sum(!ok), " units with a non-finite axis")
  Ue <- U[ok, ]
  cls_order <- order_by_mean(Ue$rate, Ue$.cl)
  p("class mean axis, fastest -> slowest : ",
    paste(sprintf("class %s (n=%d) %.4f",
                  cls_lv[cls_order],
                  as.integer(table(Ue$.cl)[as.character(cls_order)]),
                  tapply(Ue$rate, Ue$.cl, mean)[as.character(cls_order)]),
          collapse = "  >  "))

  ## --------------------------------------------------(1) rank equivalence --
  p("\n---- (1) is the class just a rank threshold on the axis? ----")
  topn <- function(rate, n) order(rate, decreasing = TRUE)[seq_len(n)]

  k1      <- min(RAPID_K, length(cls_order))
  fast_cl <- cls_order[1]
  nf      <- sum(Ue$.cl == fast_cl)
  hard    <- Ue$.cl == fast_cl
  topN    <- rep(FALSE, nrow(Ue)); topN[topn(Ue$rate, nf)] <- TRUE
  r1 <- cbind(comparison = sprintf("fastest class (%s, n=%d) vs top-%d",
                                   cls_lv[fast_cl], nf, nf), agree(hard, topN))

  rapid_cl <- cls_order[seq_len(k1)]
  nr       <- sum(Ue$.cl %in% rapid_cl)
  hard2    <- Ue$.cl %in% rapid_cl
  topN2    <- rep(FALSE, nrow(Ue)); topN2[topn(Ue$rate, nr)] <- TRUE
  r2 <- cbind(comparison = sprintf("fastest %d classes (%s, n=%d) vs top-%d",
                                   k1, paste(cls_lv[rapid_cl], collapse = "+"),
                                   nr, nr), agree(hard2, topN2))

  sizes <- as.integer(table(Ue$.cl)[as.character(cls_order)])
  qbin  <- integer(nrow(Ue)); o <- order(Ue$rate, decreasing = TRUE); st <- 0
  for (i in seq_along(sizes)) {
    qbin[o[st + seq_len(sizes[i])]] <- i; st <- st + sizes[i]
  }
  qbin_mapped <- cls_order[qbin]
  ar <- adj_rand(Ue$.cl, qbin_mapped)
  r3 <- data.frame(comparison = "whole partition vs same-size axis bins",
                   n_a = NA, n_b = NA, both = NA, a_only = NA, b_only = NA,
                   pct_agreement = 100 * mean(Ue$.cl == qbin_mapped),
                   jaccard = NA)
  res1 <- rbind(r1, r2, r3)
  print(res1, row.names = FALSE, digits = 4)
  p("adjusted Rand index, class vs axis bins : ", round(ar, 3))
  w(res1, sprintf("%s_rank_equivalence__%s.csv", TAG, e))

  ## ------------------------------------------------- (2) variance partition --
  p("\n---- (2) what does the class actually explain? ----")
  vr <- eta2(Ue$rate, Ue$.cl)
  p(sprintf("eta-squared of the class solution on the AXIS    : %5.1f %%",
            100 * vr))
  eta_rows <- data.frame(target = lab, role = "axis",
                         eta_squared = vr, n_units = nrow(Ue))
  for (lv in LEVS) {
    if (!lv %in% names(Ue)) next
    ev <- eta2(as.numeric(Ue[[lv]]), Ue$.cl)
    p(sprintf("eta-squared of the class solution on %-16s : %5.1f %%",
              paste0(lab_of(lv, LLAB), " (", lv, ")"), 100 * ev))
    eta_rows <- rbind(eta_rows,
                      data.frame(target = lab_of(lv, LLAB), role = "level",
                                 eta_squared = ev, n_units = nrow(Ue)))
  }
  w(eta_rows, sprintf("%s_eta_squared__%s.csv", TAG, e))

  if (vr >= 0.5 && all(eta_rows$eta_squared[eta_rows$role == "level"] < vr))
    p("-> the class is dominated by the AXIS: it explains the axis better than",
      " it explains the level variables")

  ## --------------------------------------------------- (3)-(5) need `out` --
  if (is.null(CFG$out)) { p("\n(no outcome table supplied: tests 3-5 skipped)"); next }

  O <- CFG$out
  O$rate  <- signed(O[[e]], e)
  O$.cl   <- cls_code(O[[CCOL]])
  if (anyNA(O$.cl)) {
    p("  note: ", sum(is.na(O$.cl)), " outcome rows dropped (no class label)")
    O <- O[!is.na(O$.cl), ]
  }
  O$rapid <- O$.cl %in% rapid_cl
  fam <- CFG$out_model
  cov <- if (!is.null(CFG$out_covars_by) && e %in% names(CFG$out_covars_by))
    CFG$out_covars_by[[e]] else
    if (!is.null(CFG$out_covars)) CFG$out_covars else character(0)
  EFF <- if (fam == "cox") "HR" else "OR"
  ev_col <- if (!is.null(CFG$out_cols$case)) CFG$out_cols$case else CFG$out_cols$event
  p("\noutcome table : ", nrow(O), " units / ", sum(O[[ev_col]]),
    " events / model ", fam,
    if (length(cov)) paste0(" / adjusted for ", paste(cov, collapse = ", ")) else "")
  O$z <- as.numeric(scale(O$rate))

  ## ------------------------------------- (3) information increment ---------
  model_sets <- list(list(nm = "crude", cov = character(0)))
  if (length(cov)) model_sets <- c(model_sets, list(list(nm = "adjusted", cov = cov)))

  tab3 <- NULL
  for (ms in model_sets) {
    m_z <- fit_one(O, "z", ms$cov, fam)
    m_c <- fit_one(O, "rapid", ms$cov, fam)
    m_b <- fit_one(O, "z + rapid", ms$cov, fam)
    a <- tidy_fit(m_z, sprintf("M1 axis only [%s]", ms$nm))
    b <- tidy_fit(m_c, sprintf("M2 rapid class only [%s]", ms$nm))
    cc<- tidy_fit(m_b, sprintf("M3 both [%s]", ms$nm))
    cc<- cc[cc$term %in% c("z", "rapidTRUE", "rapid"), ]
    tab3 <- rbind(tab3, a, b, cc)
    p(sprintf("  [%s] LRT axis | class : chi2 = %.2f  P = %s", ms$nm,
              2 * (LL(m_b) - LL(m_c)), signif(lrt_p(m_b, m_c), 3)))
    p(sprintf("  [%s] LRT class | axis : chi2 = %.2f  P = %s", ms$nm,
              2 * (LL(m_b) - LL(m_z)), signif(lrt_p(m_b, m_z), 3)))
  }
  p("\n(z-standardised axis -> the ", EFF, " is per 1 SD of the axis)")
  print(tab3, row.names = FALSE)
  w(tab3, sprintf("%s_axis_vs_class__%s.csv", TAG, e))

  ## -------------------------- (3b) the WHOLE class variable, not just the tail
  ## a binary "rapid class" indicator depends on how many classes are called
  ## rapid; this tests all K-1 dummies at once, which is the direct question
  O$.clf <- factor(O$.cl, levels = seq_along(cls_lv), labels = cls_lv)
  if (nlevels(O$.clf) > 1L) {
    p("\n---- (3b) the whole ", nlevels(O$.clf),
      "-level class variable vs the axis ----")
    p("  (class factor is entered with labels; reference class = ",
      cls_lv[1], ")", sep = "")
    df_c  <- nlevels(O$.clf) - 1L
    for (ms in model_sets) {
      m_z   <- fit_one(O, "z", ms$cov, fam)
      m_cf  <- fit_one(O, c("z", ".clf"), ms$cov, fam)
      m_cfa <- fit_one(O, ".clf", ms$cov, fam)
      q1 <- lrt_p(m_cf, m_z, df_c)
      q2 <- lrt_p(m_cf, m_cfa, 1L)
      p(sprintf("  [%s] LRT whole class (%d df) | axis : chi2 = %.2f  P = %s",
                ms$nm, df_c, 2 * (LL(m_cf) - LL(m_z)), signif(q1, 3)))
      p(sprintf("  [%s] LRT axis | whole class        : chi2 = %.2f  P = %s",
                ms$nm, 2 * (LL(m_cf) - LL(m_cfa)), signif(q2, 3)))
    }
    w(tidy_fit(m_cf, sprintf("axis + whole class factor [%s]",
                             if (length(cov)) "adjusted" else "crude")),
      sprintf("%s_axis_plus_classfactor__%s.csv", TAG, e))
  }

  ## ----------------------------------------- (4) graded below threshold ----
  p("\n---- (4) is the gradient still there below the threshold? ----")
  if (fam == "clogit") {
    keep <- ave(O$rapid, O[[CFG$out_cols$strata]],
                FUN = function(z) !any(z)) == 1
    nsub <- length(unique(O[[CFG$out_cols$strata]][keep]))
    p("  groups with NO rapid member : ", nsub, " / ", nrow(O[keep, ]), " units",
      " (of ", length(unique(O[[CFG$out_cols$strata]])), " groups)")
  } else {
    keep <- !O$rapid
    nsub <- sum(keep)
    p("  units that are NOT rapid : ", nsub, " / ", nrow(O))
  }
  tab4 <- NULL
  if (nsub >= 10) {
    Onr <- O[keep, ]
    m_nr <- fit_one(Onr, "z", cov, fam)
    tab4 <- tidy_fit(m_nr, sprintf("axis only, non-rapid subset [%s]",
                                   if (length(cov)) "adjusted" else "crude"))
    print(tab4, row.names = FALSE)
    w(tab4, sprintf("%s_below_threshold__%s.csv", TAG, e))
  } else p("  (too few units below the threshold)")

  ## --------------------------------------------- (5) quintile response -----
  NBIN <- if (!is.null(CFG$bins)) as.integer(CFG$bins) else 5L
  p(sprintf("\n---- (5) graded dose-response across %d axis bins ----", NBIN))
  br <- unique(quantile(O$rate, probs = seq(0, 1, length.out = NBIN + 1L)))
  if (length(br) == NBIN + 1L) {
    O$qb <- cut(O$rate, breaks = br, include.lowest = TRUE, labels = FALSE)
    ev_col <- if (!is.null(CFG$out_cols$case)) CFG$out_cols$case else
      CFG$out_cols$event
    qtab <- do.call(rbind, lapply(sort(unique(O$qb)), function(k) {
      s <- O[O$qb == k, ]
      data.frame(bin = k, n = nrow(s), n_event = sum(s[[ev_col]]),
                 pct_event = round(100 * mean(s[[ev_col]]), 1),
                 mean_axis = round(mean(s$rate), 4))
    }))
    print(qtab, row.names = FALSE)
    m_qt <- fit_one(O, "qb", cov, fam)
    m_qf <- fit_one(O, "factor(qb)", cov, fam)
    tab5 <- rbind(tidy_fit(m_qt, sprintf("per one %d-bin step", NBIN)),
                  tidy_fit(m_qf, sprintf("bin factor (ref = bin 1)")))
    print(tab5, row.names = FALSE)
    w(qtab, sprintf("%s_bins__%s.csv", TAG, e))
    w(tab5, sprintf("%s_bin_model__%s.csv", TAG, e))
  } else p("  (the axis has too many ties to form ", NBIN,
           " equal-size bins: skipped)")

  p("\n---- headline ----")
  p(sprintf("  (1) fastest class vs top-%d  : %.1f %% agreement, Jaccard %.3f",
            nf, r1$pct_agreement, r1$jaccard))
  p(sprintf("      whole partition vs bins  : %.1f %% agreement, adjRand %.3f",
            r3$pct_agreement, ar))
  p(sprintf("  (2) class explains %.1f %% of the axis%s", 100 * vr,
            if (length(LEVS)) paste0(" vs ",
              paste(sprintf("%.1f %% of %s", 100 * eta_rows$eta_squared[
                eta_rows$role == "level"], eta_rows$target[eta_rows$role == "level"]),
                collapse = ", ")) else ""))
  SUMMARY[[e]] <- data.frame(
    exposure = lab,
    n_units = nrow(Ue),
    fastest_class = cls_lv[fast_cl], n_fastest = nf,
    pct_agree_fastest = round(r1$pct_agreement, 1),
    jaccard_fastest = round(r1$jaccard, 3),
    pct_agree_partition = round(r3$pct_agreement, 1),
    adj_rand = round(ar, 3),
    eta2_axis = round(100 * vr, 1))
}

## ============================================================= summary ======
if (length(SUMMARY)) {
  S <- do.call(rbind, SUMMARY); rownames(S) <- NULL
  p("\n")
  p("########################################################################")
  p("###  SUMMARY  --  ", LABEL)
  p("########################################################################")
  print(S, row.names = FALSE, digits = 4)
  w(S, paste0(TAG, "_summary.csv"))
  p("\nreading : a class solution that merely thresholds the axis gives",
    " pct_agree_fastest near 100 and adj_rand high; a genuine multi-axis",
    " structure gives neither.")
}

writeLines(LOG, file.path(OUTDIR, paste0(TAG, "_log.txt")))
cat("\n[done] log ->", file.path(OUTDIR, paste0(TAG, "_log.txt")), "\n")
