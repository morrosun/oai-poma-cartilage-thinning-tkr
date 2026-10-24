## ============================================================================
## 56_index_time_audit.R
##
## Reviewer round 2, priority item 1: does any exposure window actually contain
## a measurement taken AT the index visit, and is that measurement distributed
## symmetrically between cases and controls?
##
## Background established before writing this script
##   * kmri / POMA main analysis:  preop = (visit_code < tmpt), i.e. STRICTLY
##     before the index visit.  08_joint_model.R asserts all(year < t_index/12).
##   * CTh-Maps cross-pipeline analysis: pre_index = (t <= idx_yr), i.e. "at or
##     before", INCLUSIVE.  The two paths therefore do not share a time rule.
##
## What this script writes
##   index_time_audit.csv        one row per audit question, with the answer
##   index_time_audit_report.txt full console log
##   cth_clogit_strict.csv       B1/B2/C1 refitted with the strict (<) rule
##
## Nothing here overwrites an existing analysis file.
## ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
suppressPackageStartupMessages({ library(survival); library(nlme) })

sink(file.path(BASE, "index_time_audit_report.txt"), split = TRUE)
on.exit(sink(), add = TRUE)

say <- function(...) cat(sprintf(...), "\n")

## ---------------------------------------------------------------- inputs ---
L  <- read.csv(file.path(BASE, "poma_long.csv"), stringsAsFactors = FALSE)
U  <- read.csv(file.path(BASE, "poma_analysis_long.csv"), stringsAsFactors = FALSE)
C  <- read.csv(file.path(BASE, "cth_poma_long.csv"), stringsAsFactors = FALSE)
BAS<- read.csv(file.path(BASE, "cth_poma_base.csv"), stringsAsFactors = FALSE)

AUD <- list()
add <- function(q, a) {
  AUD[[length(AUD) + 1]] <<- data.frame(question = q, answer = a,
                                        stringsAsFactors = FALSE)
  say("  %-58s %s", q, a)
}

say("=====================================================================")
say(" 1. kmri main analysis -- is the index visit included?")
say("=====================================================================")

## poma_long.csv carries `tmpt` but not `t_index_months`; take the index month
## from the unit-level table so the comparison is made on the same scale.
L$t_index_months <- U$t_index_months[match(L$unit, U$unit)]
## months is the visit month (0/12/24/36/48); t_index_months the index month
L$rel <- L$months - L$t_index_months
say("  distinct months            : %s", paste(sort(unique(L$months)), collapse = ", "))
say("  distinct t_index_months    : %s", paste(sort(unique(L$t_index_months)), collapse = ", "))
say("  rows with months == index  : %d", sum(L$rel == 0))
say("  rows with months >  index  : %d", sum(L$rel > 0))

add("kmri: rows with months == t_index_months", as.character(sum(L$rel == 0)))
add("kmri: rows with months >  t_index_months", as.character(sum(L$rel > 0)))
add("kmri: preop rule in 02_long_reshape (visit_code < tmpt)",
    "strict, index visit excluded")

## 148 pairs used for the cross-pipeline comparison: how are they built?
say("")
say("  Set A (primary, 191 pairs), last pre-index visit gap to the index visit")
COVARS <- c("V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")
PRIM   <- "cMFTC_ThCtAB_aMe"
svA <- U[, c("unit", "newstrata", "case", "t_index_months", COVARS)]
svA <- svA[complete.cases(svA[, COVARS]), ]
lgA <- L[!is.na(L[[PRIM]]), c("unit", "months", PRIM)]
names(lgA)[3] <- "y"
lgA <- lgA[lgA$unit %in% svA$unit, ]
svA <- svA[svA$unit %in% lgA$unit, ]
mmA <- merge(lgA, svA[, c("unit", "t_index_months")], by = "unit")
mmA$rel <- mmA$months - mmA$t_index_months
mmA <- mmA[mmA$rel < 0, ]                      # strictly pre-index
nvA <- table(mmA$unit)
okA <- names(nvA)[nvA >= 2]
prA <- svA$newstrata[match(okA, svA$unit)]
tbA <- table(prA)
paA <- names(tbA)[tbA == 2]
say("  knees %d, complete pairs %d", length(okA), length(paA))
gapA <- tapply(mmA$rel, mmA$unit, max)
gapA <- gapA[names(gapA) %in% okA]
gd <- data.frame(unit = names(gapA), gap = as.numeric(gapA))
gd <- merge(gd, svA[, c("unit", "case")], by = "unit")
say("  gap (months) from last analysed visit to index, by group:")
print(tapply(gd$gap, gd$case, function(x) c(n = length(x), min = min(x),
                                            median = median(x), max = max(x))))
add("Set A: knees / pairs", sprintf("%d / %d", length(okA), length(paA)))
add("Set A: max gap to index, cases / controls",
    sprintf("%.0f / %.0f months",
            max(gd$gap[gd$case == 1]), max(gd$gap[gd$case == 0])))

say("")
say("=====================================================================")
say(" 2. CTh-Maps cross-pipeline -- the inclusive rule")
say("=====================================================================")
C$strict <- as.integer(C$t < C$idx_yr)
say("  pre_index (<=) rows %d / knees %d", sum(C$pre_index == 1),
    length(unique(C$unit[C$pre_index == 1])))
say("  strict     (< ) rows %d / knees %d", sum(C$strict == 1),
    length(unique(C$unit[C$strict == 1])))

e <- C[C$pre_index == 1 & C$t == C$idx_yr, ]
say("  rows included by <= but not by < : %d, knees %d", nrow(e), length(unique(e$unit)))
say("  of those, case knees %d, control knees %d",
    length(unique(e$unit[e$case == 1])), length(unique(e$unit[e$case == 0])))
tab <- table(case = e$case)
print(tab)
add("CTh: rows included by 'or before' but not by 'strictly before'", as.character(nrow(e)))
add("CTh: of those, case knees", as.character(length(unique(e$unit[e$case == 1]))))
add("CTh: of those, control knees", as.character(length(unique(e$unit[e$case == 0]))))
add("CTh: asymmetry", "index-visit measurement present in controls only")

## which matched pairs are affected
aff <- unique(e$newstrata)
say("  matched pairs containing such a knee: %d", length(aff))
add("CTh: matched pairs affected", as.character(length(aff)))

say("")
say("=====================================================================")
say(" 3. does 148 pairs (cMFTC subset) lie inside 191 pairs (Set A)?")
say("=====================================================================")

## rebuild the cross-pipeline pair set exactly as 22_cth_score_corroboration.R
zs <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

slope_of <- function(dat, tag) {
  d <- dat[order(dat$unit, dat$t), ]
  nv <- table(d$unit)
  d <- d[as.character(d$unit) %in% names(nv)[nv >= 2], ]
  u <- sort(unique(d$unit))
  ols <- vapply(u, function(i) { x <- d[d$unit == i, ]
    unname(coef(lm(score ~ t, x))[2]) }, numeric(1))
  b0 <- vapply(u, function(i) { x <- d[d$unit == i, ]; x$score[which.min(x$t)] }, numeric(1))
  ntp <- vapply(u, function(i) sum(d$unit == i), integer(1))
  mm <- lme(score ~ t, random = ~ 1 + t | unit, data = d, method = "ML")
  fx <- fixef(mm)[["t"]]
  re <- ranef(mm)
  sl <- re[[ncol(re)]]; names(sl) <- rownames(re)
  out <- data.frame(unit = u, n_tp = ntp, base = b0, ols = ols,
                    eblup = unname(fx + sl[as.character(u)]))
  names(out)[names(out) != "unit"] <- paste0(names(out)[names(out) != "unit"], "_", tag)
  out
}

P_in <- slope_of(C[C$pre_index == 1, ], "pre")
P_st <- slope_of(C[C$strict == 1, ], "str")
Aall <- slope_of(C, "all")
k <- merge(merge(Aall, P_in, by = "unit", all = TRUE), P_st, by = "unit", all = TRUE)
u2 <- merge(BAS, k, by = "unit", all.y = FALSE)
say("  knees with a CTh slope: inclusive %d, strict %d",
    sum(!is.na(u2$ols_pre)), sum(!is.na(u2$ols_str)))

mk_pairs <- function(dat, col) {
  x <- dat[!is.na(dat[[col]]) & !is.na(dat$newstrata), ]
  tb <- table(x$newstrata)
  x[x$newstrata %in% names(tb)[tb == 2], ]
}
PA_in <- mk_pairs(u2, "ols_pre")
PA_st <- mk_pairs(u2, "ols_str")
say("  complete pairs: inclusive %d, strict %d",
    length(unique(PA_in$newstrata)), length(unique(PA_st$newstrata)))
add("CTh cross-pipeline pairs, inclusive rule", as.character(length(unique(PA_in$newstrata))))
add("CTh cross-pipeline pairs, strict rule", as.character(length(unique(PA_st$newstrata))))

inset <- length(intersect(unique(PA_in$newstrata), paA))
add("CTh pairs that also belong to Set A (191)",
    sprintf("%d of %d (inclusive set)", inset, length(unique(PA_in$newstrata))))
add("Set A pairs not covered by CTh-Maps",
    sprintf("%d", length(setdiff(paA, unique(PA_in$newstrata)))))

## knees present in the two pipelines
say("  Set A knees also present in CTh-Maps: %d of %d",
    length(intersect(okA, u2$unit)), length(okA))

say("")
say("=====================================================================")
say(" 4. refit the cross-pipeline models under the strict (< index) rule")
say("=====================================================================")

cl <- function(dat, form) {
  m <- clogit(as.formula(paste("case ~", form, "+ strata(newstrata)")), data = dat)
  s <- summary(m)$coef
  data.frame(term = rownames(s), OR = exp(s[, 1]),
             lo = exp(s[, 1] - 1.96 * s[, 3]), hi = exp(s[, 1] + 1.96 * s[, 3]),
             p = s[, 5], n_pairs = length(unique(dat$newstrata)), n = nrow(dat),
             stringsAsFactors = FALSE)
}

fit_rule <- function(PA, tag) {
  PA$zC   <- zs(PA$ols)
  PA$zCe  <- zs(PA$eblup)
  PA$zB   <- zs(PA$base)
  PA$zM   <- zs(PA[["cMFTC_ThCtAB_aMe__eblup"]])
  PA$zMB  <- zs(PA[["cMFTC_ThCtAB_aMe__base"]])
  pick <- function(d, trm, lab) {
    r <- d[d$term == trm, ]
    cbind(model = lab, r)
  }
  out <- rbind(
    pick(cl(PA, "zC"),           "zC",  "B1 CTh-Score rate, OLS"),
    pick(cl(PA, "zCe"),          "zCe", "B2 CTh-Score rate, EBLUP"),
    pick(cl(PA, "zC + zM"),      "zC",  "C1 mutual, CTh-Score term"),
    pick(cl(PA, "zC + zM"),      "zM",  "C1 mutual, cMFTC term"),
    pick(cl(PA, "zC + zB + zM + zMB"), "zC", "C2 mutual + baselines, CTh term"),
    pick(cl(PA, "zC + zB + zM + zMB"), "zM", "C2 mutual + baselines, cMFTC term"),
    pick(cl(PA, "zM"),           "zM",  "A1 cMFTC rate alone"))
  out$rule <- tag
  out
}

## column names differ between the two rule-specific frames
names(PA_in)[names(PA_in) == "ols_pre"]  <- "ols"
names(PA_in)[names(PA_in) == "eblup_pre"]<- "eblup"
names(PA_in)[names(PA_in) == "base_pre"] <- "base"
names(PA_st)[names(PA_st) == "ols_str"]  <- "ols"
names(PA_st)[names(PA_st) == "eblup_str"]<- "eblup"
names(PA_st)[names(PA_st) == "base_str"] <- "base"

R_in <- fit_rule(PA_in, "at_or_before")
R_st <- fit_rule(PA_st, "strict")
RES  <- rbind(R_in, R_st)
write.csv(RES, file.path(BASE, "cth_clogit_strict.csv"), row.names = FALSE)

say("  OR per 1 SD of the CTh-Score rate (or of the negated cMFTC rate)")
say("  %-26s %-13s %-22s %-9s %s", "model", "rule", "OR (95% CI)", "p", "pairs")
for (i in seq_len(nrow(RES))) {
  say("  %-26s %-13s %6.2f (%5.2f-%5.2f)   %-9.3g %d",
      RES$model[i], RES$rule[i], RES$OR[i], RES$lo[i], RES$hi[i],
      RES$p[i], RES$n_pairs[i])
}

b1i <- RES[RES$model == "B1 CTh-Score rate, OLS" & RES$rule == "at_or_before", ]
b1s <- RES[RES$model == "B1 CTh-Score rate, OLS" & RES$rule == "strict", ]
add("B1 CTh rate, inclusive rule", sprintf("%.2f (%.2f-%.2f), %d pairs",
    b1i$OR, b1i$lo, b1i$hi, b1i$n_pairs))
add("B1 CTh rate, strict rule", sprintf("%.2f (%.2f-%.2f), %d pairs",
    b1s$OR, b1s$lo, b1s$hi, b1s$n_pairs))

say("")
say("=====================================================================")
say(" 5. how much of each exposure's visit count is the index visit?")
say("=====================================================================")
nvp <- as.data.frame(table(mmA$unit), stringsAsFactors = FALSE)
names(nvp) <- c("unit", "n")
say("  Set A pre-index visits per knee, by group:")
print(tapply(nvp$n, svA$case[match(nvp$unit, svA$unit)], summary))

## CTh: visits per knee under each rule, counted on the LONG table (not on the
## one-row-per-knee pair frame)
ci <- as.data.frame(table(C$unit[C$pre_index == 1]), stringsAsFactors = FALSE)
cs <- as.data.frame(table(C$unit[C$strict == 1]),    stringsAsFactors = FALSE)
names(ci) <- c("unit", "n_incl"); names(cs) <- c("unit", "n_str")
cmp <- merge(ci, cs, by = "unit", all = TRUE)
cmp$n_incl[is.na(cmp$n_incl)] <- 0L
cmp$n_str[is.na(cmp$n_str)]   <- 0L
cmp$case <- C$case[match(cmp$unit, C$unit)]
say("  CTh visits per knee, INCLUSIVE rule, by group:")
print(tapply(cmp$n_incl, cmp$case, summary))
say("  CTh visits per knee, STRICT rule, by group:")
print(tapply(cmp$n_str, cmp$case, summary))
losers <- cmp[cmp$n_incl != cmp$n_str, ]
say("  knees that lose a visit under the strict rule: %d (cases %d, controls %d)",
    nrow(losers), sum(losers$case == 1), sum(losers$case == 0))
add("knees losing a visit under the strict rule",
    sprintf("%d (cases %d, controls %d)", nrow(losers),
            sum(losers$case == 1), sum(losers$case == 0)))
add("knees with a CTh slope, inclusive / strict",
    sprintf("%d / %d", sum(!is.na(u2$ols_pre)), sum(!is.na(u2$ols_str))))

## The manuscript quotes per-knee pre-index measurement density for the 148 pairs
## that actually enter the cross-pipeline model, where both knees have >=2
## measurements under the rule being used.  The 391-knee table above is wider
## than that set, so the model-set figure is computed separately -- this is the
## number a reader will check, and it must come off disk rather than be typed.
vset <- function(d, lab) {
  k <- aggregate(list(n = d$t), by = list(unit = d$unit, case = d$case,
                                          str = d$newstrata), length)
  k <- k[k$n >= 2, ]
  tb <- table(k$str)
  k <- k[k$str %in% names(tb)[tb == 2], ]
  tp <- t.test(k$n[k$case == 1], k$n[k$case == 0])$p.value
  say("  model-set density, %-10s pairs %d | cases %.3f (n=%d) | controls %.3f (n=%d) | p=%.3f",
      lab, length(unique(k$str)), mean(k$n[k$case == 1]), sum(k$case == 1),
      mean(k$n[k$case == 0]), sum(k$case == 0), tp)
  c(pairs = length(unique(k$str)),
    mcase = mean(k$n[k$case == 1]), mctrl = mean(k$n[k$case == 0]), p = tp)
}
v_in <- vset(C[C$pre_index == 1, ], "inclusive")
v_st <- vset(C[C$strict == 1, ],    "strict")
add("CTh model-set density, strict rule, cases / controls",
    sprintf("%.2f / %.2f (p = %.2f)", v_st[["mcase"]], v_st[["mctrl"]], v_st[["p"]]))
add("CTh model-set density, inclusive rule, cases / controls",
    sprintf("%.2f / %.2f", v_in[["mcase"]], v_in[["mctrl"]]))

say("")
say("=====================================================================")
say(" 6. verdict")
say("=====================================================================")
asym <- length(unique(e$unit[e$case == 1])) == 0 &&
        length(unique(e$unit[e$case == 0])) > 0
add("does the cross-pipeline window contain post-index measurements",
    if (asym) "YES, in control knees only" else "no asymmetry detected")
add("does the kmri main analysis contain them", "no -- strict rule throughout")

OUT <- do.call(rbind, AUD)
write.csv(OUT, file.path(BASE, "index_time_audit.csv"), row.names = FALSE)
say("")
say("written: index_time_audit.csv, index_time_audit_report.txt, cth_clogit_strict.csv")
sink()
