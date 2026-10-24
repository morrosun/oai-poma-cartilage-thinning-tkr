# ============================================================================
#  Which "rapid thinning" definition do the three interpretation tests use?
#
#  Reviewer's point: script 28 computes the rank-equivalence test A1 on the
#  fastest class alone (class 4, n = 25) but tests B (incremental information)
#  and C (gradient below threshold) on the two fastest classes (3+4, n = 67).
#  Those three results were then chained in the manuscript as if they all
#  referred to class 4.
#
#  This script runs B and C under BOTH definitions so the two can be reported
#  side by side, and records the number of pairs each definition drops, which
#  is what exposed the inconsistency (25 knees -> >=166 pairs must remain,
#  yet test C reported 134 pairs).
#
#  Sign convention: eblup is a thickness slope in mm/yr, NEGATIVE = thinning.
#                   rate <- -eblup   (larger = faster thinning)
# ============================================================================
options(warn = 1)
suppressMessages({ library(survival) })

B <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot/"
p <- function(...) cat(..., "\n", sep = "")

f <- readRDS(paste0(B, "lcmm_fit.rds"))
K <- f$knee
P <- f$pair_data

K$rate <- -K$eblup
P$rate <- -P$eblup
K$cl   <- as.integer(K$cl)
P$cl   <- as.integer(P$cl)

cls_by_rate <- as.integer(names(sort(tapply(K$rate, K$cl, mean), decreasing = TRUE)))
fast_cl     <- cls_by_rate[1]          # fastest class only
rapid_cl    <- cls_by_rate[1:2]        # two fastest classes

p("classes ordered fastest -> slowest : ", paste(cls_by_rate, collapse = " > "))
p("class sizes (389 classified knees) : ",
  paste(sprintf("cl%d=%d", as.integer(names(table(K$cl))), as.integer(table(K$cl))),
        collapse = "  "))
p("knees in the pair set that are class ", fast_cl, " : ", sum(P$cl == fast_cl))
p("knees in the pair set that are classes ", paste(rapid_cl, collapse = "+"),
  " : ", sum(P$cl %in% rapid_cl))
p("  of which cases / controls (class ", fast_cl, ") : ",
  sum(P$case == 1 & P$cl == fast_cl), " / ", sum(P$case == 0 & P$cl == fast_cl))
p("  of which cases / controls (classes ", paste(rapid_cl, collapse = "+"), ") : ",
  sum(P$case == 1 & P$cl %in% rapid_cl), " / ",
  sum(P$case == 0 & P$cl %in% rapid_cl))

ors_all <- function(m, nm) {
  s  <- summary(m)$coefficients
  cc <- exp(confint(m))
  data.frame(model = nm, term = rownames(s),
             OR = round(exp(coef(m)), 3),
             lo = round(cc[, 1], 3), hi = round(cc[, 2], 3),
             p = signif(s[, "Pr(>|z|)"], 3), row.names = NULL)
}

run_tests <- function(ind, lab) {
  p("\n=================== definition: ", lab, " ===================")
  P$rapid <- ind
  P$z     <- as.numeric(scale(P$rate))

  ## ---- how many pairs does this definition drop? --------------------------
  nk   <- length(unique(P$strata))
  drop <- ave(P$rapid, P$strata, FUN = function(z) any(z))
  nrem <- length(unique(P$strata[drop == 0]))
  p("pairs total = ", nk, " | pairs containing >=1 'rapid' knee = ",
    sum(tapply(drop, P$strata, function(z) z[1])), " | pairs remaining = ", nrem)

  ## ---- B : incremental information ----------------------------------------
  m_z    <- clogit(case ~ z      + strata(strata), data = P)
  m_r    <- clogit(case ~ rapid  + strata(strata), data = P)
  m_both <- clogit(case ~ z + rapid + strata(strata), data = P)
  tab <- rbind(ors_all(m_z, "B1 z-rate only"),
               ors_all(m_r, "B2 class indicator only"),
               ors_all(m_both, "B3 both"))
  print(tab, row.names = FALSE)
  lrt_z <- as.numeric(2 * (logLik(m_both) - logLik(m_r)))
  lrt_r <- as.numeric(2 * (logLik(m_both) - logLik(m_z)))
  p("LRT z-rate | class : chi2 = ", round(lrt_z, 2),
    " p = ", signif(pchisq(lrt_z, 1, lower.tail = FALSE), 3))
  p("LRT class  | z-rate: chi2 = ", round(lrt_r, 2),
    " p = ", signif(pchisq(lrt_r, 1, lower.tail = FALSE), 3))
  ## separation flag
  se <- tryCatch({
    bb <- coef(m_both); ss <- sqrt(diag(vcov(m_both)))
    max(abs(bb) / ss) > 4
  }, error = function(e) NA)
  p("any |beta|/SE > 4 in the joint model (separation warning) : ", se)

  ## ---- C : gradient below the threshold -----------------------------------
  Pnr <- P[drop == 0, ]
  p("C: pairs with no 'rapid' knee in either member = ",
    length(unique(Pnr$strata)), " pairs / ", nrow(Pnr), " knees")
  if (length(unique(Pnr$strata)) > 10) {
    m_nr <- clogit(case ~ z + strata(strata), data = Pnr)
    tb <- ors_all(m_nr, "C1 z-rate, non-rapid pairs only")
    print(tb, row.names = FALSE)
  } else {
    tb <- NULL
  }
  list(B = tab, C = tb,
       pairs_removed = nk - nrem, pairs_remaining = nrem)
}

## ---------------------------------------------------------------------------
## Separation audit (added 2026-10-24, reviewer round 2 item 7.1).
## "24 of the 25 class-4 knees were replaced" does NOT by itself establish
## complete separation: in a within-pair conditional likelihood, separation
## depends on the DIRECTION of the discordance inside each pair.  Count the
## four possible pair configurations for the fastest class before the claim is
## made in the manuscript.
## ---------------------------------------------------------------------------
sep_tab <- function(ind, lab) {
  P$flag <- ind
  s <- split(seq_len(nrow(P)), P$strata)
  cfg <- vapply(s, function(ix) {
    cc <- P$case[ix]; ff <- P$flag[ix]
    if (length(ix) != 2) return("incomplete")
    if (all(ff))            return("both_class")
    if (!any(ff))           return("neither_class")
    if (ff[cc == 1])        return("case_only")
    "control_only"
  }, character(1))
  tb <- as.data.frame(table(configuration = cfg), stringsAsFactors = FALSE)
  tb$definition <- lab
  tb
}
SEP <- rbind(sep_tab(P$cl == fast_cl,
                     sprintf("class %d only (n=%d)", fast_cl, sum(P$cl == fast_cl))),
             sep_tab(P$cl %in% rapid_cl,
                     sprintf("classes %s (n=%d)", paste(rapid_cl, collapse = "+"),
                             sum(P$cl %in% rapid_cl))))
p("\n=================== separation audit ===================")
p("pair configurations by definition")
print(SEP, row.names = FALSE)
cf <- SEP[SEP$definition == SEP$definition[1], ]
p("class-4-only : case_only = ", cf$Freq[cf$configuration == "case_only"],
  " | control_only = ", cf$Freq[cf$configuration == "control_only"],
  " | both = ", cf$Freq[cf$configuration == "both_class"],
  " | neither = ", cf$Freq[cf$configuration == "neither_class"])
write.csv(SEP, paste0(B, "lcmm_separation_pairs.csv"), row.names = FALSE)

out_fast  <- run_tests(P$cl == fast_cl,
                       sprintf("class %d only (n = %d knees)", fast_cl,
                               sum(P$cl == fast_cl)))
out_rapid <- run_tests(P$cl %in% rapid_cl,
                       sprintf("classes %s (n = %d knees)",
                               paste(rapid_cl, collapse = "+"),
                               sum(P$cl %in% rapid_cl)))

write.csv(rbind(cbind(definition = "fastest class only", out_fast$B),
                cbind(definition = "two fastest classes", out_rapid$B)),
          paste0(B, "lcmm_rate_vs_class_by_definition.csv"), row.names = FALSE)
write.csv(rbind(cbind(definition = "fastest class only", out_fast$C),
                cbind(definition = "two fastest classes", out_rapid$C)),
          paste0(B, "lcmm_rate_below_threshold_by_definition.csv"),
          row.names = FALSE)
write.csv(data.frame(definition = c("fastest class only", "two fastest classes"),
                     knees = c(sum(P$cl == fast_cl), sum(P$cl %in% rapid_cl)),
                     pairs_removed = c(out_fast$pairs_removed, out_rapid$pairs_removed),
                     pairs_remaining = c(out_fast$pairs_remaining,
                                         out_rapid$pairs_remaining)),
          paste0(B, "lcmm_definition_pair_counts.csv"), row.names = FALSE)

p("\n=============================== SUMMARY ===============================")
p("fastest class only : drops ", out_fast$pairs_removed,
  " pairs -> ", out_fast$pairs_remaining, " pairs remain for test C")
p("two fastest classes: drops ", out_rapid$pairs_removed,
  " pairs -> ", out_rapid$pairs_remaining, " pairs remain for test C")
