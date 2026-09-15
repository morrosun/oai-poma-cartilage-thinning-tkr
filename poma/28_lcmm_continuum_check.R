# ============================================================================
#  Is the POMA latent "rapid thinning" class anything more than a threshold
#  on a continuous thinning-rate ranking?
#
#  The same test that was run on CTh-Score in P2 (Scripts/oai_imaging/05),
#  now run on the manuscript's OWN cMFTC LCMM solution, because the abstract
#  claims "about 6% of knees follow a rapid-thinning trajectory, 96% of which
#  were replaced".  If the class is just the top tail of the rate, that claim
#  must be re-worded as a graded continuum rather than a discrete subtype.
#
#  Tests
#   (1) rank-threshold equivalence : fastest class  vs  top-n by observed rate
#       1a  hard class 4 (n = 25)      vs top-25 fastest
#       1b  rapid classes 3+4 (n = 67) vs top-67 fastest
#       1c  whole 4-class partition    vs same-size quantile bins of the rate
#   (2) variance partition : eta-squared of class on RATE vs on BASELINE level
#   (3) does the class add anything beyond the continuous rate?
#       within-pair conditional logistic, LRT both directions
#   (4) is the risk graded BELOW the "rapid" threshold?
#       pairs in which NEITHER knee is rapid -> is z-rate still predictive?
#   (5) graded dose-response : quintiles of the rate, OR Q5 vs Q1 + trend test
#
#  Sign convention : eblup is a thickness slope in mm/yr, NEGATIVE = thinning.
#                    rate <- -eblup   (larger = faster thinning)
# ============================================================================
options(warn = 1)
suppressMessages({ library(survival); library(dplyr) })

B <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot/"
p <- function(...) cat(..., "\n", sep = "")

f <- readRDS(paste0(B, "lcmm_fit.rds"))
K <- f$knee          # 389 knees with >= 2 visits : class, rate, baseline
P <- f$pair_data     # 382 knees = 191 complete matched pairs

K$rate <- -K$eblup   # + = faster thinning
P$rate <- -P$eblup
K$cl   <- as.integer(K$cl)
P$cl   <- as.integer(P$cl)

p("knees with a class : ", nrow(K), " | matched-pair knees : ", nrow(P),
  " | pairs : ", length(unique(P$strata)))
p("class sizes (all 389)            : ", paste(table(K$cl), collapse = "/"))
p("class mean rate (mm/yr, + = fast): ",
  paste(sprintf("cl%d=%.4f", as.integer(names(tapply(K$rate, K$cl, mean))),
                round(tapply(K$rate, K$cl, mean), 4)), collapse = "  "))

## order the classes from fastest to slowest thinning
cls_by_rate <- as.integer(names(sort(tapply(K$rate, K$cl, mean), decreasing = TRUE)))
p("classes ordered fastest -> slowest : ", paste(cls_by_rate, collapse = " > "))

## ------------------------------------------------------------ (1) equivalence
p("\n================ (1) is the class just a rank threshold? ================")

topn <- function(rate, n) order(rate, decreasing = TRUE)[seq_len(n)]

agree <- function(a, b) {
  data.frame(n_a = sum(a), n_b = sum(b), both = sum(a & b),
             a_only = sum(a & !b), b_only = sum(b & !a),
             pct_agreement = 100 * mean(a == b),
             jaccard = sum(a & b) / sum(a | b))
}

## (1a) fastest class vs the same number of fastest knees
fast_cl <- cls_by_rate[1]
nf      <- sum(K$cl == fast_cl)
hard    <- K$cl == fast_cl
topN    <- rep(FALSE, nrow(K)); topN[topn(K$rate, nf)] <- TRUE
r1 <- cbind(comparison = sprintf("class %d (n=%d) vs top-%d fastest",
                                 fast_cl, nf, nf), agree(hard, topN))

## (1b) the two fastest classes vs the same number of fastest knees
rapid_cl <- cls_by_rate[1:2]
nr       <- sum(K$cl %in% rapid_cl)
hard2    <- K$cl %in% rapid_cl
topN2    <- rep(FALSE, nrow(K)); topN2[topn(K$rate, nr)] <- TRUE
r2 <- cbind(comparison = sprintf("classes %s (n=%d) vs top-%d fastest",
                                 paste(rapid_cl, collapse = "+"), nr, nr),
            agree(hard2, topN2))

## (1c) whole partition vs quantile bins of the rate with the same sizes,
##      bins laid out in the same fastest-to-slowest order as the classes
sizes <- as.integer(table(K$cl)[as.character(cls_by_rate)])
qbin  <- integer(nrow(K))
o     <- order(K$rate, decreasing = TRUE)
st    <- 0
for (i in seq_along(sizes)) {
  qbin[o[st + seq_len(sizes[i])]] <- i
  st <- st + sizes[i]
}
qbin_mapped <- cls_by_rate[qbin]          # bin i <-> the i-th fastest class

adj_rand <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  sa <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
  si <- sum(choose(as.vector(tab), 2))
  (si - sa * sb / choose(n, 2)) / ((sa + sb) / 2 - sa * sb / choose(n, 2))
}
ar <- adj_rand(K$cl, qbin_mapped)
r3 <- data.frame(comparison = "whole 4-class partition vs same-size rate bins",
                 n_a = NA, n_b = NA, both = NA, a_only = NA, b_only = NA,
                 pct_agreement = 100 * mean(K$cl == qbin_mapped), jaccard = NA)
res1 <- rbind(r1, r2, r3)
print(res1, row.names = FALSE, digits = 4)
p("adjusted Rand index, class vs rate bins : ", round(ar, 3))
write.csv(res1, paste0(B, "lcmm_rank_equivalence.csv"), row.names = FALSE)

## ----------------------------------------------------- (2) variance partition
p("\n================ (2) what does the class actually explain? ================")
eta2 <- function(v, g) {
  a <- summary(aov(v ~ factor(g)))[[1]]
  a[["Sum Sq"]][1] / sum(a[["Sum Sq"]])
}
e_rate <- eta2(K$rate, K$cl)
e_base <- eta2(K$b0,   K$cl)
p("eta-squared of the 4-class solution")
p("   on the RATE     : ", round(100 * e_rate, 1), "%")
p("   on the BASELINE : ", round(100 * e_base, 1), "%")
write.csv(data.frame(target = c("rate", "baseline"),
                     eta_squared = c(e_rate, e_base)),
          paste0(B, "lcmm_eta_squared.csv"), row.names = FALSE)

## ------------------------- (3) does the class add anything beyond the rate? --
p("\n=========== (3) class vs continuous rate, within-pair conditional ========")
P$z     <- as.numeric(scale(P$rate))
P$rapid <- P$cl %in% rapid_cl

ors_all <- function(m, nm) {
  s  <- summary(m)$coefficients
  cc <- exp(confint(m))
  data.frame(model = nm, term = rownames(s),
             OR = round(exp(coef(m)), 3),
             lo = round(cc[, 1], 3), hi = round(cc[, 2], 3),
             p = signif(s[, "Pr(>|z|)"], 3), row.names = NULL)
}

m_z    <- clogit(case ~ z         + strata(strata), data = P)
m_r    <- clogit(case ~ rapid     + strata(strata), data = P)
m_both <- clogit(case ~ z + rapid + strata(strata), data = P)

tab3 <- rbind(ors_all(m_z,    "M1  z-rate only"),
              ors_all(m_r,    "M2  rapid class only"),
              ors_all(m_both, "M3  both"))
print(tab3, row.names = FALSE)

lrt_z_given_r <- as.numeric(2 * (logLik(m_both) - logLik(m_r)))
lrt_r_given_z <- as.numeric(2 * (logLik(m_both) - logLik(m_z)))
pz <- pchisq(lrt_z_given_r, 1, lower.tail = FALSE)
pr <- pchisq(lrt_r_given_z, 1, lower.tail = FALSE)
p("LRT  z-rate  | rapid class : chi2 = ", round(lrt_z_given_r, 2),
  "  df = 1  p = ", signif(pz, 3))
p("LRT  rapid class | z-rate  : chi2 = ", round(lrt_r_given_z, 2),
  "  df = 1  p = ", signif(pr, 3))
write.csv(tab3, paste0(B, "lcmm_rate_vs_class.csv"), row.names = FALSE)

## ------------------- (4) graded risk BELOW the rapid threshold ---------------
p("\n=========== (4) is the gradient still there below the threshold? =========")
keep <- ave(P$rapid, P$strata, FUN = function(z) !any(z))
Pnr  <- P[keep == 1, ]
p("pairs with NO rapid knee in either member : ", length(unique(Pnr$strata)),
  " pairs / ", nrow(Pnr), " knees")
if (length(unique(Pnr$strata)) > 10) {
  m_nr <- clogit(case ~ z + strata(strata), data = Pnr)
  tab4 <- ors_all(m_nr, "z-rate, non-rapid pairs only")
  print(tab4, row.names = FALSE)
  write.csv(tab4, paste0(B, "lcmm_rate_below_threshold.csv"), row.names = FALSE)
  tab4_z <- tab4
} else {
  tab4_z <- NULL
}

## --------------------------- (5) quintile dose-response ----------------------
p("\n============= (5) graded dose-response across rate quintiles =============")
P$q5 <- cut(P$rate, breaks = quantile(P$rate, probs = seq(0, 1, .2)),
            include.lowest = TRUE, labels = FALSE)
qtab <- P %>% group_by(q5) %>%
  summarise(n = n(), n_TKR = sum(case), pct_TKR = round(100 * mean(case), 1),
            mean_rate = round(mean(rate), 4), .groups = "drop")
print(as.data.frame(qtab), row.names = FALSE)

m_qt <- clogit(case ~ q5                 + strata(strata), data = P)
m_qf <- clogit(case ~ factor(q5)         + strata(strata), data = P)
tab5 <- rbind(ors_all(m_qt, "OR per one-quintile step"),
              ors_all(m_qf, "quintile factor (ref = Q1)"))
print(tab5, row.names = FALSE)
write.csv(qtab, paste0(B, "lcmm_rate_quintiles.csv"), row.names = FALSE)
write.csv(tab5, paste0(B, "lcmm_rate_quintile_OR.csv"), row.names = FALSE)

## ------------------------------------------------------------------ summary --
p("\n=============================== SUMMARY ===============================")
p("(1) fastest class (n=", nf, ") vs top-", nf, " by rate : ",
  round(r1$pct_agreement, 1), "% agreement, Jaccard ", round(r1$jaccard, 3))
p("    rapid classes (n=", nr, ") vs top-", nr, " by rate : ",
  round(r2$pct_agreement, 1), "% agreement, Jaccard ", round(r2$jaccard, 3))
p("    whole partition vs rate bins  : ", round(r3$pct_agreement, 1),
  "% agreement, adjusted Rand ", round(ar, 3))
p("(2) class explains ", round(100 * e_rate, 1), "% of RATE variance vs ",
  round(100 * e_base, 1), "% of BASELINE variance")
p("(3) z-rate given the class : p = ", signif(pz, 3),
  " ; class given the z-rate : p = ", signif(pr, 3))
if (!is.null(tab4_z))
  p("(4) z-rate within non-rapid pairs : OR = ", tab4_z$OR, " (",
    tab4_z$lo, "-", tab4_z$hi, "), p = ", tab4_z$p)
p("(5) OR per quintile step = ", tab5$OR[1], " (", tab5$lo[1], "-",
  tab5$hi[1], "), p = ", tab5$p[1])
p("\n-> ", if (r1$pct_agreement > 80 || e_base > e_rate)
  "CONTINUUM READING : the latent class carries little beyond the rate ranking"
  else "the class is NOT reducible to a rank threshold")
