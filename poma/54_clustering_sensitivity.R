# ============================================================================
#  Clustering in the POMA matched design, and a cluster-robust refit.
#
#  The released POMA file gives 450 knee units from 415 participants.  35
#  participants contribute BOTH knees, each to a DIFFERENT matched pair, and
#  no pair contains two knees from the same participant.  Matched pairs are
#  therefore not all independent: those 35 participants each link two pairs.
#
#  This script (1) documents the clustering from the source release and
#  (2) refits the primary within-pair conditional logistic models with a
#  participant-level cluster-robust variance, so the reported P values can be
#  checked against the independence assumption.
# ============================================================================
options(warn = 1)
suppressMessages({ library(survival); library(dplyr) })

B <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot/"
p <- function(...) cat(..., "\n", sep = "")

d <- read.csv(paste0(B, "poma_analysis_long.csv"), stringsAsFactors = FALSE)
u <- d %>% distinct(unit, .keep_all = TRUE)

p("=================== clustering structure ===================")
p("knee units : ", nrow(u), " | participants : ", n_unique <- length(unique(u$id)),
  " | pairs : ", length(unique(u$newstrata)))
pairs_per_person <- u %>% group_by(id) %>% summarise(np = n_distinct(newstrata),
                                                     nk = n(), .groups = "drop")
p("participants contributing 2 pairs : ", sum(pairs_per_person$np == 2),
  " (their knees: ", sum(pairs_per_person$nk[pairs_per_person$np == 2]), ")")
p("pairs containing 2 knees from ONE participant : ",
  sum(tapply(u$id, u$newstrata, function(z) length(unique(z)))) == 0)
p("  -> pairs with >1 distinct participant : ",
  sum(tapply(u$id, u$newstrata, function(z) length(unique(z))) > 1))

## primary analysis set = complete matched pairs (Set A), built exactly as in
## 10c_clogit_seq_primary.R so the cluster-robust refit is the same model
PRIM <- "cMFTC_ThCtAB_aMe"
d$slope <- -d[[paste0(PRIM, "__eblup")]]
d$base  <- d[[paste0(PRIM, "__base")]]
need <- c("slope", "base", "V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")
A <- d %>%
  group_by(newstrata) %>%
  mutate(n_ok = sum(complete.cases(across(all_of(need)))), n_tot = n()) %>%
  ungroup() %>%
  filter(n_ok == n_tot, n_tot == 2, complete.cases(across(all_of(need))))
ua <- A %>% distinct(unit, .keep_all = TRUE)
pp2 <- ua %>% group_by(id) %>% summarise(np = n_distinct(newstrata), nk = n(),
                                         .groups = "drop")
p("\n--- primary Set A ---")
p("knees : ", nrow(ua), " | participants : ", length(unique(ua$id)),
  " | pairs : ", length(unique(ua$newstrata)))
p("participants contributing 2 pairs : ", sum(pp2$np == 2),
  " (knees: ", sum(pp2$nk[pp2$np == 2]), ")")
p("pairs linked by a shared participant : ", 2 * sum(pp2$np == 2))

## --------------------------- cluster-robust refit ---------------------------
zscale <- function(x) as.numeric(scale(x))
P <- ua %>% mutate(rate = slope, z = zscale(rate)) %>%
  mutate(kl = V00XRKL, strata = newstrata)

p("\n=================== cluster-robust refit (Set A) ===================")
ors <- function(m, nm) {
  s  <- summary(m)$coefficients
  cc <- exp(confint(m))
  data.frame(model = nm, term = rownames(s),
             OR = round(exp(coef(m)), 3),
             lo = round(cc[, 1], 3), hi = round(cc[, 2], 3),
             p = signif(s[, "Pr(>|z|)"], 3), row.names = NULL)
}

m0  <- clogit(case ~ z + strata(strata), data = P, method = "efron")
m0c <- clogit(case ~ z + strata(strata) + cluster(id), data = P,
              method = "efron")
m6  <- clogit(case ~ z + V00AGE + P01BMI + kl + WOMAC_pain + strata(strata),
              data = P, method = "efron")
m6c <- clogit(case ~ z + V00AGE + P01BMI + kl + WOMAC_pain + strata(strata) +
                cluster(id), data = P, method = "efron")

tab <- rbind(cbind(variance = "model-based", ors(m0, "M0 crude")),
             cbind(variance = "cluster-robust (participant)", ors(m0c, "M0 crude")),
             cbind(variance = "model-based", ors(m6, "M6 adjusted")),
             cbind(variance = "cluster-robust (participant)", ors(m6c, "M6 adjusted")))
print(tab, row.names = FALSE)
write.csv(tab, paste0(B, "clogit_cluster_robust.csv"), row.names = FALSE)

p("\nSE ratio (cluster-robust / model-based), exposure term:")
for (nm in c("M0 crude", "M6 adjusted")) {
  a <- tab %>% filter(model == nm, term == "z", variance == "model-based")
  b <- tab %>% filter(model == nm, term == "z", variance == "cluster-robust (participant)")
  p("  ", nm, " : (log-OR SE) ", round(log(a$hi / a$lo) / 3.92, 4), " -> ",
    round(log(b$hi / b$lo) / 3.92, 4))
}
