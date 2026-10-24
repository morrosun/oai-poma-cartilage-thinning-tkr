# ============================================================================
#  57_pair_connectivity_cluster.R
#
#  Reviewer round 3, item 7.  The primary conditional logistic model is fitted
#  with `cluster(id)` (participant).  That groups the variance by participant,
#  which is a defensible grouping for a knee-level exposure, but the dependence
#  induced by the matched design runs BETWEEN PAIRS: a participant who supplies
#  both knees to two different pairs makes those two pairs dependent, and the
#  pair is the stratum, i.e. the unit that the conditional likelihood treats as
#  independent.  Clustering on the participant does not necessarily coincide
#  with the connected components of that pair graph, because a component can
#  chain three or more pairs through successive participants.
#
#  This script builds the pair graph, takes its connected components, and
#  refits the same two models with components as the clustering unit.  It does
#  not introduce a new estimator: same data, same covariates, same strata, only
#  the grouping vector of the robust variance changes.
#
#  Writes  clogit_connectivity.csv  and  connectivity_report.txt
# ============================================================================
options(warn = 1)
suppressMessages({ library(survival); library(dplyr) })

B <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot/"
p <- function(...) cat(..., "\n", sep = "")
con <- file(paste0(B, "connectivity_report.txt"), open = "wt")
say <- function(...) { s <- paste0(..., collapse = ""); cat(s, "\n"); cat(s, "\n", file = con) }

d <- read.csv(paste0(B, "poma_analysis_long.csv"), stringsAsFactors = FALSE)
u <- d %>% distinct(unit, .keep_all = TRUE)

## ---- primary analysis set, built exactly as in 54_clustering_sensitivity.R
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

## ---- pair graph: nodes = pairs, an edge joins two pairs sharing a participant
pairs_id <- ua %>% group_by(newstrata) %>% summarise(ids = list(unique(id)),
                                                     .groups = "drop")
np <- nrow(pairs_id)
parent <- seq_len(np)
find <- function(i) { while (parent[i] != i) { parent[i] <<- parent[parent[i]]; i <- parent[i] }; i }
union <- function(i, j) { ri <- find(i); rj <- find(j); if (ri != rj) parent[rj] <<- ri }

owner <- list()   # participant -> first pair index seen
for (i in seq_len(np)) {
  for (pid in pairs_id$ids[[i]]) {
    k <- as.character(pid)
    if (!is.null(owner[[k]])) union(owner[[k]], i) else owner[[k]] <- i
  }
}
comp_root <- vapply(seq_len(np), find, integer(1))
comp <- match(comp_root, unique(comp_root))
pairs_id$comp <- comp

sz <- as.data.frame(table(comp))
say("=================== pair connectivity ===================")
say("pairs: ", np, " | participants in the primary set: ",
    length(unique(ua$id)))
say("connected components: ", nrow(sz),
    " | of size >1: ", sum(sz$Freq > 1),
    " | largest: ", max(sz$Freq))
say("pairs sitting in a component of size >1: ", sum(sz$Freq[sz$Freq > 1]))
say("component size distribution: ",
    paste(names(table(sz$Freq)), table(sz$Freq), sep = "->", collapse = ", "))
say("identical to participant clustering? ",
    ifelse(nrow(sz) == length(unique(ua$id)), "yes", "no"))

## ---- refit the two models with each grouping vector
zscale <- function(x) as.numeric(scale(x))
P <- ua %>% mutate(rate = slope, z = zscale(rate)) %>%
  mutate(kl = V00XRKL, strata = newstrata)
P$comp <- pairs_id$comp[match(P$strata, pairs_id$newstrata)]

ors <- function(m, nm, var) {
  s  <- summary(m)$coefficients
  cc <- exp(confint(m))
  data.frame(variance = var, model = nm, term = rownames(s),
             OR = round(exp(coef(m)), 3),
             lo = round(cc[, 1], 3), hi = round(cc[, 2], 3),
             p = signif(s[, "Pr(>|z|)"], 3), row.names = NULL)
}

f0 <- case ~ z + strata(strata)
f6 <- case ~ z + V00AGE + P01BMI + kl + WOMAC_pain + strata(strata)

tab <- rbind(
  ors(clogit(f0, data = P, method = "efron"), "M0 crude", "model-based"),
  ors(clogit(update(f0, . ~ . + cluster(id)), data = P, method = "efron"),
      "M0 crude", "cluster-robust (participant)"),
  ors(clogit(update(f0, . ~ . + cluster(comp)), data = P, method = "efron"),
      "M0 crude", "cluster-robust (connected component)"),
  ors(clogit(f6, data = P, method = "efron"), "M6 adjusted", "model-based"),
  ors(clogit(update(f6, . ~ . + cluster(id)), data = P, method = "efron"),
      "M6 adjusted", "cluster-robust (participant)"),
  ors(clogit(update(f6, . ~ . + cluster(comp)), data = P, method = "efron"),
      "M6 adjusted", "cluster-robust (connected component)"))

say("")
say("=================== exposure term, z ===================")
zz <- tab[tab$term == "z", c("variance", "model", "OR", "lo", "hi", "p")]
print(zz, row.names = FALSE)

say("")
say("SE ratio against model-based, exposure term:")
for (nm in c("M0 crude", "M6 adjusted")) {
  a <- zz[zz$model == nm & zz$variance == "model-based", ]
  for (v in c("cluster-robust (participant)", "cluster-robust (connected component)")) {
    b <- zz[zz$model == nm & zz$variance == v, ]
    say("  ", nm, " | ", v, " : ",
        sprintf("%.3f", (log(b$hi / b$lo) / 3.92) / (log(a$hi / a$lo) / 3.92)))
  }
}

write.csv(tab, paste0(B, "clogit_connectivity.csv"), row.names = FALSE)

## the counts quoted in Section 3.5 come off disk too, not out of the log
write.csv(data.frame(n_pairs = np,
                     n_components = nrow(sz),
                     n_multi_components = sum(sz$Freq > 1),
                     n_pairs_in_multi = sum(sz$Freq[sz$Freq > 1]),
                     largest_component = max(sz$Freq)),
          paste0(B, "connectivity_summary.csv"), row.names = FALSE)
close(con)
