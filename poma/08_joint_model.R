## 08 : joint model (JMbayes2), stratified by matched pair
suppressPackageStartupMessages({
  library(JMbayes2); library(nlme); library(survival); library(splines)
})
BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
d  <- readRDS(file.path(BASE, "jm_input.rds")); lg <- d$lg; sv <- d$sv
sv$strf <- factor(sv$newstrata); lg$strf <- factor(lg$newstrata)

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
if (quick) {
  sub <- unique(sv$newstrata)[1:40]
  sv <- droplevels(sv[sv$newstrata %in% sub, ]); lg <- droplevels(lg[lg$newstrata %in% sub, ])
  cat("[QUICK MODE] 40 strata\n")
}
cat("knees:", nrow(sv), " long rows:", nrow(lg), " strata:", length(unique(sv$newstrata)), "\n")
cat("events:", sum(sv$case), "\n")

## long rows must be pre-index
cat("all visits pre-index?", all(lg$year < lg$t_index_months/12), "\n\n")

## ---- longitudinal submodel: linear thinning ----
mf <- lme(y ~ year, random = ~ year | unit, data = lg,
          control = lmeControl(opt = "optim", maxIter = 200))
cat("--- longitudinal fixed effects ---\n"); print(round(fixef(mf), 4))
cat("VarCorr:\n"); print(VarCorr(mf))

## ---- event submodel: stratified by matched pair ----
sf <- coxph(Surv(t_idx, case) ~ strata(strf), data = sv)
cat("\n--- event submodel (strata only) ---\n"); print(sf)

cat("\nfunctional forms available:  value(), slope(), Delta()\n")
t0 <- Sys.time()
set.seed(2026)
jm1 <- try(jm(sf, mf, time_var = "year",
              functional_forms = ~ value(y) + slope(y),
              n_iter = if (quick) 2000 else 10000,
              n_burnin = if (quick) 500 else 2000,
              n_thin = if (quick) 5 else 10,
              cores = 1))
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs")),1), "s\n")
cat("class:", class(jm1)[1], "\n")
if (inherits(jm1, "try-error")) { print(jm1) } else {
  cat("\n=== Survival (event submodel) ===\n"); print(summary(jm1)$Survival)
  cat("\n=== Associations ===\n"); print(summary(jm1)$Associations)
  cat("\n=== Longitudinal ===\n"); print(summary(jm1)$Longitudinal)
  if (!quick) saveRDS(jm1, file.path(BASE, "jm_fit_primary.rds"))
}
