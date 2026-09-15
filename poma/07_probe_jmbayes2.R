## probe v2: strata support + speed in JMbayes2
suppressPackageStartupMessages({
  library(JMbayes2); library(nlme); library(survival); library(splines)
})
BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
unit <- read.csv(file.path(BASE, "poma_analysis_long.csv"))   # 1 row / knee
vst  <- read.csv(file.path(BASE, "poma_visits_long.csv"))     # 1 row / knee / visit

## ---- long table ----
lg <- vst[!is.na(vst$MFTC_ThCtAB_aMe),
          c("unit","years","MFTC_ThCtAB_aMe","newstrata","case",
            "V00AGE","female","P01BMI","V00XRKL")]
lg <- merge(lg, unit[, c("unit","t_index_months")], by="unit", all.x=TRUE)
names(lg)[names(lg)=="years"] <- "year"
names(lg)[names(lg)=="MFTC_ThCtAB_aMe"] <- "y"
lg$strf  <- factor(lg$newstrata)
lg$age5  <- lg$V00AGE/5
lg$bmi5  <- lg$P01BMI/5
lg       <- lg[!is.na(lg$y) & !is.na(lg$t_index_months), ]

## ---- survival table (1 row / knee) ----
sv <- lg[!duplicated(lg$unit), c("unit","t_index_months","case","newstrata",
                                 "age5","bmi5","female","V00XRKL","strf")]
sv$t_idx <- sv$t_index_months/12

cat("long rows:", nrow(lg), " units:", nrow(sv),
    " strata:", length(unique(sv$newstrata)), "\n")

## keep strata with exactly 2 units
tb   <- table(sv$newstrata)
keep <- as.integer(names(tb)[tb == 2])
sv   <- sv[sv$newstrata %in% keep, ]
lg   <- lg[lg$newstrata %in% keep, ]
cat("balanced strata:", length(keep), "units:", nrow(sv), "long rows:", nrow(lg), "\n\n")

saveRDS(list(lg=lg, sv=sv), file.path(BASE, "jm_input.rds"))

## ---- Test 1: strata in coxph ----
cat("--- Test 1: coxph strata(strf) on full data ---\n")
cf1 <- try(coxph(Surv(t_idx, case) ~ strata(strf) + age5 + bmi5 + V00XRKL, data=sv))
cat("class:", class(cf1)[1], "\n")
if (!inherits(cf1, "try-error")) print(round(summary(cf1)$coef, 4)) else print(cf1)

cat("\n--- Test 1b: same but unstratified ---\n")
cf1b <- coxph(Surv(t_idx, case) ~ age5 + bmi5 + V00XRKL, data=sv)
print(round(summary(cf1b)$coef, 4))

## ---- longitudinal fit ----
lm2 <- lme(y ~ ns(year, 2), random = ~ ns(year, 2) | unit, data=lg,
           control=lmeControl(opt="optim"))
cat("\n--- lme ---\n"); print(round(fixef(lm2), 4))

## ---- Test 2: jm() with stratified survFit, small subset, tiny MCMC ----
sub <- keep[1:40]
svs <- droplevels(sv[sv$newstrata %in% sub, ])
lvs <- droplevels(lg[lg$newstrata %in% sub, ])
cat("\nsubset units:", nrow(svs), " long rows:", nrow(lvs), "\n")
cfs <- coxph(Surv(t_idx, case) ~ strata(strf) + age5 + bmi5 + V00XRKL, data=svs)
lms <- lme(y ~ ns(year, 2), random = ~ ns(year, 2) | unit, data=lvs,
           control=lmeControl(opt="optim"))

t0 <- Sys.time()
jm1 <- try(jm(cfs, lms, time_var="year", n_iter=2000, n_burnin=500, n_thin=5,
              cores=1, verbose=FALSE))
cat("Test 2 (strata) elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs")),1), "s  class:", class(jm1)[1], "\n")
if (!inherits(jm1,"try-error")) { print(summary(jm1)$Survival); print(summary(jm1)$Associations) } else print(jm1)

## ---- Test 3: same subset, unstratified ----
cfs2 <- coxph(Surv(t_idx, case) ~ age5 + bmi5 + V00XRKL, data=svs)
t0 <- Sys.time()
jm2 <- try(jm(cfs2, lms, time_var="year", n_iter=2000, n_burnin=500, n_thin=5,
              cores=1, verbose=FALSE))
cat("\nTest 3 (no strata) elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs")),1), "s  class:", class(jm2)[1], "\n")
if (!inherits(jm2,"try-error")) { print(summary(jm2)$Survival); print(summary(jm2)$Associations) } else print(jm2)
