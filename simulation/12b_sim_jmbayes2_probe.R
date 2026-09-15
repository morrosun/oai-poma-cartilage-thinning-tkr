# ============================================================================
# 12b : controlled grid -- when does a stratified joint model fail?
#
#   Competing explanations for the real-data failure
#       (a) baseline-hazard degeneracy from the shared index visit
#           -> a failure of the SAMPLER covariance, message "chol()"
#       (b) within-pair-constant covariates in a stratified Cox submodel,
#           which are unidentifiable and return NA coefficients
#           -> a failure at coxph()/initialisation, different message
#
#   The grid separates them: REAL data and SIMULATED data (known truth), each
#   with and without the three matched (pair-constant) covariates.  Every fit is
#   capped at TIME_LIM s so that a non-converging chain cannot block the run.
#
#   Env: SIM_NPAIR (40) SIM_NITER (300) SIM_TL (150)
# ============================================================================
source("D:/BaiduSyncdisk/OAI/Scripts/sim/sim_core.R")
suppressPackageStartupMessages({library(JMbayes2); library(splines)})

BASE <- "D:/BaiduSyncdisk/OAI"
OUT  <- file.path(BASE, "Analysis/sim")
design <- read_design(file.path(OUT, "sim_visit_schedule.csv"))
CAL <- readRDS(file.path(OUT, "sim_as.rds")); aS <- CAL$a; aS <- CAL$aS; CONST <- CAL$const

NPAIR <- as.integer(Sys.getenv("SIM_NPAIR", "40"))
NITER <- as.integer(Sys.getenv("SIM_NITER", "300"))
TL    <- as.numeric(Sys.getenv("SIM_TL", "150"))

LOG <- list(); P <- function(...) { s <- sprintf(...); cat(s, "\n"); LOG[[length(LOG)+1]] <<- s }
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

get_assoc <- function(f) {
  a <- try(summary(f)$Associations, silent = TRUE)
  if (inherits(a, "try-error")) return(NA_real_)
  v <- suppressWarnings(as.numeric(a[[1]][1, 1])); if (length(v)) v[1] else NA_real_
}

## ---- assemble long + survival tables from a long data frame ---------------
## covtype: "none" | "pair" (constant within a pair, as age/BMI/KL after
##          matching) | "knee" (varies within a pair)
tables_from <- function(lg, covtype) {
  lg$year <- lg$t; lg$strf <- factor(lg$pair)
  if (covtype == "pair") {
    set.seed(4242)
    zv <- setNames(rnorm(nlevels(lg$strf)), levels(lg$strf))
    lg$age5 <- 12 + 2 * zv[as.character(lg$strf)]
    lg$bmi5 <- 5.6 + 0.5 * zv[as.character(lg$strf)]^2
    lg$kl   <- round(1.5 + abs(zv[as.character(lg$strf)]))
  } else if (covtype == "knee") {
    set.seed(4243)
    u <- unique(lg$unit)
    zu <- setNames(rnorm(length(u)), as.character(u))
    lg$age5 <- 12 + 2 * zu[as.character(lg$unit)]
    lg$bmi5 <- 5.6 + 0.5 * zu[as.character(lg$unit)]^2
    lg$kl   <- round(1.5 + abs(zu[as.character(lg$unit)]))
  }
  cols <- c("unit", "tstar", "case", "pair", "strf",
            if (covtype != "none") c("age5", "bmi5", "kl"))
  sv <- lg[!duplicated(lg$unit), cols]
  names(sv)[names(sv) == "tstar"] <- "t_idx"
  list(lg = lg, sv = sv)
}

attempt <- function(lg, use_strata, covtype) {
  tb <- tables_from(lg, covtype); lgl <- tb$lg; sv <- tb$sv
  lm2 <- try(lme(y ~ ns(year, 2), random = ~ ns(year, 2) | unit, data = lgl,
                 control = lmeControl(opt = "optim")), silent = TRUE)
  if (inherits(lm2, "try-error")) return(list(stage = "lme", msg = "lme failed", sec = NA))
  rhs <- c(if (use_strata) "strata(strf)", if (covtype != "none") c("age5", "bmi5", "kl"))
  if (!length(rhs)) rhs <- "1"
  cf <- try(coxph(as.formula(paste("Surv(t_idx, case) ~", paste(rhs, collapse = " + "))),
                  data = sv), silent = TRUE)
  if (inherits(cf, "try-error")) return(list(stage = "coxph", msg = as.character(cf)[1], sec = NA))
  nas <- sum(is.na(coef(cf)))
  t0 <- Sys.time()
  setTimeLimit(elapsed = TL, transient = TRUE)
  fit <- try(jm(cf, lm2, time_var = "year", n_iter = NITER, n_burnin = 100,
                n_thin = 2, cores = 1), silent = TRUE)
  setTimeLimit(elapsed = Inf, transient = TRUE)
  sec <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  if (inherits(fit, "try-error"))
    return(list(stage = "jm-FAIL", msg = gsub("\\s+", " ", as.character(fit))[1],
                sec = sec, na_coef = nas))
  list(stage = "jm-ok", msg = "", sec = sec, na_coef = nas, assoc = get_assoc(fit))
}

## ---- REAL data ------------------------------------------------------------
P("=== 12b  controlled probe:  %d pairs, %d iterations, cap %.0f s ===", NPAIR, NITER, TL)
suppressPackageStartupMessages(library(survival))
U <- read.csv(file.path(BASE, "Analysis/poma_pilot/poma_analysis_long.csv"))
V <- read.csv(file.path(BASE, "Analysis/poma_pilot/poma_visits_long.csv"))
L <- read.csv(file.path(BASE, "Analysis/poma_pilot/poma_long.csv"))
PRIM <- "cMFTC_ThCtAB_aMe"
nv <- table(L$unit[!is.na(L[[PRIM]])]); ok <- names(nv)[nv >= 2]
U <- U[U$unit %in% ok, ]; tb <- table(U$newstrata)
U <- U[U$newstrata %in% names(tb)[tb == 2], ]
keep <- unique(U$newstrata)[seq_len(min(NPAIR, length(unique(U$newstrata))))]
U <- U[U$newstrata %in% keep, ]
V <- V[V$unit %in% U$unit & !is.na(V[[PRIM]]), ]
lgR <- data.frame(unit = V$unit, pair = V$newstrata, t = V$years,
                  y = V[[PRIM]], case = V$case)
lgR <- merge(lgR, U[, c("unit", "t_index_months")], by = "unit", all.x = TRUE)
lgR$tstar <- lgR$t_index_months / 12
P("\n--- REAL POMA data, %d pairs ---", length(unique(lgR$pair)))
for (wc in c("none", "pair", "knee")) for (us in c(TRUE, FALSE)) {
  z <- attempt(lgR, us, wc)
  P("  strata=%-5s covariates=%-5s -> %-9s [%s s] NA-coef=%s %s",
    us, wc, z$stage, format(z$sec), format(z$na_coef %||% NA), substr(z$msg, 1, 80))
  if (z$stage == "jm-ok" && is.finite(z$assoc %||% NA)) P("        association = %.3f", z$assoc)
}

## ---- SIMULATED data ------------------------------------------------------
res <- list()
for (r in 1:2) {
  d <- gen_data(design, NPAIR, CONST$tau, CONST$mu1, CONST$b0mu, CONST$b0sd,
                CONST$sigma_e, aS = aS, seed = 8000 + r)
  lgS <- data.frame(unit = d$unit, pair = d$pair, t = d$t, y = d$y,
                    case = d$case, tstar = d$tstar)
  cond <- clogit_fit(knee_slopes(d), "slope_eblup")[["beta"]]
  P("\n--- SIMULATED data, replicate %d, %d pairs (truth aS = %.3f, conditional = %.3f) ---",
    r, NPAIR, aS, cond)
  for (wc in c("none", "pair", "knee")) for (us in c(TRUE, FALSE)) {
    z <- attempt(lgS, us, wc)
    P("  strata=%-5s covariates=%-5s -> %-9s [%s s] NA-coef=%s %s",
      us, wc, z$stage, format(z$sec), format(z$na_coef %||% NA), substr(z$msg, 1, 80))
    if (z$stage == "jm-ok" && is.finite(z$assoc %||% NA)) P("        association = %.3f", z$assoc)
  }
  res[[r]] <- list(cond = cond)
}
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

saveRDS(list(log = LOG, real_pairs = length(unique(lgR$pair))),
        file.path(OUT, "sim1b_jmbayes2_probe.rds"))
writeLines(unlist(LOG), file.path(OUT, "sim1b_jmbayes2_probe_log.txt"))
P("\nwritten -> Analysis/sim/sim1b_jmbayes2_probe.rds")
