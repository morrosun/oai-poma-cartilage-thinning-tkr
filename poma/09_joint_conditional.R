# ============================================================================
# 09 : Bayesian joint model for the POMA nested case-control design
#
#   longitudinal : y_ij = b0 + b1*t_ij + u0_i + u1_i*t_ij + e_ij ,  e ~ N(0,s2)
#                  (u0_i,u1_i)' ~ N(0, D)
#   event        : WITHIN-PAIR CONDITIONAL likelihood
#                  P(case in pair s) = plogis( aV*[d0_s + d1_s*t*_s] + aS*d1_s )
#                  where d0_s, d1_s are the case-minus-control latent contrasts
#
#   KEY POINT: both knees of a pair share the same index time t*, so the
#   baseline hazard AND every within-pair-constant term (b0, b1, and the
#   matching covariates) cancel exactly.  The event submodel therefore needs
#   NO baseline hazard -> this is the correct joint model for this design,
#   and it is what a 224-stratum Cox joint model reduces to.
#
#   sampler: Gibbs for (beta, sigma2, D) + independence MH for the random
#            effects (proposal = longitudinal conditional) + adaptive RW MH
#            for (aV, aS).  Fully vectorised.
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"

## ---------------------------------------------------------------- data ----
prep <- function(min_visits = 1) {
  d  <- readRDS(file.path(BASE, "jm_input.rds"))
  lg <- d$lg; sv <- d$sv
  lg$y <- as.numeric(lg$y); lg$year <- as.numeric(lg$year)
  lg <- lg[!is.na(lg$y) & !is.na(lg$year), ]
  sv$pair <- as.integer(factor(sv$newstrata))
  lg$pair <- sv$pair[match(lg$unit, sv$unit)]
  lg <- lg[!is.na(lg$pair), ]
  for (iter in 1:6) {
    tb <- table(sv$pair); sv <- sv[sv$pair %in% as.integer(names(tb)[tb == 2]), ]
    lg <- lg[lg$pair %in% sv$pair, ]
    if (min_visits > 1) {
      nv <- table(lg$unit)
      lg <- lg[lg$unit %in% names(nv)[nv >= min_visits], ]
      sv <- sv[sv$unit %in% lg$unit, ]
    }
    tb <- table(sv$pair); sv <- sv[sv$pair %in% as.integer(names(tb)[tb == 2]), ]
    lg <- lg[lg$pair %in% sv$pair, ]
    if (all(table(sv$pair) == 2) && all(table(lg$unit) >= min_visits)) break
  }
  sv <- sv[order(sv$pair, -sv$case), ]
  list(lg = lg, sv = sv)
}

build <- function(lg, sv) {
  n <- nrow(sv); u <- sv$unit
  idx <- match(lg$unit, u)
  M <- cbind(y = lg$y, t = lg$year, ty = lg$year * lg$y,
             tt = lg$year^2, yy = lg$y^2, one = 1)
  R <- rowsum(M, group = idx)
  stopifnot(nrow(R) == n)
  rownames(R) <- NULL
  cas <- which(sv$case == 1); ctl <- which(sv$case == 0)
  stopifnot(length(cas) == length(ctl), all(sv$pair[cas] == sv$pair[ctl]))
  list(R = R, n = n, cas = cas, ctl = ctl, tstar = sv$t_idx,
       pair = sv$pair, case = sv$case, unit = u, newstrata = sv$newstrata)
}

## ------------------------------------------------------------- sampler ----
log1pexp  <- function(x) ifelse(x > 30, x + log1p(exp(-x)), log1p(exp(pmin(x, 30))))
logplogis <- function(x) -log1pexp(-x)

fit_jm <- function(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 1,
                   aS_fixed = NULL, prior_sd = 20, verbose = TRUE) {
  set.seed(seed)
  R <- dt$R; n <- dt$n; cas <- dt$cas; ctl <- dt$ctl
  tstar_i <- dt$tstar
  Sy <- R[, "y"]; Sty <- R[, "ty"]; St <- R[, "t"]; Stt <- R[, "tt"]
  Syy <- R[, "yy"]; n_i <- R[, "one"]
  partner <- integer(n); partner[cas] <- ctl; partner[ctl] <- cas
  sgn <- ifelse(dt$case == 1, 1, -1)
  tstar_c <- tstar_i[cas]

  b0 <- rep(0, n); b1 <- rep(0, n)
  beta0 <- 2.9; beta1 <- -0.08
  s2 <- 0.02
  D  <- matrix(c(1.2, 0.02, 0.02, 0.01), 2, 2)
  D0 <- matrix(c(1.2, 0, 0, 0.01), 2, 2); df0 <- 4
  aV <- 0; aS <- if (is.null(aS_fixed)) 0 else aS_fixed
  cstep <- 1.1

  keep <- seq(BURN + 1, NITER, by = THIN)
  nk <- length(keep)
  out <- list(aV = numeric(nk), aS = numeric(nk), beta0 = numeric(nk),
              beta1 = numeric(nk), s2 = numeric(nk), D = matrix(0, nk, 3),
              SDval = numeric(nk), SDslo = numeric(nk))
  eta_store <- matrix(0, nk, length(cas))
  acc_b <- 0; acc_b_tot <- 0; acc_a <- 0; acc_a_win <- 0; k <- 0

  for (it in 1:NITER) {
    ## ---- 1. random effects: independence MH, proposal = longitudinal cond.
    Dinv <- solve(D)
    a11 <- n_i / s2 + Dinv[1, 1]; a12 <- St / s2 + Dinv[1, 2]
    a22 <- Stt / s2 + Dinv[2, 2]
    det <- a11 * a22 - a12^2
    V11 <- a22 / det; V22 <- a11 / det; V12 <- -a12 / det
    r0 <- (Sy  - (n_i * beta0 + St  * beta1)) / s2
    r1 <- (Sty - (St  * beta0 + Stt * beta1)) / s2
    m0 <- V11 * r0 + V12 * r1
    m1 <- V12 * r0 + V22 * r1
    L11 <- sqrt(pmax(V11, 1e-12)); L21 <- V12 / L11
    L22 <- sqrt(pmax(V22 - L21^2, 1e-12))
    z0 <- rnorm(n); z1 <- rnorm(n)
    b0p <- m0 + L11 * z0
    b1p <- m1 + L21 * z0 + L22 * z1

    d0n <- sgn * (b0p - b0[partner]); d1n <- sgn * (b1p - b1[partner])
    d0o <- sgn * (b0  - b0[partner]); d1o <- sgn * (b1  - b1[partner])
    e_n <- aV * (d0n + d1n * tstar_i) + aS * d1n
    e_o <- aV * (d0o + d1o * tstar_i) + aS * d1o
    acc_i <- runif(n) < exp(logplogis(e_n) - logplogis(e_o))
    acc_b <- acc_b + sum(acc_i); acc_b_tot <- acc_b_tot + n
    b0[acc_i] <- b0p[acc_i]; b1[acc_i] <- b1p[acc_i]

    ## ---- 2. fixed effects (Gibbs; they cancel in the event term) ----
    Ab0 <- n_i * b0 + St * b1
    Ab1 <- St * b0 + Stt * b1
    P  <- matrix(c(sum(n_i), sum(St), sum(St), sum(Stt)), 2, 2)
    q  <- c(sum(Sy) - sum(Ab0), sum(Sty) - sum(Ab1))
    bb <- solve(P, q); beta0 <- bb[1]; beta1 <- bb[2]

    ## ---- 3. residual variance (Gibbs, InvGamma(0.01, 0.01)) ----
    g0 <- beta0 + b0; g1 <- beta1 + b1
    SSE <- sum(Syy - 2 * (g0 * Sy + g1 * Sty) + g0^2 * n_i +
                 2 * g0 * g1 * St + g1^2 * Stt)
    s2 <- 1 / rgamma(1, shape = 0.01 + sum(n_i) / 2, rate = 0.01 + SSE / 2)

    ## ---- 4. D (Gibbs, inverse-Wishart) ----
    B <- cbind(b0, b1)
    W <- rWishart(1, df0 + n, solve(D0 + t(B) %*% B))[, , 1]
    D <- solve(W)

    ## ---- 5. association parameters: Laplace-proposal MH (flat prior) ----
    ##   the posterior of (aV,aS) is strongly correlated; a plain random walk
    ##   mixes badly (verified against a known-truth simulation), so we build a
    ##   proposal that targets the local mode with the local curvature and
    ##   correct it with the full Hastings ratio.
    d0 <- b0[cas] - b0[ctl]; d1 <- b1[cas] - b1[ctl]
    u <- d0 + d1 * tstar_c; v <- d1
    Xa <- if (is.null(aS_fixed)) cbind(u, v) else cbind(u)
    a_cur <- if (is.null(aS_fixed)) c(aV, aS) else aV
    ll_a <- function(a) sum(logplogis(drop(Xa %*% a)))
    ll_pr <- function(a) -0.5 * sum(a^2) / prior_sd^2     # N(0, prior_sd^2)
    lp_prop <- function(a) {
      e <- drop(Xa %*% a); p <- plogis(e); w <- p * (1 - p)
      g <- drop(crossprod(Xa, 1 - p))
      H <- -crossprod(Xa * w, Xa)
      Si <- try(solve(-H), silent = TRUE)
      if (inherits(Si, "try-error")) Si <- diag(length(a)) * 0.05
      det <- determinant(Si, logarithm = TRUE)
      if (!is.finite(det$modulus) || det$modulus <= -700) Si <- diag(length(a)) * 0.05
      list(mu = a + drop(Si %*% g), Sig = Si * cstep^2)
    }
    dmv <- function(x, mu, S) {
      d <- length(x); Si <- solve(S)
      -0.5 * (d * log(2 * pi) + as.numeric(determinant(S, logarithm = TRUE)$modulus) +
                drop(t(x - mu) %*% Si %*% (x - mu)))
    }
    q1 <- lp_prop(a_cur)
    Ch <- try(chol(q1$Sig), silent = TRUE)
    if (inherits(Ch, "try-error")) Ch <- chol(diag(length(a_cur)) * 0.05)
    a_prop <- q1$mu + drop(t(Ch) %*% rnorm(length(a_cur)))
    q2 <- lp_prop(a_prop)
    lr <- (ll_a(a_prop) + ll_pr(a_prop)) - (ll_a(a_cur) + ll_pr(a_cur)) +
          dmv(a_cur, q2$mu, q2$Sig) - dmv(a_prop, q1$mu, q1$Sig)
    if (log(runif(1)) < lr) {
      aV <- a_prop[1]
      if (is.null(aS_fixed)) aS <- a_prop[2] else aS <- aS_fixed
      acc_a <- acc_a + 1
    }
    acc_a_win <- acc_a_win + 1
    if (acc_a_win == 500) { acc_a <- 0; acc_a_win <- 0 }

    ## ---- store ----
    if (k < nk && it == keep[k + 1]) {
      k <- k + 1
      out$aV[k] <- aV; out$aS[k] <- aS
      out$beta0[k] <- beta0; out$beta1[k] <- beta1; out$s2[k] <- s2
      out$D[k, ] <- c(D[1, 1], D[1, 2], D[2, 2])
      out$SDval[k] <- sd(b0 + b1 * tstar_i)
      out$SDslo[k] <- sd(b1)
      eta_store[k, ] <- aV * u + aS * v
    }
    if (verbose && it %% 10000 == 0)
      cat(sprintf("  iter %6d  aV=%+.3f aS=%+.3f  beta1=%.4f  s2=%.4f\n",
                  it, aV, aS, beta1, s2))
  }
  list(draws = out, eta = eta_store, acc_b = acc_b / acc_b_tot,
       n = n, npair = length(cas), tstar = tstar_i, keep = keep)
}

waic_pairs <- function(fit) {
  E <- fit$eta
  lppd <- sum(log(colMeans(exp(E))))
  pw   <- sum(apply(E, 2, var))
  c(waic = -2 * (lppd - pw), lppd = lppd, p_waic = pw)
}

summ <- function(x) {
  q <- quantile(x, c(0.025, 0.5, 0.975))
  sprintf("%.3f (%.3f, %.3f)", q[2], q[1], q[3])
}

report <- function(fit, label) {
  d <- fit$draws
  cat(sprintf("\n================ %s ================\n", label))
  cat(sprintf("pairs %d | knees %d | draws %d | MH acc (random effects) %.2f\n",
              fit$npair, fit$n, length(d$aV), fit$acc_b))
  cat("\n-- longitudinal submodel --\n")
  cat("  beta0 (mm, at t=0)       :", summ(d$beta0), "\n")
  cat("  beta1 (mm/yr)            :", summ(d$beta1), "\n")
  cat("  residual SD              :", summ(sqrt(d$s2)), "\n")
  cat("  SD(u0)                   :", summ(sqrt(d$D[, 1])), "\n")
  cat("  SD(u1) (mm/yr)           :", summ(sqrt(d$D[, 3])), "\n")
  cat("  corr(u0,u1)              :", summ(d$D[, 2] / sqrt(d$D[, 1] * d$D[, 3])), "\n")
  cat("\n-- event submodel (within-pair conditional) --\n")
  cat("  aV (per mm latent level) :", summ(d$aV), "\n")
  cat("  aS (per mm/yr latent)    :", summ(d$aS), "\n")
  cat(sprintf("  P(aV<0)=%.3f   P(aS<0)=%.3f\n", mean(d$aV < 0), mean(d$aS < 0)))
  cat("  SD latent level (mm)     :", summ(d$SDval), "\n")
  cat("  SD latent slope (mm/yr)  :", summ(d$SDslo), "\n")
  orV <- exp(-d$aV * d$SDval); orS <- exp(-d$aS * d$SDslo)
  orV1 <- exp(-d$aV); orS01 <- exp(-0.1 * d$aS)
  cat("\n-- OR per 1 SD (positive = FASTER thinning / THINNER level) --\n")
  cat("  OR per 1 SD thinner latent level :", summ(orV), "\n")
  cat("  OR per 1 SD faster latent slope  :", summ(orS), "\n")
  cat("-- OR on natural dosing --\n")
  cat("  OR per 1 mm thinner latent level :", summ(orV1), "\n")
  cat("  OR per 0.1 mm/yr faster thinning :", summ(orS01), "\n")
  invisible(list(orV = orV, orS = orS, orV1 = orV1, orS01 = orS01))
}

## ------------------------------------------------------ simulation check --
sim_check <- function(npair = 224, seed = 7) {
  set.seed(seed)
  D <- matrix(c(1.2, 0.02, 0.02, 0.01), 2, 2)
  true <- c(beta0 = 2.9, beta1 = -0.08, s2 = 0.021, aV = -0.45, aS = -6)
  cat(sprintf("\n[SIM] true aV=%.2f  aS=%.2f\n", true["aV"], true["aS"]))
  ## generate 2*npair knees
  n <- 2 * npair
  B <- matrix(rnorm(n * 2), n, 2) %*% chol(D)
  nv <- sample(1:5, n, TRUE)
  tstar <- sample(c(1, 2, 3, 4, 5), npair, TRUE)
  rows <- do.call(rbind, lapply(seq_len(n), function(i) {
    tt <- sort(sample(0:4, nv[i]))
    data.frame(unit = i, year = tt,
               y = true["beta0"] + true["beta1"] * tt + B[i, 1] + B[i, 2] * tt +
                   rnorm(length(tt), 0, sqrt(true["s2"])))
  }))
  ## assign case/control so that the conditional likelihood is the target
  cas <- seq(1, n, by = 2); ctl <- cas + 1
  pair <- rep(seq_len(npair), each = 2)
  ts <- rep(tstar, each = 2)
  eta <- true["aV"] * ((B[cas, 1] - B[ctl, 1]) + (B[cas, 2] - B[ctl, 2]) * tstar) +
         true["aS"] * (B[cas, 2] - B[ctl, 2])
  ## assign which member of each pair is the case -> conditional likelihood
  ## is exactly plogis(eta) and every pair is retained
  is_case <- rbinom(npair, 1, plogis(eta))
  sv <- data.frame(unit = seq_len(n), case = 0, pair = pair, t_idx = ts,
                   newstrata = pair)
  sv$case[cas] <- is_case
  sv$case[ctl] <- 1 - is_case
  lg <- merge(rows, sv[, c("unit", "pair", "t_idx")], by = "unit")
  sv <- sv[order(sv$pair, -sv$case), ]
  M <- cbind(y = lg$y, t = lg$year, ty = lg$year * lg$y,
             tt = lg$year^2, yy = lg$y^2, one = 1)
  idx <- match(lg$unit, sv$unit)
  R <- rowsum(M, group = idx); rownames(R) <- NULL
  dt <- list(R = R, n = nrow(sv), cas = which(sv$case == 1),
             ctl = which(sv$case == 0), tstar = sv$t_idx, pair = sv$pair,
             case = sv$case, unit = sv$unit)
  f <- fit_jm(dt, NITER = 20000, BURN = 4000, THIN = 5, seed = 3, verbose = FALSE)
  d <- f$draws
  cat(sprintf("  recovered aV = %s   (true %.2f)\n", summ(d$aV), true["aV"]))
  cat(sprintf("  recovered aS = %s   (true %.2f)\n", summ(d$aS), true["aS"]))
  cat(sprintf("  recovered b1 = %s   (true %.3f)\n", summ(d$beta1), true["beta1"]))
  invisible(f)
}

## ================================================================= main ====
cat("========================================================================\n")
cat("POMA joint model, within-pair conditional event likelihood\n")
cat("========================================================================\n")

p  <- prep(1)
dt <- build(p$lg, p$sv)
cat(sprintf("data: %d knees, %d complete pairs, %d longitudinal visits\n",
            dt$n, length(dt$cas), nrow(dt$R)))
cat(sprintf("mean visits/knee: %.2f\n", nrow(p$lg) / dt$n))

fit_full <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 2026)
report(fit_full, "MODEL 1 : latent level + latent slope")
w1 <- waic_pairs(fit_full)
cat(sprintf("  WAIC (pair level) = %.1f   (lppd %.1f, p_waic %.1f)\n",
            w1["waic"], w1["lppd"], w1["p_waic"]))

fit_val <- fit_jm(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 2027, aS_fixed = 0)
report(fit_val, "MODEL 2 : latent level only (aS = 0)")
w2 <- waic_pairs(fit_val)
cat(sprintf("  WAIC (pair level) = %.1f   (lppd %.1f, p_waic %.1f)\n",
            w2["waic"], w2["lppd"], w2["p_waic"]))
cat(sprintf("\nDelta WAIC (level-only minus level+slope) = %.1f  (>0 favours level+slope)\n",
            w2["waic"] - w1["waic"]))
cat("NOTE: this WAIC is computed on the pair-conditional likelihood in which the\n")
cat("      random effects act as parameters, so p_waic approaches the number of\n")
cat("      pairs and the criterion is unstable.  The formal comparison of\n")
cat("      'level vs rate' is the conditional-logistic LR test in step 01/06.\n")

## ---- Savage-Dickey Bayes factors from the joint model (prior N(0,20^2)) ----
sd_bf01 <- function(x, prior_sd = 20) {
  d <- density(x, n = 4096)
  approx(d$x, d$y, xout = 0, rule = 2)$y / dnorm(0, 0, prior_sd)
}
bfS <- sd_bf01(fit_full$draws$aS)
bfV <- sd_bf01(fit_full$draws$aV)
cat("\n-- Savage-Dickey Bayes factors (prior aV,aS ~ N(0,20^2), flat otherwise) --\n")
cat(sprintf("  H0: aS = 0 (rate adds nothing) : BF01 = %.3g   -> BF10 = %.3g\n", bfS, 1 / bfS))
cat(sprintf("  H0: aV = 0 (level adds nothing): BF01 = %.3g   -> BF10 = %.3g\n", bfV, 1 / bfV))
cat("  interpret BF10 > 10 as strong, > 100 as decisive evidence.\n")

## sensitivity: only knees with >= 3 visits
p3 <- prep(3); dt3 <- build(p3$lg, p3$sv)
cat(sprintf("\n[>=3 visits] %d knees, %d pairs\n", dt3$n, length(dt3$cas)))
fit_s <- fit_jm(dt3, NITER = 30000, BURN = 5000, THIN = 5, seed = 2028, verbose = FALSE)
report(fit_s, "MODEL 3 (sensitivity) : >= 3 visits per knee")

## validate the longitudinal part against nlme
suppressPackageStartupMessages(library(nlme))
lme_fit <- lme(y ~ year, random = ~ year | unit, data = p$lg,
               control = lmeControl(opt = "optim"))
cat("\n-- validation: longitudinal submodel vs nlme::lme --\n")
vc <- suppressWarnings(as.numeric(VarCorr(lme_fit)))   # col-wise: [4:6] = SDs
cat("  nlme  beta  :", sprintf("%+.4f %+.4f", fixef(lme_fit)[1], fixef(lme_fit)[2]), "\n")
cat("  Bayes beta  :", sprintf("%+.4f %+.4f", median(fit_full$draws$beta0),
                                 median(fit_full$draws$beta1)), "\n")
cat("  nlme  SD(u0, u1, resid):",
    sprintf("%.4f %.4f %.4f", vc[4], vc[5], vc[6]), "\n")
cat("  Bayes SD(u0, u1, resid):",
    sprintf("%.4f %.4f %.4f", median(sqrt(fit_full$draws$D[, 1])),
            median(sqrt(fit_full$draws$D[, 3])),
            median(sqrt(fit_full$draws$s2))), "\n")

sim_check()

## ------------------------------------------------------------- figure -----
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
df <- data.frame(
  or01 = exp(-0.1 * fit_full$draws$aS),        # per 0.1 mm/yr faster thinning
  orSD = exp(-fit_full$draws$aS * fit_full$draws$SDslo)
)
p1 <- ggplot(df, aes(or01)) +
  geom_density(fill = "#C44E52", alpha = .35, colour = "#8C2D30") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = median(df$or01), colour = "#8C2D30") +
  scale_x_log10() + theme_bw(base_size = 11) +
  labs(x = "OR per 0.1 mm/yr faster latent thinning", y = "posterior density",
       title = "Joint model: rate of change")
p2 <- ggplot(data.frame(x = exp(-fit_full$draws$aV)),
             aes(x)) +
  geom_density(fill = "#4C72B0", alpha = .35, colour = "#2F4A75") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_x_log10() + theme_bw(base_size = 11) +
  labs(x = "OR per 1 mm thinner latent level", y = "posterior density",
       title = "Joint model: level")
ggsave(file.path(BASE, "fig_jm_posterior.png"),
       p1 + p2, width = 9, height = 3.6, dpi = 200)

est <- data.frame(
  quantity = c("OR per 0.1 mm/yr faster thinning (latent)",
               "OR per 1 SD faster thinning (latent)",
               "OR per 1 mm thinner level (latent)",
               "OR per 1 SD thinner level (latent)",
               "residual SD (mm)", "SD of latent slope (mm/yr)"),
  median = c(median(exp(-0.1 * fit_full$draws$aS)),
             median(exp(-fit_full$draws$aS * fit_full$draws$SDslo)),
             median(exp(-fit_full$draws$aV)),
             median(exp(-fit_full$draws$aV * fit_full$draws$SDval)),
             median(sqrt(fit_full$draws$s2)),
             median(fit_full$draws$SDslo)),
  lo = c(quantile(exp(-0.1 * fit_full$draws$aS), .025),
         quantile(exp(-fit_full$draws$aS * fit_full$draws$SDslo), .025),
         quantile(exp(-fit_full$draws$aV), .025),
         quantile(exp(-fit_full$draws$aV * fit_full$draws$SDval), .025),
         quantile(sqrt(fit_full$draws$s2), .025),
         quantile(fit_full$draws$SDslo, .025)),
  hi = c(quantile(exp(-0.1 * fit_full$draws$aS), .975),
         quantile(exp(-fit_full$draws$aS * fit_full$draws$SDslo), .975),
         quantile(exp(-fit_full$draws$aV), .975),
         quantile(exp(-fit_full$draws$aV * fit_full$draws$SDval), .975),
         quantile(sqrt(fit_full$draws$s2), .975),
         quantile(fit_full$draws$SDslo, .975)))
write.csv(est, file.path(BASE, "jm_estimates.csv"), row.names = FALSE)
cat("\n-> jm_estimates.csv / fig_jm_posterior.png\n")

saveRDS(list(full = fit_full, value = fit_val, sens3 = fit_s,
             waic = rbind(full = w1, value = w2)),
        file.path(BASE, "jm_conditional_fit.rds"))
cat("-> jm_conditional_fit.rds\n")
