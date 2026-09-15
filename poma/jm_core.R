# ============================================================================
# jm_core.R : shared sampler for the POMA within-pair conditional joint model
#
#   longitudinal : y_ij = b0 + b1*t_ij + u0_i + u1_i*t_ij + e_ij , e ~ N(0,s2)
#   event        : WITHIN-PAIR CONDITIONAL likelihood
#                  P(case in pair s) = plogis( aV*[d0_s + d1_s*t*_s] + aS*d1_s )
#
#   Both knees of a pair share the index visit t*, so the baseline hazard and
#   every within-pair-constant term (b0, b1, matching covariates) cancel.  This
#   is why no baseline hazard is required and why a 224-stratum Cox joint model
#   is not identifiable (verified: JMbayes2 chol() fails / never converges).
#
#   uses only: dt$R dt$n dt$cas dt$ctl dt$tstar dt$case
# ============================================================================

## numerically stable log(1+exp(x)) and log(p / (1-p))
log1pexp  <- function(x) ifelse(x > 30, x + log1p(exp(-x)), log1p(exp(pmin(x, 30))))
logplogis <- function(x) -log1pexp(-x)

## ------------------------------------------------------------------ build --
## lg : long table with columns unit, year, y ;  sv : one row per unit
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

## prune to complete pairs (two knees each) and optionally >= min_visits
prune_pairs <- function(lg, sv, min_visits = 1) {
  sv$pair <- as.integer(factor(sv$newstrata))
  lg$pair <- sv$pair[match(lg$unit, sv$unit)]
  lg <- lg[!is.na(lg$pair), ]
  for (iter in 1:8) {
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

## ---------------------------------------------------------------- sampler --
## Zd : optional (npair x m) matrix of WITHIN-PAIR covariate differences
##      (case minus control).  Only covariates that vary within a matched pair
##      are identifiable; within-pair-constant ones cancel like b0 and b1.
fit_jm <- function(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 1,
                   aS_fixed = NULL, prior_sd = 20, Zd = NULL, verbose = TRUE) {
  set.seed(seed)
  R <- dt$R; n <- dt$n; cas <- dt$cas; ctl <- dt$ctl
  tstar_i <- dt$tstar
  Sy <- R[, "y"]; Sty <- R[, "ty"]; St <- R[, "t"]; Stt <- R[, "tt"]
  Syy <- R[, "yy"]; n_i <- R[, "one"]
  partner <- integer(n); partner[cas] <- ctl; partner[ctl] <- cas
  sgn <- ifelse(dt$case == 1, 1, -1)
  tstar_c <- tstar_i[cas]

  b0 <- rep(0, n); b1 <- rep(0, n)
  beta0 <- mean(Sy / pmax(n_i, 1)); beta1 <- -0.08
  s2 <- 0.02
  D  <- matrix(c(1.2, 0.02, 0.02, 0.01), 2, 2)
  D0 <- matrix(c(1.2, 0, 0, 0.01), 2, 2); df0 <- 4
  aV <- 0; aS <- if (is.null(aS_fixed)) 0 else aS_fixed
  nG <- if (is.null(Zd)) 0 else ncol(Zd)
  aG <- numeric(nG)
  cstep <- 1.1

  keep <- seq(BURN + 1, NITER, by = THIN)
  nk <- length(keep)
  out <- list(aV = numeric(nk), aS = numeric(nk), beta0 = numeric(nk),
              beta1 = numeric(nk), s2 = numeric(nk), D = matrix(0, nk, 3),
              SDval = numeric(nk), SDslo = numeric(nk),
              gamma = if (nG > 0) matrix(0, nk, nG) else NULL)
  eta_store <- matrix(0, nk, length(cas))
  acc_b <- 0; acc_b_tot <- 0; acc_a <- 0; k <- 0

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

    ## ---- 5. association parameters: Laplace-proposal MH ----
    ##   the posterior of (aV,aS) is strongly correlated (banana shaped); a
    ##   plain random walk mixes badly -> use local mode + curvature, corrected
    ##   with the full Hastings ratio.
    d0 <- b0[cas] - b0[ctl]; d1 <- b1[cas] - b1[ctl]
    u <- d0 + d1 * tstar_c; v <- d1
    Xa <- if (is.null(aS_fixed)) cbind(u, v) else cbind(u)
    if (nG > 0) Xa <- cbind(Xa, Zd)
    a_cur <- if (is.null(aS_fixed)) c(aV, aS, aG) else c(aV, aG)
    ll_a <- function(a) sum(logplogis(drop(Xa %*% a)))
    ll_pr <- function(a) -0.5 * sum(a^2) / prior_sd^2
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
      if (is.null(aS_fixed)) {
        aS <- a_prop[2]; if (nG > 0) aG <- a_prop[-(1:2)]
      } else {
        if (nG > 0) aG <- a_prop[-1]
      }
      acc_a <- acc_a + 1
    }

    if (k < nk && it == keep[k + 1]) {
      k <- k + 1
      out$aV[k] <- aV; out$aS[k] <- aS
      out$beta0[k] <- beta0; out$beta1[k] <- beta1; out$s2[k] <- s2
      out$D[k, ] <- c(D[1, 1], D[1, 2], D[2, 2])
      out$SDval[k] <- sd(b0 + b1 * tstar_i)
      out$SDslo[k] <- sd(b1)
      if (nG > 0) out$gamma[k, ] <- aG
      eta_store[k, ] <- drop(Xa %*% (if (is.null(aS_fixed)) c(aV, aS, aG) else c(aV, aG)))
    }
    if (verbose && it %% 10000 == 0)
      cat(sprintf("  iter %6d  aV=%+.3f aS=%+.3f  beta1=%.4f  s2=%.4f\n",
                  it, aV, aS, beta1, s2))
  }
  list(draws = out, eta = eta_store, acc_b = acc_b / acc_b_tot,
       n = n, npair = length(cas), tstar = tstar_i, keep = keep)
}

## ------------------------------------------------------------- summaries --
summ <- function(x) {
  q <- quantile(x, c(0.025, 0.5, 0.975))
  sprintf("%.3f (%.3f, %.3f)", q[2], q[1], q[3])
}

or_tab <- function(fit, label = "") {
  d <- fit$draws
  v <- function(x) list(med = median(x), lo = quantile(x, .025), hi = quantile(x, .975))
  a <- v(exp(-0.1 * d$aS)); b <- v(exp(-d$aS * d$SDslo))
  c1 <- v(exp(-d$aV));      d1 <- v(exp(-d$aV * d$SDval))
  data.frame(label = label,
             pairs = fit$npair, knees = fit$n,
             beta1 = median(d$beta1),
             resid_sd = median(sqrt(d$s2)),
             sd_slope = median(d$SDslo),
             aS = median(d$aS), aS_lo = quantile(d$aS, .025), aS_hi = quantile(d$aS, .975),
             p_aS_neg = mean(d$aS < 0),
             OR_01 = a$med, OR_01_lo = a$lo, OR_01_hi = a$hi,
             OR_SD = b$med, OR_SD_lo = b$lo, OR_SD_hi = b$hi,
             OR_1mm = c1$med, OR_1mm_lo = c1$lo, OR_1mm_hi = c1$hi,
             OR_SDlev = d1$med, OR_SDlev_lo = d1$lo, OR_SDlev_hi = d1$hi)
}

sd_bf01 <- function(x, prior_sd = 20) {
  d <- density(x, n = 4096)
  approx(d$x, d$y, xout = 0, rule = 2)$y / dnorm(0, 0, prior_sd)
}

report <- function(fit, label) {
  d <- fit$draws
  cat(sprintf("\n================ %s ================\n", label))
  cat(sprintf("pairs %d | knees %d | draws %d | MH acc (random effects) %.2f\n",
              fit$npair, fit$n, length(d$aV), fit$acc_b))
  cat("  beta1 (mm/yr)            :", summ(d$beta1), "\n")
  cat("  residual SD              :", summ(sqrt(d$s2)), "\n")
  cat("  SD(u1) (mm/yr)           :", summ(sqrt(d$D[, 3])), "\n")
  cat("  aV (per mm latent level) :", summ(d$aV), "\n")
  cat("  aS (per mm/yr latent)    :", summ(d$aS), "\n")
  cat(sprintf("  P(aV<0)=%.3f   P(aS<0)=%.3f\n", mean(d$aV < 0), mean(d$aS < 0)))
  cat("  OR / 0.1 mm/yr faster    :", summ(exp(-0.1 * d$aS)), "\n")
  cat("  OR / 1 SD faster         :", summ(exp(-d$aS * d$SDslo)), "\n")
  cat("  OR / 1 mm thinner (level):", summ(exp(-d$aV)), "\n")
  if (!is.null(d$gamma)) {
    cat("  within-pair covariate adjustments (OR per unit difference):\n")
    for (j in seq_len(ncol(d$gamma)))
      cat(sprintf("    %-14s %s\n", colnames(d$gamma)[j] %||% paste0("Z", j),
                  summ(exp(d$gamma[, j]))))
  }
  invisible(or_tab(fit, label))
}
`%||%` <- function(a, b) if (is.null(a)) b else a
