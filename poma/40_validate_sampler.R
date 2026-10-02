# ============================================================================
# 40_validate_sampler.R : external validation of the POMA joint-model sampler
#
#   WHY THIS FILE EXISTS
#   Up to manuscript v4 the sampler in jm_core.R contained three defects
#   (see the header of jm_core.R).  A genealogy change on its own is not
#   evidence that the corrected kernel targets the right posterior, so here the
#   corrected sampler, the legacy sampler and a completely INDEPENDENT
#   brute-force sampler are run on the same data and compared.
#
#   THREE ROUTES TO THE SAME POSTERIOR
#     ROUTE 1  gold / importance sampling.  With every parameter except the
#              random effects frozen, independence importance sampling with
#              proposal product(g_i) and weights product(plogis(eta_s)) gives an
#              EXACT (self-normalised) estimate of p(b | rest) plus an honest
#              Monte-Carlo standard error.  This isolates the random-effect step,
#              where defect (1) lived.
#     ROUTE 2  brute force componentwise random-walk MH written from scratch:
#              no sufficient statistics, no conjugacy, the joint log posterior
#              is re-evaluated from the RAW ROWS at every single proposal.
#              Computationally stupid, statistically unimpeachable.
#     ROUTE 3  the two versions of the production sampler (legacy / corrected).
#
#   Anything that survives all three routes is very unlikely to be an artefact
#   of one piece of algebra.
#
#   Usage:  Rscript 40_validate_sampler.R [stage]     stage = A | B | C | all
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
SDIR <- "D:/BaiduSyncdisk/OAI/Scripts/poma"

stage <- if (length(commandArgs(trailingOnly = TRUE)) > 0)
  commandArgs(trailingOnly = TRUE)[1] else "all"

source(file.path(SDIR, "jm_core.R"))                       # ROUTE 3 (corrected)
LEG  <- new.env()
sys.source(file.path(SDIR, "jm_core_legacy_sampler.R"), envir = LEG)   # legacy

logplogis <- function(x) {
  ifelse(x > 30, -log1p(exp(-x)), -log1p(exp(pmin(-x, 30))))
}

## truth used to simulate.  Deliberately in the same regime as the real POMA
## estimates (OR per SD of the latent slope ~ 2-3) so that any bias the kernels
## introduce is visible at the scale at which the results are reported.
##   aS = -5 with sd(u1) = 0.173 gives OR/SD = exp(0.87) = 2.4.
TRUE_PAR <- list(beta0 = 2.9, beta1 = -0.08, s2 = 0.02,
                 D = matrix(c(0.60, 0.02, 0.02, 0.030), 2, 2),
                 aV = -0.50, aS = -5.0, aG = 0.50)
IW_PAR   <- list(D0 = matrix(c(1.2, 0, 0, 0.01), 2, 2), df0 = 4)

## ------------------------------------------------------------- simulator --
## every row of the long table is drawn from exactly the model the sampler
## assumes, and the case label inside each pair is drawn from the event model
## itself, so the simulated likelihood IS the assumed likelihood.
simulate_pairs <- function(npair = 20, visits = 0:3, withZ = FALSE, seed = 1,
                           tp = TRUE_PAR) {
  set.seed(seed)
  n  <- 2 * npair
  Lc <- chol(tp$D)
  U  <- t(replicate(n, drop(t(Lc) %*% rnorm(2))))          # unit random effects
  nv <- length(visits)
  tstar <- runif(npair, 2.5, 5)

  ## within-pair covariate (e.g. KL grade): one value per unit -> difference
  xu <- rnorm(n, 0, 1)
  p1 <- seq(1, n, by = 2); p2 <- p1 + 1
  Zd <- if (withZ) cbind(kl_diff = xu[p1] - xu[p2]) else NULL

  ## SIGN CONVENTION (fixed 2026-10-02).
  ## The model's conditional likelihood for a pair is  P(case = u) =
  ## plogis( aV*[(b0_u - b0_other) + (b1_u - b1_other)*t*] + aS*(b1_u - b1_other) ),
  ## i.e. it is written in the CASE-MINUS-CONTROL direction.  An earlier version
  ## of this simulator drew the label from plogis(eta) with eta built in the
  ## p1-minus-p2 direction, which is exactly the WRONG orientation: the simulated
  ## labels then behaved as if aV and aS had the OPPOSITE sign, so the posterior
  ## came back at (+0.4, +8) when the truth was (-0.4, -8).  Everything below is
  ## therefore built in the p2-minus-p1 direction, matching the model.
  eta_p2 <- tp$aV * ((U[p2, 1] - U[p1, 1]) + (U[p2, 2] - U[p1, 2]) * tstar) +
            tp$aS * (U[p2, 2] - U[p1, 2])
  if (withZ) eta_p2 <- eta_p2 + tp$aG * (xu[p2] - xu[p1])

  ## which member becomes the case is decided BY THE EVENT MODEL
  flip <- rbinom(npair, 1, plogis(eta_p2))
  case_unit  <- rep(0L, n); ctl_unit <- rep(0L, n)
  case_unit[ifelse(flip == 1, p2, p1)] <- 1L
  ctl_unit[ifelse(flip == 0, p2, p1)]  <- 1L

  ord      <- as.vector(rbind(which(case_unit == 1L), which(ctl_unit == 1L)))
  unit_id  <- ord                                     # new unit id, cases first
  ## rows of sv/lg are ordered case(pair 1..S) then control(pair 1..S), so the
  ## per-pair index time repeats WITH THE SAME ORDER TWICE, not interleaved
  tstar_u  <- rep(tstar, times = 2)

  lg <- do.call(rbind, lapply(seq_len(n), function(v) {
    b <- U[ord[v], ]
    mu <- tp$beta0 + tp$beta1 * visits + b[1] + b[2] * visits
    data.frame(unit = unit_id[v], year = visits,
               y = mu + rnorm(length(visits), 0, sqrt(tp$s2)))
  }))
  sv <- data.frame(unit = unit_id,
                   newstrata = rep(seq_len(npair), each = 2),
                   pair      = rep(seq_len(npair), each = 2),
                   case      = rep(c(1, 0), npair),
                   t_idx     = tstar_u / 12,
                   stringsAsFactors = FALSE)
  Zd_out <- if (withZ) cbind(kl_diff = Zd[, 1] * ifelse(flip == 1, -1, 1)) else NULL
  list(lg = lg, sv = sv, dt = build(lg, sv), Zd = Zd_out, U = U[ord, ], tstar = tstar_u)
}

## ------------------------------------------ ROUTE 2 : brute force sampler --
## joint log posterior written directly from the model statement
make_ll <- function(dt, lg, Zd, prior_sd = 20, P0 = IW_PAR) {
  y <- lg$y; tv <- lg$year; ui <- match(lg$unit, dt$unit)
  N <- length(y); n <- dt$n
  cas <- dt$cas; ctl <- dt$ctl; tstar_c <- dt$tstar[cas]
  D0 <- P0$D0; df0 <- P0$df0
  function(b0, b1, beta, s2, D, a) {
    mu  <- beta[1] + beta[2] * tv + b0[ui] + b1[ui] * tv
    res <- y - mu
    ldat <- -0.5 * N * log(2 * pi * s2) - 0.5 * sum(res^2) / s2
    Di   <- solve(D)
    ldD  <- as.numeric(determinant(D, logarithm = TRUE)$modulus)
    quad <- sum(Di[1, 1] * b0^2 + 2 * Di[1, 2] * b0 * b1 + Di[2, 2] * b1^2)
    lb   <- -n * log(2 * pi) - 0.5 * n * ldD - 0.5 * quad
    d0 <- b0[cas] - b0[ctl]; d1 <- b1[cas] - b1[ctl]
    eta <- a[1] * (d0 + d1 * tstar_c) + a[2] * d1
    if (!is.null(Zd)) eta <- eta + drop(Zd %*% a[-(1:2)])
    lev  <- sum(logplogis(eta))
    ls2  <- -1.01 * log(s2) - 0.01 / s2 + log(s2)      # IG(.01,.01) + log-Jacobian
    lDp  <- -(df0 + 2 + 1) / 2 * ldD - 0.5 * sum(diag(D0 %*% Di))
    la   <- -0.5 * sum(a^2) / prior_sd^2
    if (!is.finite(ldat + lb + lev + ls2 + lDp + la)) -1e12 else
      ldat + lb + lev + ls2 + lDp + la
  }
}

run_bruteforce <- function(dt, lg, Zd = NULL, NITER = 100000, BURN = 25000,
                           THIN = 10, seed = 101, prior_sd = 20, verbose = TRUE) {
  set.seed(seed)
  ll <- make_ll(dt, lg, Zd, prior_sd)
  n <- dt$n; cas <- dt$cas; ctl <- dt$ctl
  b0 <- rep(0, n); b1 <- rep(0, n)
  beta <- c(mean(lg$y), 0); s2 <- 0.05
  D <- matrix(c(0.6, 0.02, 0.02, 0.02), 2, 2)
  na <- if (is.null(Zd)) 2 else 2 + ncol(Zd)
  a <- rep(0, na)
  sc <- list(b = c(0.20, 0.12), beta = c(0.06, 0.06), ls2 = 0.20,
             D = c(0.25, 0.06, 0.02), a = rep(0.45, na))
  cnt <- list(b = 0L, beta = 0L, ls2 = 0L, D = 0L, a = 0L)
  tot <- list(b = 0L, beta = 0L, ls2 = 0L, D = 0L, a = 0L)
  keep <- seq(BURN + 1, NITER, by = THIN); nk <- length(keep); k <- 0
  out <- list(aV = numeric(nk), aS = numeric(nk), beta0 = numeric(nk),
              beta1 = numeric(nk), s2 = numeric(nk), D = matrix(0, nk, 3),
              SDval = numeric(nk), SDslo = numeric(nk),
              gamma = if (!is.null(Zd)) matrix(0, nk, ncol(Zd)) else NULL)
  cur <- ll(b0, b1, beta, s2, D, a)
  stopifnot(is.finite(cur))                 # nothing is allowed to start at -Inf

  for (it in 1:NITER) {
    ## --- each knee separately, random walk, acceptance from the RAW density
    for (i in seq_len(n)) {
      bb0 <- b0[i] + rnorm(1, 0, sc$b[1]); bb1 <- b1[i] + rnorm(1, 0, sc$b[2])
      if (is.infinite(cur)) {
        cand <- ll(replace(b0, i, bb0), replace(b1, i, bb1), beta, s2, D, a)
        if (cand > -1e11) { b0[i] <- bb0; b1[i] <- bb1; cur <- cand }
        next
      }
      cand <- ll(replace(b0, i, bb0), replace(b1, i, bb1), beta, s2, D, a)
      tot$b <- tot$b + 1L
      if (log(runif(1)) < cand - cur) {
        b0[i] <- bb0; b1[i] <- bb1; cur <- cand; cnt$b <- cnt$b + 1L
      }
    }
    ## --- fixed effects
    for (j in 1:2) {
      bp <- beta; bp[j] <- bp[j] + rnorm(1, 0, sc$beta[j]); tot$beta <- tot$beta + 1L
      cand <- ll(b0, b1, bp, s2, D, a)
      if (log(runif(1)) < cand - cur) { beta <- bp; cur <- cand; cnt$beta <- cnt$beta + 1L }
    }
    ## --- residual variance (log scale)
    sp <- s2 * exp(rnorm(1, 0, sc$ls2)); tot$ls2 <- tot$ls2 + 1L
    cand <- ll(b0, b1, beta, sp, D, a)
    if (log(runif(1)) < cand - cur) { s2 <- sp; cur <- cand; cnt$ls2 <- cnt$ls2 + 1L }
    ## --- D : additive symmetric proposal on the three free entries
    ep <- rnorm(3, 0, sc$D)
    Dp <- D; Dp[1, 1] <- Dp[1, 1] + ep[1]
    Dp[1, 2] <- Dp[2, 1] <- Dp[1, 2] + ep[2]; Dp[2, 2] <- Dp[2, 2] + ep[3]
    pd <- try(chol(Dp), silent = TRUE)
    if (!inherits(pd, "try-error")) {
      tot$D <- tot$D + 1L
      cand <- ll(b0, b1, beta, s2, Dp, a)
      if (log(runif(1)) < cand - cur) { D <- Dp; cur <- cand; cnt$D <- cnt$D + 1L }
    }
    ## --- association parameters
    for (j in seq_len(na)) {
      ap <- a; ap[j] <- ap[j] + rnorm(1, 0, sc$a[j]); tot$a <- tot$a + 1L
      cand <- ll(b0, b1, beta, s2, D, ap)
      if (log(runif(1)) < cand - cur) { a <- ap; cur <- cand; cnt$a <- cnt$a + 1L }
    }
    ## --- Robbins-Monro step-size tuning during the first half of the burn-in
    if (it < BURN / 2 && it %% 250 == 0) {
      adj <- function(key, target = 0.30) {
        if (tot[[key]] > 0) {
          r <- cnt[[key]] / tot[[key]]
          sc[[key]] <<- pmin(pmax(sc[[key]] * exp(1.2 * (r - target)), 1e-4), 5)
          cnt[[key]] <<- 0L; tot[[key]] <<- 0L
        }
      }
      adj("b", 0.35); adj("beta", 0.35); adj("ls2", 0.35); adj("D", 0.35); adj("a", 0.30)
    }
    if (k < nk && it == keep[k + 1]) {
      k <- k + 1
      out$aV[k] <- a[1]; out$aS[k] <- a[2]
      out$beta0[k] <- beta[1]; out$beta1[k] <- beta[2]; out$s2[k] <- s2
      out$D[k, ]  <- c(D[1, 1], D[1, 2], D[2, 2])
      out$SDval[k] <- sd(b0 + b1 * dt$tstar)
      out$SDslo[k] <- sd(b1)
      if (!is.null(Zd)) out$gamma[k, ] <- a[-(1:2)]
    }
    if (verbose && it %% 20000 == 0)
      cat(sprintf("    brute force iter %6d  aV=%+.3f aS=%+.3f  s2=%.4f  beta1=%.4f\n",
                  it, a[1], a[2], s2, beta[2]))
  }
  list(draws = out, acc_b = NA, n = n, npair = length(cas))
}

## ---------------------------------------------------------------- helpers --
qs <- function(x) sprintf("%+.3f (%.3f, %.3f)",
                          quantile(x, .5), quantile(x, .025), quantile(x, .975))

summ_draws <- function(d) {
  SDslo <- d$SDslo
  c(beta1 = sprintf("%.4f", median(d$beta1)),
    resid_sd = sprintf("%.4f", median(sqrt(d$s2))),
    sd_slope = sprintf("%.4f", median(SDslo)),
    aV = qs(d$aV), aS = qs(d$aS),
    OR_01 = qs(exp(-0.1 * d$aS)),
    OR_SD = qs(exp(-d$aS * SDslo)),
    D11 = sprintf("%.3f", median(d$D[, 1])),
    D22 = sprintf("%.4f", median(d$D[, 3])))
}

out <- character(0)
say <- function(...) {
  s <- paste0(...)
  cat(s, "\n")
  out <<- c(out, s)
  writeLines(out, file.path(BASE, "sampler_validation.txt"))
}

## ==========================================================================
## STAGE A : the two candidate kernels for the random-effect step against
##   (i)  the EXACT 1-D quadrature reference below, and
##   (ii) self-normalised importance sampling,
##   with every parameter except the random effects frozen at the truth.
## ==========================================================================
## ---------------------------------------------- exact 1-D quadrature gold --
##   p(b | rest) factorises over pairs, and inside a pair it is
##        N(d ; mu_cas - mu_ctl , V_cas + V_ctl) x plogis(eta(d)) ,  d = (d0,d1)
##   plogis depends on d only through eta = c'd, so integrating out the
##   orthogonal direction leaves the ONE-DIMENSIONAL density
##        N(eta ; c'mu , c'(V_cas+V_ctl)c) x plogis(eta),
##   which a fine grid integrates to machine precision.  No Monte Carlo, no
##   tail problems: this is the reference the two kernels have to match.
stage_a3 <- function() {
  S <- simulate_pairs(npair = 20, seed = 7)
  dt <- S$dt; tp <- TRUE_PAR
  R <- dt$R; n_i <- R[, "one"]; Sy <- R[, "y"]; St <- R[, "t"]
  Sty <- R[, "ty"]; Stt <- R[, "tt"]
  Dinv <- solve(tp$D)
  a11 <- n_i / tp$s2 + Dinv[1, 1]; a12 <- St / tp$s2 + Dinv[1, 2]
  a22 <- Stt / tp$s2 + Dinv[2, 2]; det <- a11 * a22 - a12^2
  V11 <- a22 / det; V22 <- a11 / det; V12 <- -a12 / det
  r0 <- (Sy - (n_i * tp$beta0 + St * tp$beta1)) / tp$s2
  r1 <- (Sty - (St * tp$beta0 + Stt * tp$beta1)) / tp$s2
  m0 <- V11 * r0 + V12 * r1; m1 <- V12 * r0 + V22 * r1
  aV <- tp$aV; aS <- tp$aS
  gr <- seq(-14, 14, length.out = 40001)
  eta_q <- numeric(length(dt$cas)); pq <- eta_q
  for (s in seq_along(dt$cas)) {
    i <- dt$cas[s]; j <- dt$ctl[s]
    Sig <- matrix(c(V11[i] + V11[j], V12[i] + V12[j],
                    V12[i] + V12[j], V22[i] + V22[j]), 2, 2)
    mud <- c(m0[i] - m0[j], m1[i] - m1[j])
    cv  <- c(aV, aV * dt$tstar[i] + aS)
    me  <- sum(cv * mud); se <- sqrt(drop(t(cv) %*% Sig %*% cv))
    g   <- me + gr * se
    w   <- dnorm(g, me, se) * plogis(g)
    eta_q[s] <- sum(g * w) / sum(w)
    pq[s]    <- sum(plogis(g) * w) / sum(w)
  }
  list(ETAmean = mean(eta_q), Pmean = mean(pq))
}

if (stage %in% c("A", "all")) {
  say("\n################ STAGE A : random-effect step vs exact gold ################")
  S <- simulate_pairs(npair = 20, seed = 7)
  dt <- S$dt; tp <- TRUE_PAR
  b0 <- rep(0, dt$n); b1 <- rep(0, dt$n)
  beta <- c(tp$beta0, tp$beta1); s2 <- tp$s2; D <- tp$D; a <- c(tp$aV, tp$aS)
  ## everything except the random effects is FROZEN at the truth here

  ## ---- proposal g_i built per unit straight from the raw rows ------------
  proposal_unit <- function(i, M) {
    rows <- which(S$lg$unit == dt$unit[i])
    tv <- S$lg$year[rows]; yv <- S$lg$y[rows]
    X <- cbind(1, tv); XtX <- crossprod(X); Xtr <- c(sum(yv - beta[1] - beta[2] * tv),
                                                     sum(tv * (yv - beta[1] - beta[2] * tv)))
    pr <- solve(D) + XtX / s2; mu <- solve(pr, Xtr / s2); Si <- solve(pr)
    Ch <- chol(Si)
    z  <- matrix(rnorm(2 * M), ncol = 2)
    res <- t(t(z %*% Ch) + mu); colnames(res) <- c("g0", "g1")
    res
  }
  set.seed(4242); M <- 60000
  P <- lapply(seq_len(dt$n), proposal_unit, M = M)

  d0 <- P[[1]]  # placeholder
  E <- function(f, lw) {                    # self-normalised importance estimate
    ## the weights MUST be passed in.  An earlier version read `lw` lexically
    ## from the global frame (all zeros) and silently returned UNWEIGHTED
    ## proposal expectations - a textbook scoping trap.
    w <- exp(lw - max(lw)); w <- w / sum(w)
    list(est = sum(w * f), se = sqrt(sum((w * (f - sum(w * f)))^2)))
  }
  is_est <- function() {
    Q <- list()
    ## weights: product over pairs of plogis(eta)
    lw <- numeric(M)
    Slog <- matrix(0, M, length(dt$cas))
    for (s in seq_along(dt$cas)) {
      d0 <- P[[dt$cas[s]]][, 1] - P[[dt$ctl[s]]][, 1]
      d1 <- P[[dt$cas[s]]][, 2] - P[[dt$ctl[s]]][, 2]
      eta <- a[1] * (d0 + d1 * dt$tstar[dt$cas[s]]) + a[2] * d1
      lw <- lw + logplogis(eta); Slog[, s] <- eta
    }
    Q$eta_mean <- E(rowMeans(Slog), lw)
    Q$p_mean   <- E(rowMeans(1 / (1 + exp(-Slog))), lw)
    Q$sd_b1    <- E(apply(sapply(P, function(x) x[, 2]), 1, sd), lw)
    ww <- exp(lw - max(lw)); Q$ess <- sum(ww)^2 / sum(ww^2)
    Q
  }
  G <- is_est()
  say(sprintf("gold (importance sampling, M=%d, ESS=%.0f):", M, G$ess))
  nmv <- setdiff(names(G), "ess")
  for (nm in nmv)
    say(sprintf("   %-9s = %+.4f  (MC se %.4f)", nm, G[[nm]]$est, G[[nm]]$se))

  ## ---- the two competing kernels, everything else frozen ------------------
  kern <- function(kind, NITER = 1000000, BURN = 50000, seed = 9) {
    set.seed(seed)
    b0 <- rep(0, dt$n); b1 <- rep(0, dt$n)
    n <- dt$n; cas <- dt$cas; ctl <- dt$ctl; npair <- length(cas)
    tstar_c <- dt$tstar[cas]
    keep <- seq(BURN + 1, NITER, by = 5); nk <- length(keep); k <- 0
    A <- numeric(nk); Bv <- numeric(nk); Cs <- numeric(nk)
    tstar_i <- dt$tstar
    partner <- integer(n); partner[cas] <- ctl; partner[ctl] <- cas
    sgn <- ifelse(dt$case == 1, 1, -1)
    Dinv <- solve(D)
    ## sufficient statistics identical to production code
    R <- dt$R
    Sy <- R[, "y"]; Sty <- R[, "ty"]; St <- R[, "t"]; Stt <- R[, "tt"]; n_i <- R[, "one"]
    for (it in 1:NITER) {
      a11 <- n_i / s2 + Dinv[1, 1]; a12 <- St / s2 + Dinv[1, 2]
      a22 <- Stt / s2 + Dinv[2, 2]; det <- a11 * a22 - a12^2
      V11 <- a22 / det; V22 <- a11 / det; V12 <- -a12 / det
      r0 <- (Sy - (n_i * beta[1] + St * beta[2])) / s2
      r1 <- (Sty - (St * beta[1] + Stt * beta[2])) / s2
      m0 <- V11 * r0 + V12 * r1; m1 <- V12 * r0 + V22 * r1
      L11 <- sqrt(pmax(V11, 1e-12)); L21 <- V12 / L11
      L22 <- sqrt(pmax(V22 - L21^2, 1e-12))
      ## NOTE the SAME z0 must feed both components: that is what makes the
      ## proposal covariance equal V instead of its diagonal.  Writing three
      ## separate rnorm(n) calls here (an earlier version of this harness did)
      ## silently sets cov(b0,b1) = 0 and moves the chain off target.
      z0 <- rnorm(n); z1 <- rnorm(n)
      b0p <- m0 + L11 * z0
      b1p <- m1 + L21 * z0 + L22 * z1
      if (kind == "jacobi") {
        d0n <- sgn * (b0p - b0[partner]); d1n <- sgn * (b1p - b1[partner])
        d0o <- sgn * (b0 - b0[partner]);  d1o <- sgn * (b1 - b1[partner])
        e_n <- a[1] * (d0n + d1n * tstar_i) + a[2] * d1n
        e_o <- a[1] * (d0o + d1o * tstar_i) + a[2] * d1o
        acc <- runif(n) < exp(logplogis(e_n) - logplogis(e_o))
        b0[acc] <- b0p[acc]; b1[acc] <- b1p[acc]
      } else {
        d0n <- b0p[cas] - b0p[ctl]; d1n <- b1p[cas] - b1p[ctl]
        d0o <- b0[cas]  - b0[ctl];  d1o <- b1[cas]  - b1[ctl]
        e_n <- a[1] * (d0n + d1n * tstar_c) + a[2] * d1n
        e_o <- a[1] * (d0o + d1o * tstar_c) + a[2] * d1o
        acc <- runif(npair) < exp(logplogis(e_n) - logplogis(e_o))
        if (any(acc)) {
          i1 <- cas[acc]; i2 <- ctl[acc]
          b0[i1] <- b0p[i1]; b1[i1] <- b1p[i1]; b0[i2] <- b0p[i2]; b1[i2] <- b1p[i2]
        }
      }
      if (k < nk && it == keep[k + 1]) {
        k <- k + 1
        d0 <- b0[cas] - b0[ctl]; d1 <- b1[cas] - b1[ctl]
        ee <- a[1] * (d0 + d1 * tstar_c) + a[2] * d1
        A[k] <- mean(ee); Bv[k] <- mean(plogis(ee)); Cs[k] <- sd(b1)
      }
    }
    list(eta_mean = mean(A), p_mean = mean(Bv), sd_b1 = mean(Cs),
         series = list(eta = A, p = Bv, sdb1 = Cs))
  }
  kj <- kern("jacobi");  kp <- kern("pairblock")
  GD <- stage_a3()                       # exact, no Monte Carlo at all
  mcse <- function(x) {                  # batch-means Monte Carlo standard error
    nb <- 25; sz <- floor(length(x) / nb)
    m <- colMeans(matrix(x[1:(nb * sz)], nrow = sz))
    sqrt(var(m) / nb)
  }
  say("")
  say(sprintf("deterministic gold (1-D quadrature, EXACT) : eta_mean %+.6f | p_mean %+.6f",
              GD$ETAmean, GD$Pmean))
  say(sprintf("importance sampling (M=%d, ESS=%.0f)       : eta_mean %+.4f (se %.4f) | p_mean %+.4f (se %.4f) | sd_b1 %+.4f (se %.4f)",
              M, G$ess, G$eta_mean$est, G$eta_mean$se,
              G$p_mean$est, G$p_mean$se, G$sd_b1$est, G$sd_b1$se))
  gv <- c(eta_mean = GD$ETAmean, p_mean = GD$Pmean)
  for (nm in c("eta_mean", "p_mean")) {
    se_j <- mcse(kj$series[[if (nm == "eta_mean") "eta" else "p"]])
    se_p <- mcse(kp$series[[if (nm == "eta_mean") "eta" else "p"]])
    say(sprintf("  %-9s exact %+.6f | jacobi %+.4f (mcse %.4f, %s) | pair-block %+.4f (mcse %.4f, %s)",
                nm, gv[[nm]], kj[[nm]], se_j,
                ifelse(abs(kj[[nm]] - gv[[nm]]) < 4 * max(se_j, 1e-4), "PASS", "FAIL"),
                kp[[nm]], se_p,
                ifelse(abs(kp[[nm]] - gv[[nm]]) < 4 * max(se_p, 1e-4), "PASS", "FAIL")))
  }
  say(sprintf("  sd_b1     (no exact value) | jacobi %+.4f | pair-block %+.4f | I.S. %+.4f",
              kj$sd_b1, kp$sd_b1, G$sd_b1$est))
}

## ==========================================================================
## STAGE A2 : ARBITER.  Stage A produced a contradiction (the two production
## kernels agreed with each other but not with the importance-sample gold), so
## the b-step target is recomputed here by a chain that never uses the g_i
## algebra at all: every proposal is a plain random walk and every acceptance
## uses the joint log density re-evaluated from the raw rows.
##   Superseded by STAGE A3 (the exact deterministic gold): kept for provenance
##   only, and deliberately NOT part of "all" because it costs ~3 minutes.
## ==========================================================================
if (stage %in% c("A2")) {
  say("\n################ STAGE A2 : random-walk arbiter for the b step ################")
  S <- simulate_pairs(npair = 20, seed = 7)
  dt <- S$dt; tp <- TRUE_PAR
  n <- dt$n; cas <- dt$cas; ctl <- dt$ctl
  ll <- make_ll(dt, S$lg, NULL, prior_sd = 20)
  beta <- c(tp$beta0, tp$beta1); s2 <- tp$s2; D <- tp$D; a <- c(tp$aV, tp$aS)
  set.seed(999)
  b0 <- rep(0, n); b1 <- rep(0, n)
  cur <- ll(b0, b1, beta, s2, D, a)
  stopifnot(is.finite(cur))
  sc <- c(0.10, 0.06); cnt <- 0L; tot <- 0L
  NITER <- 120000; BURN <- 20000; keep <- seq(BURN + 1, NITER, by = 10)
  nk <- length(keep); k <- 0
  A <- numeric(nk); Bv <- numeric(nk); Cs <- numeric(nk)
  tstar_c <- dt$tstar[cas]
  for (it in 1:NITER) {
    for (i in seq_len(n)) {
      nb0 <- b0[i] + rnorm(1, 0, sc[1]); nb1 <- b1[i] + rnorm(1, 0, sc[2])
      cand <- ll(replace(b0, i, nb0), replace(b1, i, nb1), beta, s2, D, a)
      tot <- tot + 1L
      if (log(runif(1)) < cand - cur) {
        b0[i] <- nb0; b1[i] <- nb1; cur <- cand; cnt <- cnt + 1L
      }
    }
    if (it < BURN / 2 && it %% 500 == 0 && tot > 0) {
      sc <- pmin(pmax(sc * exp(1.5 * (cnt / tot - 0.35)), 1e-4), 2); cnt <- 0L; tot <- 0L
    }
    if (k < nk && it == keep[k + 1]) {
      k <- k + 1
      ee <- a[1] * ((b0[cas] - b0[ctl]) + (b1[cas] - b1[ctl]) * tstar_c) +
            a[2] * (b1[cas] - b1[ctl])
      A[k] <- mean(ee); Bv[k] <- mean(plogis(ee)); Cs[k] <- sd(b1)
    }
    if (it %% 40000 == 0) cat(sprintf("    arbiter iter %6d  eta_mean=%+.4f\n", it, mean(A[1:max(k,1)])))
  }
  say(sprintf("  arbiter (random walk) : eta_mean %+.4f | p_mean %+.4f | sd_b1 %+.4f",
              mean(A), mean(Bv), mean(Cs)))
  say("  (compare with the deterministic gold and the two production kernels")
  say("   printed by STAGE A; the arbiter shares no algebra with either.)")
}

## ==========================================================================
## STAGE A3 : DETERMINISTIC GOLD for the b step (replaces the importance sample)
##   p(b | rest) factorises over pairs, and inside a pair it is
##        N(d ; mu_cas - mu_ctl , V_cas + V_ctl) x plogis(eta(d)) ,  d = (d0,d1)
##   Because plogis depends on d only through eta = c'd, integrating out the
##   orthogonal direction leaves a ONE-DIMENSIONAL density
##        N(eta ; c'mu , c'(V_cas+V_ctl)c) x plogis(eta)
##   so every expectation of eta (and of plogis(eta)) is a 1-D integral that a
##   fine grid evaluates to machine precision.  No Monte Carlo, no tail problems
##   -- this is the reference the two kernels have to match.
## ==========================================================================
if (stage %in% c("A3", "all")) {
  say("\n################ STAGE A3 : deterministic (1-D quadrature) gold ################")
  G <- stage_a3()
  say(sprintf("  gold exact : eta_mean %+.6f | p_mean %+.6f", G$ETAmean, G$Pmean))
  say("  (the same figures are used as the reference line inside STAGE A)")
}

## ==========================================================================
## STAGE B : full model, three routes side by side (no covariates)
## ==========================================================================
if (stage %in% c("B", "all")) {
  say("\n################ STAGE B : full model, simulated truth ################")
  S <- simulate_pairs(npair = 30, seed = 11)
  say(sprintf("30 pairs / %d knees / %d rows | truth aV=%.2f aS=%.2f beta1=%.3f s2=%.3f D22=%.3f",
              S$dt$n, nrow(S$lg), TRUE_PAR$aV, TRUE_PAR$aS, TRUE_PAR$beta1,
              TRUE_PAR$s2, TRUE_PAR$D[2, 2]))
  cat("  running brute force reference ...\n")
  t0 <- Sys.time()
  BF <- run_bruteforce(S$dt, S$lg, NITER = 80000, BURN = 20000, THIN = 8, seed = 101)
  say(sprintf("  brute force done in %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  FL <- fit_jm(S$dt, NITER = 60000, BURN = 10000, THIN = 5, seed = 202, verbose = FALSE)
  LX <- LEG$fit_jm(S$dt, NITER = 60000, BURN = 10000, THIN = 5, seed = 202, verbose = FALSE)
  tab <- rbind(brute_force = summ_draws(BF$draws),
               corrected    = summ_draws(FL$draws),
               legacy       = summ_draws(LX$draws),
               truth        = c(beta1 = sprintf("%.4f", TRUE_PAR$beta1),
                                resid_sd = sprintf("%.4f", sqrt(TRUE_PAR$s2)),
                                sd_slope = sprintf("%.4f", sqrt(TRUE_PAR$D[2, 2])),
                                aV = "---", aS = "---", OR_01 = "---", OR_SD = "---",
                                D11 = sprintf("%.3f", TRUE_PAR$D[1, 1]),
                                D22 = sprintf("%.4f", TRUE_PAR$D[2, 2])))
  say(paste(capture.output(print(tab, quote = FALSE)), collapse = "\n"))
}

## ==========================================================================
## STAGE C : same thing WITH the within-pair covariate term (defect 3)
## ==========================================================================
if (stage %in% c("C", "all")) {
  say("\n################ STAGE C : with within-pair covariate Zd ################")
  S <- simulate_pairs(npair = 30, withZ = TRUE, seed = 13)
  say(sprintf("30 pairs | truth aV=%.2f aS=%.2f aG=%.2f",
              TRUE_PAR$aV, TRUE_PAR$aS, TRUE_PAR$aG))
  cat("  running brute force reference ...\n")
  BF <- run_bruteforce(S$dt, S$lg, Zd = S$Zd, NITER = 80000, BURN = 20000,
                       THIN = 8, seed = 303)
  FL <- fit_jm(S$dt, NITER = 60000, BURN = 10000, THIN = 5, seed = 404,
               Zd = S$Zd, verbose = FALSE)
  LX <- LEG$fit_jm(S$dt, NITER = 60000, BURN = 10000, THIN = 5, seed = 404,
                   Zd = S$Zd, verbose = FALSE)
  gl <- function(d, nm) sprintf("%s = %+.3f (%.3f, %.3f)", nm,
                                median(d$gamma[, 1]), quantile(d$gamma[, 1], .025),
                                quantile(d$gamma[, 1], .975))
  say(sprintf("  truth gamma = %+.2f", TRUE_PAR$aG))
  say(paste0("  brute force  ", gl(BF$draws, "gamma")))
  say(paste0("  corrected    ", gl(FL$draws, "gamma")))
  say(paste0("  legacy       ", gl(LX$draws, "gamma")))
  tab <- rbind(brute_force = summ_draws(BF$draws),
               corrected    = summ_draws(FL$draws),
               legacy       = summ_draws(LX$draws))
  say(paste(capture.output(print(tab, quote = FALSE)), collapse = "\n"))
}

say("\n-> sampler_validation.txt")
