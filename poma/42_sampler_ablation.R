# ============================================================================
# 42 : ABLATION.  The corrected sampler differs from the legacy one in exactly
#      three places (verified by diff).  On the real 191 matched pairs the two
#      samplers disagree far beyond Monte-Carlo error, so this script switches
#      each correction on and off separately and measures which one moves the
#      answer:
#
#        fix1 : pair-blocked random-effect update (legacy = Jacobi sweep)
#        fix2 : beta drawn from its full conditional (legacy = conditional mean)
#        fix3 : Zd %*% aG retained in the random-effect acceptance ratio
#
#      Only the sampler body is copied; everything else is jm_core.R.
#      -> sampler_ablation.csv
# ============================================================================

BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
source("D:/BaiduSyncdisk/OAI/Scripts/poma/jm_core.R")

fit_abl <- function(dt, NITER = 30000, BURN = 5000, THIN = 5, seed = 1,
                    prior_sd = 20, Zd = NULL,
                    fix1 = TRUE, fix2 = TRUE, fix3 = TRUE) {
  set.seed(seed)
  R <- dt$R; n <- dt$n; cas <- dt$cas; ctl <- dt$ctl; tstar_i <- dt$tstar
  Sy <- R[, "y"]; Sty <- R[, "ty"]; St <- R[, "t"]; Stt <- R[, "tt"]
  Syy <- R[, "yy"]; n_i <- R[, "one"]
  partner <- integer(n); partner[cas] <- ctl; partner[ctl] <- cas
  sgn <- ifelse(dt$case == 1, 1, -1)
  tstar_c <- tstar_i[cas]; npair <- length(cas)
  b0 <- rep(0, n); b1 <- rep(0, n)
  beta0 <- mean(Sy / pmax(n_i, 1)); beta1 <- -0.08
  s2 <- 0.02
  D <- matrix(c(1.2, 0.02, 0.02, 0.01), 2, 2)
  D0 <- matrix(c(1.2, 0, 0, 0.01), 2, 2); df0 <- 4
  aV <- 0; aS <- 0
  nG <- if (is.null(Zd)) 0 else ncol(Zd); aG <- numeric(nG)
  cstep <- 1.1
  keep <- seq(BURN + 1, NITER, by = THIN); nk <- length(keep)
  out <- list(aV = numeric(nk), aS = numeric(nk), beta0 = numeric(nk),
              beta1 = numeric(nk), s2 = numeric(nk), D = matrix(0, nk, 3),
              SDval = numeric(nk), SDslo = numeric(nk),
              gamma = if (nG > 0) matrix(0, nk, nG) else NULL)
  acc_b <- 0; acc_b_tot <- 0; k <- 0
  for (it in 1:NITER) {
    ## ---- 1. random effects --------------------------------------------------
    Dinv <- solve(D)
    a11 <- n_i / s2 + Dinv[1, 1]; a12 <- St / s2 + Dinv[1, 2]
    a22 <- Stt / s2 + Dinv[2, 2]; det <- a11 * a22 - a12^2
    V11 <- a22 / det; V22 <- a11 / det; V12 <- -a12 / det
    r0 <- (Sy  - (n_i * beta0 + St  * beta1)) / s2
    r1 <- (Sty - (St  * beta0 + Stt * beta1)) / s2
    m0 <- V11 * r0 + V12 * r1; m1 <- V12 * r0 + V22 * r1
    L11 <- sqrt(pmax(V11, 1e-12)); L21 <- V12 / L11
    L22 <- sqrt(pmax(V22 - L21^2, 1e-12))
    z0 <- rnorm(n); z1 <- rnorm(n)
    b0p <- m0 + L11 * z0; b1p <- m1 + L21 * z0 + L22 * z1
    zt <- if (nG > 0) drop(Zd %*% aG) else 0
    if (fix1) {
      d0n <- b0p[cas] - b0p[ctl]; d1n <- b1p[cas] - b1p[ctl]
      d0o <- b0[cas]  - b0[ctl];  d1o <- b1[cas]  - b1[ctl]
      zt_b <- if (fix3) zt else 0
      e_n <- aV * (d0n + d1n * tstar_c) + aS * d1n + zt_b
      e_o <- aV * (d0o + d1o * tstar_c) + aS * d1o + zt_b
      acc_i <- runif(npair) < exp(logplogis(e_n) - logplogis(e_o))
      acc_b <- acc_b + sum(acc_i); acc_b_tot <- acc_b_tot + npair
      if (any(acc_i)) {
        i1 <- cas[acc_i]; i2 <- ctl[acc_i]
        b0[i1] <- b0p[i1]; b1[i1] <- b1p[i1]
        b0[i2] <- b0p[i2]; b1[i2] <- b1p[i2]
      }
    } else {
      d0n <- sgn * (b0p - b0[partner]); d1n <- sgn * (b1p - b1[partner])
      d0o <- sgn * (b0  - b0[partner]); d1o <- sgn * (b1  - b1[partner])
      e_n <- aV * (d0n + d1n * tstar_i) + aS * d1n
      e_o <- aV * (d0o + d1o * tstar_i) + aS * d1o
      acc_i <- runif(n) < exp(logplogis(e_n) - logplogis(e_o))
      acc_b <- acc_b + sum(acc_i); acc_b_tot <- acc_b_tot + n
      b0[acc_i] <- b0p[acc_i]; b1[acc_i] <- b1p[acc_i]
    }
    ## ---- 2. fixed effects ---------------------------------------------------
    Ab0 <- n_i * b0 + St * b1; Ab1 <- St * b0 + Stt * b1
    P  <- matrix(c(sum(n_i), sum(St), sum(St), sum(Stt)), 2, 2)
    q  <- c(sum(Sy) - sum(Ab0), sum(Sty) - sum(Ab1))
    if (fix2) {
      Pinv <- try(solve(P), silent = TRUE)
      if (inherits(Pinv, "try-error")) Pinv <- diag(2) * 1e-8
      Chb <- try(chol(s2 * Pinv), silent = TRUE)
      if (inherits(Chb, "try-error")) Chb <- chol(diag(2) * 1e-8)
      bb <- solve(P, q) + drop(t(Chb) %*% rnorm(2))
    } else bb <- solve(P, q)
    beta0 <- bb[1]; beta1 <- bb[2]
    ## ---- 3. residual variance ----------------------------------------------
    g0 <- beta0 + b0; g1 <- beta1 + b1
    SSE <- sum(Syy - 2 * (g0 * Sy + g1 * Sty) + g0^2 * n_i +
                 2 * g0 * g1 * St + g1^2 * Stt)
    s2 <- 1 / rgamma(1, shape = 0.01 + sum(n_i) / 2, rate = 0.01 + SSE / 2)
    ## ---- 4. D ---------------------------------------------------------------
    B <- cbind(b0, b1)
    D <- solve(rWishart(1, df0 + n, solve(D0 + t(B) %*% B))[, , 1])
    ## ---- 5. association parameters -----------------------------------------
    d0 <- b0[cas] - b0[ctl]; d1 <- b1[cas] - b1[ctl]
    u <- d0 + d1 * tstar_c; v <- d1
    Xa <- cbind(u, v); if (nG > 0) Xa <- cbind(Xa, Zd)
    a_cur <- c(aV, aS, aG)
    ll_a <- function(a) sum(logplogis(drop(Xa %*% a)))
    ll_pr <- function(a) -0.5 * sum(a^2) / prior_sd^2
    lp_prop <- function(a) {
      e <- drop(Xa %*% a); p <- plogis(e); w <- p * (1 - p)
      g <- drop(crossprod(Xa, 1 - p)); H <- -crossprod(Xa * w, Xa)
      Si <- try(solve(-H), silent = TRUE)
      if (inherits(Si, "try-error")) Si <- diag(length(a)) * 0.05
      dt <- determinant(Si, logarithm = TRUE)
      if (!is.finite(dt$modulus) || dt$modulus <= -700) Si <- diag(length(a)) * 0.05
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
      aV <- a_prop[1]; aS <- a_prop[2]; if (nG > 0) aG <- a_prop[-(1:2)]
    }
    if (k < nk && it == keep[k + 1]) {
      k <- k + 1
      out$aV[k] <- aV; out$aS[k] <- aS
      out$beta0[k] <- beta0; out$beta1[k] <- beta1; out$s2[k] <- s2
      out$D[k, ] <- c(D[1, 1], D[1, 2], D[2, 2])
      out$SDval[k] <- sd(b0 + b1 * tstar_i); out$SDslo[k] <- sd(b1)
      if (nG > 0) out$gamma[k, ] <- aG
    }
  }
  list(draws = out, acc_b = acc_b / acc_b_tot, n = n, npair = npair)
}

## --------------------------------------------------------------- real data --
L <- read.csv(file.path(BASE, "poma_long.csv"))
U <- read.csv(file.path(BASE, "poma_analysis_long.csv"))
COVARS <- c("V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain")
metric <- "cMFTC_ThCtAB_aMe"
sub <- L[!is.na(L[[metric]]), ]
lg  <- data.frame(unit = sub$unit, year = sub$years, y = sub[[metric]],
                  newstrata = sub$newstrata, case = sub$case)
sv  <- U[, c("unit", "newstrata", "case", "t_index_months", COVARS, "V00PASE")]
sv  <- sv[sv$unit %in% lg$unit, ]
sv  <- sv[complete.cases(sv[, COVARS]), ]
lg  <- lg[lg$unit %in% sv$unit, ]
sv$t_idx <- sv$t_index_months / 12
A <- prune_pairs(lg, sv, 2); dtA <- build(A$lg, A$sv)
uA <- U[match(A$sv$unit, U$unit), ]
Zd <- cbind(kl_diff    = uA$V00XRKL[dtA$cas]    - uA$V00XRKL[dtA$ctl],
            womac_diff = uA$WOMAC_pain[dtA$cas] - uA$WOMAC_pain[dtA$ctl])
cat(sprintf("Set A : %d pairs / %d knees | visits per knee %.2f | median %.0f\n",
            length(dtA$cas), dtA$n, nrow(A$lg) / dtA$n,
            median(table(A$lg$unit))))

one <- function(lab, Z, f1, f2, f3, seeds = c(2026, 2027, 2028)) {
  r <- lapply(seeds, function(s) {
    f <- fit_abl(dtA, seed = s, Zd = Z, fix1 = f1, fix2 = f2, fix3 = f3)
    d <- f$draws
    data.frame(label = lab, seed = s,
               ORSD = median(exp(-d$aS * d$SDslo)),
               OR01 = median(exp(-0.1 * d$aS)),
               aS = median(d$aS), aV = median(d$aV),
               sd_slope = median(d$SDslo),
               ORSDlev = median(exp(-d$aV * d$SDval)),
               g_kl = if (is.null(d$gamma)) NA_real_ else median(d$gamma[, 1]),
               acc_re = f$acc_b)
  })
  do.call(rbind, r)
}

res <- rbind(
  one("crude | legacy      (f1=F f2=F)",  NULL, FALSE, FALSE, TRUE),
  one("crude | fix1 only   (f1=T f2=F)",  NULL, TRUE,  FALSE, TRUE),
  one("crude | fix2 only   (f1=F f2=T)",  NULL, FALSE, TRUE,  TRUE),
  one("crude | corrected   (f1=T f2=T)",  NULL, TRUE,  TRUE,  TRUE),
  one("adj   | legacy      (f1=F f2=F)",  Zd,   FALSE, FALSE, TRUE),
  one("adj   | corrected   (f1=T f2=T)",  Zd,   TRUE,  TRUE,  TRUE),
  one("adj   | corrected, zt dropped",    Zd,   TRUE,  TRUE,  FALSE)
)
write.csv(res, file.path(BASE, "sampler_ablation.csv"), row.names = FALSE)

cat("\n============ ABLATION on real Set A (mean over 3 seeds) ============\n")
ag <- aggregate(cbind(ORSD, OR01, aS, aV, sd_slope, ORSDlev, g_kl, acc_re) ~ label,
                data = res, FUN = mean)
ag$ORSD_sd <- aggregate(ORSD ~ label, res, sd)$ORSD
print(ag, row.names = FALSE, digits = 4)
cat("\n-> sampler_ablation.csv\n")
