# -*- coding: utf-8 -*-
r"""
03 : Is CTh-Score progression a continuum or a set of discrete classes?  (P2, decisive)

Where this sits
---------------
`01_cth_trajectory_lcmm.R`  fits the authors' specification (lcmm::hlme).  The BIC
                            keeps falling with K, so "the number of classes" is not
                            pinned down by fit alone.
`02_cth_blup_gmm.py`        gives the distribution-light route (two-step BLUP + GMM).
                            Its classes turn out to be ordered almost entirely by
                            BASELINE level, not by slope.

This script answers the underlying question directly:

    Is the per-knee rate of CTh-Score change a continuously distributed quantity,
    or the mixture of a small number of well-separated rate classes?

Three rate/level variables are examined, because each carries a different bias and
the three together bracket the truth:

  slope_blup   BLUP per-knee slope from the mixed model.  Best predictor of the
               individual rate, but shrinkage compresses the tails -> if anything it
               UNDERSTATES heterogeneity.
  slope_ols_ge4
               unshrunk per-knee OLS slope, restricted to knees with >= 4 visits.
               Unbiased per knee, but each knee carries estimation noise -> the
               observed density is the true density convolved with noise, which BLURS
               modes -> if anything it OVERSTATES smoothness.  The restriction matters:
               on all 8,278 knees the OLS rate has excess kurtosis 39 and a maximum of
               +53 pts/yr, produced by two-visit knees with a short interval, and those
               outliers dominate any modality diagnostic.
  baseline     baseline CTh-Score (severity level), carried along as a contrast: the
               level is the part of the data that does look lumpy.

If the rate is unimodal under BOTH the shrunk and the noisy estimator, the conclusion
is not an artefact of the estimator.

Evidence
--------
(a) modality of the density as a function of bandwidth.
(b) Silverman (1981, JRSS-B) critical-bandwidth bootstrap test of H0: unimodal.
(c) one-dimensional Gaussian mixture, K = 1..8, with BIC/AIC plus the geometry of the
    fitted components: pairwise separation  |mu_i-mu_j| / sqrt((s_i^2+s_j^2)/2)  and
    the overlap of the best-separated pair.  Real classes => separation >= 2 SD and
    overlap near zero.
(d) class-membership certainty for the best 1-D solution on the rate.

Inputs   Analysis/oai_imaging/cth_blup_fit.npz   (from 02_cth_blup_gmm.py)
Outputs  cth_slope_modality.csv          cth_silverman_test.csv
         cth_slope_gmm_selection.csv     cth_slope_gmm_components.csv
         fig_cth_slope_kde.png           fig_cth_modality.png
         fig_cth_slope_gmm.png           p2_continuum_run.log
"""
import os
import time
import numpy as np
import pandas as pd
from scipy import stats as sps
from sklearn.mixture import GaussianMixture
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"D:\BaiduSyncdisk\OAI\Analysis\oai_imaging"
NPZ = os.path.join(BASE, "cth_blup_fit.npz")
LONG = os.path.join(BASE, "cth_long.csv")
RNG = np.random.default_rng(20260914)

log = []
def p(*a):
    s = " ".join(str(x) for x in a)
    log.append(s)
    print(s, flush=True)

# ------------------------------------------------------------------ data ----
z = np.load(NPZ, allow_pickle=True)
uid = z["uid"].astype(str)
slope = z["eblup_slope"].astype(float)      # BLUP per-knee slope, points / year
oslope = z["ols_slope"].astype(float)       # unshrunk OLS slope
base = z["score00"].astype(float)           # baseline CTh-Score (level)

ok = np.isfinite(slope) & np.isfinite(oslope) & np.isfinite(base)
uid, slope, oslope, base = uid[ok], slope[ok], oslope[ok], base[ok]
n = len(slope)

# number of visits per knee -- used to build a precision-restricted OLS variable
_long = pd.read_csv(LONG, usecols=["uid", "t"])
ntp = _long.groupby("uid")["t"].nunique()
nv = pd.Series(ntp.reindex(uid).to_numpy(), index=range(len(uid)))
p("visits per knee: %s" % dict(sorted(pd.Series(nv).value_counts().items())))

VARS = [("slope_blup", slope),                 # shrunk, model-consistent
        ("slope_ols_ge4", oslope[nv.to_numpy() >= 4]),   # unshrunk, >= 4 visits
        ("baseline", base)]                    # severity level, as a contrast
V = dict(VARS)

p("knees = %d" % n)
for nm, v in VARS:
    p("%-11s mean %+9.3f   sd %8.3f   median %+8.3f   IQR %+7.3f..%+8.3f   "
      "range %+7.2f..%+8.2f   skew %+6.2f   exkurt %+6.2f"
      % (nm, v.mean(), v.std(ddof=1), np.median(v), *np.percentile(v, [25, 75]),
         v.min(), v.max(), sps.skew(v, bias=False), sps.kurtosis(v, bias=False)))
p("worsening (slope > 0): BLUP %.1f%%, OLS %.1f%%"
  % (100 * (slope > 0).mean(), 100 * (oslope > 0).mean()))

# ----------------------------------------------------------- KDE helpers ----
def kde_eval(x, grid, h):
    u = (grid[:, None] - x[None, :]) / h
    return np.exp(-0.5 * u * u).sum(axis=1) / (len(x) * h * np.sqrt(2 * np.pi))

def n_modes(x, h, grid):
    d = kde_eval(x, grid, h)
    m = (d[1:-1] > d[:-2]) & (d[1:-1] >= d[2:])
    return int(m.sum()), d

def rule_of_thumb(x):
    return 1.06 * x.std(ddof=1) * len(x) ** (-1 / 5)

def h_crit(x, grid, hi_mult=3.0, n_grid=60):
    """Smallest bandwidth at which the KDE has exactly one mode (Silverman 1981).

    Search ascends in h over [0.05 sd, hi_mult sd]; the mode count is monotone
    non-increasing in h, so the loop normally exits after a few evaluations.
    """
    sd = x.std(ddof=1)
    hs = np.exp(np.linspace(np.log(0.05 * sd), np.log(hi_mult * sd), n_grid))
    for h in hs:
        k, _ = n_modes(x, h, grid)
        if k <= 1:
            return h, hs
    return np.nan, hs

# ------------------------------------------------- (a) modality vs bandwidth -
# Each variable is swept on multiples of ITS OWN rule-of-thumb bandwidth; a shared
# absolute bandwidth would be meaningless because the three variables differ by more
# than an order of magnitude in sd.  h/sd is reported alongside for comparability.
p("\n== (a) modality as a function of bandwidth ==")
FACT = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0]
grids, rules = {}, {}
for nm, v in VARS:
    grids[nm] = np.linspace(v.min() - 3 * v.std(ddof=1) / 10,
                            v.max() + 3 * v.std(ddof=1) / 10, 512)
    rules[nm] = rule_of_thumb(v)
    p("%-11s rule-of-thumb h = %.4f  (sd %.3f)" % (nm, rules[nm], v.std(ddof=1)))

mod = []
for f in FACT:
    for nm, v in VARS:
        h = f * rules[nm]
        k, _ = n_modes(v, h, grids[nm])
        mod.append(dict(variable=nm, h_factor=f, bandwidth=h,
                        h_over_sd=h / v.std(ddof=1), n_modes=k))
moddf = pd.DataFrame(mod)
moddf.to_csv(os.path.join(BASE, "cth_slope_modality.csv"), index=False)
print(moddf.pivot(index="h_factor", columns="variable", values="n_modes").to_string())

first_uni = {}
for nm, v in VARS:
    s = moddf[(moddf.variable == nm) & (moddf.n_modes <= 1)].sort_values("h_factor")
    first_uni[nm] = (float(s.iloc[0].h_factor), float(s.iloc[0].h_over_sd)) if len(s) else None
    p("%-11s first unimodal at h = %s x its own rule of thumb (h/sd = %s)"
      % (nm, "%.1f" % first_uni[nm][0] if first_uni[nm] else "never",
         "%.3f" % first_uni[nm][1] if first_uni[nm] else "-"))
p("  the rule of thumb uses sd only, so it undersmooths heavy-tailed data; a rate that"
  " needs more smoothing than the level is showing tail weight, not extra modes.")

# ------------------------------------------------- (b) Silverman's test -----
# Smoothing bootstrap: under H0 the bootstrap distribution of h_crit is obtained by
# resampling with replacement from a kernel density estimate of the observed data at
# h_crit, then adding N(0, h_crit^2) jitter (Silverman 1981).  n is large, so the test
# runs on independent subsamples and the p-values are summarised by their median.
p("\n== (b) Silverman (1981) critical-bandwidth test, H0: unimodal ==")
N_SUB, B, REPS = 1500, 200, 4
SILCACHE = os.path.join(BASE, "cth_silverman_test.csv")

def silverman_once(xs, grid_):
    hc, _ = h_crit(xs, grid_)
    if not np.isfinite(hc):
        return np.nan, np.nan
    cnt = 0
    for _ in range(B):
        idx = RNG.integers(0, len(xs), len(xs))
        xb = xs[idx] + hc * RNG.standard_normal(len(xs))
        hb, _ = h_crit(xb, grid_)
        if np.isfinite(hb) and hb > hc:
            cnt += 1
    return hc, cnt / B

LEGACY = False
if os.path.exists(SILCACHE):
    _tmp = pd.read_csv(SILCACHE)
    LEGACY = "variable" not in _tmp.columns
if os.path.exists(SILCACHE) and not LEGACY:
    sildf = _tmp
    have = set(sildf.variable.unique())
    p("  [cache: %s]" % ", ".join(sorted(have)))
else:
    if LEGACY:
        p("  [cache dropped: written by an earlier version of the search grid]")
    sildf = pd.DataFrame(columns=["variable", "rep", "n_sub", "B", "h_crit",
                                  "p_value", "elapsed_s"])
    have = set()

t0 = time.time()
new = []
for nm, v in VARS:
    if nm in have:
        continue
    for r in range(REPS):
        xs = RNG.choice(v, N_SUB, replace=False)
        hc, pv = silverman_once(xs, np.linspace(v.min() - .1, v.max() + .1, 224))
        new.append(dict(variable=nm, rep=r, n_sub=N_SUB, B=B, h_crit=hc,
                        p_value=pv, elapsed_s=time.time() - t0))
        p("  %-11s rep %d : h_crit=%.4f (h/sd %.2f)  p=%.3f  (%.0fs)"
          % (nm, r, hc, hc / v.std(ddof=1), pv, time.time() - t0))
if new:
    sildf = pd.concat([sildf, pd.DataFrame(new)], ignore_index=True)
    sildf.to_csv(SILCACHE, index=False)

silk = sildf.groupby("variable").agg(h_crit=("h_crit", "median"),
                                     p_value=("p_value", "median"),
                                     p_min=("p_value", "min"),
                                     p_max=("p_value", "max"),
                                     reps=("rep", "size")).reset_index()
print(silk.to_string(index=False))
for _, r in silk.iterrows():
    p("  %-11s median p = %.3f (range %.3f-%.3f) -> %s"
      % (r.variable, r.p_value, r.p_min, r.p_max,
         "do NOT reject unimodality" if r.p_value > 0.05 else "reject unimodality"))

# ------------------------------------------------- (c) 1-D Gaussian mixture -
p("\n== (c) one-dimensional Gaussian mixture ==")

def gmm1d(v, Kmax=8):
    X = ((v - v.mean()) / v.std(ddof=1)).reshape(-1, 1)
    sd_raw, mu_raw = v.std(ddof=1), v.mean()
    rows, comps, fits = [], [], {}
    for k in range(1, Kmax + 1):
        g = GaussianMixture(n_components=k, covariance_type="full", n_init=20,
                            random_state=0, max_iter=1000, tol=1e-6).fit(X)
        resp = g.predict_proba(X)
        ent = np.nan if k == 1 else 1 - (-(resp * np.log(np.clip(resp, 1e-12, None))).sum()
                                         / (len(X) * np.log(k)))
        sz = np.bincount(g.predict(X), minlength=k)
        mu = g.means_.ravel() * sd_raw + mu_raw
        sd = np.sqrt(g.covariances_.ravel()) * sd_raw
        o = np.argsort(mu)
        sep_max = sep_min = ovl_best = np.nan
        if k > 1:
            sep = np.abs(np.diff(mu[o])) / np.sqrt((sd[o][:-1] ** 2 + sd[o][1:] ** 2) / 2)
            sep_max, sep_min = float(sep.max()), float(sep.min())
            ovl = np.array([2 * sps.norm.cdf(-abs(mu[o][i] - mu[o][i + 1]) /
                                             (2 * np.sqrt((sd[o][i] ** 2 + sd[o][i + 1] ** 2) / 2)))
                            for i in range(k - 1)])
            ovl_best = float(ovl.min())          # the least-overlapping pair
        rows.append(dict(K=k, BIC=g.bic(X), AIC=g.aic(X), entropy=ent,
                         min_pct=100 * sz.min() / len(X),
                         sizes="/".join(map(str, sz))))
        comps.append(dict(K=k, means=";".join("%.3f" % x for x in mu[o]),
                          sds=";".join("%.3f" % x for x in sd[o]),
                          max_separation=sep_max, min_separation=sep_min,
                          overlap_best_pair=ovl_best))
        fits[k] = g
    return pd.DataFrame(rows), pd.DataFrame(comps), fits

SEL, COMP, FITS = {}, {}, {}
for nm, v in VARS:
    SEL[nm], COMP[nm], FITS[nm] = gmm1d(v, Kmax=8)
    kb = int(SEL[nm].loc[SEL[nm].BIC.idxmin(), "K"])
    p("\n[%s]  BIC-selected K = %d" % (nm, kb))
    print(SEL[nm].to_string(index=False))
    p("[%s] component geometry, raw units (separation wants >= 2, overlap wants ~ 0):"
      % nm)
    print(COMP[nm][["K", "means", "sds", "max_separation", "min_separation",
                    "overlap_best_pair"]].to_string(index=False))
    r = COMP[nm][COMP[nm].K == kb].iloc[0]
    p("  -> K=%d : best-pair separation %.2f SD, remaining overlap %.0f%%, entropy %.3f,"
      " smallest class %.1f%%"
      % (kb, r.max_separation, 100 * r.overlap_best_pair,
         float(SEL[nm].loc[SEL[nm].K == kb, "entropy"].iloc[0]),
         float(SEL[nm].loc[SEL[nm].K == kb, "min_pct"].iloc[0])))
    r2 = COMP[nm][COMP[nm].K == 2].iloc[0]
    p("     the plain 2-component split of this variable has separation %.2f SD and "
      "overlap %.0f%%" % (r2.max_separation, 100 * r2.overlap_best_pair))
    if nm in ("slope_blup", "slope_ols_ge4"):
        tag = "" if nm == "slope_blup" else "_ols_ge4"
        SEL[nm].to_csv(os.path.join(BASE, "cth_slope%s_gmm_selection.csv" % tag), index=False)
        COMP[nm].to_csv(os.path.join(BASE, "cth_slope%s_gmm_components.csv" % tag), index=False)
p("  reference: a clinically usable 2-class split needs separation >= 2 SD and overlap"
  " near zero; ~45% overlap means the two components share nearly half their mass.")

# ------------------------------------------------- (d) classification check -
p("\n== (d) certainty of class membership for the progression rate ==")
for nm in ("slope_blup", "slope_ols_ge4"):
    v = V[nm]
    X = ((v - v.mean()) / v.std(ddof=1)).reshape(-1, 1)
    kb = int(SEL[nm].loc[SEL[nm].BIC.idxmin(), "K"])
    for k in sorted({2, 3, kb}):
        pm = FITS[nm][k].predict_proba(X).max(axis=1)
        p("  %-11s K=%d : mean max-posterior %.3f ; %.1f%% of knees below 0.80 certainty"
          " ; entropy %.3f"
          % (nm, k, pm.mean(), 100 * (pm < 0.80).mean(),
             float(SEL[nm].loc[SEL[nm].K == k, "entropy"].iloc[0])))

# ------------------------------------------------------------- figures -----
COL = {"slope_blup": "#1f4e79", "slope_ols_ge4": "#2E6DA4", "baseline": "#B03A2E"}
LBL = {"slope_blup": "progression rate, BLUP (shrunk)",
       "slope_ols_ge4": "progression rate, OLS (unshrunk, >=4 visits)",
       "baseline": "severity level, baseline score"}

v, nm = V["slope_blup"], "slope_blup"
g = grids[nm]
h_sel = [0.2 * rules[nm], 0.6 * rules[nm], 1.0 * rules[nm], 3.0 * rules[nm]]
fig, ax = plt.subplots(figsize=(8.2, 4.6), dpi=200)
ax.hist(v, bins=90, density=True, color="0.85", edgecolor="0.7", linewidth=.3,
        label="per-knee BLUP slope (histogram)")
for h, c in zip(h_sel, ["#1f4e79", "#2E6DA4", "#7BA7D7", "#B03A2E"]):
    k, d = n_modes(v, h, g)
    ax.plot(g, d, color=c, linewidth=1.5,
            label="KDE  h=%.2f (%.1f x rule)  modes=%d" % (h, h / rules[nm], k))
    loc = g[1:-1][(d[1:-1] > d[:-2]) & (d[1:-1] >= d[2:])]
    ax.plot(loc, np.interp(loc, g, d), "o", color=c, ms=5, mec="white", mew=.8)
ax.axvline(0, color="0.4", linestyle=":", linewidth=1)
ax.set_xlim(np.percentile(v, .2), np.percentile(v, 99.8))
ax.set_xlabel("CTh-Score progression slope (points / year)")
ax.set_ylabel("density")
ax.set_title("Modality of the per-knee CTh-Score progression rate", fontsize=10)
ax.legend(fontsize=7, frameon=False)
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "fig_cth_slope_kde.png"))
plt.close(fig)

fig, ax = plt.subplots(figsize=(7.4, 4.3), dpi=200)
for nm, _ in VARS:
    s = moddf[moddf.variable == nm]
    ax.plot(s.h_over_sd, s.n_modes, "-o", color=COL[nm], ms=4, label=LBL[nm])
ax.axhline(1, color="0.5", linestyle="--", linewidth=1)
ax.text(0.0065, 1.25, "unimodal", fontsize=8, color="0.35")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("bandwidth as a fraction of the variable's own standard deviation (log)")
ax.set_ylabel("number of modes of the KDE (log)")
ax.set_title("No stable multi-modal structure at any smoothing level", fontsize=10)
ax.legend(fontsize=7.5, frameon=False)
ax.grid(alpha=.25, linewidth=.4, which="both")
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "fig_cth_modality.png"))
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.2), dpi=200, sharey=False)
for ax, nm in zip(axes, ("slope_blup", "slope_ols_ge4")):
    v = V[nm]
    ax.hist(v, bins=90, density=True, color="0.88", edgecolor="0.72", linewidth=.3)
    xs = np.linspace(np.percentile(v, .05), np.percentile(v, 99.95), 800)
    for k, c in zip([1, 2, 3], ["#1f4e79", "#2E6DA4", "#B03A2E"]):
        gg = FITS[nm][k]
        mu = gg.means_.ravel() * v.std(ddof=1) + v.mean()
        sd = np.sqrt(gg.covariances_.ravel()) * v.std(ddof=1)
        w = gg.weights_.ravel()
        d = sum(wi * np.exp(-0.5 * ((xs - mi) / si) ** 2) / (si * np.sqrt(2 * np.pi))
                for wi, mi, si in zip(w, mu, sd))
        ax.plot(xs, d, color=c, linewidth=1.5,
                label="K=%d  BIC %.0f  entropy %.3f"
                      % (k, float(SEL[nm].loc[SEL[nm].K == k, "BIC"].iloc[0]),
                         float(SEL[nm].loc[SEL[nm].K == k, "entropy"].iloc[0])))
        for mi in mu:
            ax.axvline(mi, color=c, linestyle="--", linewidth=.7, alpha=.6)
    ax.set_xlim(np.percentile(v, .2), np.percentile(v, 99.8))
    ax.set_xlabel("progression slope (points / year)")
    ax.set_title(LBL[nm], fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=.2, linewidth=.4)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
axes[0].set_ylabel("density")
fig.suptitle("Mixture components fitted to the progression rate overlap rather than "
             "separate", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, .95])
fig.savefig(os.path.join(BASE, "fig_cth_slope_gmm.png"))
plt.close(fig)

open(os.path.join(BASE, "p2_continuum_run.log"), "w", encoding="utf-8").write("\n".join(log))
p("\n-> cth_slope_modality.csv / cth_silverman_test.csv / cth_slope_gmm_selection.csv")
p(" / cth_slope_gmm_components.csv / cth_slope_ols_gmm_selection.csv")
p(" / cth_slope_ols_gmm_components.csv")
p(" / fig_cth_slope_kde.png / fig_cth_modality.png / fig_cth_slope_gmm.png")
