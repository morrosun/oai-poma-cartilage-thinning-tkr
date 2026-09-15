# -*- coding: utf-8 -*-
r"""
02 : CTh-Score progression classes -- BLUP-based latent class analysis (P2 robustness)

Why this script exists
----------------------
`01_cth_trajectory_lcmm.R` fits the full latent-class mixed model (lcmm::hlme),
which is the specification used in the authors' MCI / POMA manuscripts.
On this particular metric hlme lands in degenerate local optima (one class
absorbing all knees), which is a known failure mode when the outcome has a
floor at zero and very large between-knee level differences.

This script therefore provides an independent, distribution-light route to the
same question: two-step latent class analysis
    1. per-knee BLUP of (intercept, slope) from a linear mixed model
    2. Gaussian mixture / k-means on the standardised BLUP pair
Model selection by BIC over K = 2..6 plus an entropy criterion.

It is deliberately labelled as an approximation.  If the two routes agree on
the number and ordering of classes, the class structure is not an artefact of
the estimation strategy.

Inputs   Analysis/oai_imaging/cth_long.csv      (ID, SIDE, uid, t, score)
         Analysis/oai_imaging/zib_klinfo.csv    (base 00m, ... KLGrade)
Outputs  cth_blup_model_selection.csv
         cth_blup_class_profile.csv
         cth_blup_class_vs_KL.csv
         fig_cth_blup_scatter.png / fig_cth_blup_trajectories.png
         cth_blup_fit.npz
"""
import os
import json
import time
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = r"D:\BaiduSyncdisk\OAI\Analysis\oai_imaging"
LONG = os.path.join(BASE, "cth_long.csv")
KLF = os.path.join(BASE, "zib_klinfo.csv")

log = []
def p(*a):
    s = " ".join(str(x) for x in a)
    log.append(s)
    print(s, flush=True)

# ----------------------------------------------------------------- data -----
d = pd.read_csv(LONG)
d = d.dropna(subset=["t", "score"])
nv = d.groupby("uid")["t"].nunique()
keep = nv[nv >= 2].index
d = d[d["uid"].isin(keep)].sort_values(["uid", "t"])
p("knees with >= 2 visits : %d ; rows : %d" % (d["uid"].nunique(), len(d)))

# --------------------------------------- step 1 : mixed model -> BLUPs ------
t0 = time.time()
md = smf.mixedlm("score ~ t", d, groups=d["uid"], re_formula="~ t")
mf = md.fit(method="lbfgs", maxiter=200)
p("mixed model fitted in %.1fs ; converged=%s" % (time.time() - t0, mf.converged))
p("fixed effects : intercept %.3f, t %.4f" % (mf.fe_params["Intercept"], mf.fe_params["t"]))

re = mf.random_effects
uids = list(d["uid"].unique())
rows = []
for u in uids:
    r = re.get(u)
    if r is None:
        continue
    b = float(r.iloc[0])                 # random intercept
    s = float(r.iloc[1]) if len(r) > 1 else 0.0
    rows.append((u, b, s))
bl = pd.DataFrame(rows, columns=["uid", "re_int", "re_slope"])
bl["eblup_int"] = mf.fe_params["Intercept"] + bl["re_int"]
bl["eblup_slope"] = mf.fe_params["t"] + bl["re_slope"]
p("BLUPs extracted for %d knees" % len(bl))

# descriptive per-knee summaries
agg = d.groupby("uid").agg(score00=("score", "first"), n_tp=("score", "size"),
                           t_max=("t", "max"), mean_score=("score", "mean")).reset_index()
ols = (d.groupby("uid").apply(lambda x: np.polyfit(x["t"], x["score"], 1)[0])
         .rename("ols_slope").reset_index())
K = bl.merge(agg, on="uid").merge(ols, on="uid")
p("overall slope : mean %+.3f, sd %.3f ; worsening %.1f%%"
  % (K["ols_slope"].mean(), K["ols_slope"].std(), 100 * (K["ols_slope"] > 0).mean()))

# --------------------------------- step 2 : Gaussian mixture / k-means ------
X = StandardScaler().fit_transform(K[["eblup_int", "eblup_slope"]].to_numpy())
sel = []
fits = {}
for k in range(2, 9):
    gm = GaussianMixture(n_components=k, covariance_type="full", n_init=10,
                         random_state=0, max_iter=500).fit(X)
    resp = gm.predict_proba(X)
    ent = 1 - (-(resp * np.log(np.clip(resp, 1e-12, None))).sum() /
               (len(X) * np.log(k)))
    sizes = np.bincount(gm.predict(X), minlength=k)
    sel.append(dict(model="GMM", K=k, BIC=gm.bic(X), AIC=gm.aic(X), entropy=ent,
                    min_pct=100 * sizes.min() / len(X),
                    sizes="/".join(map(str, sizes))))
    fits[k] = gm
km = KMeans(n_clusters=4, n_init=25, random_state=0).fit(X)
sel.append(dict(model="kmeans", K=4, BIC=np.nan, AIC=np.nan, entropy=np.nan,
                min_pct=100 * np.bincount(km.labels_).min() / len(X),
                sizes="/".join(map(str, np.bincount(km.labels_)))))
seldf = pd.DataFrame(sel)
seldf.to_csv(os.path.join(BASE, "cth_blup_model_selection.csv"), index=False)
p("\n== model selection (BLUP-based) ==")
p(seldf.to_string(index=False))

kbest = int(seldf[seldf.model == "GMM"].sort_values("BIC").iloc[0]["K"])
gm = fits[kbest]
K["cl"] = gm.predict(X)
K["pmax"] = gm.predict_proba(X).max(axis=1)
p("\nBIC-selected K = %d (entropy %.3f)" % (kbest, seldf.loc[seldf.K == kbest, "entropy"].iloc[0]))

prof = (K.groupby("cl")
          .agg(n=("uid", "size"), slope_pts_yr=("eblup_slope", "mean"),
               slope_ols=("ols_slope", "mean"), baseline=("score00", "mean"),
               visits=("n_tp", "mean"), t_max=("t_max", "mean"))
          .reset_index())
prof["pct"] = 100 * prof["n"] / len(K)
prof = prof.sort_values("slope_pts_yr").reset_index(drop=True)
prof["label"] = ["C%d" % (i + 1) for i in range(len(prof))]
prof.to_csv(os.path.join(BASE, "cth_blup_class_profile.csv"), index=False)
p("\n== class profile (slowest to fastest) ==")
p(prof.to_string(index=False))

# ------------------------------------------------ external criterion: KL ---
kl = pd.read_csv(KLF)
K["ID"] = K["uid"].str.rsplit("_", n=1).str[0].astype(int)
K["SIDE"] = K["uid"].str.rsplit("_", n=1).str[1]
m = K[K["SIDE"] == "RIGHT"].merge(kl[["SubjectID", "KLGrade", "Age", "BMI"]],
                                 left_on="ID", right_on="SubjectID")
m = m.dropna(subset=["KLGrade"])
if len(m):
    byc = m.groupby("cl").agg(n=("uid", "size"), mean_KL=("KLGrade", "mean"),
                              mean_slope=("eblup_slope", "mean")).reset_index()
    byc.to_csv(os.path.join(BASE, "cth_blup_class_vs_KL.csv"), index=False)
    p("\n== baseline KL by class (right knees, n = %d) ==" % len(m))
    p(byc.to_string(index=False))
    from scipy.stats import spearmanr
    p("  Spearman (class mean KL vs class slope) = %+.3f"
      % spearmanr(byc["mean_KL"], byc["mean_slope"]).statistic)
    p("  Spearman (per-knee KL  vs BLUP slope)   = %+.3f"
      % spearmanr(m["KLGrade"], m["eblup_slope"]).statistic)

# ------------------------------------------------------------- figures -----
labmap = dict(zip(prof["cl"], prof["label"]))
K["lab"] = K["cl"].map(labmap)
cols = plt.cm.tab10(np.linspace(0, 1, 10))
order = list(prof["cl"])
fig, ax = plt.subplots(figsize=(7.2, 5.2), dpi=200)
for i, c in enumerate(order):
    s = K[K["cl"] == c]
    ax.scatter(s["score00"], s["eblup_slope"], s=6, alpha=.55,
               color=cols[i], label="%s (n=%d, %.2f pts/yr)" % (labmap[c], len(s),
                                                                prof.loc[i, "slope_pts_yr"]))
ax.set_xlabel("baseline CTh-Score (0-100, higher = worse)")
ax.set_ylabel("BLUP progression slope (points / year)")
ax.set_title("CTh-Score progression classes\n"
             "two-step latent class analysis on BLUPs (K = %d, BIC-selected)" % kbest,
             fontsize=10)
ax.legend(fontsize=7, frameon=False, loc="upper left")
ax.grid(alpha=.25, linewidth=.4)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "fig_cth_blup_scatter.png"))
plt.close(fig)

K["pred0"] = K["eblup_int"]
fig, ax = plt.subplots(figsize=(7.4, 4.4), dpi=200)
for i, c in enumerate(order):
    s = K[K["cl"] == c]
    tt = np.array([0, 8])
    yy = s["eblup_int"].mean() + s["eblup_slope"].mean() * tt
    ax.plot(tt, yy, color=cols[i], linewidth=2,
            label="%s  %+.2f pts/yr, base %.0f (n=%d)" %
                  (labmap[c], prof.loc[i, "slope_pts_yr"], prof.loc[i, "baseline"], len(s)))
    ax.fill_between(tt, yy - 1.96 * s["eblup_int"].std() * 0,
                    yy + 1.96 * s["eblup_int"].std() * 0, color=cols[i], alpha=0)
ax.set_xlabel("years from baseline")
ax.set_ylabel("CTh-Score (0-100)")
ax.set_title("Mean CTh-Score trajectory by class\n"
             "lines = class-mean BLUP intercept + BLUP slope", fontsize=10)
ax.legend(fontsize=7, frameon=False)
ax.grid(alpha=.25, linewidth=.4)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(BASE, "fig_cth_blup_trajectories.png"))
plt.close(fig)

np.savez(os.path.join(BASE, "cth_blup_fit.npz"),
         uid=K["uid"].to_numpy(), cl=K["cl"].to_numpy(),
         eblup_int=K["eblup_int"].to_numpy(), eblup_slope=K["eblup_slope"].to_numpy(),
         ols_slope=K["ols_slope"].to_numpy(), score00=K["score00"].to_numpy(),
         K=kbest)
open(os.path.join(BASE, "p2_blup_run.log"), "w", encoding="utf-8").write("\n".join(log))
p("\n-> cth_blup_model_selection.csv / cth_blup_class_profile.csv / cth_blup_class_vs_KL.csv")
p(" / fig_cth_blup_scatter.png / fig_cth_blup_trajectories.png / cth_blup_fit.npz")
