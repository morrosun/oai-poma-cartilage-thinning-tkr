# -*- coding: utf-8 -*-
"""
POMA TKR nested case-control — reshape to long + pilot analysis (v2)
====================================================================
Input : Dataset/OAIZIB/subject_info_source/xlsx/kmri_poma_tkr_chondrometrics.xlsx
Output: Analysis/poma_pilot/
        poma_long.csv / poma_wide.csv / poma_unit_slopes.csv
        pilot_report.txt / fig_slope_distribution.png

Fixes vs v1
  * conditional logistic: correct within-stratum Hessian (A - B), Newton step
    solved with +H, covariance = inv(H)  (v1 diverged)
  * matched pairs now required to be COMPLETE (both members present with a
    finite slope); v1 let single-member strata through
  * 2-point OLS slopes are extremely noisy -> added mixed-model (random
    intercept + random slope) EBLUP shrinkage slopes as the primary exposure
"""

import os
import warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

BASE = r"D:\BaiduSyncdisk\OAI"
SRC = os.path.join(BASE, "Dataset", "OAIZIB", "subject_info_source", "xlsx",
                   "kmri_poma_tkr_chondrometrics.xlsx")
OUT = os.path.join(BASE, "Analysis", "poma_pilot")
os.makedirs(OUT, exist_ok=True)

rep = []
def P(s=""):
    print(s); rep.append(str(s))

# ============================================================ 1. long format
df = pd.read_excel(SRC)
P("=" * 80)
P("POMA TKR nested case-control — long reshape + pilot  (v2)")
P("=" * 80)
P(f"raw shape = {df.shape}")

VISIT_MONTHS = {0: 0, 1: 12, 3: 24, 5: 36, 6: 48}
IDVARS = ["id", "side", "version", "visit", "class", "tmpt", "newstrata"]

long = df.copy()
long["visit_code"] = long["visit"]
long["months"] = long["visit"].map(VISIT_MONTHS)
long["years"] = long["months"] / 12.0
long["unit"] = long["id"].astype(str) + "_" + long["side"].astype(str)
long["case"] = (long["class"] == "TKR").astype(int)
long["preop"] = (long["visit_code"] < long["tmpt"]).astype(int)

P("\n[check] every row is pre-index (visit < tmpt)? "
  f"{bool(long['preop'].all())}   violations = {int((long['preop']==0).sum())}")
P("[check] month mapping complete ? "
  f"{bool(long['months'].notna().all())}")

metric_cols = [c for c in df.columns if c not in IDVARS]
meta = IDVARS + ["visit_code", "months", "years", "unit", "case", "preop"]
long = long[meta + metric_cols].sort_values(["unit", "months"]).reset_index(drop=True)
long.to_csv(os.path.join(OUT, "poma_long.csv"), index=False)
P(f"\n-> poma_long.csv        {long.shape[0]} rows x {long.shape[1]} cols")

KEEP = ["MFTC_ThCtAB_aMe", "cMFTC_ThCtAB_aMe", "MFTC_VC",
        "MT_ThCtAB_aMe", "cMF_ThCtAB_aMe"]
w = long.pivot_table(index="unit", columns="visit_code", values=KEEP, aggfunc="first")
w.columns = [f"{m}_v{v}" for m, v in w.columns]
w.reset_index().to_csv(os.path.join(OUT, "poma_wide.csv"), index=False)
P(f"-> poma_wide.csv        {w.shape[0]} knees")

# ============================================================ 2. per-knee slopes
METRICS = [
    ("MFTC_ThCtAB_aMe",   "PRIMARY  MFTC mean thickness (medial femorotibial)"),
    ("cMFTC_ThCtAB_aMe",  "central weight-bearing MFTC, mean thickness"),
    ("cMFTC_ThCtAB_aMiv", "central weight-bearing MFTC, minimum thickness"),
    ("MT_ThCtAB_aMe",     "medial tibia, mean thickness"),
    ("cMF_ThCtAB_aMe",    "central medial femur, mean thickness"),
    ("LT_ThCtAB_aMe",     "lateral tibia, mean thickness"),
    ("cLF_ThCtAB_aMe",    "central lateral femur, mean thickness"),
    ("cLFTC_ThCtAB_aMe",  "central weight-bearing LFTC, mean thickness"),
    ("MFTC_VC",           "MFTC cartilage volume"),
]

def ols_slope(t, y):
    m = np.isfinite(t) & np.isfinite(y)
    t, y = np.asarray(t)[m], np.asarray(y)[m]
    if len(t) < 2 or np.ptp(t) == 0:
        return np.nan, np.nan
    b, a = np.polyfit(t, y, 1)
    return b, float(y[np.argmin(t)])

rows = []
for unit, g in long.groupby("unit", sort=False):
    g = g.sort_values("months")
    t = g["months"].to_numpy() / 12.0
    r = dict(unit=unit, id=g["id"].iloc[0], side=g["side"].iloc[0],
             cls=g["class"].iloc[0], case=int(g["case"].iloc[0]),
             newstrata=g["newstrata"].iloc[0], tmpt=g["tmpt"].iloc[0],
             n_visits=len(g),
             span_yr=float(g["months"].max() - g["months"].min()) / 12)
    for m, _ in METRICS:
        b, base = ols_slope(t, g[m].to_numpy(float))
        r[f"{m}__slope"] = b
        r[f"{m}__base"] = base
    rows.append(r)

sl = pd.DataFrame(rows)
PRIM = "MFTC_ThCtAB_aMe"

# ---------- mixed-model EBLUP slopes (random intercept + random slope) ------
P("\n" + "-" * 80)
P("2. per-knee slopes: OLS  vs  mixed-model EBLUP (shrunk)")
P("-" * 80)
import statsmodels.formula.api as smf

# y is standardised before fitting (volume metrics are ~2900 units and make
# the likelihood surface badly conditioned); BLUPs are back-transformed.
diag = []
for m, _ in METRICS:
    sub = long[["unit", "years", m]].dropna().rename(columns={m: "y"})
    mu, sg = sub["y"].mean(), sub["y"].std(ddof=1)
    sub["yz"] = (sub["y"] - mu) / sg
    md, used = None, None
    for meth in ("powell", "lbfgs", "bfgs"):
        try:
            cand = smf.mixedlm("yz ~ years", sub, groups=sub["unit"],
                               re_formula="~years").fit(reml=True,
                                                        method=meth,
                                                        maxiter=800)
            if cand.converged:
                md, used = cand, meth
                break
            md = md or cand
        except Exception:
            continue
    if md is None:
        P(f"   [mixedlm failed for {m}]")
        sl[f"{m}__eblup"] = np.nan
        continue
    fe = md.fe_params["years"] * sg
    re = pd.DataFrame(md.random_effects).T
    blup = (fe / sg + re["years"]) * sg
    sl[f"{m}__eblup"] = sl["unit"].map(blup)
    sl[f"{m}__fx"] = fe
    vc = md.cov_re.to_numpy() * sg ** 2
    diag.append(dict(metric=m, optimizer=used, converged=bool(md.converged),
                     sd_slope_true=np.sqrt(max(vc[1, 1], 0)),
                     sd_resid=np.sqrt(md.scale * sg ** 2),
                     fx_slope=fe))
    P(f"   {m:<20s} {used:<7s} conv={bool(md.converged)}  "
      f"SD(true slope)={np.sqrt(max(vc[1,1],0)):.4f} mm/yr  "
      f"SD(resid)={np.sqrt(md.scale*sg**2):.4f}")

sl.to_csv(os.path.join(OUT, "poma_unit_slopes.csv"), index=False)
P(f"-> poma_unit_slopes.csv {sl.shape[0]} knees x {sl.shape[1]} cols")
P("\nvariance components (mm/yr): SD(true between-knee slope) vs SD(residual)")
P(pd.DataFrame(diag).to_string(index=False,
                               float_format=lambda v: f"{v:.4g}"))

o = sl[f"{PRIM}__slope"]; e = sl[f"{PRIM}__eblup"]
P(f"\nfixed-effect mean slope (whole cohort) : {sl[f'{PRIM}__fx'].iloc[0]:+.4f} "
  f"mm/yr  ->  {abs(sl[f'{PRIM}__fx'].iloc[0])/sl[f'{PRIM}__base'].mean()*100:.2f}%/yr")
P(f"OLS  slope: SD = {o.std(ddof=1):.4f} mm/yr, IQR = "
  f"[{o.quantile(.25):.4f}, {o.quantile(.75):.4f}]")
P(f"EBLUP slope: SD = {e.std(ddof=1):.4f} mm/yr, IQR = "
  f"[{e.quantile(.25):.4f}, {e.quantile(.75):.4f}]")
P("   -> EBLUP SD is far smaller: the OLS spread is dominated by "
  "measurement noise, so ORs on OLS slopes are attenuated toward 1.")

P(f"\nknees with >=2 visits: {int((sl['n_visits']>=2).sum())}/{len(sl)}   "
  f">=3 visits: {int((sl['n_visits']>=3).sum())}/{len(sl)}")

# ============================================================ 3. matching
P("\n" + "-" * 80)
P("3. matching structure (knee level)")
P("-" * 80)
P(f"strata: {sl['newstrata'].nunique()};  all exactly 1 case + 1 control: "
  f"{bool((sl.groupby('newstrata')['case'].agg(['sum','size']) == [1,2]).all().all())}")
dup = sl.groupby("id")["unit"].nunique()
P(f"persons contributing 2 knees: {int((dup>1).sum())} / {sl['id'].nunique()}")

# ============================================================ 4. clogit engine
def clogit(d, xcols, groupcol="newstrata", case_col="case"):
    """Exact conditional logistic (one case per stratum).
    Newton-Raphson with the correct within-stratum information."""
    g = pd.factorize(d[groupcol])[0]
    y = d[case_col].to_numpy(float)
    X = d[xcols].to_numpy(float)
    n, p = X.shape
    ng = g.max() + 1
    M = np.zeros((ng, n)); M[g, np.arange(n)] = 1.0
    beta = np.zeros(p)
    H = None
    for _ in range(200):
        eta = X @ beta
        gm = np.full(ng, -np.inf); np.maximum.at(gm, g, eta)
        s = gm[g]
        ee = np.exp(eta - s)
        denom = M.T @ (M @ ee)
        pr = ee / denom                       # within-stratum softmax
        grad = X.T @ (y - pr)
        WX = X * pr[:, None]
        A = X.T @ WX                          # sum_j p_j x x'
        U = M @ WX                            # per stratum sum_j p_j x
        H = A - U.T @ U                       # within-stratum covariance
        try:
            step = np.linalg.solve(H, grad)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(H, grad, rcond=None)[0]
        step = np.clip(step, -5, 5)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-10:
            break
    cov = np.linalg.pinv(H)
    se = np.sqrt(np.abs(np.diag(cov)))
    return beta, se, int(d[groupcol].nunique()), int(n)

def build(slope_col, min_visits=2, one_knee_per_person=False):
    """analysis set = COMPLETE matched pairs with finite slope + baseline"""
    m = slope_col.replace("__slope", "").replace("__eblup", "")
    s = sl[sl["n_visits"] >= min_visits].copy()
    if one_knee_per_person:
        s = s[~s["id"].isin(dup[dup > 1].index)]
    gg = s.groupby("newstrata")
    ok = (gg["case"].transform("nunique") == 2) & (gg["case"].transform("size") == 2)
    s = s[ok]
    s = s[s.groupby("newstrata")[slope_col].transform(lambda v: v.notna().all())]
    s = s.dropna(subset=[slope_col, f"{m}__base"])
    s = s[s.groupby("newstrata")[slope_col].transform("size") == 2]
    # thinning rate: positive = losing cartilage faster
    s["loss"] = -s[slope_col]
    sd = s["loss"].std(ddof=1)
    s["z_loss"] = (s["loss"] - s["loss"].mean()) / sd
    step = 10.0 if "VC" in m else 0.1      # 10 mm3/yr or 0.1 mm/yr
    s["loss01"] = s["loss"] / step
    bs = f"{m}__base"
    s["z_base"] = (s[bs] - s[bs].mean()) / s[bs].std(ddof=1)
    return s, sd

def fit(s, xcols):
    b, se, nk, n = clogit(s, xcols)
    z = b / se
    p = 2 * stats.norm.sf(np.abs(z))
    return dict(beta=b[0], se=se[0], or_=np.exp(b[0]),
                lo=np.exp(b[0] - 1.96 * se[0]), hi=np.exp(b[0] + 1.96 * se[0]),
                p=p[0], n_strata=nk, n_knees=n)

P("\n" + "-" * 80)
P("4. conditional logistic — cartilage thinning RATE vs KR (matched strata)")
P("-" * 80)
P("exposure z-standardised as THINNING RATE (positive = faster loss)")
P("=> OR > 1  means  faster thinning raises the odds of knee replacement\n")

res = []
for m, lab in METRICS:
    for kind, col in [("OLS", f"{m}__slope"), ("EBLUP", f"{m}__eblup")]:
        if col not in sl.columns or sl[col].notna().sum() < 20:
            continue
        s, sd = build(col)
        if s["newstrata"].nunique() < 20:
            continue
        r = fit(s, ["z_loss"])
        s2, _ = build(col)
        r2 = fit(s2, ["z_loss", "z_base"])
        s3, _ = build(col)
        r3 = fit(s3, ["loss01"])
        step = 10.0 if "VC" in m else 0.1
        res.append(dict(metric=m, label=lab, est=kind, sd_loss=sd,
                        step=step, strata=r["n_strata"],
                        OR_SD=r["or_"], lo=r["lo"], hi=r["hi"], p=r["p"],
                        OR_SD_adj=r2["or_"], p_adj=r2["p"],
                        OR_01=r3["or_"], lo01=r3["lo"], hi01=r3["hi"]))

tab = pd.DataFrame(res)
# BH across the 9 EBLUP (primary estimator) p-values
mask = tab["est"] == "EBLUP"
p = tab.loc[mask, "p"].to_numpy()
o_ = np.argsort(p); q = np.empty_like(p); k = len(p)
q[o_] = np.minimum.accumulate((p[o_] * k / (np.arange(k) + 1))[::-1])[::-1]
tab.loc[mask, "q_BH"] = np.clip(q, 0, 1)

tab["unitlab"] = np.where(tab["step"] == 0.1, "OR/0.1mm.yr",
                          "OR/10mm3.yr")
sh = tab[["metric", "est", "sd_loss", "strata", "OR_SD", "lo", "hi", "p",
          "q_BH", "OR_SD_adj", "p_adj", "unitlab", "OR_01", "lo01", "hi01"]].copy()
sh.columns = ["metric", "est", "SD(loss)", "pairs", "OR/SD", "lo", "hi", "P",
              "q_BH", "OR/SD adj*", "P adj*", "unit", "OR/unit", "lo", "hi"]
P(sh.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
P("\n* additionally adjusted for the baseline value of the same metric "
  "(per 1 SD).\n  SD(loss) is in mm/yr (thickness) or mm3/yr (volume).")
tab.to_csv(os.path.join(OUT, "pilot_clogit_results.csv"), index=False)
P("-> pilot_clogit_results.csv")

# ============================================================ 5. distributions
P("\n" + "-" * 80)
P("5. slope distribution — primary metric (EBLUP, complete pairs)")
P("-" * 80)
s, sd_e = build(f"{PRIM}__eblup")
s_ols, sd_o = build(f"{PRIM}__slope")

for nm, ss, col in [("EBLUP", s, "loss"), ("OLS", s_ols, "loss")]:
    piv = ss.pivot_table(index="newstrata", columns="case", values=col)
    piv.columns = ["control", "tkr"]; piv = piv.dropna()
    d = piv["tkr"] - piv["control"]
    tt, pt = stats.ttest_rel(piv["tkr"], piv["control"])
    ww, pw = stats.wilcoxon(piv["tkr"], piv["control"])
    conc = float((piv["tkr"] > piv["control"]).mean())
    P(f"\n-- {nm}  (complete matched pairs n = {len(piv)}) --")
    P(f"   thinning rate  control {piv['control'].mean():+.4f} "
      f"(SD {piv['control'].std(ddof=1):.4f}) mm/yr")
    P(f"                  TKR     {piv['tkr'].mean():+.4f} "
      f"(SD {piv['tkr'].std(ddof=1):.4f}) mm/yr")
    P(f"   paired difference {d.mean():+.4f} "
      f"(95% CI {d.mean()-1.96*d.sem():+.4f} to {d.mean()+1.96*d.sem():+.4f})")
    P(f"   paired t = {tt:.3f}, P = {pt:.3g}    Wilcoxon P = {pw:.3g}")
    P(f"   KR knee thinning FASTER in {conc*100:.1f}% of pairs")
    if nm == "EBLUP":
        piv_e, d_e, pw_e, conc_e = piv, d, pw, conc

P(f"\ncohort baseline MFTC thickness = "
  f"{sl[f'{PRIM}__base'].mean():.2f} mm")
P(f"   => control loses {s[s.case==0]['loss'].mean()/sl[f'{PRIM}__base'].mean()*100:.2f}%/yr, "
  f"TKR loses {s[s.case==1]['loss'].mean()/sl[f'{PRIM}__base'].mean()*100:.2f}%/yr")
P(f"   SRM (mean/SD of within-person change) is the right benchmark for "
  f"measurement error; SD of paired difference = {d_e.std(ddof=1):.4f} mm/yr")

# ============================================================ 6. sensitivity
P("\n" + "-" * 80)
P("6. sensitivity — primary metric")
P("-" * 80)
sens = []
for lab, kw in [("EBLUP, >=2 visits (primary)", dict(c=f"{PRIM}__eblup")),
                ("EBLUP, >=3 visits", dict(c=f"{PRIM}__eblup", min_visits=3)),
                ("EBLUP, one knee/person",
                 dict(c=f"{PRIM}__eblup", one_knee_per_person=True)),
                ("OLS,   >=2 visits", dict(c=f"{PRIM}__slope")),
                ("OLS,   >=3 visits", dict(c=f"{PRIM}__slope", min_visits=3))]:
    c = kw.pop("c")
    ss, _ = build(c, **kw)
    r = fit(ss, ["z_loss"])
    sens.append(dict(analysis=lab, pairs=r["n_strata"],
                     OR_SD=r["or_"], lo=r["lo"], hi=r["hi"], P=r["p"]))
P(pd.DataFrame(sens).to_string(index=False, float_format=lambda v: f"{v:.4g}"))

# cross-check vs statsmodels
try:
    from statsmodels.discrete.conditional_models import ConditionalLogit
    ss, _ = build(f"{PRIM}__eblup")
    mod = ConditionalLogit(ss["case"].to_numpy(float),
                           ss[["z_loss"]].to_numpy(float),
                           groups=ss["newstrata"].to_numpy()).fit(disp=0)
    own = fit(ss, ["z_loss"])
    P(f"\n[cross-check] statsmodels ConditionalLogit beta = {mod.params[0]:+.5f}"
      f" (SE {mod.bse[0]:.5f}); own implementation beta = {own['beta']:+.5f}"
      f" (SE {own['se']:.5f})")
except Exception as ex:
    P(f"\n[cross-check] statsmodels failed: {ex}")

# ============================================================ 7. figure
try:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
    bins = np.histogram_bin_edges(
        pd.concat([piv_e["control"], piv_e["tkr"]]), bins=26)
    ax[0].hist(piv_e["control"], bins=bins, color="#4C72B0", alpha=.75,
               label="control")
    ax[0].hist(piv_e["tkr"], bins=bins, color="#C44E52", alpha=.75,
               label="TKR")
    ax[0].axvline(0, color="grey", lw=.8, ls="--")
    ax[0].set_xlabel("MFTC thinning rate (mm/yr)"); ax[0].set_ylabel("knees")
    ax[0].set_title("EBLUP slope distribution", fontsize=11)
    ax[0].legend(frameon=False, fontsize=9)

    ax[1].hist(d_e, bins=26, color="#55A868", alpha=.8)
    ax[1].axvline(0, color="grey", lw=.8, ls="--")
    ax[1].axvline(d_e.mean(), color="black", lw=1.4)
    ax[1].set_xlabel("within-pair difference (mm/yr)")
    ax[1].set_title(f"Paired difference\nmean {d_e.mean():+.3f}, "
                    f"P = {pw_e:.2g}", fontsize=11)

    ax[2].scatter(piv_e["control"], piv_e["tkr"], s=18, alpha=.6,
                  color="#8172B2", edgecolor="none")
    lim = [min(piv_e.min()), max(piv_e.max())]
    ax[2].plot(lim, lim, color="grey", lw=.8, ls="--")
    ax[2].set_xlabel("control thinning rate (mm/yr)")
    ax[2].set_ylabel("KR knee thinning rate (mm/yr)")
    ax[2].set_title(f"Matched pairs\n{conc_e*100:.1f}% above the diagonal",
                    fontsize=11)
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.suptitle("POMA nested case-control: pre-operative medial femorotibial "
                 "cartilage thinning (EBLUP slopes)", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_slope_distribution.png"), dpi=200)
    P("\n-> fig_slope_distribution.png")
except Exception as ex:
    P(f"\n[figure skipped] {ex}")

with open(os.path.join(OUT, "pilot_report.txt"), "w", encoding="utf-8") as f:
    f.write("\n".join(rep))
P("-> pilot_report.txt")
P(f"\nALL OUTPUTS: {OUT}")
