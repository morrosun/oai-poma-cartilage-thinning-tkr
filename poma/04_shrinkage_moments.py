# -*- coding: utf-8 -*-
"""
Moment-based cross-check of the EBLUP shrinkage.
Decomposition:  Var(OLS slope across knees) = Var(true slope) + E[Var(slope_hat)]
  Var(slope_hat_k) = sigma2_resid / Sxx_k
  sigma2_resid estimated by pooling within-knee residuals from knees with >=3 visits.
"""
import numpy as np
import pandas as pd

BASE = r"D:\BaiduSyncdisk\OAI"
long = pd.read_csv(BASE + r"\Analysis\poma_pilot\poma_long.csv")
sl = pd.read_csv(BASE + r"\Analysis\poma_pilot\poma_unit_slopes.csv")

METRICS = ["MFTC_ThCtAB_aMe", "cMFTC_ThCtAB_aMe", "cMFTC_ThCtAB_aMiv",
           "MT_ThCtAB_aMe", "cMF_ThCtAB_aMe", "LT_ThCtAB_aMe",
           "cLF_ThCtAB_aMe", "cLFTC_ThCtAB_aMe", "MFTC_VC"]
EBLUP_SD_OBS = {"MFTC_ThCtAB_aMe": 0.0984, "cMFTC_ThCtAB_aMe": None}

rows = []
for m in METRICS:
    ssr = 0.0; dof = 0
    for u, g in long.groupby("unit"):
        y = g[m].to_numpy(float); t = g["years"].to_numpy(float)
        mm = np.isfinite(y) & np.isfinite(t)
        y, t = y[mm], t[mm]
        if len(t) < 3:
            continue
        b, a = np.polyfit(t, y, 1)
        ssr += float(((y - (a + b * t)) ** 2).sum()); dof += len(t) - 2
    s2 = ssr / dof

    s = sl[sl["n_visits"] >= 2].copy()
    b, se, nk, n = 0, 0, 0, 0
    # Sxx per knee
    sxx = (long[long["unit"].isin(s["unit"])]
           .groupby("unit")["years"]
           .apply(lambda t: float(((t - t.mean()) ** 2).sum())))
    s["sxx"] = s["unit"].map(sxx)
    ok = s[f"{m}__slope"].notna() & (s["sxx"] > 0)
    ols = s.loc[ok, f"{m}__slope"]
    var_est = s2 / s.loc[ok, "sxx"]
    var_obs = ols.var(ddof=1)
    var_true = var_obs - var_est.mean()
    lam = var_true / (var_true + var_est)

    # EBLUP under the moment estimate
    mu = ols.mean()
    eb = mu + lam * (ols - mu)
    s.loc[ok, "eb"] = eb

    rows.append(dict(
        metric=m, sd_resid=np.sqrt(s2), dof=dof,
        sd_OLS=ols.std(ddof=1), mean_var_est=np.sqrt(var_est.mean()),
        sd_true_mom=np.sqrt(max(var_true, 0)),
        lambda_mean=lam.mean(), lambda_med=lam.median(),
        sd_EBLUP_mom=eb.std(ddof=1),
        sd_EBLUP_mixedlm=sl[f"{m}__eblup"].std(ddof=1),
    ))

d = pd.DataFrame(rows)
pd.set_option("display.width", 220)
print("moment-based variance decomposition of the per-knee slope\n")
print(d.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
print("\nlambda = Var_true/(Var_true+Var_est) = per-knee shrinkage weight")
print("sd_EBLUP_mom vs sd_EBLUP_mixedlm: if similar, the mixed-model EBLUPs "
      "are trustworthy")
