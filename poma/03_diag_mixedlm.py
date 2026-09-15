# -*- coding: utf-8 -*-
"""Diagnose MixedLM convergence + variance components for the POMA slope model."""
import warnings
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

warnings.simplefilter("always")

BASE = r"D:\BaiduSyncdisk\OAI"
long = pd.read_csv(BASE + r"\Analysis\poma_pilot\poma_long.csv")

METRICS = ["MFTC_ThCtAB_aMe", "cMFTC_ThCtAB_aMe", "cMFTC_ThCtAB_aMiv",
           "MT_ThCtAB_aMe", "cMF_ThCtAB_aMe", "LT_ThCtAB_aMe",
           "cLF_ThCtAB_aMe", "cLFTC_ThCtAB_aMe", "MFTC_VC"]

rows = []
for m in METRICS:
    sub = long[["unit", "years", m]].dropna().rename(columns={m: "y"})
    mu, sg = sub["y"].mean(), sub["y"].std(ddof=1)
    sub["yz"] = (sub["y"] - mu) / sg
    out = {"metric": m, "sd_y": sg}
    for meth in ["lbfgs", "powell", "bfgs", "cg"]:
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                md = smf.mixedlm("yz ~ years", sub, groups=sub["unit"],
                                 re_formula="~years").fit(reml=True,
                                                          method=meth,
                                                          maxiter=500)
                conv = md.converged
            vc = md.cov_re.to_numpy()
            out[f"{meth}_conv"] = conv
            out[f"{meth}_var_ri"] = float(vc[0, 0]) * sg ** 2
            out[f"{meth}_var_rs"] = float(vc[1, 1]) * sg ** 2
            out[f"{meth}_cov"] = float(vc[0, 1]) * sg ** 2
            out[f"{meth}_resid"] = float(md.scale) * sg ** 2
            out[f"{meth}_nwarn"] = len(w)
            if conv:
                out["best"] = meth
                out["beta_years"] = float(md.fe_params["years"]) * sg
                break
        except Exception as e:
            out[f"{meth}_conv"] = f"ERR {e}"
    rows.append(out)

d = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(d[["metric", "sd_y", "lbfgs_conv", "powell_conv", "best",
         "beta_years", "lbfgs_var_ri", "lbfgs_var_rs", "lbfgs_resid",
         "lbfgs_nwarn"]].to_string(index=False,
                                   float_format=lambda v: f"{v:.4g}"))
print("\n--- full variance components ---")
print(d[["metric", "best", "beta_years"] +
        [c for c in d.columns if "_var_" in c or "_resid" in c]]
      .to_string(index=False, float_format=lambda v: f"{v:.4g}"))
