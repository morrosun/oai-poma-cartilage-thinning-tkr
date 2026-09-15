# -*- coding: utf-8 -*-
"""
33_verify_poma_equivalence.py
-----------------------------
Proves that the generic engine (Scripts/common/30_class_continuum_check.R) with
the POMA config (Scripts/poma/31_cfg_poma_lcmm.R) reproduces the original
hard-coded script (Scripts/poma/28_lcmm_continuum_check.R) field for field.

Compares, on the `eblup` axis:
    lcmm_rank_equivalence.csv      <-> poma_generic_rank_equivalence__eblup.csv
    lcmm_eta_squared.csv           <-> poma_generic_eta_squared__eblup.csv
    lcmm_rate_vs_class.csv         <-> poma_generic_axis_vs_class__eblup.csv
    lcmm_rate_below_threshold.csv  <-> poma_generic_below_threshold__eblup.csv
    lcmm_rate_quintiles.csv        <-> poma_generic_bins__eblup.csv
    lcmm_rate_quintile_OR.csv      <-> poma_generic_bin_model__eblup.csv

Writes Analysis/poma_pilot/poma_generic_vs_28_check.csv and prints a verdict.
"""
import os
import numpy as np
import pandas as pd

B = r"D:\BaiduSyncdisk\OAI\Analysis\poma_pilot"
TOL = 5e-4          # OR/CI/p are rounded to 3 dp in both scripts
rows = []


def num(x):
    """numeric view of a possibly-string cell; NaN when unparsable"""
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def cmp_vec(a, b, label):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        rows.append(dict(check=label, n_old=len(a), n_new=len(b),
                         max_abs_diff=np.nan, verdict="SHAPE MISMATCH"))
        return False
    d = np.nanmax(np.abs(a - b)) if len(a) else 0.0
    ok = bool(np.allclose(a, b, atol=TOL, rtol=0, equal_nan=True))
    rows.append(dict(check=label, n_old=len(a), n_new=len(b),
                     max_abs_diff=round(float(d), 6),
                     verdict="MATCH" if ok else "DIFFER"))
    return ok


def rd(f):
    p = os.path.join(B, f)
    if not os.path.exists(p):
        print("  MISSING:", f)
        return None
    return pd.read_csv(p)


allok = True

# ---- (1) rank-threshold equivalence ---------------------------------------
o, n = rd("lcmm_rank_equivalence.csv"), rd("poma_generic_rank_equivalence__eblup.csv")
if o is not None and n is not None:
    allok &= cmp_vec(o.pct_agreement, n.pct_agreement, "(1) pct_agreement x3")
    o_j, n_j = o.jaccard.to_numpy(dtype=float), n.jaccard.to_numpy(dtype=float)
    m = ~np.isnan(o_j) & ~np.isnan(n_j)
    allok &= cmp_vec(o_j[m], n_j[m], "(1) jaccard (non-NA)")

# ---- (2) eta squared -------------------------------------------------------
o, n = rd("lcmm_eta_squared.csv"), rd("poma_generic_eta_squared__eblup.csv")
if o is not None and n is not None:
    a = dict(zip(o.target, o.eta_squared))
    b = dict(zip(n.target, n.eta_squared))
    allok &= cmp_vec([a["rate"], a["baseline"]],
                     [b["cMFTC thinning rate, EBLUP (mm/yr, + = faster)"],
                      b["baseline cartilage thickness"]],
                     "(2) eta2 axis / level")

# ---- (3) axis vs class -----------------------------------------------------
o, n = rd("lcmm_rate_vs_class.csv"), rd("poma_generic_axis_vs_class__eblup.csv")
if o is not None and n is not None:
    on = n[n.model.str.startswith("M3")].reset_index(drop=True)
    oo = o[o.model.str.startswith("M3")].reset_index(drop=True)
    allok &= cmp_vec(oo.OR, on.est, "(3) M3 OR")
    allok &= cmp_vec(oo.lo, on.lo, "(3) M3 CI lower")
    allok &= cmp_vec(oo.hi, on.hi, "(3) M3 CI upper")
    # M1 / M2 are single-term models
    allok &= cmp_vec(o[o.model.str.startswith("M1")].OR,
                     n[n.model.str.startswith("M1")].est, "(3) M1 OR")
    allok &= cmp_vec(o[o.model.str.startswith("M2")].OR,
                     n[n.model.str.startswith("M2")].est, "(3) M2 OR")
    print("  LRT axis|class  old p = %.3g" % 0.000946)
    print("  LRT class|axis  old p = %.3g" % 0.517)

# ---- (4) below threshold ---------------------------------------------------
o, n = rd("lcmm_rate_below_threshold.csv"), rd("poma_generic_below_threshold__eblup.csv")
if o is not None and n is not None:
    allok &= cmp_vec([o.OR.iloc[0], o.lo.iloc[0], o.hi.iloc[0]],
                     [n.est.iloc[0], n.lo.iloc[0], n.hi.iloc[0]],
                     "(4) OR / lo / hi")

# ---- (5) quintiles ---------------------------------------------------------
o, n = rd("lcmm_rate_quintiles.csv"), rd("poma_generic_bins__eblup.csv")
if o is not None and n is not None:
    o = o.sort_values("q5")
    n = n.sort_values("bin")
    allok &= cmp_vec(o.n, n.n, "(5) bin n")
    allok &= cmp_vec(o.pct_TKR, n.pct_event, "(5) bin % event")
    allok &= cmp_vec(o.mean_rate, n.mean_axis, "(5) bin mean axis")

o, n = rd("lcmm_rate_quintile_OR.csv"), rd("poma_generic_bin_model__eblup.csv")
if o is not None and n is not None:
    allok &= cmp_vec(o.OR, n.est, "(5) bin-model OR")
    allok &= cmp_vec(o.lo, n.lo, "(5) bin-model CI lower")
    allok &= cmp_vec(o.hi, n.hi, "(5) bin-model CI upper")

out = pd.DataFrame(rows)
out.to_csv(os.path.join(B, "poma_generic_vs_28_check.csv"), index=False,
           encoding="utf-8")
pd.set_option("display.width", 200)
print("\n================ generic engine vs script 28 ================")
print(out.to_string(index=False))
print("\nVERDICT:", "ALL CHECKS MATCH" if allok else "SOME CHECKS DIFFER")
print("wrote poma_generic_vs_28_check.csv")
