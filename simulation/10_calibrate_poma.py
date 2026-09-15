# -*- coding: utf-8 -*-
"""
10 : calibrate the simulation to the real POMA design
=====================================================
Everything the simulation needs is taken from the observed data, not assumed:

  * the per-knee visit schedule of the 191-pair primary set (months),
    kept as (case, control) PAIRS so that any asymmetry in visit density
    between cases and controls is preserved;
  * the exposure distribution: mean / SD of the EBLUP and OLS cMFTC slopes;
  * the target log-OR per SD of the rate;
  * the observed latent-class fit (cMFTC LCMM K=1..6) for the continuum study.

Output -> Analysis/sim/
    sim_calibration.json
    sim_visit_schedule.csv      (pair, member, months 0/12/24/36/48 ...)
"""
import json, os
import numpy as np
import pandas as pd

BASE = r"D:\BaiduSyncdisk\OAI"
PP   = os.path.join(BASE, "Analysis", "poma_pilot")
OUT  = os.path.join(BASE, "Analysis", "sim")
os.makedirs(OUT, exist_ok=True)

PRIM = "cMFTC_ThCtAB_aMe"
L = pd.read_csv(os.path.join(PP, "poma_long.csv"))
U = pd.read_csv(os.path.join(PP, "poma_analysis_long.csv"))

# ---------------------------------------------------------------- primary set
sub = L[~L[PRIM].isna()].copy()
nv = sub.groupby("unit").size()
ok_knee = set(nv[nv >= 2].index)                      # >= 2 usable visits
U = U[U.unit.isin(ok_knee)].copy()                    # prune A
tb = U.groupby("newstrata").size()
U = U[U.newstrata.isin(tb[tb == 2].index)].copy()     # keep complete pairs
sub = sub[sub.unit.isin(set(U.unit))].copy()
print("primary set : knees=%d  pairs=%d  obs=%d"
      % (len(U), U.newstrata.nunique(), len(sub)))

# --------------------------------------------------- visit schedule per knee
sub = sub.sort_values(["unit", "months"])
sched = sub.groupby("unit")["months"].apply(lambda s: sorted(s.tolist())).to_dict()

cls = U.set_index("unit")["case"].to_dict()
strata = U.set_index("unit")["newstrata"].to_dict()
tmpt = U.set_index("unit")["t_index_months"].to_dict()

rows = []
for pair, g in U.groupby("newstrata"):
    mm = {int(r.case): r.unit for r in g.itertuples()}
    if len(mm) != 2:
        continue
    for role in (1, 0):                      # 1 = case, 0 = control
        u = mm[role]
        rows.append(dict(pair=int(pair), role=role, unit=u,
                         months=";".join(str(int(x)) for x in sched[u]),
                         t_index=tmpt[u], n_visits=len(sched[u])))
S = pd.DataFrame(rows)
S.to_csv(os.path.join(OUT, "sim_visit_schedule.csv"), index=False)
print("schedule rows = %d pairs = %d" % (len(S), S.pair.nunique()))

# case/control visit-density asymmetry (drives differential shrinkage)
nc = S[S.role == 1].n_visits
nk = S[S.role == 0].n_visits
print("visits  case: mean %.2f  control: mean %.2f  (paired mean diff %+.3f)"
      % (nc.mean(), nk.mean(),
         (nc.values - nk.values).mean()))

# ------------------------------------------------------------- exposure moments
sl = U["cMFTC_ThCtAB_aMe__slope"].dropna()
eb = U["cMFTC_ThCtAB_aMe__eblup"].dropna()
print("OLS  slope: mean %+.4f sd %.4f  %%negative %.1f" % (sl.mean(), sl.std(), (sl < 0).mean() * 100))
print("EBLUP slope: mean %+.4f sd %.4f  %%negative %.1f" % (eb.mean(), eb.std(), (eb < 0).mean() * 100))

# ------------------------------------------------------- target log-OR per SD
# R-independent conditional logistic (Newton), orientation: LOSS = -(slope)
def clogit_logOR(df, col, sign=-1.0):
    a = df.pivot_table(index="newstrata", columns="case", values=col, aggfunc="first").dropna()
    x = sign * (a[1] - a[0]).values              # case minus control, in LOSS units
    b = 0.0
    for _ in range(300):
        p = 1.0 / (1.0 + np.exp(-b * x))
        g = (x * (1 - p)).sum()
        h = (x ** 2 * p * (1 - p)).sum()
        b = b + g / h
    return b, 1.0 / np.sqrt(h), x.std(ddof=1)

b_eb, se_eb, sd_eb = clogit_logOR(U, "cMFTC_ThCtAB_aMe__eblup")
b_sl, se_sl, sd_sl = clogit_logOR(U, "cMFTC_ThCtAB_aMe__slope")
print("clogit  EBLUP : logOR/SD %+.4f (se %.4f) -> OR %.3f" % (b_eb, se_eb, np.exp(b_eb)))
print("clogit  OLS   : logOR/SD %+.4f (se %.4f) -> OR %.3f" % (b_sl, se_sl, np.exp(b_sl)))

# paired differences (LOSS units, positive = cases thinner)
def paired(df, col, sign=-1.0):
    a = df.pivot_table(index="newstrata", columns="case", values=col, aggfunc="first").dropna()
    d = sign * (a[1] - a[0])
    return d.mean(), d.std(), len(d)
print("paired diff  EBLUP loss %+.4f (sd %.4f)" % paired(U, "cMFTC_ThCtAB_aMe__eblup")[:2])
print("paired diff  OLS   loss %+.4f (sd %.4f)" % paired(U, "cMFTC_ThCtAB_aMe__slope")[:2])

CAL = dict(
    n_pairs=int(U.newstrata.nunique()),
    n_knees=int(len(U)),
    n_obs=int(len(sub)),
    visits_per_knee={str(k): int(v) for k, v in U.n_visits.value_counts().sort_index().items()},
    slope_mean=float(sl.mean()), slope_sd=float(sl.std(ddof=1)),
    eblup_mean=float(eb.mean()), eblup_sd=float(eb.std(ddof=1)),
    frac_worsening_ols=float((sl < 0).mean()),
    frac_worsening_eblup=float((eb < 0).mean()),
    base_mean=float(U["cMFTC_ThCtAB_aMe__base"].mean()),
    base_sd=float(U["cMFTC_ThCtAB_aMe__base"].std(ddof=1)),
    logOR_per_sd_eblup=float(b_eb), se_logOR_per_sd_eblup=float(se_eb),
    logOR_per_sd_ols=float(b_sl),   se_logOR_per_sd_ols=float(se_sl),
    paired_diff_eblup=float(paired(U, "cMFTC_ThCtAB_aMe__eblup")[0]),
    paired_diff_ols=float(paired(U, "cMFTC_ThCtAB_aMe__slope")[0]),
    case_mean_visits=float(nc.mean()), control_mean_visits=float(nk.mean()),
)
with open(os.path.join(OUT, "sim_calibration.json"), "w", encoding="utf-8") as f:
    json.dump(CAL, f, indent=2, ensure_ascii=False)
print("\nwritten -> sim_calibration.json")
for k, v in CAL.items():
    print("   %-26s %s" % (k, v))
