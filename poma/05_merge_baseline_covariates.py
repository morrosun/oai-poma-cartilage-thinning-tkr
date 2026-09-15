# -*- coding: utf-8 -*-
"""
Merge OAI baseline covariates onto the POMA knee-level slope table.

sources (all local, no NDA request needed)
  subjectchar00.xlsx : V00AGE
  enrollees.xlsx     : P02SEX, V00COHORT
  allclinical00.xlsx : P01BMI, P01HEIGHT, P01WEIGHT,
                       V00WOMKPL/R (WOMAC pain L/R), V00WOMTSL/R (WOMAC total),
                       V00PASE
  kxr_sq_bu00.xlsx   : V00XRKL  (baseline KL grade, per ID x SIDE)

output: Analysis/poma_pilot/poma_analysis_long.csv  (knee-level)
        Analysis/poma_pilot/poma_visits_long.csv    (visit-level, for joint model)
"""
import os
import numpy as np
import pandas as pd

B = r"D:\BaiduSyncdisk\OAI\Dataset\OAIZIB\subject_info_source\xlsx\\"
OUT = r"D:\BaiduSyncdisk\OAI\Analysis\poma_pilot"
rep = []


def P(s=""):
    print(s); rep.append(str(s))


sl = pd.read_csv(os.path.join(OUT, "poma_unit_slopes.csv"))
long = pd.read_csv(os.path.join(OUT, "poma_long.csv"))

P("=" * 78)
P("merge baseline covariates onto POMA knee-level table")
P("=" * 78)
P(f"knee-level input : {sl.shape}")

# ---------------------------------------------------------------- subjectchar
sc = pd.read_excel(B + "subjectchar00.xlsx",
                   usecols=["ID", "V00AGE"])
P(f"\nsubjectchar00 {sc.shape}  V00AGE missing = {sc['V00AGE'].isna().sum()}")
P(f"  V00AGE range {sc['V00AGE'].min()}-{sc['V00AGE'].max()}")

# ---------------------------------------------------------------- enrollees
en = pd.read_excel(B + "enrollees.xlsx",
                   usecols=["ID", "P02SEX", "V00COHORT"])
P(f"enrollees {en.shape}  P02SEX: {en['P02SEX'].value_counts().to_dict()}")
P(f"  V00COHORT: {en['V00COHORT'].value_counts().to_dict()}")

# ---------------------------------------------------------------- allclinical00
ac = pd.read_excel(B + "allclinical00.xlsx",
                   usecols=["ID", "P01BMI", "P01HEIGHT", "P01WEIGHT",
                            "V00WOMKPL", "V00WOMKPR",
                            "V00WOMTSL", "V00WOMTSR", "V00PASE"])
P(f"allclinical00 {ac.shape}")

# ---------------------------------------------------------------- KL
kl = pd.read_excel(B + "kxr_sq_bu00.xlsx", usecols=["ID", "SIDE", "V00XRKL"])
P(f"\nkxr_sq_bu00 {kl.shape}")
P(f"  SIDE codes: {sorted(kl['SIDE'].unique())}")
P(f"  V00XRKL distribution:\n{kl['V00XRKL'].value_counts().sort_index().to_string()}")
dup = kl.duplicated(["ID", "SIDE"]).sum()
P(f"  duplicated ID x SIDE rows: {dup}")
if dup:
    kl = kl.sort_values(["ID", "SIDE", "V00XRKL"]).drop_duplicates(
        ["ID", "SIDE"], keep="first")
    P("  -> kept first after sorting")

# ---------------------------------------------------------------- merge
d = sl.copy()
d = d.merge(sc, left_on="id", right_on="ID", how="left").drop(columns=["ID"])
d = d.merge(en, left_on="id", right_on="ID", how="left").drop(columns=["ID"])
d = d.merge(ac, left_on="id", right_on="ID", how="left").drop(columns=["ID"])
d = d.merge(kl, left_on=["id", "side"], right_on=["ID", "SIDE"],
            how="left").drop(columns=["ID", "SIDE"])
P(f"\nafter merge: {d.shape}")

# ---------------------------------------------------------------- derive
# side 1 = RIGHT, 2 = LEFT (OAI convention) -> side-specific WOMAC
d["WOMAC_pain"] = np.where(d["side"] == 1, d["V00WOMKPR"], d["V00WOMKPL"])
d["WOMAC_total"] = np.where(d["side"] == 1, d["V00WOMTSR"], d["V00WOMTSL"])
d["female"] = (d["P02SEX"] == 2).astype(int)
d["kl_ge2"] = (d["V00XRKL"] >= 2).astype("float")
d.loc[d["V00XRKL"].isna(), "kl_ge2"] = np.nan
d["cohort"] = d["V00COHORT"].map({1: "progression", 2: "incidence",
                                  3: "control"})
# event "time" for the design (months to index visit); index_months by tmpt code
d["t_index_months"] = d["tmpt"].map({1: 12, 3: 24, 5: 36, 6: 48, 7: 60})

P("\n-- coverage of covariates (knee level, n=%d) --" % len(d))
for v in ["V00AGE", "P02SEX", "P01BMI", "V00XRKL", "WOMAC_pain", "V00PASE"]:
    P(f"  {v:<12s} non-missing {d[v].notna().sum():>4d} "
      f"({d[v].notna().mean()*100:5.1f}%)")

P("\n-- baseline characteristics by case status (knee level) --")
for v in ["V00AGE", "P01BMI", "V00XRKL", "WOMAC_pain", "V00PASE"]:
    g = d.groupby("case")[v].agg(["count", "mean", "std", "median"])
    P(f"\n  {v}:"); P("   " + g.to_string().replace("\n", "\n   "))
P("\n  female (%):")
P("   " + (d.groupby("case")["female"].mean() * 100).round(1).to_string())
P("\n  cohort distribution:")
P("   " + pd.crosstab(d["case"], d["cohort"]).to_string())

# ---------------------------------------------------------------- key design check
P("\n" + "=" * 78)
P("DESIGN CHECK: is tmpt identical within a matched stratum?")
P("=" * 78)
gg = d.groupby("newstrata")["tmpt"].nunique()
P(f"strata where BOTH members share the same tmpt : {int((gg == 1).sum())} / {len(gg)}")
P(f"strata where the two members differ           : {int((gg > 1).sum())}")
if (gg > 1).any():
    ex = gg[gg > 1].index[:5]
    P("example strata with differing tmpt:")
    P("   " + d[d["newstrata"].isin(ex)][["newstrata", "cls", "tmpt", "unit"]]
      .sort_values(["newstrata", "cls"]).to_string(index=False).replace("\n", "\n   "))

P("\nwithin-stratum difference in t_index_months (case - control):")
piv = d.pivot_table(index="newstrata", columns="case",
                    values="t_index_months")
dd = (piv[1] - piv[0]).dropna()
P(f"  mean {dd.mean():+.2f} months, SD {dd.std():.2f}, "
  f"identical in {(dd == 0).mean()*100:.1f}% of strata")

# ---------------------------------------------------------------- write
keep = ["unit", "id", "side", "cls", "case", "newstrata", "tmpt",
        "t_index_months", "n_visits", "span_yr",
        "V00AGE", "P02SEX", "female", "P01BMI", "P01HEIGHT", "P01WEIGHT",
        "V00XRKL", "kl_ge2", "WOMAC_pain", "WOMAC_total", "V00PASE",
        "V00COHORT", "cohort"] + \
       [c for c in d.columns if c.endswith("__slope") or c.endswith("__base")
        or c.endswith("__eblup")]
d[keep].to_csv(os.path.join(OUT, "poma_analysis_long.csv"), index=False)
P(f"\n-> poma_analysis_long.csv  {d[keep].shape}")

# visit-level file (for the joint model / longitudinal submodel)
vl = long[["unit", "id", "side", "visit_code", "months", "years", "case",
           "newstrata", "tmpt", "class", "MFTC_ThCtAB_aMe", "cMFTC_ThCtAB_aMe",
           "MT_ThCtAB_aMe", "cMF_ThCtAB_aMe", "MFTC_VC"]].copy()
vl = vl.merge(d[["unit", "V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain"]],
              on="unit", how="left")
vl.to_csv(os.path.join(OUT, "poma_visits_long.csv"), index=False)
P(f"-> poma_visits_long.csv    {vl.shape}")

with open(os.path.join(OUT, "covariate_merge_report.txt"), "w",
          encoding="utf-8") as f:
    f.write("\n".join(rep))
print("\n-> covariate_merge_report.txt")
