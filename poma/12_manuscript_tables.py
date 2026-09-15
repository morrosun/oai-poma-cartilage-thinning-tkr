#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
12_manuscript_tables.py -- build the REAL numeric tables (1-4) for the POMA
cMFTC manuscript, because until now the manuscript carried only legends with no
table bodies.

Table 1 (baseline characteristics) is RECOMPUTED here from
`poma_analysis_long.csv` on the primary analysis set (Set A), and is validated
against the golden values already printed in `clogit_cMFTC_seq_report.txt`.
Tables 2-4 are pure REFORMATTING of the already-validated result files
(nothing is re-estimated), so no number can drift.

Outputs (Analysis/poma_pilot/):
  Table1_baseline.csv / Table2_primary.csv / Table3_gradient.csv / Table4_lcmm.csv
  manuscript_tables.md      <- ready-to-paste markdown blocks

Usage:
  python "Scripts/poma/12_manuscript_tables.py"
"""
import io, os, re, sys
import numpy as np
import pandas as pd
from scipy import stats

ANA = "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot/"
TOL = 0.06          # absolute tolerance for golden-value check on means

# ---------------------------------------------------------------- helpers
SUP = str.maketrans("0123456789-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b")


def fmt_p(p):
    """P value in the house style: 7.2 × 10⁻⁶ / 0.069 / 0.0018"""
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "—"
    if p < 1e-4:
        e = int(np.floor(np.log10(p)))
        m = p / 10 ** e
        return "%.1f \u00d7 10%s" % (m, str(e).translate(SUP))
    if p < 0.001:
        return "%.4f" % p
    if p < 0.01:
        return "%.4f" % p
    return "%.3f" % p


def ci(x):
    """mean and 95% CI of a sample."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = len(x)
    m = x.mean()
    se = x.std(ddof=1) / np.sqrt(n)
    t = stats.t.ppf(0.975, n - 1)
    return m, m - t * se, m + t * se


def paired(b, a):
    """paired difference b-a, dropping pairs with a missing member (as R does):
    mean, 95% CI, paired-t P."""
    bb = np.asarray(b, float)
    aa = np.asarray(a, float)
    ok = np.isfinite(bb) & np.isfinite(aa)
    ff = ok.sum()
    xd = bb[ok] - aa[ok]
    m, lo, hi = ci(xd)
    t = stats.t.ppf(0.975, ff - 1)
    return m, m - t * xd.std(ddof=1) / np.sqrt(ff), m + t * xd.std(ddof=1) / np.sqrt(ff), \
        stats.ttest_rel(bb[ok], aa[ok]).pvalue


def ms(x, nd=1):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return "%.*f (%.*f)" % (nd, x.mean(), nd, x.std(ddof=1))


def orci(o, lo, hi, nd=2):
    return "%.*f (%.*f\u2013%.*f)" % (nd, o, nd, lo, nd, hi)


issues = []


def check(tag, got, want):
    ok = abs(got - want) <= TOL
    print("   %-26s got %9.4f   expected %9.4f   %s"
          % (tag, got, want, "OK" if ok else "*** MISMATCH ***"))
    if not ok:
        issues.append(tag)
    return ok


# ================================================================ Set A
print("=" * 78)
print("rebuilding Set A (primary analysis set) from poma_analysis_long.csv")
print("=" * 78)
d = pd.read_csv(ANA + "poma_analysis_long.csv", low_memory=False)
print("rows in file: %d" % len(d))

# NOTE: Set A is defined by `need <- c("slope","base","V00AGE","female","P01BMI",
# "V00XRKL","WOMAC_pain")` in 10c_clogit_seq_primary.R -- PASE is deliberately
# NOT part of the completeness rule (one control knee has a missing PASE), which
# is why the analysis set is 191 pairs and not 190. PASE is handled with na.rm
# below, exactly as the R script does.
REQ = ["V00AGE", "female", "P01BMI", "V00XRKL", "WOMAC_pain"]
d["thin"] = -d["cMFTC_ThCtAB_aMe__eblup"]          # positive = faster thinning
# Unshrunk per-knee OLS slope -- used ONLY for the differential-shrinkage
# sensitivity (3.6 (v)): the EBLUP shrinks each knee's slope toward the
# population mean with a weight that depends on the number of visits, so if
# cases happen to be followed more densely a reviewer will ask whether the
# contrast is an artefact of that shrinkage.
d["thin_ols"] = -d["cMFTC_ThCtAB_aMe__slope"]
keep = (d["n_visits"] >= 2) & d["thin"].notna() & d["cMFTC_ThCtAB_aMe__base"].notna()
for c in REQ:
    keep &= d[c].notna()
s = d[keep].copy()
cnt = s.groupby("newstrata").size()
s = s[s["newstrata"].isin(cnt[cnt == 2].index)].copy()
n_pairs, n_knees = s["newstrata"].nunique(), len(s)
print("Set A  ->  %d pairs / %d knees   (expect 191 / 382)" % (n_pairs, n_knees))
if (n_pairs, n_knees) != (191, 382):
    issues.append("Set A size")

a = s[s["case"] == 0].set_index("newstrata")
b = s[s["case"] == 1].set_index("newstrata")
idx = a.index.intersection(b.index)
print("matched pairs with both members: %d" % len(idx))

# ---- golden-value gate (all from clogit_cMFTC_seq_report.txt)
print("\n--- golden-value validation (vs clogit_cMFTC_seq_report.txt) ---")
check("age control", a.loc[idx, "V00AGE"].mean(), 63.9)
check("age case", b.loc[idx, "V00AGE"].mean(), 64.3)
check("female% both", 100 * a.loc[idx, "female"].mean(), 55.5)
check("BMI control", a.loc[idx, "P01BMI"].mean(), 29.5)
check("BMI case", b.loc[idx, "P01BMI"].mean(), 29.3)
check("KL control", a.loc[idx, "V00XRKL"].mean(), 2.80)
check("KL case", b.loc[idx, "V00XRKL"].mean(), 2.90)
check("WOMAC control", a.loc[idx, "WOMAC_pain"].mean(), 4.6)
check("WOMAC case", b.loc[idx, "WOMAC_pain"].mean(), 4.7)
check("thin rate control", a.loc[idx, "thin"].mean(), 0.1074)
check("thin rate case", b.loc[idx, "thin"].mean(), 0.1817)
check("SD thin rate (pooled)", s["thin"].std(ddof=1), 0.1411)
dm, dlo, dhi, dp = paired(b.loc[idx, "thin"], a.loc[idx, "thin"])
check("paired d_thin", dm, 0.074)
check("paired d_thin lo", dlo, 0.047)
check("paired d_thin hi", dhi, 0.101)
check("paired d_thin log10P", np.log10(dp), np.log10(1.43e-07))
# unshrunk OLS slope -- the differential-shrinkage rebuttal (sensitivity (v))
n_ols = int(s["thin_ols"].notna().sum())
print("   Set A knees with an OLS slope: %d / %d" % (n_ols, len(s)))
if n_ols != len(s):
    issues.append("OLS slope missing for some Set A knees")
dmo, dlo, dhi, dpo = paired(b.loc[idx, "thin_ols"], a.loc[idx, "thin_ols"])
check("OLS paired d_thin", dmo, 0.119)
check("OLS paired d_thin lo", dlo, 0.072)
check("OLS paired d_thin hi", dhi, 0.166)
print("   OLS thin rate control %.4f (%.4f)  case %.4f (%.4f)  poolSD %.4f"
      % (a.loc[idx, "thin_ols"].mean(), a.loc[idx, "thin_ols"].std(ddof=1),
         b.loc[idx, "thin_ols"].mean(), b.loc[idx, "thin_ols"].std(ddof=1),
         s["thin_ols"].std(ddof=1)))
r_vis_ols = np.corrcoef(s["n_visits"], s["thin_ols"])[0, 1]
r_vis_ebl = np.corrcoef(s["n_visits"], s["thin"])[0, 1]
print("   corr(n_visits, OLS thin) = %+.3f ; corr(n_visits, EBLUP thin) = %+.3f"
      % (r_vis_ols, r_vis_ebl))
check("visits/knee control", a.loc[idx, "n_visits"].mean(), 3.04)
check("visits/knee case", b.loc[idx, "n_visits"].mean(), 3.50)
# baseline thickness difference (the level signal that the rate later absorbs)
dm, dlo, dhi, dp = paired(b.loc[idx, "cMFTC_ThCtAB_aMe__base"],
                          a.loc[idx, "cMFTC_ThCtAB_aMe__base"])
check("paired d_base thickness", dm, -0.268)
check("paired d_base thickness P", dp, 0.030)
dm, _, _, dp = paired(b.loc[idx, "V00XRKL"], a.loc[idx, "V00XRKL"])
check("paired d_KL", dm, 0.094)
check("paired d_KL P", dp, 0.00184)
# PASE has one missing control value in Set A -> the R script's paired test uses
# the finite differences only (n = 190). Reproduce that exactly.
pm = np.asarray(b.loc[idx, "V00PASE"], float) - np.asarray(a.loc[idx, "V00PASE"], float)
pm = pm[np.isfinite(pm)]
print("   PASE paired differences available: %d / %d" % (len(pm), len(idx)))
check("paired d_PASE", pm.mean(), 4.342)
check("paired d_PASE lo", pm.mean() - stats.t.ppf(0.975, len(pm) - 1) * pm.std(ddof=1) / np.sqrt(len(pm)), -10.656)
check("paired d_PASE hi", pm.mean() + stats.t.ppf(0.975, len(pm) - 1) * pm.std(ddof=1) / np.sqrt(len(pm)), 19.340)
check("paired d_PASE P", stats.ttest_1samp(pm, 0).pvalue, 0.569)

if issues:
    print("\n!! validation failures: %s -- ABORTING, do not paste tables" % issues)
    sys.exit(2)
print("\nall golden values reproduced -> Set A reconstruction is trustworthy\n")

# ================================================================ Table 1
ROWS = [                       # label, column, decimals, footnote
    ("年龄（岁）", "V00AGE", 1, ""),
    ("女性，n (%)", "female", None, ""),
    ("体质指数（kg/m\u00b2）", "P01BMI", 1, ""),
    ("Kellgren\u2013Lawrence 分级", "V00XRKL", 2, ""),
    ("WOMAC 疼痛（0\u201320）", "WOMAC_pain", 1, ""),
    ("PASE 体力活动评分", "V00PASE", 1, "\u2020"),
    ("cMFTC 基线软骨厚度（mm）", "cMFTC_ThCtAB_aMe__base", 2, ""),
    ("cMFTC 变薄速率（mm/yr）", "thin", 3, ""),
    ("　其中：未收缩逐膝 OLS 速率（mm/yr）", "thin_ols", 3, "\u2021"),
    ("术前访视次数/膝", "n_visits", 2, ""),
]
t1 = []
for lab, col, nd, fn in ROWS:
    ac, bc = a.loc[idx, col], b.loc[idx, col]
    if nd is None:                                   # binary / percentage row
        n0 = int(np.nansum(ac)); n1 = int(np.nansum(bc))
        c0 = "%d (%.1f)" % (n0, 100 * n0 / len(idx))
        c1 = "%d (%.1f)" % (n1, 100 * n1 / len(idx))
        t1.append([lab + fn, c0, c1, "0（配对内恒定）", "\u2014"])
        continue
    c0, c1 = ms(ac, nd), ms(bc, nd)
    m, lo, hi, p = paired(bc, ac)
    t1.append([lab + fn, c0, c1, "%+.3f (%.3f\u2013%+.3f)" % (m, lo, hi), fmt_p(p)])
# KL grade distribution (categorical, n (%))
kl0 = a.loc[idx, "V00XRKL"].value_counts().sort_index()
kl1 = b.loc[idx, "V00XRKL"].value_counts().sort_index()
for g in sorted(set(kl0.index) | set(kl1.index)):
    n0, n1 = int(kl0.get(g, 0)), int(kl1.get(g, 0))
    t1.append(["　KL 分级 %d，n (%%)" % int(g),
               "%d (%.1f)" % (n0, 100 * n0 / len(idx)),
               "%d (%.1f)" % (n1, 100 * n1 / len(idx)), "—", "—"])

df1 = pd.DataFrame(t1, columns=["变量", "对照（n = %d）" % len(idx),
                                "病例（n = %d）" % len(idx), "配对内差值（病例\u2212对照）", "P"])
df1.to_csv(ANA + "Table1_baseline.csv", index=False, encoding="utf-8-sig")

md1 = ["| " + " | ".join(df1.columns) + " |",
       "|" + "---|" * len(df1.columns)]
for _, r in df1.iterrows():
    md1.append("| " + " | ".join(str(v) for v in r.values) + " |")
T1 = "\n".join(md1) + (
    "\n\n> \u2020 PASE 在 Set A 中有 1 例对照膝缺失（配对 118），该行均数按可用值计算，"
    "配对差值与 P 值基于 190 个配对；其余各行均为 191 个配对。\n"
    ">\n"
    "> \u2021 未收缩的逐膝最小二乘（OLS）斜率仅用于评估**经验贝叶斯收缩**是否制造组间差异"
    "（见 3.6 节 (v)）；其分析集内标准差大于 EBLUP 速率，不作为 OR 的「每 1 SD」换算基准。")

# ================================================================ Table 2
seq = pd.read_csv(ANA + "clogit_cMFTC_seq_OR.csv")
jm = pd.read_csv(ANA + "jm_multimetric.csv")
jm = jm[jm["metric"] == "cMFTC_ThCtAB_aMe"].iloc[0]
LAB = {"M0 crude": "M0　粗模型", "M1 + age(/5y)": "M1　+ 年龄（每 5 岁）",
       "M2 + sex": "M2　+ 性别", "M3 + BMI(/5)": "M3　+ BMI（每 5）",
       "M4 + KL grade": "M4　+ KL 分级", "M5 + WOMAC pain(/5)": "M5　+ WOMAC 疼痛（每 5）",
       "M6 full adjusted": "M6　完整校正", "M7 full + PASE": "M7　M6 + PASE"}
t2 = []
for _, r in seq.iterrows():
    t2.append([LAB.get(r["model"], r["model"]), False,
               orci(r["OR"], r["lo"], r["hi"]), orci(r["OR_01"], r["OR_01_lo"], r["OR_01_hi"]),
               fmt_p(r["p"])])
t2.append(["联合模型　未校正（配对内条件似然）", False,
           orci(jm["JM_A_ORSD"], jm["JM_A_ORSD_lo"], jm["JM_A_ORSD_hi"]),
           orci(jm["JM_A_OR01"], jm["JM_A_OR01_lo"], jm["JM_A_OR01_hi"]), "—"])
t2.append(["联合模型　+ 配对内 KL/疼痛差值", False,
           "%.2f" % jm["JM_Aadj_ORSD"],
           orci(jm["JM_Aadj_OR01"], jm["JM_Aadj_OR01_lo"], jm["JM_Aadj_OR01_hi"]), "—"])
df2 = pd.DataFrame(t2, columns=["模型", "panel", "OR／每快 1 SD（95% CI）",
                                "OR／每快 0.1 mm/yr（95% CI）", "P"])
df2[["模型", "OR／每快 1 SD（95% CI）", "OR／每快 0.1 mm/yr（95% CI）", "P"]] \
    .to_csv(ANA + "Table2_primary.csv", index=False, encoding="utf-8-sig")

md2 = ["| 模型 | OR／每快 1 SD（95% CI） | OR／每快 0.1 mm/yr（95% CI） | P |",
       "|---|---|---|---|",
       "| **A. 配对内条件 logistic 递进校正（n = 191 配对）** | | | |"]
for _, r in df2[df2["panel"] == False].iterrows():
    if str(r["模型"]).startswith("联合模型"):
        continue
    md2.append("| %s | %s | %s | %s |" % (r["模型"], r["OR／每快 1 SD（95% CI）"],
                                          r["OR／每快 0.1 mm/yr（95% CI）"], r["P"]))
md2.append("| **B. 配对内条件似然联合模型（同一 191 配对）** | | | |")
for _, r in df2[df2["panel"] == False].iterrows():
    if not str(r["模型"]).startswith("联合模型"):
        continue
    md2.append("| %s | %s | %s | %s |" % (r["模型"], r["OR／每快 1 SD（95% CI）"],
                                          r["OR／每快 0.1 mm/yr（95% CI）"], r["P"]))
T2 = "\n".join(md2)

# ================================================================ Table 3
REG = {"cMFTC_ThCtAB_aMe": "cMFTC　中央内侧胫-股（**负重**，主指标）",
       "MFTC_ThCtAB_aMe": "MFTC　内侧胫-股（整体）",
       "MT_ThCtAB_aMe": "MT　内侧胫骨",
       "cMF_ThCtAB_aMe": "cMF　中央内侧股骨",
       "cLFTC_ThCtAB_aMe": "cLFTC　中央外侧胫-股（**负重**）",
       "cLF_ThCtAB_aMe": "cLF　中央外侧股骨"}
order = ["cMFTC_ThCtAB_aMe", "MFTC_ThCtAB_aMe", "MT_ThCtAB_aMe", "cMF_ThCtAB_aMe",
         "cLFTC_ThCtAB_aMe", "cLF_ThCtAB_aMe"]
grd = pd.read_csv(ANA + "gradient_contrast_test.csv")
log = io.open(ANA + "gradient_contrast.log", encoding="utf-8", errors="replace").read()


def diff_block(name):
    # anchor on the "JM on DIFF_xxx" banner, otherwise the first bare mention of
    # "DIFF_cMFTC" (the one-line "DIFF_cMFTC : 190 pairs / 380 knees" note) would
    # terminate the capture before any estimate is reached.
    m = re.search(r"JM on " + re.escape(name) + r"[^\n]*\n(.*?)(?:={10,}|\Z)", log, re.S)
    if not m:
        raise RuntimeError("cannot parse DIFF block for " + name)
    t = m.group(1)
    out = {}
    for k, pat in [("or01", r"OR / 0\.1 mm/yr faster\s*:\s*([-\d.]+) \(([-\d.]+), ([-\d.]+)\)"),
                   ("orsd", r"OR / 1 SD faster\s*:\s*([-\d.]+) \(([-\d.]+), ([-\d.]+)\)"),
                   ("clogit", r"\[clogit on the same difference\] OR/SD = ([\d.]+) \(([\d.]+)-([\d.]+)\)"),
                   ("pas", r"P\(aV<0\)=[\d.]+\s+P\(aS<0\)=([\d.]+)"),
                   ("pairs", r"pairs (\d+) \| knees (\d+)")]:
        mm = re.search(pat, t)
        if mm is None:
            raise RuntimeError("pattern %r not found in %s block" % (k, name))
        out[k] = [float(x) if "." in x else int(x) for x in mm.groups()]
    return out


Dc, Dm = diff_block("DIFF_cMFTC"), diff_block("DIFF_MFTC")
t3 = []
for m in order:
    r = jm_all = pd.read_csv(ANA + "jm_multimetric.csv")
    r = r[r["metric"] == m].iloc[0]
    t3.append([REG[m], "内侧" if r["region"] == "Medial" else "外侧",
               "%d" % r["n_pairs_A"],
               orci(r["JM_A_ORSD"], r["JM_A_ORSD_lo"], r["JM_A_ORSD_hi"]),
               orci(r["JM_A_OR01"], r["JM_A_OR01_lo"], r["JM_A_OR01_hi"])])
df3a = pd.DataFrame(t3, columns=["区室指标", "分组", "配对", "OR／每快 1 SD（95% CrI）",
                                 "OR／每快 0.1 mm/yr（95% CrI）"])

g = {r["contrast"]: r for _, r in grd.iterrows()}
comp = g["composite medial vs lateral"]
t3b = []
for key in ["cMFTC vs cLF", "cMFTC vs cLFTC", "MFTC vs cLF", "MFTC vs cLFTC"]:
    r = g[key]
    # the contrast file stores the CI only for the medial term and for the ratio;
    # the lateral term is carried as a point estimate (its SE enters the ratio
    # variance V11 + V22 - 2*V12, which is what the test uses).
    t3b.append(["(1) 双变量条件 logistic：%s" % key,
                "%.3f (%.3f\u2013%.3f)" % (r["OR_medial"], r["OR_medial_lo"], r["OR_medial_hi"]),
                "%.3f" % r["OR_lateral"],
                "**%.3f (%.3f\u2013%.3f)**" % (r["OR_ratio"], r["OR_ratio_lo"], r["OR_ratio_hi"]),
                fmt_p(r["p"])])
t3b.append(["(2) 复合评分：内侧族 vs 外侧族",
            "%.3f (%.3f\u2013%.3f)" % (comp["OR_medial"], comp["OR_medial_lo"], comp["OR_medial_hi"]),
            "%.3f" % comp["OR_lateral"],
            "**%.3f (%.3f\u2013%.3f)**" % (comp["OR_ratio"], comp["OR_ratio_lo"], comp["OR_ratio_hi"]),
            fmt_p(comp["p"])])
df3b = pd.DataFrame(t3b, columns=["正式对照", "内侧 OR／SD（95% CI）",
                                  "外侧 OR／SD", "内侧/外侧 OR 比值（95% CI）", "P"])
t3c = []
for nm, dd in [("cMFTC \u2212 cLF（主指标）", Dc), ("MFTC \u2212 cLF", Dm)]:
    t3c.append(["膝内 %s" % nm,
                orci(dd["orsd"][0], dd["orsd"][1], dd["orsd"][2]),
                orci(dd["or01"][0], dd["or01"][1], dd["or01"][2]),
                "%.3f" % dd["pas"][0],
                orci(dd["clogit"][0], dd["clogit"][1], dd["clogit"][2])])
df3c = pd.DataFrame(t3c, columns=["膝内内侧\u2212外侧差值（联合模型，%d 配对）" % Dc["pairs"][0],
                                  "OR／每快 1 SD（95% CrI）",
                                  "OR／每快 0.1 mm/yr（95% CrI）",
                                  "P(aS < 0)", "同批条件 logistic OR／SD（95% CI）"])
df3a.to_csv(ANA + "Table3_gradient.csv", index=False, encoding="utf-8-sig")


def mdtbl(df):
    out = ["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(v) for v in r.values) + " |")
    return "\n".join(out)


T3 = (mdtbl(df3a)
      + "\n\n**下栏　三种正式梯度检验**\n\n"
      + "**(1)(2) 逐指标与复合评分的对比**\n\n" + mdtbl(df3b)
      + "\n\n**(3) 膝内差值（消除个体层面的共同缩放）**\n\n" + mdtbl(df3c)
      + ("\n\n> 外侧一列的 95% 区间未在结果文件中单独保存（比值检验已使用完整协方差 "
         "V\u2081\u2081 + V\u2082\u2082 \u2212 2V\u2081\u2082，故内侧与外侧的相关性已被考虑）；"
         "第 (3) 组以 190 对可用配对拟合，P(aS < 0) 为后验方向概率。"))

# ================================================================ Table 4
prof = pd.read_csv(ANA + "lcmm_class_profile.csv")
cor = pd.read_csv(ANA + "lcmm_class_OR.csv")
cormap = {int(r["class"]): r for _, r in cor.iterrows()}
REF = 2
prof = prof.sort_values("slope_mm_yr", ascending=False).reset_index(drop=True)  # slow->fast
t4 = []
for _, r in prof.iterrows():
    cl = int(r["class"])
    if cl == REF:
        oc, oratio, pgt = "参考", "\u2014", "\u2014"
    else:
        cr = cormap[cl]
        oratio = orci(cr["OR"], cr["lo"], cr["hi"])
        pgt = "%.2f" % cr["P_OR_gt1"]
    t4.append(["class %d%s" % (cl, "（最慢，参考）" if cl == REF else
                               ("（**快速变薄型**）" if cl == 4 else "")),
               "%d (%.1f)" % (r["n"], r["pct"]),
               "%.3f" % r["slope_mm_yr"], "%.2f" % r["baseline_mm"],
               "**%.1f**" % r["pct_KR"], "%.2f" % r["mean_KL"], "%.2f" % r["visits"],
               oratio, pgt])
df4 = pd.DataFrame(t4, columns=["轨迹类", "n (%)", "平均变薄速率 (mm/yr)", "平均基线厚度 (mm)",
                                "置换率 (%)", "平均 KL", "平均访视次数", "OR（vs 第 2 类）",
                                "P(OR>1)"])
df4.to_csv(ANA + "Table4_lcmm.csv", index=False, encoding="utf-8-sig")
T4 = mdtbl(df4)

# ================================================================ write blocks
blocks = [("# 表 1　主分析集基线特征（191 个匹配配对 / 382 膝）", T1),
          ("# 表 2　主指标（cMFTC）变薄速率与全膝置换的配对内关联：递进校正与两法对照", T2),
          ("# 表 3　内侧与外侧区室变薄率关联的正式对照（191 个完整配对）", T3),
          ("# 表 4　术前 cMFTC 轨迹的潜类别解、类特征与置换结局关联", T4)]
with io.open(ANA + "manuscript_tables.md", "w", encoding="utf-8") as f:
    f.write("<!-- generated by Scripts/poma/12_manuscript_tables.py ; do not edit by hand -->\n\n")
    for h, b in blocks:
        f.write(h + "\n\n" + b + "\n\n")

print("=" * 78)
print("Table 1")
print(df1.to_string(index=False))
print("\nTable 2")
print(df2[["模型", "OR／每快 1 SD（95% CI）", "OR／每快 0.1 mm/yr（95% CI）", "P"]].to_string(index=False))
print("\nTable 3 (upper)")
print(df3a.to_string(index=False))
print("\nTable 3 (lower, tests 1-2)")
print(df3b.to_string(index=False))
print("\nTable 3 (lower, test 3)")
print(df3c.to_string(index=False))
print("\nTable 4")
print(df4.to_string(index=False))
print("\nwrote: Table1_baseline.csv Table2_primary.csv Table3_gradient.csv "
      "Table4_lcmm.csv manuscript_tables.md")
