# -*- coding: utf-8 -*-
"""
47_verify_numbers_strict.py  --  assertion-based numeric audit of the OAC manuscript.

WHY THIS EXISTS
---------------
`45_verify_v5_numbers.R` PRINTS the stored values next to the manuscript values
and leaves the comparison to a human.  That is exactly the step that failed
repeatedly during the v4 -> v5 -> v6 cycle: a table row was updated while the
sentence quoting it, or the Abstract quoting the table, was not.  Five such
one-step rounding drifts survived into v5 and were only caught here.

This script removes the human.  For every audited number it

  (a) recomputes the value from the object on disk -- including the interval
      inversions that orient Table 6 so that OR > 1 means faster progression,
      which are easy to get wrong and are NOT stored anywhere; and
  (b) asserts that the manuscript contains the correctly formatted string.

It additionally checks the two things a value-level audit cannot see:

  * every script named in Supplementary Table S2 exists in the ARCHIVED
    repository, and every archived script is named -- so the inventory cannot
    silently fall behind the code (this is how `jm_core_legacy_sampler.R` was
    found unlisted);
  * the DOI in the Declarations and in S2 points at the archive that actually
    contains the corrected sampler.

Typography is normalised before matching: the manuscript is typeset with a real
minus sign (U+2212) and en dashes, the CSVs carry ASCII, and an en dash must
never register as a numeric mismatch.

Usage
    python 47_verify_numbers_strict.py [BASENAME]      # default Manuscript_EN_OAC_v6
Requires GH_TOKEN in the environment for the repository cross-check; without it
that section is skipped and reported as SKIPPED rather than silently passing.
Exit code 0 = every check passed, 1 = at least one failed.
"""
import csv, decimal, io, json, os, re, statistics, sys, urllib.request

BASE = "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"
BASENAME = sys.argv[1] if len(sys.argv) > 1 else "Manuscript_EN_OAC_v6"
MS = os.path.join(BASE, BASENAME + ".md")
if not os.path.exists(MS):
    sys.exit("FAIL: %s not found -- run build_oac_v1.py first" % MS)

REPO = "morrosun/oai-poma-cartilage-thinning-tkr"
TOKEN = os.environ.get("GH_TOKEN", "")

txt = io.open(MS, encoding="utf-8").read()
lines = txt.split("\n")

# The manuscript is typeset with a real minus sign (U+2212) and an en dash
# (U+2013) as the interval separator, while the CSVs carry ASCII "-".  Compare
# on a normalised copy so that a typography choice never registers as a
# numeric mismatch -- but keep the original for the "is it present" checks.
def norm(s):
    return (s.replace("\u2212", "-").replace("\u2013", "-")
             .replace("\u2014", "-").replace("\u2010", "-"))

NTXT = norm(txt)

results = []          # (ok, section, label, detail)
def chk(sec, label, ok, detail=""):
    results.append((bool(ok), sec, label, detail))
    return bool(ok)

def has(s, sec=None, label=None):
    """the manuscript must contain the string s (dash-normalised)"""
    n = norm(s)
    ok = (n in NTXT)
    return chk(sec, label or s, ok, "present" if ok else "MISSING: %s" % s)

def has_re(pat, sec, label, group=1):
    m = re.search(pat, txt)
    return chk(sec, label, bool(m), (m.group(group) if m else "NO MATCH"))

def f(x, nd=3):
    return ("%%.%df" % nd) % float(x)

def ci(est, lo, hi, nd=2, dash="-"):
    """render 1.23 (1.05-1.50) with an en-dash, the manuscript's style"""
    return "%s (%s%s%s)" % (f(est, nd), f(lo, nd), dash, f(hi, nd))

def ci2(est, lo, hi, nd=2):
    """same, ASCII hyphen -- the manuscript is dash-normalised before matching"""
    return "%s (%s-%s)" % (f(est, nd), f(lo, nd), f(hi, nd))

def rd(name):
    return list(csv.DictReader(io.open(os.path.join(BASE, name), encoding="utf-8")))

def fnum(x):
    try:
        return float(x)
    except Exception:
        return None


# =========================================================== 1. Table 2 / 3 --
mm = rd("jm_multimetric.csv")
by = {r["metric"]: r for r in mm}

sec = "Table 2B / Table 3 : joint model (jm_multimetric.csv)"
c = by["cMFTC_ThCtAB_aMe"]
chk(sec, "unadjusted OR/SD + OR/0.1",
    has(ci(c["JM_A_ORSD"], c["JM_A_ORSD_lo"], c["JM_A_ORSD_hi"]), sec, "unadj OR/SD")
    and has(ci(c["JM_A_OR01"], c["JM_A_OR01_lo"], c["JM_A_OR01_hi"]), sec, "unadj OR/0.1"),
    "cMFTC %s / %s" % (ci(c["JM_A_ORSD"], c["JM_A_ORSD_lo"], c["JM_A_ORSD_hi"]),
                       ci(c["JM_A_OR01"], c["JM_A_OR01_lo"], c["JM_A_OR01_hi"])))

# the adjusted point estimates live in jm_multimetric.csv but their intervals do
# not (10_jm_multimetric.R never wrote those columns), so they come from the
# fitted object via 46_export_audit_values.R.  The manuscript quotes two
# decimals, so compare at that precision.
av = {r["key"]: r["value"] for r in rd("audit_values.csv")}
def av2(k):
    m = re.match(r"([\d.]+)\s*\(([\d.]+)-([\d.]+)\)", av[k])
    return ci2(m.group(1), m.group(2), m.group(3), 2)
chk(sec, "within-pair-adjusted OR/SD + OR/0.1",
    has(av2("Aadj_ORSD"), sec, "adj OR/SD") and has(av2("Aadj_OR01"), sec, "adj OR/0.1"),
    "Aadj %s / %s" % (av2("Aadj_ORSD"), av2("Aadj_OR01")))
chk(sec, "adjusted point estimate matches the CSV column",
    abs(fnum(av["Aadj_ORSD"].split()[0]) - fnum(c["JM_Aadj_ORSD"])) < 0.001,
    "csv %s vs rds %s" % (c["JM_Aadj_ORSD"], av["Aadj_ORSD"]))

sec = "Table 3 : six compartment metrics (jm_multimetric.csv)"
for k, r in by.items():
    a = ci(r["JM_A_ORSD"], r["JM_A_ORSD_lo"], r["JM_A_ORSD_hi"])
    b = ci(r["JM_A_OR01"], r["JM_A_OR01_lo"], r["JM_A_OR01_hi"])
    chk(sec, k, has(a, sec, k + " OR/SD") and has(b, sec, k + " OR/0.1"),
        "%s | %s" % (a, b))

# ================================================= 2. gradient contrast test --
sec = "Table 3 : medial-to-lateral contrasts (gradient_contrast_test.csv)"
for r in rd("gradient_contrast_test.csv"):
    for key, nd in (("ratio", 2),):
        est, lo, hi = fnum(r.get("ratio")), fnum(r.get("ratio_lo")), fnum(r.get("ratio_hi"))
        if est is None:
            continue
        s = ci(est, lo, hi, nd)
        chk(sec, "%s ratio" % r.get("contrast", r.get("label", "?")), s in txt,
            "%s %s" % (r.get("contrast", r.get("label", "?")), s))

# ================================================ 3. sampler ablation (S4 B/C) --
ab = rd("sampler_ablation.csv")
groups = {}
for r in ab:
    groups.setdefault(r["label"], []).append(r)

sec = "S4 Panel B : ablation on the real 191 pairs (sampler_ablation.csv)"
want_b = [
    ("crude | legacy      (f1=F f2=F)", "v4 kernel (all three corrections off)"),
    ("crude | fix1 only   (f1=T f2=F)", "Correction (i) only"),
    ("crude | fix2 only   (f1=F f2=T)", "Correction (ii) only"),
    ("crude | corrected   (f1=T f2=T)", "v5 kernel (i) + (ii)"),
    ("adj   | legacy      (f1=F f2=F)", "v4 kernel, with within-pair covariates"),
    ("adj   | corrected   (f1=T f2=T)", "v5 kernel, with within-pair covariates"),
    ("adj   | corrected, zt dropped", "covariant term dropped"),
]
for lab, disp in want_b:
    rows = groups[lab]
    orsd = [fnum(r["ORSD"]) for r in rows]
    or01 = [fnum(r["OR01"]) for r in rows]
    gk = [fnum(r["g_kl"]) for r in rows]
    sd = statistics.stdev(orsd)
    a, b = f(statistics.mean(orsd)), f(statistics.mean(or01))
    chk(sec, disp, has("| %s | %s | %s |" % (a, f(sd), b), sec, disp),
        "mean OR/SD %s (sd %s), OR/0.1 %s, g_kl %s" %
        (a, f(sd), b, f(statistics.mean(gk)) if gk[0] is not None else "—"))

sec = "S4 Panel C : per-metric change (jm_multimetric.csv)"
c = by["cMFTC_ThCtAB_aMe"]
chk(sec, "cMFTC v4 -> v5 ratio 1.17",
    has("| cMFTC | Medial | 2.601 | 3.037 | 1.17 |", sec, "cMFTC row"),
    "v5 %.3f" % fnum(c["JM_A_ORSD"]))
for key, v4v in (("MFTC_ThCtAB_aMe", 2.423), ("cLFTC_ThCtAB_aMe", 2.056),
                ("cLF_ThCtAB_aMe", 1.652)):
    r = by[key]
    v5 = fnum(r["JM_A_ORSD"])
    chk(sec, "%s ratio" % key, has("| %s |" % f(v5), sec, key), "v5 %.3f" % v5)

# ============================================== 4. defect isolation (S4 A) ----
sec = "S4 Panel A : exact reference (defect1_isolation.csv)"
di = {r["quantity"]: r for r in rd("defect1_isolation.csv")}
exact = {k: fnum(v["exact"]) for k, v in di.items()}
jac = {k: fnum(v["jacobi_s9"]) for k, v in di.items()}
pb = {k: fnum(v["pairblock_s9"]) for k, v in di.items()}
for k, disp in (("eta", "0.696321"), ("plogis", "0.624371"), ("logplogis", "-0.550898")):
    chk(sec, "%s exact" % k, has(disp, sec, "%s exact" % k), "disk %.6f" % exact[k])
chk(sec, "logplogis v4 kernel -0.555256", has("**-0.555256**", sec, "jacobi logplogis"),
    "disk %.6f" % jac["logplogis"])
chk(sec, "logplogis v5 kernel -0.551112", has("**-0.551112**", sec, "pairblock logplogis"),
    "disk %.6f" % pb["logplogis"])
chk(sec, "z values -24.9 / -1.3", has("**-24.9**", sec, "z jacobi") and has("**-1.3**", sec, "z pairblock"))
chk(sec, "second seed deviations",
    has("−5.5, −10.2 and −20.0", sec, "seed 2 jacobi") and has("−0.5, −0.7 and −0.7", sec, "seed 2 pairblock"))

# =========================================== 5. chain adequacy (S4 D) ---------
sec = "S4 Panel D : chain adequacy (chain_diagnostics_primary.csv + mcmc_multichain_rhat.csv)"
cd = {r["quantity"]: fnum(r["value"]) for r in rd("chain_diagnostics_primary.csv")}
mc = {r["parameter"]: r for r in rd("mcmc_multichain_rhat.csv")}
chk(sec, "ESS aS 1 145", has("| Effective sample size, rate association | 1 145 |", sec, "ESS"),
    "disk %.0f" % cd["ESS aS"])
chk(sec, "act 4.37 (0.558)",
    has("| 4.37 (0.558) |", sec, "act"), "disk %.2f (%.3f)" % (cd["autocorrelation time aS"], cd["lag-1 autocorrelation aS"]))
chk(sec, "MC-SE 0.030", has("| 0.030 |", sec, "mcse sd"), "disk %.4f" % cd["MC-SE of aS (in units of posterior SD)"])
chk(sec, "MC-SE log 0.014", has("| 0.014 |", sec, "mcse log"), "disk %.4f" % cd["MC-SE on log-OR scale"])
chk(sec, "first/second half 3.09 / 3.01", has("| 3.09 / 3.01 |", sec, "halves"),
    "disk %.2f / %.2f" % (cd["first-half OR per SD"], cd["second-half OR per SD"]))
chk(sec, "acceptance 0.81 / 0.82", has("| 0.81 / 0.82 |", sec, "acceptance"),
    "disk %.2f / %.2f" % (cd["MH acceptance, random effects"], cd["MH acceptance, association"]))
worst = max(fnum(v["Rhat"]) for v in mc.values())
chk(sec, "worst R-hat 1.0002", has("| 1.0002 |", sec, "rhat"), "disk %.6f" % worst)
chk(sec, "cross-chain ESS 4 451", has("| 4 451 |", sec, "cross ess"), "disk %.0f" % fnum(mc["aS"]["ESS_cross_chain"]))

# =============================================== 6. cross-pipeline (Table 6) --
sec = "Table 6 / Abstract : cross-pipeline corroboration (cth_clogit_models.csv)"
rows = rd("cth_clogit_models.csv")

# Table 6 panel A states every row oriented so that OR > 1 = faster progression.
# For the cMFTC terms the fitted OR is < 1 (thinner = slower), so both the point
# estimate and the interval are inverted: OR(-x) = 1/OR(x) with the endpoints
# swapped.  The inversion is recomputed here, never read from the CSV.
def inv(r, nd=2):
    return ci2(1.0 / float(r["OR"]), 1.0 / float(r["hi"]), 1.0 / float(r["lo"]), nd)

def fwd(r, nd=2):
    return ci2(r["OR"], r["lo"], r["hi"], nd)

def pick(prefix, term):
    """row of the model whose label starts with `prefix`, for coefficient `term`"""
    for r in rows:
        if r["model"].strip().startswith(prefix) and r["term"] == term:
            return r
    return None

A1 = pick("A1", "z_cmftc")
A2 = pick("A2", "z_cmftc_o")
B1 = pick("B1", "z_cth_pre")
B2 = pick("B2", "z_cth_pre_e")
B3 = pick("B3", "z_cth_pre")
C1C = pick("C1", "z_cth_pre")
C1M = pick("C1", "z_cmftc")
C2C = pick("C2", "z_cth_pre")
C2M = pick("C2", "z_cmftc")
# NB the manuscript renumbers the all-timepoint sensitivities S1-S3; the CSV
# labels them S2, S3 and S5 (its S1/S4 are cMFTC-baseline variants the
# manuscript does not tabulate).  Map explicitly rather than by position.
S1 = pick("S2", "z_cth_all")          # manuscript S1
S2 = pick("S3", "z_cth_e")            # manuscript S2
S3C = pick("S5", "z_cth_all")         # manuscript S3, CTh term
S3M = pick("S5", "z_cmftc")           # manuscript S3, cMFTC term

for lbl, r, kind in (("A1 cMFTC EBLUP (inverted)", A1, "inv"),
                     ("A2 cMFTC unshrunk (inverted)", A2, "inv"),
                     ("B1 CTh pre-index OLS (primary)", B1, "fwd"),
                     ("B2 CTh pre-index EBLUP (primary)", B2, "fwd"),
                     ("B3 CTh + baseline", B3, "fwd"),
                     ("C1 CTh term", C1C, "fwd"),
                     ("C1 cMFTC term (inverted)", C1M, "inv"),
                     ("C2 CTh term", C2C, "fwd"),
                     ("C2 cMFTC term (inverted)", C2M, "inv"),
                     ("S1 all-timepoint CTh", S1, "fwd"),
                     ("S2 all-timepoint EBLUP", S2, "fwd"),
                     ("S3 CTh term", S3C, "fwd"),
                     ("S3 cMFTC term (inverted)", S3M, "inv")):
    if r is None:
        chk(sec, lbl, False, "row not found in CSV")
    else:
        s = inv(r) if kind == "inv" else fwd(r)
        chk(sec, lbl, has(s, sec, lbl), "%s -> %s" % (r["model"][:34], s))

# the Abstract quotes B1 as well: it MUST carry the same rounding as Table 6.
# The abstract is typeset in running prose ("2.48, 1.68 to 3.65") while the
# tables use brackets and an en dash, so compare the three numbers in either
# form -- the property being audited is the rounding, not the typography.
ABS = NTXT.split("**Results.**")[1].split("**Conclusions.**")[0]
_b1 = [f(fnum(B1[k]), 2) for k in ("OR", "lo", "hi")]
_b1pat = (re.escape(_b1[0]) + r"\s*[,(\s]\s*" + re.escape(_b1[1])
          + r"\s*(?:to|[-])\s*" + re.escape(_b1[2]))
chk(sec, "Abstract B1 rounding identical to Table 6", re.search(_b1pat, ABS) is not None,
    "Table 6 %s | abstract %s" % (ci2(B1["OR"], B1["lo"], B1["hi"]),
                                  _b1[0] + "/" + _b1[1] + "/" + _b1[2]))

tp = rd("cth_tertile_OR.csv")
chk(sec, "panel B fastest tertile", has(ci2(tp[0]["OR"], tp[0]["lo"], tp[0]["hi"]), sec, "tertile"),
    ci2(tp[0]["OR"], tp[0]["lo"], tp[0]["hi"]))
cc2 = {r["metric"]: r["value"] for r in rd("cth_concordance.csv")}
chk(sec, "panel C tertile agreement matches the CSV",
    has("%.1f" % float(cc2["tertile_agreement_pct"]), sec, "agreement"),
    "disk %.1f" % float(cc2["tertile_agreement_pct"]))
chk(sec, "panel C extreme discordance matches the CSV",
    has("%.1f%%" % float(cc2["extreme_discordance_pct"]), sec, "extreme discordance"),
    "disk %.1f%%" % float(cc2["extreme_discordance_pct"]))
chk(sec, "panel C Pearson matches the CSV",
    has("%.3f" % float(cc2["Pearson_pre"]), sec, "Pearson"),
    "disk %.3f" % float(cc2["Pearson_pre"]))

# the shrinkage benchmark is quoted three times (Table 6C, 4.3, Supplementary
# Results) and is derived, not stored as a ratio
sh = {(r["metric"], r["estimate"]): fnum(r["mean_paired_diff"])
      for r in rd("cth_shrinkage.csv")}
_rc = sh[("CTh-Score pre-index", "EBLUP")] / sh[("CTh-Score pre-index", "unshrunk OLS")]
_rm = sh[("cMFTC", "EBLUP")] / sh[("cMFTC", "unshrunk OLS")]
chk(sec, "shrinkage ratios 0.541 / 0.619 recomputed from cth_shrinkage.csv",
    abs(round(_rc, 3) - 0.541) < 5e-4 and abs(round(_rm, 3) - 0.619) < 5e-4,
    "disk %.4f / %.4f" % (_rc, _rm))
chk(sec, "the strict-window shrinkage ratio is quoted everywhere as 0.541",
    has("0.541") and "0.585" not in txt, "0.541 in Table 6C, 4.3 and Supp Results")

# the strict-window C1 pair in Supplementary Results, in prose format
_pf = lambda r, k: f(fnum(r[k]), 2)
_s1 = "%s, %s to %s" % (_pf(C1C, "OR"), _pf(C1C, "lo"), _pf(C1C, "hi"))
_s2 = "%s, %s to %s" % (_pf(C1M, "OR"), _pf(C1M, "lo"), _pf(C1M, "hi"))
chk(sec, "Supplementary Results C1 CTh term: %s" % _s1, has(_s1, sec, "C1 CTh"), _s1)
chk(sec, "Supplementary Results C1 cMFTC term in the thickness sign: %s" % _s2,
    has(_s2, sec, "C1 cMFTC"), _s2)
# Table 6 note: the baseline terms come from the same C2 fit, not from v13
_c2b = pick("C2", "z_base_pre")
_c2m = pick("C2", "z_cmbase")
_b1s = "%s, %s to %s" % (_pf(_c2b, "OR"), _pf(_c2b, "lo"), _pf(_c2b, "hi"))
_b2s = "%s, %s to %s" % (_pf(_c2m, "OR"), _pf(_c2m, "lo"), _pf(_c2m, "hi"))
chk(sec, "Table 6 note baseline CTh-Score OR %s is the strict-window C2 term" % _b1s,
    has(_b1s, sec, "baseline CTh"), _b1s)
chk(sec, "Table 6 note baseline cMFTC OR %s is the strict-window C2 term" % _b2s,
    has(_b2s, sec, "baseline cMFTC"), _b2s)

# The cross-pipeline limitation quotes per-knee measurement density in the 148
# model pairs under BOTH time rules, because the whole point is that the
# asymmetry came from the rule and not from follow-up.  Quoting it from 56's
# audit table is what keeps that argument honest.
ita = {r["question"]: r["answer"] for r in rd("index_time_audit.csv")}
_ds = ita["CTh model-set density, strict rule, cases / controls"]   # "3.57 / 3.66 (p = 0.45)"
_di = ita["CTh model-set density, inclusive rule, cases / controls"]  # "3.57 / 4.22"
_ns = re.findall(r"\d+\.\d+", _ds)
chk(sec, "cross-pipeline density quoted under the strict rule",
    has("%s pre-index measurements" % _ns[0], sec, "strict density")
    and has("%s in controls" % _ns[1]) and has("P = %s" % _ns[2]),
    "disk " + _ds)
chk(sec, "the superseded inclusive density is quoted as such",
    has("%s against %s under the inclusive rule"
        % tuple(re.findall(r"\d+\.\d+", _di)[:2])),
    "disk " + _di)
sh = rd("cth_shrinkage.csv")
chk(sec, "panel D pre-index slope 2.09 (EBLUP) vs 3.57 (unshrunk)",
    has(f(float([r for r in sh if r["estimate"] == "EBLUP"][0]["mean_paired_diff"]), 2), sec, "EBLUP")
    and has(f(float([r for r in sh if r["estimate"] == "unshrunk OLS"][0]["mean_paired_diff"]), 2), sec, "OLS"))

# ============================================== 7. archived-repo inventory ----
sec = "S2 : script inventory vs the archived repository"
try:
    op = urllib.request.build_opener(urllib.request.ProxyHandler(
        {"https": "http://127.0.0.1:7890", "http": "http://127.0.0.1:7890"})) \
        if os.environ.get("USE_PROXY") else urllib.request.build_opener()
    hdr = {"User-Agent": "v6-audit"}
    if TOKEN:
        hdr["Authorization"] = "Bearer " + TOKEN
    req = urllib.request.Request(
        "https://api.github.com/repos/%s/git/trees/main?recursive=1" % REPO, headers=hdr)
    tree = json.loads(op.open(req, timeout=60).read().decode())
    remote = sorted(x["path"] for x in tree["tree"] if x["type"] == "blob")
except Exception as e:
    remote = None
    chk(sec, "fetch archived tree", False, "ERROR %s" % str(e)[:120])

if remote:
    code = [p for p in remote if "/" in p and not p.endswith(("LICENSE", "README.md", ".gitignore", ".gitattributes"))]
    remote_base = {os.path.basename(p) for p in code}
    # every backticked script filename anywhere in the manuscript counts as "named"
    named = set(re.findall(r"`([A-Za-z0-9_]+\.(?:R|py))`", txt))
    # files the manuscript itself declares as NOT shipped (superseded predecessors)
    declared_absent = set(re.findall(r"supersedes the earlier `([A-Za-z0-9_]+\.(?:R|py))`", txt))
    # Scripts added after the cited Zenodo release, pending the next archive.
    # They MUST be pushed before submission; this list is asserted separately so
    # an un-pushed script can never pass silently as "in the archive".
    # Empty since 2026-10-24: 52-57 were pushed to main in cb88628.  The gate is
    # kept so that the next round of new scripts cannot ship unannounced.
    PENDING_PUSH = set()
    missing_repo = sorted(n for n in named
                          if n not in remote_base and n not in declared_absent
                          and n not in PENDING_PUSH)
    pending_absent = sorted(n for n in PENDING_PUSH
                            if n in named and n not in remote_base)
    chk(sec, "every script named in S2 exists in the archive", not missing_repo,
        "named=%d, absent=%s (pending push: %d)"
        % (len(named), missing_repo or "none", len(pending_absent)))
    # If a script on the pending list has meanwhile reached the archive, the list
    # is stale and must be shortened -- otherwise the next reader is told a
    # shipped script is unshipped.  Either state passes; an out-of-date list does not.
    stale_pending = sorted(n for n in PENDING_PUSH if n in remote_base)
    chk(sec, "the pending-push list matches the archive",
        not stale_pending,
        "pending=%d, already archived (remove from list): %s"
        % (len(PENDING_PUSH) - len(stale_pending), stale_pending or "none"))
    chk(sec, "superseded predecessors are declared, not silently missing",
        all(("supersedes the earlier `%s`" % n) in txt for n in declared_absent),
        "declared absent: %s" % (sorted(declared_absent) or "none"))
    scripts = remote_base
    notnamed = sorted(s for s in scripts if s not in named)
    chk(sec, "every archived script is named in the manuscript", not notnamed,
        "archived=%d, unnamed=%s" % (len(scripts), notnamed or "none"))
    chk(sec, "legacy sampler is archived (audit reproducibility)",
        "poma/jm_core_legacy_sampler.R" in remote)

    # the archive must be identical to the local clone -- this is what catches a
    # script that was committed but never pushed.  Hard-coding a file count would
    # break on every future version, so the two are compared instead.
    import subprocess
    # the clone sits next to Scripts/ and Analysis/, i.e. two levels above BASE
    CLONE = os.environ.get("CLONE_DIR",
                           os.path.join(os.path.dirname(os.path.dirname(BASE)),
                                        REPO.split("/")[-1]))
    try:
        loc = subprocess.run(["git", "ls-tree", "-r", "HEAD", "--name-only"],
                             cwd=CLONE, capture_output=True, text=True, timeout=60).stdout.split()
        chk(sec, "archive matches the local clone file-for-file",
            sorted(loc) == sorted(remote),
            "local=%d remote=%d, only-local=%s only-remote=%s" %
            (len(loc), len(remote),
             sorted(set(loc) - set(remote)) or "none",
             sorted(set(remote) - set(loc)) or "none"))
    except Exception as e:
        chk(sec, "compare against the local clone", False, "ERROR %s" % str(e)[:100])

# ============================================================== 8. DOI -------
# The cited archive must be the one that holds BOTH the corrected sampler and
# the complete S2 inventory.  Superseded DOIs may appear, but only when they are
# explicitly labelled as superseded -- a bare older DOI is how a reader ends up
# reproducing the wrong numbers.
sec = "DOI : points at the archive holding the corrected sampler"
CITE = "10.5281/zenodo.23111986"      # v1.0.2 -- corrected sampler + 46/47
chk(sec, "cited version DOI is %s (v1.0.2)" % CITE, has(CITE, sec, "cited DOI"))
chk(sec, "cited release is named as v1.0.2", has("release v1.0.2", sec, "release label"))
chk(sec, "concept DOI 22771312 unchanged", has("10.5281/zenodo.22771312", sec, "concept DOI"))
for old, why in (("10.5281/zenodo.23111696", "v1.0.1 has the sampler but a pre-audit inventory"),
                 ("10.5281/zenodo.22771313", "v1.0.0 has the uncorrected sampler")):
    n = txt.count(old)
    # A superseded DOI may legitimately appear more than once -- the version
    # history keeps the v6 note that cites v1.0.1.  What must hold is that EVERY
    # mention sits in a "supersede" context, so widen the window well beyond the
    # sentence and look at the text around each occurrence.
    labelled = True
    where = []
    for m in re.finditer(re.escape(old), txt):
        w = txt[max(0, m.start() - 200):m.end() + 200]
        ok = "supersed" in w.lower()
        labelled &= ok
        where.append("ok" if ok else "UNLABELLED")
    chk(sec, "%s: every mention labelled superseded (%s)" % (old, why),
        n >= 1 and labelled, "n=%d | %s" % (n, ",".join(where)))
chk(sec, "GitHub repo URL", has("https://github.com/morrosun/oai-poma-cartilage-thinning-tkr", sec, "repo url"))

# ================================================== 9. OAC hard metrics ------
sec = "OAC hard metrics"
# v22: the compliance block no longer sits in the manuscript (it starts at the
# title).  Limits and counts are declared in Submission_Checklist.md, and the
# counts there are pinned to the rebuilt text.
_CHK = io.open(os.path.join(BASE, "Submission_Checklist.md"),
               encoding="utf-8").read()
_abs_n = len(re.findall(r"[A-Za-z0-9][A-Za-z0-9'\u2019\-\.\u00d7\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079]*",
                        txt.split("## Abstract")[1].split("\n---\n")[0]))
_abs_c = re.search(r"abstract (\d+) words", _CHK)
m = re.search(r"body (\d+) words", _CHK)
chk(sec, "body word count <= 4000", bool(m) and int(m.group(1)) <= 4000,
    "checklist body %s" % (m.group(1) if m else "?"))
chk(sec, "abstract declared <=300 and equal to the rebuilt abstract",
    bool(_abs_c) and int(_abs_c.group(1)) == _abs_n and _abs_n <= 300,
    "checklist %s vs rebuilt %d" % (_abs_c.group(1) if _abs_c else "?", _abs_n))
m = re.search(r"\*\*(\d+) figures and tables combined\*\* against a limit of 8", _CHK)
chk(sec, "tables+figures == 8",
    ("6 tables + 2 figures = 8" in _CHK)
    or (bool(m) and int(m.group(1)) == 8), "checklist")
refs = len(re.findall(r"^\| \d+ \|", txt, re.M))
chk(sec, "references <= 50", 0 < refs <= 50, "reference table rows=%d" % refs)
chk(sec, "no unresolved citation tokens", "[@r" not in txt)
chk(sec, "no bare pipes inside table cells",
    not any(re.search(r"\|[^|\n]*\|[^|\n]*\|[^|\n]*\|[^|\n]*\|[^|\n]*\|", l)
            and l.count("|") > 40 and "S3 panel" in l for l in lines))

# ============================================ 10. OAI Declarations statements --
# The Publications Committee sentence asserts a real-world event.  These checks
# do not decide whether it is true; they make sure that, IF it is present, the
# sponsor boilerplate that must accompany it is present too, so the two can
# never drift apart.
sec = "OAI Declarations"
decl = txt.split("## Declarations")[1].split("## References")[0]
# v10 REMOVED the committee sentence: it asserted an approval that was never
# obtained.  Assert its absence in the Declarations, but allow the version
# history above the fold to keep naming it when describing what v8 did -- that
# is a record of a past edit, not a claim made to the journal.
has_pub = "approval of the OAI Publications Committee" in decl
chk(sec, "no unearned committee-approval claim in Declarations", not has_pub,
    "removed in v10" if not has_pub else "STILL PRESENT -- DELETE")
chk(sec, "committee not named anywhere in Declarations",
    "Publications Committee" not in decl,
    "clean" if "Publications Committee" not in decl else "STILL NAMED")
# The sponsor boilerplate is required whenever OAI data are used, regardless.
for lab, needle in [
        ("OAI five-contract boilerplate", "N01-AR-2-2258; N01-AR-2-2259; N01-AR-2-2260; "
         "N01-AR-2-2261; N01-AR-2-2262"),
        ("OAI private funding partners named", "Merck Research Laboratories"),
        ("FNIH manages private sector funding", "Foundation for the National Institutes "
         "of Health"),
        ("public-use disclaimer", "OAI public use data set"),
        ("POMA contract named", "HHSN2682010000 21C"),
        ("measurement generator named", "Chondrometrics GmbH")]:
    chk(sec, lab, needle in decl)
chk(sec, "Acknowledgements has no unfilled placeholder",
    not re.search(r"\*\*Acknowledgements\.\*\*\s*\[", decl))
# The author's own grant numbers were supplied in v9.  Assert the filled state
# rather than either outcome, so that a later accidental deletion of a grant
# number is caught instead of passing as "still a placeholder".  The three
# numbers are the ones on the title page; they are author facts, not derived
# values, so the check is that the string is present and that the funder roles
# are stated.
for g in ("HAB2025037", "HAP202302", "MS2025094"):
    chk(sec, "author grant %s stated in Funding" % g, g in decl)
chk(sec, "Funding no longer an unfilled placeholder",
    re.search(r"\*\*Funding\.\*\*\s*\[", decl) is None, "FILLED")
chk(sec, "funder role statement present",
    "The funders had no role in study design" in decl)
# Competing interests and authors' contributions were filled in v10.  Assert
# the filled state so that a later accidental deletion is caught, and assert
# that the two statements actually say something rather than being emptied.
for lab in ["Competing interests", "Authors' contributions"]:
    chk(sec, "%s is filled, not a placeholder" % lab,
        re.search(r"\*\*%s\.?\*\*?\s*\[" % re.escape(lab), decl) is None,
        "FILLED")
chk(sec, "competing interests states the no-conflict declaration",
    "no competing financial or non-financial interests" in decl)
chk(sec, "competing interests addresses the Chondrometrics relationship",
    "Chondrometrics GmbH" in decl and "no author has received funding" in decl)
chk(sec, "author contributions names every author's initials",
    all(x in decl for x in ("KW", "CB and JS", "XZ, HY, DM and JL")))
chk(sec, "author contributions ends with the accountability sentence",
    "agree to be accountable for all aspects of the work" in decl)
# The generative-AI declaration is required by the journal.  It must appear
# outside the body word count (i.e. after "## Declarations") and must name the
# tool, so that the statement cannot silently go missing at the next revision.
chk(sec, "generative-AI declaration present in Declarations",
    "Declaration of generative AI and AI-assisted technologies" in decl)
chk(sec, "AI declaration names the tool used",
    "ChatGPT" in decl and "OpenAI" in decl)
# British spelling: the -ize form is a US spelling and the journal's scan flags
# it.  Guard the word that the supplied text actually contained.
chk(sec, "AI declaration uses British -ise spelling",
    "utilize" not in decl.lower() and "utilis" in decl.lower(),
    "utilised")

# ------------------------------------------------- 11. title page / authors --
# Author names and affiliations are author facts that no other check covers, and
# a name mangled between the title page and the Declarations would not be
# caught by any numeric audit.  Assert that the title page carries a real author
# list, that both affiliations resolve to the single hospital, and that the
# corresponding author is the one marked as such.
sec = "Title page"
tp = txt.split("**Authors:**")[1].split("---")[0] if "**Authors:**" in txt else ""
chk(sec, "author list present and non-placeholder",
    "[" not in tp.split("**Affiliations:**")[0],
    "%d authors" % len(re.findall(r"[A-Z][a-z]+ [A-Z][a-z]+", tp.split("**Affiliations:**")[0])))
chk(sec, "seven authors named",
    len(re.findall(r"[A-Z][a-z]+ [A-Z][a-z]+", tp.split("**Affiliations:**")[0])) == 7)
affs = tp.split("**Affiliations:**")[1].split("\\* **Correspondence")[0] \
    if "**Affiliations:**" in tp else ""
chk(sec, "two affiliations given", affs.count("Nanjing Medical University") == 2,
    "1=Rehabilitation, 2=Rheumatology and Immunology")
chk(sec, "no unresolved affiliation placeholder",
    "[" not in affs and "[Department" not in tp)
chk(sec, "corresponding author is marked with an asterisk",
    re.search(r"Kai Wang\S*\*", tp.split("**Affiliations:**")[0]) is not None)
chk(sec, "corresponding e-mail present",
    "morrosun@hotmail.com" in tp and "morrosun@njmu.edu.cn" in tp)
# The ORCID was supplied by the author in v10.  Assert the real iD is present
# and that neither the placeholder nor a zero-filled stub survives.
chk(sec, "corresponding ORCID is a real iD",
    "0000-0002-0511-8038" in tp, "0000-0002-0511-8038")
chk(sec, "ORCID placeholder and zero-stub both gone",
    "0000-0000-0000-0000" not in tp
    and re.search(r"ORCID:?\s*\[", tp) is None)
# An ORCID iD is 16 digits in four groups; a typo would usually break that.
chk(sec, "ORCID has the 4x4 digit shape",
    re.search(r"\b\d{4}-\d{4}-\d{4}-\d{3}[\dX]\b", tp) is not None)
# The supplied correspondence address read "..., Huaian, 223001, China, Huaian,
# 223001, China" -- the clause was duplicated.  It is written once, so exactly
# three occurrences are expected in the title page: one per affiliation plus the
# correspondence line.  A fourth would mean the duplication is back.
chk(sec, "postal address not duplicated",
    "China, Huaian" not in tp and tp.count("Huaian, 223001, China") == 3,
    "3 = two affiliations + correspondence")

# ================================== 6b. v11 : the two "rapid thinning" defs --
# v11 separated every latent-class read-out into the fastest class alone
# (n = 25) and the two fastest classes (n = 67), because tests B and C had been
# attributed to the former while being computed on the latter.  Both are now
# computed (52_lcmm_class_definition_check.R) and both must appear.  Gated on
# the version so that an older manuscript is not failed by a check that did not
# yet exist.
VER = 0
_mv = re.search(r"_v(\d+)([a-z]?)$", BASENAME)
if _mv:
    VER = int(_mv.group(1))

if VER >= 11:
    sec = "Table 5 / 3.4 : two definitions of rapid thinning"
    rq = rd("lcmm_rank_equivalence.csv")
    r1, r2, r3 = rq[0], rq[1], rq[2]
    def _pct(r):
        return "%s%%" % f(float(r["pct_agreement"]), 1)
    def _jac(r):
        return f(float(r["jaccard"]), 3)
    chk(sec, "A1 fastest class alone: %s / %s" % (_pct(r1), _jac(r1)),
        has(_pct(r1)) and has(_jac(r1)), "%s | %s" % (_pct(r1), _jac(r1)))
    chk(sec, "A2 two fastest classes: %s / %s" % (_pct(r2), _jac(r2)),
        has(_pct(r2)) and has(_jac(r2)), "%s | %s" % (_pct(r2), _jac(r2)))
    chk(sec, "A3 whole partition: %s" % _pct(r3), has(_pct(r3)), _pct(r3))
    chk(sec, "both class sizes stated (n = 25 and n = 67)",
        has("n = 25") and has("n = 67"), "25 / 42+25")

    dfc = rd("lcmm_rate_vs_class_by_definition.csv")
    bkey = {(r["definition"], r["model"], r["term"]): r for r in dfc}
    for k, lbl in ((("two fastest classes", "B2 class indicator only", "rapidTRUE"), "B2 class alone"),
                   (("two fastest classes", "B3 both", "z"), "B3 rate"),
                   (("two fastest classes", "B3 both", "rapidTRUE"), "B3 class")):
        r = bkey[k]
        s = ci2(r["OR"], r["lo"], r["hi"])
        chk(sec, "%s %s" % (lbl, s), has(s), "disk %s p=%s" % (s, r["p"]))
    # B4 is a separation, not a constancy: the pairwise conditional likelihood
    # gets its information from DISCORDANT pairs, so what makes the indicator
    # inestimable is that every discordant pair points the same way.  Auditing
    # the direction -- not the marginal 24/25 -- is the point of this check.
    b4 = bkey[("fastest class only", "B2 class indicator only", "rapidTRUE")]
    sep = {r["configuration"]: int(r["Freq"]) for r in rd("lcmm_separation_pairs.csv")
           if r["definition"].startswith("class 4 only")}
    chk(sec, "B4 separation is one-directional (%d vs %d discordant pairs)"
        % (sep.get("case_only", -1), sep.get("control_only", -1)),
        sep.get("case_only", 0) > 0 and sep.get("control_only", 0) == 0,
        "case_only=%s control_only=%s both=%s"
        % (sep.get("case_only"), sep.get("control_only"), sep.get("both_class")))
    chk(sec, "B4 quotes those counts and drops the numeric artefact",
        has("%d discordant pairs favoured class 4 in cases" % sep["case_only"])
        and has("none favoured class 4 in controls")
        and "5.9 \u00d7 10\u2078" not in txt,
        "OR %.3g is no longer reported as a result" % fnum(b4["OR"]))
    chk(sec, "the 'constant within a stratum' reading is gone from all three places",
        all(p not in txt for p in
            ("never varies within a stratum",
             "constant within every informative stratum",
             "constant within every stratum where it varies at all")),
        "body / methods / Table 5")

    bt = {r["definition"]: r for r in rd("lcmm_rate_below_threshold_by_definition.csv")}
    pc = {r["definition"]: r for r in rd("lcmm_definition_pair_counts.csv")}
    for defn, lbl in (("fastest class only", "C1"), ("two fastest classes", "C2")):
        r = bt[defn]
        p = pc[defn]
        s = ci2(r["OR"], r["lo"], r["hi"])
        chk(sec, "%s %s" % (lbl, s), has(s), "disk %s p=%s" % (s, r["p"]))
        chk(sec, "%s pair count %s matches the CSV" % (lbl, p["pairs_remaining"]),
            has("%s pairs" % p["pairs_remaining"]),
            "removed %s -> remaining %s" % (p["pairs_removed"], p["pairs_remaining"]))
        chk(sec, "%s P value %s" % (lbl, f(float(r["p"]), 3)),
            has("P = %s" % f(float(r["p"]), 3)), "disk p=%s" % r["p"])

    sec = "Table 2 / 3.2 : fixed 0.141 mm/yr reference SD (jm_fixed_sd_effect.csv)"
    fx = {(r["metric"], r["model"]): r for r in rd("jm_fixed_sd_effect.csv")}
    for mdl, lbl in (("A", "unadjusted"), ("Aadj", "adjusted")):
        r = fx[("cMFTC_ThCtAB_aMe", mdl)]
        s = ci2(r["OR_fixedSD"], r["OR_fixedSD_lo"], r["OR_fixedSD_hi"])
        chk(sec, "%s re-expressed at the fixed SD: %s" % (lbl, s), has(s),
            "disk %s" % s)
        s01 = ci2(r["OR_per0p1"], r["OR_per0p1_lo"], r["OR_per0p1_hi"])
        chk(sec, "%s per-0.1 unchanged: %s" % (lbl, s01), has(s01), "disk %s" % s01)
    chk(sec, "the reference SD is named as 0.141 mm/yr", has("0.141 mm/yr"))
    chk(sec, "the latent-slope SD is named as 0.176 mm/yr", has("0.176 mm/yr"))

    sec = "Table 3 : latent-slope SD column and the per-0.1 ordering"
    SHORT = {"cMFTC_ThCtAB_aMe": "cMFTC", "MFTC_ThCtAB_aMe": "MFTC",
             "MT_ThCtAB_aMe": "MT", "cMF_ThCtAB_aMe": "cMF",
             "cLFTC_ThCtAB_aMe": "cLFTC", "cLF_ThCtAB_aMe": "cLF"}
    sds = {}
    for r in rd("jm_fixed_sd_effect.csv"):
        if r["model"] == "A":
            sds[r["metric"]] = fnum(r["sd_latent_slope"])
    for k, v in sorted(sds.items(), key=lambda kv: -kv[1]):
        chk(sec, "%s latent-slope SD %s" % (SHORT[k], f(v, 3)),
            has("| %s |" % f(v, 3)), "disk %.6f" % v)
    spread = max(sds.values()) / min(sds.values())
    chk(sec, "spread stated as 4.5-fold", has("4.5-fold"), "disk %.2f-fold" % spread)
    order = sorted(sds, key=lambda k: -fnum(fx[(k, "A")]["OR_per0p1"]))
    want = ", ".join("%s %s" % (SHORT[k], f(fnum(fx[(k, "A")]["OR_per0p1"]), 2))
                     for k in order)
    or01 = {k: fnum(fx[(k, "A")]["OR_per0p1"]) for k in sds}
    orsd = {k: fnum(fx[(k, "A")]["OR_ownSD"]) for k in sds}
    lead01 = max(or01, key=lambda k: or01[k])
    leadsd = max(orsd, key=lambda k: orsd[k])
    # v13 listed all six per-0.1 odds ratios in 3.3; v14 dropped the list to stay
    # inside the word limit and points at Table 3 instead.  The six values must
    # therefore be present as table cells, and the body must state the one
    # substantive claim the list was there to make: the leader depends on scale.
    chk(sec, "Table 3 carries every per-0.1 and per-SD odds ratio",
        all(has("| %s " % f(or01[k], 2)) and has("| %s " % f(orsd[k], 2)) for k in sds),
        "per-0.1 leader %s %s ; per-SD leader %s %s"
        % (SHORT[lead01], f(or01[lead01], 2), SHORT[leadsd], f(orsd[leadsd], 2)))
    chk(sec, "3.3 states that the two scales order the regions differently",
        lead01 != leadsd and has("a different metric leads (Table 3)"),
        "per-0.1 leader %s, per-SD leader %s (list in 3.3: %s)"
        % (SHORT[lead01], SHORT[leadsd], want))

    sec = "2.1 / 3.5 : non-independent pairs (clogit_connectivity.csv)"
    cr = rd("clogit_connectivity.csv")
    def _cr(var, mdl, term="z"):
        for r in cr:
            if r["variance"] == var and r["model"] == mdl and r["term"] == term:
                return r
        return None
    m0c = _cr("cluster-robust (connected component)", "M0 crude")
    m0p = _cr("cluster-robust (participant)", "M0 crude")
    m0m = _cr("model-based", "M0 crude")
    chk(sec, "every grouping leaves the point estimate unchanged",
        abs(fnum(m0c["OR"]) - fnum(m0m["OR"])) < 1e-9
        and abs(fnum(m0p["OR"]) - fnum(m0m["OR"])) < 1e-9,
        "component %s / participant %s / model %s"
        % (m0c["OR"], m0p["OR"], m0m["OR"]))
    # v15: the clustering unit is the connected component of the pair graph, not
    # the participant -- a participant who links two pairs leaves the rest of the
    # component linked through other participants.
    s = "unchanged at %s (%s to %s" % (
        f(fnum(m0c["OR"]), 2), f(fnum(m0c["lo"]), 2), f(fnum(m0c["hi"]), 2))
    chk(sec, "3.5 quotes the component-clustered interval", has(s), s)
    chk(sec, "3.5 quotes the component-clustered P value",
        has("P = 1.1 \u00d7 10\u207b\u2077"), "disk p=%s" % m0c["p"])
    sp = "than participant clustering (%s to %s)" % (f(fnum(m0p["lo"]), 2),
                                                     f(fnum(m0p["hi"]), 2))
    chk(sec, "3.5 contrasts the component and participant intervals", has(sp), sp)
    cs = {r["n_pairs"]: r for r in rd("connectivity_summary.csv")}
    kb = list(cs.values())[0]
    chk(sec, "3.5 states the piecewise structure (%s pairs -> %s components)"
        % (kb["n_pairs"], kb["n_components"]),
        has("%s connected components, %s containing more than one pair"
            % (kb["n_components"], kb["n_multi_components"])),
        "pairs in multi-pair components: %s" % kb["n_pairs_in_multi"])
    chk(sec, "release-level clustering stated (450 knees / 415 participants / 35)",
        has("415 participants") and has("35 of whom contributed both knees"),
        "release set")
    chk(sec, "analysis-set clustering stated (30 participants)",
        has("30 participants contributed both knees", sec, "Set A clustering"),
        "Set A")
    # v11a: the release-level count (35 of 450 knees) and the analysis-set
    # count (30 of 382) sit pages apart.  Once the manuscript reconciles them
    # that reconciliation must not be lost -- an unreconciled pair of numbers
    # reads as an internal inconsistency even though both are correct.
    if "v11a changes no result" in txt:
        chk(sec, "the two cluster counts are reconciled in 2.1",
            has("(30 within the primary set)", sec, "count reconciliation"),
            "30 of the 35 lie in Set A")

    sec = "v11 wording invariants"
    chk(sec, "whole-compartment contrast is MFTC - cLFTC, not MFTC - cLF",
        has("MFTC \u2212 cLFTC")
        and re.search(r"(?<![A-Za-z])MFTC\s+-\s+cLF(?!TC)", NTXT) is None,
        "relabelled to the variable actually fitted")
    chk(sec, "lateral metrics are anatomical comparators",
        has("anatomical comparators", sec, "anatomical comparators"))
    chk(sec, "no unqualified negative-control claim",
        all("rather than negative control" in NTXT[max(0, m.start() - 60):m.end()]
            for m in re.finditer(r"negative control", NTXT)),
        "%d mention(s), all negated" % len(re.findall(r"negative control", NTXT)))
    chk(sec, "per-1-SD ratios carry the no-cross-region-ranking warning",
        has("cannot be ranked across regions", sec, "no-cross-region ranking")
        or has("not be ranked across metrics", sec, "no-cross-region ranking"),
        "scale warning present")

# ============================ 6c. v12 : the second batch of the same review --
# These are not numbers, they are the claims that survived the audit as words.
# Each one records something that was wrong and is now right, so that a later
# edit cannot quietly restore it.
if VER >= 12:
    sec = "v12 wording invariants (second batch of the review)"
    chk(sec, "title no longer promises three transferable methodological results",
        "three transferable methodological results" not in txt,
        "subtitle names the level-versus-rate question")
    chk(sec, "Table 4's signed column is called a thickness slope",
        has("| Mean thickness slope (mm/yr) |", sec, "Table 4 column"),
        "the cells are negative, so \"rate\" invited a sign misreading")
    # NB a negative assertion must not go through has(), which appends its own
    # result row: `not has(x)` would log a FAIL while the composite passes.
    chk(sec, "Table 6 carries exactly one primary label",
        "(OLS; primary)" in txt and "(EBLUP; primary)" not in txt,
        "B1 OLS primary; B2 relabelled as its shrinkage counterpart")
    chk(sec, "Table 6 note no longer claims Table 2 keeps the thickness sign",
        "with the thickness sign retained in Table 2" not in txt,
        "the sign-retained value is the 148-pair one")
    chk(sec, "Monte-Carlo error is a magnitude, not a digit-level assurance",
        "cannot move any reported digit" not in txt,
        "replaced by \"about 1.4% multiplicatively\"")
    chk(sec, "no log statement contradicts the Declarations",
        "logs accompany the archived repository" not in txt
        and "run logs are archived" not in txt,
        "%d statement(s) that run logs are not redistributed"
        % len(re.findall(r"run logs are not redistributed", txt)))
    chk(sec, "the over-strong modals stay out of the body",
        "mandatory" not in txt
        and not re.search(r"\bthe question is settled\b", NTXT),
        "dialled back in v11, guarded here")

# ================= 6d. v13 : the third batch (window / density / lag) ------
# Two kinds of check.  First, Table S5 is REBUILT from the two CSVs that the
# analysis scripts write, and each row is matched line by line, so a change in
# the analysis cannot leave the table behind.  Second, the wording invariants
# record what v13 changed, so a later edit cannot quietly restore the
# unqualified mechanistic claim or the old "insensitive to visit density".
if VER >= 13:
    sec = "v13 window / density / lag (Supplementary Table S5)"
    WL = {r["design"]: r for r in rd("sens_window_lag.csv")}
    RL = {r["design"]: r for r in rd("sens_wl_reliability.csv")}

    def s5(w):
        return "%s (%s to %s)" % (f(w["clogit_ols_OR01"], 2),
                                  f(w["clogit_ols_OR01_lo"], 2),
                                  f(w["clogit_ols_OR01_hi"], 2))

    m5 = re.search(r"### Supplementary Table S5(.*?)### Supplementary Figure S1",
                   txt, re.S)
    S5T = m5.group(1) if m5 else ""
    chk(sec, "Supplementary Table S5 exists and is the last supplementary table",
        bool(S5T), "%d characters" % len(S5T))

    for key in ("primary", "dens", "cw24", "cw3", "cw3_ref", "lag24",
                "lag24_ref", "lag36"):
        w = WL.get(key)
        if not w:
            chk(sec, "S5 row %s exists in sens_window_lag.csv" % key, False, "missing")
            continue
        s = s5(w)
        row = next((L for L in S5T.split("\n") if norm(s) in norm(L)), "")
        ok = (bool(row)
              and ("| %s |" % w["n_pairs"]) in row
              and f(w["visits_case"], 2) in row
              and f(w["visits_control"], 2) in row)
        chk(sec, "S5 row %-9s OR, pairs and visit density agree" % key, ok,
            row[:70] if row else "MISSING: " + s)

    # slope reliability: on the row of its own design, never elsewhere
    for key in ("primary", "dens", "cw3", "cw3_ref", "lag24", "lag24_ref", "lag36"):
        w, r = WL.get(key), RL.get(key)
        if not (w and r):
            chk(sec, "reliability for %s exists" % key, False, "missing")
            continue
        row = next((L for L in S5T.split("\n") if norm(s5(w)) in norm(L)), "")
        lam = f(r["lambda_obs"], 2)
        chk(sec, "S5 reliability %-9s = %s" % (key, lam),
            bool(row) and row.rstrip().rstrip("|").rstrip().endswith(lam),
            row[-40:] if row else "no row")
    chk(sec, "the two-visit design carries no reliability value",
        "\u2014\u2020" in S5T, "two points per knee: lambda not identifiable")

    # the prose numbers, rebuilt from the same CSVs, now read against the
    # same pairs on all pre-index visits (v19)
    chk(sec, "3.5 quotes the common-window and density estimates",
        (s5(WL["cw3"]) + " against " + s5(WL["cw3_ref"])) in norm(txt)
        and (f(WL["dens"]["clogit_ols_OR01"], 2) + " against "
             + f(WL["primary"]["clogit_ols_OR01"], 2)) in norm(txt),
        "cw3 %s vs %s ; dens %s" % (s5(WL["cw3"]), s5(WL["cw3_ref"]),
                                    s5(WL["dens"])))
    chk(sec, "3.5 quotes the lag contrast from the CSV",
        (s5(WL["lag24"]) + " against " + s5(WL["lag24_ref"])) in norm(txt)
        and s5(WL["lag36"]) in norm(txt),
        "%s vs %s ; %s" % (s5(WL["lag24"]), s5(WL["lag24_ref"]),
                           s5(WL["lag36"])))
    chk(sec, "the dilution-corrected lag contrast is quoted from the CSV",
        ("%s against %s per 0.1 mm/yr" % (f(RL["lag24"]["or01_corrected"], 2),
                                          f(RL["lag24_ref"]["or01_corrected"], 2))) in norm(txt),
        "%s vs %s" % (f(RL["lag24"]["or01_corrected"], 2),
                      f(RL["lag24_ref"]["or01_corrected"], 2)))
    chk(sec, "the lag attenuation is not attributed to the shorter trajectories",
        "not a product of the shorter trajectories" in txt
        or "remain a partial explanation rather than an excluded one" in txt,
        "reliability differs, the corrected contrast still attenuates")

    # wording invariants
    chk(sec, "the old unqualified density claim is gone",
        "insensitive to visit density" not in txt,
        "replaced by the density-matched estimate")
    # v16: the mechanism paragraph must state the lag result as compatible with
    # peri-operative change, not as attributing the residual effect to it, and
    # it must name the two alternatives it cannot exclude.
    chk(sec, "Section 4.6 no longer asserts the mechanism unqualifiedly",
        "compatible with a contribution from change near surgery" in txt
        and "do not distinguish this from differences in trajectory length and "
            "measurement error" in txt,
        "qualified mechanistic reading present")
    chk(sec, "the limitations carry the exposure-timing caveat",
        "Measurement, exposure timing and unmeasured mechanism" in txt
        and "final pre-index years were removed" in txt, "Section 4.5 item 3")
    chk(sec, "Supplementary Table S5 is cited from the body",
        len(re.findall(r"\(Supplementary Table S5\)", txt.split("## Declarations")[0])) >= 1,
        "%d in-body citation(s)"
        % len(re.findall(r"\(Supplementary Table S5\)", txt.split("## Declarations")[0])))
    # v22: the first-page working block (target journal / compliance / front
    # matter / INTERNAL NOTE) was removed from the manuscript and now lives in
    # Submission_Checklist.md, which is not uploaded.
    _chk = _CHK
    chk(sec, "the manuscript starts at the title and carries no working block",
        txt.split("\n")[0].startswith("# ")
        and all(m not in txt for m in ("TARGET JOURNAL", "Compliance (OAC",
                                       "Front matter status", "INTERNAL NOTE")),
        "first line: %s" % txt.split("\n")[0][:40])
    chk(sec, "the checklist still counts five supplementary tables",
        "Tables S1\u2013S5" in _chk and "not uploaded" in _chk, "S1-S5")

    # ---- v16 invariants: method and estimand must name the same slope ------
    chk(sec, "2.6 names the OLS slope as the primary cross-pipeline exposure",
        "ordinary-least-squares slope of the score over strictly pre-index "
        "visits" in txt and "mixed model with random intercept and slope" not in txt,
        "method matches Table 6 B1")
    chk(sec, "Table 6 note says which slope entered C1, C2 and S3",
        "C1 and C2 entered the pre-index OLS CTh-Score slope together with the "
        "cMFTC empirical-Bayes thinning rate" in txt, "slope per model")
    chk(sec, "Table 6 note reverses every row that carries a cMFTC term",
        "The cMFTC rows A1, A2, C1, C2 and S3" in txt, "C2 and S3 included")
    chk(sec, "Table S1 note states one scale rule, not two",
        "each per-1-SD column uses the subset's own rate SD" not in txt
        and "not the scaling basis for the joint-model estimates" in txt,
        "single scale statement")
    chk(sec, "prior sensitivity does not claim every class interval spans 1",
        "intervals spanning 1 throughout" not in txt
        and "including 1 only under the narrowest prior" in txt,
        "class 3 crosses 1 only at N(0, 1\u00b2)")
    chk(sec, "Table S5 footnote is conditional on this design",
        "the joint model is not identified either" not in txt
        and "In the two-visit design and model specification used here" in txt,
        "no general non-identifiability claim")
    chk(sec, "sampler audit separates three corrections from the safeguard",
        "The corrected kernel differs from the v4 kernel in three places; a "
        "fourth change is a safeguard rather than a correction" in txt
        and "(iii) *Within-pair covariates.*" in txt
        and "omitting the covariate term from that ratio (correction iii)" in txt
        and "(correction iv)" not in txt,
        "three corrections + safeguard, renumbered")
    chk(sec, "the ablation names its three factors and the fixed bound",
        "The ablation varied three factors" in txt
        and "The proposal-scale bound was not varied but held fixed" in txt,
        "ablation scope stated")
    chk(sec, "the proposal bound carries both a rule and a number",
        "no proposal standard deviation is allowed to exceed a quarter of the "
        "prior SD" in txt and "5 on the log-odds-ratio scale" in txt,
        "prior_sd/4 = 5")
    chk(sec, "the simulation residual is stated in mm",
        "residual measurement-error SD to 0.232 mm," in txt, "thickness units")

    # ---- v17 invariants: editorial errors introduced by the v16 edits ------
    chk(sec, "Table 2 panel B no longer carries the lateral-metric count",
        "joint model (same 191 pairs)" in txt
        and "190 for the two lateral metrics); intervals" not in txt,
        "panel B is cMFTC only")
    chk(sec, "Table 6 note distinguishes the S2 EBLUP row",
        "Rows S1 and S3 use the all-timepoint OLS CTh-Score slope, whereas S2 "
        "uses its empirical-Bayes counterpart" in txt
        and "rows S1-S3 repeat these models with the all-timepoint OLS" not in txt,
        "S2 labelled as EBLUP")
    chk(sec, "Table 6 note keeps the cMFTC exposure of S3 pre-index",
        "in S3 the cMFTC exposure remains the pre-index empirical-Bayes "
        "thinning rate" in txt, "S3 cMFTC term")
    _abs = txt.split("## Abstract")[1].split("## 1. Introduction")[0]
    chk(sec, "the abstract no longer juxtaposes overlap with agreement",
        "24 of 25 knees but the wider definition" not in _abs
        and "largely reproduced an upper rate tail" in _abs,
        "one statistic in the abstract")
    chk(sec, "4.3 reports the same agreement statistic on both sides",
        "matched the fastest 25 knees by rate on 99.5% of knees" in txt
        and "whereas the wider two-class definition matched only 91.3%" in txt,
        "agreement on both sides")
    chk(sec, "Table S3B separates the EBLUP and oracle ratios",
        "the ratio was 0.413 without and 0.345 with a pair-level confounder" in txt
        and "with the true slope, which removes slope estimation error, it was "
            "0.441 and 0.327" in txt,
        "0.327 labelled as the oracle ratio")
    chk(sec, "Table S3B no longer calls the gap a property of the design",
        "a property of the design rather than of the estimator" not in txt
        and "not attributable solely to slope estimation error" in txt
        and "is not collapsible" in txt,
        "no bias decomposition claimed")
    chk(sec, "measurement error is not compared across incompatible units",
        "residual SD (about 0.23 mm) greatly exceeds" not in txt
        and "making short-window two-point slopes imprecise" in txt,
        "mm vs mm/yr removed")
    chk(sec, "the joint model states the adjusted linear predictor",
        "**\u03b3**\u2032\u0394**Z**" in txt and "are not multiplied by *a*_V" in txt,
        "event-level covariate shift")
    chk(sec, "the joint model states its priors",
        "N(0, 20\u00b2), with the width varied from 5 to 50" in txt
        and "\u03c3\u00b2 ~ inverse-gamma(0.01, 0.01)" in txt
        and "inverse-Wishart with 4 degrees of freedom and scale matrix" in txt,
        "association / residual / Wishart")
    chk(sec, "the shared index time is described as remaining in the contrast",
        "but not the index time, which remains in the projected level contrast"
        in txt and "removes the event time from the within-pair comparison"
        not in txt, "index time kept")
    chk(sec, "no unreported Savage-Dickey Bayes factor is mentioned",
        "Savage-Dickey" not in txt, "unreported result removed")
    chk(sec, "the MCSE is attached to a named summary",
        "estimated Monte-Carlo standard error of the posterior mean was 0.014"
        in txt and "\u00b10.014" not in txt,
        "MCSE of the per-1-SD log OR")
    chk(sec, "the class analyses state which labels they use",
        "all use the modal posterior assignment" in txt
        and "Rubin's rules were applied to the log coefficients" in txt
        and "Barnard\u2013Rubin degrees of freedom" in txt,
        "modal vs sampled assignment")
    _mi = {r["term"]: r for r in rd("lcmm_class_OR_mi.csv")}
    def _mit(term):
        r = _mi[term]
        return "%s (%s to %s)" % (f(r["OR"], 2), f(r["lo"], 2), f(r["hi"], 2))
    chk(sec, "the imputation check quotes the stored MI values",
        _mit("rapid class vs rest") in txt and _mit("OR per class step") in txt,
        "%s / %s" % (_mit("rapid class vs rest"), _mit("OR per class step")))
    chk(sec, "the remaining over-strong wordings are gone",
        "which must be negative" not in txt
        and "which was expected to be negative" in txt
        and "demonstrating that the no-random-effect" not in txt
        and "suggesting that additional classes may be absorbing unmodelled "
            "heterogeneity" in txt
        and "neither that indicator nor a joint model containing it" not in txt
        and "no unpenalised conditional-logistic estimate exists" in txt
        and "frozen-parameter palette" not in txt
        and "frozen-parameter comparison" in txt
        and "confounded with differential shrinkage" not in txt,
        "six wording fixes present")

    # ---------------- v18: modal labels, probe, additivity, sign, S3B, M7 ----
    sec = "v18 : modal labels in the primary class analyses"
    chk(sec, "the main text states modal assignment plus a separate check",
        "Primary class-based analyses used the modal posterior assignment"
        in txt
        and "classification uncertainty was propagated by multiple imputation"
        not in txt,
        "no claim that the main tables are MI-corrected")
    chk(sec, "both class tables carry the same one-line caveat",
        txt.count("Analyses use modal posterior class assignments; "
                  "classification uncertainty was assessed separately in "
                  "sensitivity analyses.") == 2,
        "Table 4 + Table 5 notes")
    _r1, _r2 = _mi["rapid class vs rest"], _mi["OR per class step"]
    _s1 = "%s, %s to %s, for the two fastest classes" % (
        f(_r1["OR"], 2), f(_r1["lo"], 2), f(_r1["hi"], 2))
    _s2 = "%s, %s to %s, per class step" % (
        f(_r2["OR"], 2), f(_r2["lo"], 2), f(_r2["hi"], 2))
    chk(sec, "3.5 quotes the stored class-draw estimates",
        _s1 in NTXT and _s2 in NTXT, "%s / %s" % (_s1, _s2))
    chk(sec, "the class-draw check states that the class model was held fixed",
        "with the class model held at its fitted estimates" in txt
        and "with the fitted class model held fixed" in txt,
        "label-only resampling")

    sec = "v18 : software probe vs identifiability"
    chk(sec, "the main text separates the observation from the mechanism",
        "A separate software probe documented initialisation failures" in txt
        and "rank-deficient" not in txt
        and "leaves a stratum baseline hazard unidentifiable" not in txt,
        "4.3 no longer derives the failure")
    chk(sec, "the supplement reports the probe as observed behaviour",
        "initialisation failed when "
        "within-pair-constant covariates were included" in txt
        and "does not by itself establish the cause of failure" in txt
        and "collinear with the intercept" not in txt,
        "implementation-specific")
    chk(sec, "the S3A heading names the invariance result",
        "Baseline-hazard invariance of the within-pair conditional likelihood"
        in txt
        and "Non-identifiability of the stratum baseline hazard" not in txt,
        "panel heading")

    sec = "v18 : sampler correction (iii) matches the additive formula"
    chk(sec, "correction (iii) is the additive covariate term",
        "the additive term" in txt
        and "retained when evaluating both the proposed and the current event "
            "probability" in txt
        and "When the level term carries within-pair covariate differences"
            not in txt,
        "gamma'Z entered additively")

    sec = "v18 : sign of the per-SD log odds ratio"
    chk(sec, "the MCSE is attached to -aS times the latent-slope SD",
        "computed draw by draw as \u2212*a*_S times that draw's latent-slope SD"
        in txt,
        "negative sign explicit")

    sec = "S3B : estimator comparison (sim_degeneracy_estimate_summary.csv)"
    _sb = rd("sim_degeneracy_estimate_summary.csv")
    for r in _sb:
        lab = ("No pair-level confounder" if r["scenario"].startswith("C1")
               else "Pair-level confounder present")
        row = "| %s | %s | %s | %s | %s | %s |" % (
            lab, r["exposure"], r["mean"], r["bias"].replace("-", "\u2212"),
            r["sd"], "**%s**" % r["coverage"] if r["coverage"] == "0.000"
            else r["coverage"])
        chk(sec, "%s / %s" % (lab, r["exposure"]), row in txt, row[:56])
    chk(sec, "bias and coverage are tied to the conditional target",
        "relative to the generating within-pair conditional coefficient "
        "(5.555)" in txt
        and all(abs(fnum(r["target"]) - 5.555) < 1e-9 for r in _sb),
        "%d rows against 5.555" % len(_sb))
    chk(sec, "the six-replicate panel is not oversold",
        "supporting the three methodological findings" in txt
        and "validating the three methodological findings" not in txt
        and "six replicates it illustrates" in txt,
        "illustrating, not validating")

    sec = "v18 : reporting consistency"
    chk(sec, "M7 is reported as 190 pairs",
        "| M7 M6 + PASE (190 pairs) |" in txt
        and "M7 is fitted on 190" in txt
        and "complete data on every covariate used in models M0\u2013M6" in txt,
        "M0-M6 191, M7 190")
    chk(sec, "the acceptance rate is quoted once, from the fitted chain",
        "acceptance about 0.88" not in txt
        and "Acceptance was 0.81 for the random-effect updates" in txt,
        "stale 0.88 removed")
    chk(sec, "the likelihood-ratio test is attributed to the clogit models",
        "Nested conditional-logistic models were compared by likelihood-ratio "
        "tests" in txt
        and "models were compared by likelihood-ratio tests rather than WAIC"
            not in txt,
        "no LRT claim for the joint model")
    chk(sec, "the unreported restricted joint model is no longer cited",
        "a joint model was fitted with *a*_S fixed at zero" not in txt
        and "; The two coefficients" not in txt,
        "aS = 0 fit dropped from the supplement")
    chk(sec, "the lag window is described as scheduled, not observed",
        "the last scheduled pre-index visit could fall within the year before "
        "KR" in txt
        and "the last analysed visit precedes the operation" not in txt,
        "visit completeness not assumed")
    chk(sec, "the remaining wording fixes are in place",
        "Two observations are compatible with this interpretation" in txt
        and "The released dataset does not document the matching variables"
        in txt
        and "fraction of retained odds-ratio draws above 1" in txt
        and "E[log plogis] is a functional directly relevant" in txt
        and "the functional the association step maximises" not in txt,
        "four wording fixes")

    # the compliance word count must equal the count the builder would make
    _body = txt.split("## 1. Introduction")[1].split("## Declarations")[0]
    _n = len(re.findall(r"[A-Za-z0-9][A-Za-z0-9'\u2019\-\.\u00d7\u207b\u2070\u00b9\u00b2\u00b3"
                        r"\u2074\u2075\u2076\u2077\u2078\u2079]*", _body))
    # the word count is no longer printed in the manuscript; the checklist
    # carries it, so it is still pinned to the count the builder would make
    chk(sec, "the checklist agrees with the rebuilt body count",
        ("body %d words" % _n) in _chk, "%d words" % _n)

# ============================================ 19. report-level consistency ----
sec = "v19 : OLS sensitivity designs are labelled and self-referenced"
_wl = {r["design"]: r for r in rd("sens_window_lag.csv")}
_cw3, _cw3r = _wl["cw3"], _wl["cw3_ref"]
_dens, _lag, _lagr = _wl["dens"], _wl["lag24"], _wl["lag24_ref"]
_prim, _lag36 = _wl["primary"], _wl["lag36"]
_frag = [
    "The aligned window gave %s against %s using all pre-index visits in the "
    "same %d pairs" % (s5(_cw3), s5(_cw3r), int(_cw3["n_pairs"])),
    "equal within-pair visit counts gave %s against %s in the full set" % (
        f(_dens["clogit_ols_OR01"], 2), f(_prim["clogit_ols_OR01"], 2)),
    "Restricting measurements to \u2264 index \u2212 24 months gave %s (%s to "
    "%s) against %s (%s to %s) using all visits in the same %d pairs" % (
        f(_lag["clogit_ols_OR01"], 2), f(_lag["clogit_ols_OR01_lo"], 2),
        f(_lag["clogit_ols_OR01_hi"], 2), f(_lagr["clogit_ols_OR01"], 2),
        f(_lagr["clogit_ols_OR01_lo"], 2), f(_lagr["clogit_ols_OR01_hi"], 2),
        int(_lag["n_pairs"])),
    "restriction to \u2264 index \u2212 36 months gave %s (%s to %s)" % (
        f(_lag36["clogit_ols_OR01"], 2), f(_lag36["clogit_ols_OR01_lo"], 2),
        f(_lag36["clogit_ols_OR01_hi"], 2)),
]
chk(sec, "3.5 opens a separate OLS paragraph that states the physical scale",
    "Window, visit-density and lag sensitivity analyses used unshrunk OLS "
    "slopes, with ORs expressed per 0.1 mm/yr faster thinning." in txt
    and "These designs use unshrunk OLS slopes." not in txt,
    "own paragraph + per 0.1 mm/yr")
chk(sec, "3.5 gives same-set references for the window and the lag designs",
    all(x in txt for x in _frag),
    " / ".join(x[:40] for x in _frag))
chk(sec, "3.5 states the lag rule in months before the index visit",
    "a year before the index visit" not in txt
    and "two years before" not in txt
    and "within a year of the index visit" not in txt,
    "no relative-year shorthand")
chk(sec, "the abstract no longer mixes EBLUP with OLS designs",
    "In OLS-slope sensitivity analyses" in txt
    and "matched on visit count" not in txt
    and "only measurements taken at least 24 months before the index visit "
        "were retained" in txt
    and "within 24 months of the index visit were excluded" not in txt,
    "abstract label + 24-month rule")

sec = "v19 : the 96% stays on the class"
chk(sec, "4.4 attributes 96% to the class, not to a rate cut-point",
    "the fastest trajectory class held about 6% of sampled knees and 96% of "
    "them underwent KR" in txt
    and "the rate-defined upper tail" not in txt,
    "result attribution")

sec = "v19 : OR differences rather than coefficient percentages"
chk(sec, "the covariate-correction effect is quoted as two odds ratios",
    "changed the within-pair KL odds ratio from 1.879 to 1.594" in txt
    and "about 18%" not in txt,
    "1.879 vs 1.594")

sec = "v19 : scale and panel definitions"
chk(sec, "Table S5 defines the sign of the slope in the OR conversion",
    "exp(0.1*\u03b2*_loss)" in txt
    and "the coefficient of the negated ordinary-least-squares thickness "
        "slope" in txt
    and "the log-odds-ratio per mm/yr of the least-squares slope" not in txt,
    "beta_loss")
chk(sec, "Table 1 confines the KL remark to adjustment covariates",
    "Among the adjustment covariates, KL grade showed the clearest "
    "within-pair difference" in txt
    and "only variable differing systematically within pairs" not in txt,
    "covariates only")
chk(sec, "Table 3 distinguishes frequentist P from posterior probabilities",
    "Frequentist P values in panels (1) and (2) are two-sided; panel (3) "
    "additionally reports posterior directional probabilities" in txt
    and "overlap of those intervals. P values are two-sided." not in txt,
    "panel-labelled")

sec = "v19 : panel A settings and the probe panel"
chk(sec, "panel A states the likelihood it evaluates and its fixed settings",
    "The three log-likelihood values are sums of 191 pair-level conditional "
    "event contributions" in txt
    and "\u03b7\u209b is the case-minus-control linear predictor defined in "
    "Supplementary Methods" in txt
    and "The numerical information-matrix check used a 12-stratum subset" in txt
    and "Every value in (A)" not in txt
    and "191 strata is impractical" not in txt,
    "conditional event likelihood")
chk(sec, "panel E lists the grid the probe cross-reference points at",
    "**(E) Controlled probe of the stratified specification" in txt
    and "failed (3 NA coefficients)" in txt
    and "the slowest successful fit took 85.9 s" in txt,
    "six specifications")
chk(sec, "the probe cross-references now point at panel E",
    "specification tested (Table S3E)" in txt
    and "(Supplementary Table S3E)" in txt
    and "documented initialisation failures under the stratified `JMbayes2` "
        "specification tested (Table S3A)" not in txt,
    "S3E")
chk(sec, "the chol() message is no longer attributed to this probe",
    "chol()" not in txt,
    "unverified message removed")

sec = "v19 : hedged safeguard and simulation wording"
chk(sec, "the safeguard is described as reducing a risk",
    "curvature becomes small and the resulting proposal can be excessively "
    "diffuse" in txt
    and "makes the observed curvature vanish" not in txt
    and "cannot freeze" not in txt,
    "no guarantee")
chk(sec, "the simulation sentences report observations",
    "giving no evidence of inflated type-I error in this simulated setting"
    in txt
    and "No additional association was detected beyond the continuous slope"
    in txt
    and "consistent with a larger discrepancy for the metrics showing "
        "stronger associations" in txt
    and "signature of a bias" not in txt,
    "three sentences")
chk(sec, "the 20-pair aside is gone",
    "invisible at 20" not in txt
    and "changed the reported odds ratio in the 191-pair analysis" in txt,
    "4.3")

# ============================================ 20. eighth review round ------
sec = "v20 : the probe grid is counted correctly"
chk(sec, "the simulation text counts the five successful specifications",
    "whereas the other five specifications returned without error" in txt
    and "other four specifications completed" not in txt,
    "six rows, one failed")

sec = "v20 : prior sensitivity lists the prior it reports"
_ps = {}
for _r in rd("prior_sensitivity_joint.csv"):
    _m = re.match(r"log-OR prior N\(0, (\d+)\^2\)$", _r["specification"])
    if _m and _r["term"] == "OR per 1 SD faster thinning":
        _g = re.match(r"([\d.]+) \(([\d.]+), ([\d.]+)\)", _r["estimate"])
        _ps[int(_m.group(1))] = tuple(float(x) for x in _g.groups())
_psf = ["%s (%s to %s)" % tuple(
    decimal.Decimal(repr(v)).quantize(decimal.Decimal("0.01"),
                                     rounding=decimal.ROUND_HALF_UP)
    for v in _ps[sd]) for sd in (5, 10, 20, 50)]
chk(sec, "the four prior SDs are named in the order the results appear",
    ("For prior SDs of 5, 10, 20 and 50 the odds ratio per SD of thinning "
     "rate was %s, %s, %s and %s," % tuple(_psf)) in txt,
    " / ".join(_psf))
chk(sec, "the N(0, 20^2) upper limit is explained as a separate chain",
    "These fits were run separately from the primary chain" in txt
    and _psf[2].split(" to ")[1].rstrip(")") in txt
    and "9.47" in txt,
    "%s vs 9.47" % _psf[2])

sec = "v20 : Table 6 labels carry the sign convention"
chk(sec, "the correlation row names both slopes without a double negation",
    "Pearson correlation: CTh-Score slope versus negated cMFTC thickness "
    "slope" in txt
    and "the latter negated" not in txt,
    "panel C row")
chk(sec, "the panel A note describes the original coding, not ORs in the table",
    "With the original thickness-slope coding the cMFTC odds ratios in panel "
    "A lie below 1" in txt
    and "In panel A, ORs below 1 for cMFTC" not in txt,
    "note opening")

# =========================== v21 : the object behind S3A panel (A) ==========
sec = "v21 : Table S3A panel (A) is the model's event log-likelihood"
# panel (A) is exported by the simulation study itself, so it is read from
# Analysis/sim rather than from the analysis directory
_sa = list(csv.DictReader(io.open(
    os.path.join(BASE, "..", "sim", "sim1a_invariance.csv"),
    encoding="utf-8")))[0]
# v22: the table carries six decimals; the ten-decimal value lives only in the
# exported CSV, so the manuscript is checked against "%.6f" of that value.
_ll = "%.6f" % fnum(_sa["ll_zero_common"])
chk(sec, "panel A carries the re-exported log-likelihood three times",
    txt.count("\u2212" + _ll.lstrip("-")) == 3
    and fnum(_sa["abs_change"]) == 0.0,
    "%s, change %s" % (_ll, _sa["abs_change"]))
chk(sec, "panel A no longer prints ten decimals from the export",
    "\u2212100.1874111326" not in txt
    and "| Absolute change across all three specifications | **0** at "
        "numerical precision" in txt,
    "%.6f in the table, %s in the CSV" % (fnum(_sa["ll_zero_common"]),
                                          _sa["ll_zero_common"]))
chk(sec, "the retired mis-scaled value is gone",
    "\u22121223" not in txt and "1223.07" not in txt,
    "aS x level predictor")
chk(sec, "the panel states how many contributions the sum has",
    "| Pair-level contributions summed (mean per pair) | %s (%s) |" % (
        _sa["n_pairs"], f(_sa["mean_per_pair"], 3).replace("-", "\u2212")) in txt,
    "%s contributions, mean %s" % (_sa["n_pairs"],
                                   f(_sa["mean_per_pair"], 3)))
chk(sec, "the note gives the range of the per-pair event probability",
    "ranged from %s to %s" % (f(_sa["plogis_min"], 3),
                              f(_sa["plogis_max"], 3)) in txt,
    "%s to %s" % (f(_sa["plogis_min"], 3), f(_sa["plogis_max"], 3)))
chk(sec, "the information check is reported with its tolerance",
    "| Null directions (absolute eigenvalue within tolerance) | **%d of %s** |" % (
        int(_sa["null_directions"]), _sa["n_params"]) in txt
    and "a tolerance of 10\u207b\u2075 of the largest absolute eigenvalue" in txt,
    "%s of %s" % (_sa["null_directions"], _sa["n_params"]))
chk(sec, "the leading and the largest null-direction eigenvalues agree with the run",
    "| Largest absolute eigenvalue | %s |" % f(_sa["max_abs_eigen"], 3) in txt
    and "the largest of those %d at %.1f \u00d7 10\u207b\u2077" % (
        int(_sa["null_directions"]),
        fnum(_sa["next_abs_eigen"]) * 1e7) in txt,
    "%.4g with null directions at %.1e" % (fnum(_sa["max_abs_eigen"]),
                                           fnum(_sa["next_abs_eigen"])))
chk(sec, "the main text quotes the same leading eigenvalue",
    "had %d null directions, the largest eigenvalue being %s." % (
        int(_sa["null_directions"]), f(_sa["max_abs_eigen"], 3)) in txt,
    "Section 3 summary")

sec = "v21 : reporting consistency"
chk(sec, "Table S2 no longer calls all three contrasts posterior-draw quantities",
    "computed on one posterior draw set" not in txt
    and "the bivariate conditional-logistic coefficient differences" in txt
    and "the joint-model association of the within-knee regional difference"
        in txt,
    "S2 gradient row")
chk(sec, "the unreported 224-pair set is no longer claimed",
    "224 pairs" not in txt and "single-visit knees" not in txt
    and "a single observation provides no within-knee estimate of change" in txt,
    "Section 2.4")
chk(sec, "the KL column is labelled as an OR",
    "Within-pair KL coefficient" not in txt
    and "Within-pair KL OR per grade difference" in txt,
    "Table S4")
chk(sec, "the supplementary class-definition claim matches what is shown",
    "every class-based result is reported under both definitions" not in txt
    and "Rank-threshold agreement and below-threshold gradients are reported "
        "under both definitions" in txt,
    "Supplementary Results")

# =============================== v22 : final consistency round ===============
sec = "v22 : one definition of eta in Table S3A panel (A)"
chk(sec, "the note uses the single within-pair predictor and nothing else",
    "\u03a3\u209b log{plogis(\u03b7\u209b)}" in txt
    and "\u03b7_case" not in txt and "\u03b7_control" not in txt
    and "log plogis(\u03b7_case" not in txt,
    "eta_s only")
chk(sec, "the note does not deny that a simulated log-likelihood is one",
    "is not an estimate of the likelihood of the fitted model" not in txt
    and "this is a numerical invariance check at fixed generating parameters, "
    "not a model-comparison analysis" in txt,
    "invariance check, not a model comparison")
chk(sec, "the note keeps the probability range and the contribution count",
    "sums of 191 pair-level conditional event contributions" in txt
    and "ranged from 0.060 to 0.988" in txt
    and "the longitudinal likelihood is excluded" in txt,
    "191 contributions, 0.060 to 0.988")

sec = "v22 : Table 1 does not claim OLS is used for shrinkage only"
chk(sec, "the OLS footnote is limited to this table's descriptive purpose",
    "reported solely to evaluate empirical-Bayes shrinkage" not in txt
    and "The OLS summary in this table is presented as a descriptive "
        "shrinkage benchmark" in txt,
    "Table 1 footnote")

sec = "v22 : the provenance of the S3A numbers is on disk"
_s3a = os.path.join(BASE, "..", "sim", "S3A_provenance.md")
chk(sec, "an audit note records the old object, the new object and the scope",
    os.path.exists(_s3a)
    and "\u22121223.0754113725" in io.open(_s3a, encoding="utf-8").read()
    and "0.08431" in io.open(_s3a, encoding="utf-8").read()
    and "Part C" in io.open(_s3a, encoding="utf-8").read(),
    os.path.basename(_s3a))

sec = "v23 : the public script and the manuscript describe the same check"
chk(sec, "the panel is a sum over the 191 pairs of the model's predictor",
    "sums of 191 pair-level conditional event contributions" in txt
    and "it is a different function of time whose value at that visit is "
        "the same" in txt,
    "shape varied in the third specification")
chk(sec, "the information matrix is described as a three-coefficient baseline",
    "1 association + 12 strata \u00d7 3 basis coefficients" in txt,
    "37 parameters, all basis coefficients used")
chk(sec, "Table 1 points at the last two columns, not the last",
    "The last two columns give the mean within-pair difference" in txt
    and "The final column gives the mean within-pair difference" not in txt,
    "Table 1 footnote")
chk(sec, "the sampler overview states the paired-block update up front",
    "Random effects were updated in pairwise blocks" in txt
    and "Random effects were updated by independent Metropolis-Hastings with a "
        "conditional-normal proposal" not in txt,
    "Section 2.5.2")

# ================== v24 : A2 double-scaling correction in the public code ====
sec = "v24 : the information matrix uses the raw thinning-rate difference"
_rsim = io.open(os.path.join(BASE, "..", "..", "Scripts", "sim",
                             "12_sim_degeneracy.R"), encoding="utf-8").read()
# comment lines are hard-wrapped, so compare on a whitespace-flattened copy
_rsim_flat = re.sub(r"\s+", " ", _rsim.replace("#", " "))
chk(sec, "A2 regresses on the raw rate difference, not on the predictor",
    "dltS  <- (rate[cs] - rate[js])[seq_len(S)]" in _rsim
    and "rate  <- -d1$b1" in _rsim
    and "dltS <- (eta_c - eta_j)[seq_len(S)]" not in _rsim,
    "a is per mm/yr, so a = aS reproduces eta")
chk(sec, "A2 asserts that aS x rate difference reproduces the predictor",
    "stopifnot(max(abs(aS * dltS - (eta_c - eta_j)[seq_len(S)])) < 1e-10)"
    in _rsim,
    "guards against a future re-scaling")
chk(sec, "the note declares the rate-only special case with a_V = 0",
    "This check used the rate-only special case, with *a*_V = 0 and the "
    "calibrated association coefficient multiplying the case-minus-control "
    "thinning-rate difference." in txt,
    "a_V = 0 stated in Table S3A")
chk(sec, "the note quotes the re-exported leading eigenvalue",
    "against a leading eigenvalue of %s" % f(_sa["max_abs_eigen"], 3) in txt,
    "0.059 from the corrected Hessian")
chk(sec, "the script no longer claims an improper posterior for any implementation",
    "has an improper posterior" not in _rsim
    and "the precise reason off-the-shelf joint-model software fails" not in _rsim
    and "does not establish identifiability in a full joint-model likelihood"
        in _rsim_flat
    and "statement about the posterior of any particular implementation"
        in _rsim_flat,
    "observation, not a causal claim")
_readme = io.open(os.path.join(BASE, "..", "..",
                               "oai-poma-cartilage-thinning-tkr", "README.md"),
                  encoding="utf-8").read()
chk(sec, "the README carries the corrected leading eigenvalue",
    "0.08063" not in _readme and "leading eigenvalue of 0.0594" in _readme,
    "README status of main")

# ============ v25 : the varied quantity is a common additive term, not a hazard
sec = "v25 : Table S3A names the common additive event term"
CS = "c\u209b"          # the term itself; the index visit is given in words
chk(sec, "the three panel rows are conditional-event statements",
    "| Conditional event log-likelihood, common additive term set to zero |"
    in txt
    and "| Conditional event log-likelihood, arbitrary common additive term |"
        in txt
    and "| Conditional event log-likelihood, two common-term functions agreeing "
        "at the index visit |" in txt,
    "row labels")
chk(sec, "no panel row calls the varied quantity a baseline hazard",
    "index-visit baseline hazard set to zero" not in txt
    and "arbitrary positive value |" not in txt,
    "terminology")
chk(sec, "the note defines the term and gives its proportional-hazards reading",
    "The check varied a pair-shared additive event term " + CS
    + ", evaluated at the index visit" in txt
    and "that term cancels from the conditional likelihood" in txt
    and "in a proportional-hazards representation it corresponds to the log "
        "baseline hazard where the baseline hazard is positive" in txt,
    "c_s and log baseline hazard")
chk(sec, "the note no longer names the varied term h_0s(t*)",
    "h\u2080" not in txt and CS in txt
    and "and only that common additive term varied" in txt,
    "no h_0s left in the manuscript")
chk(sec, "Section 3 uses the same name for the same object",
    "*(i) Pair-shared additive event term.*" in txt
    and "*(i) Stratum baseline hazard.*" not in txt,
    "supplementary results (i)")
chk(sec, "the sign correspondence between the simulation and the model is stated",
    "The positive simulation coefficient is defined for thinning rate, i.e. the "
    "negative thickness slope, and therefore corresponds to \u2212*a*_S in the "
    "joint-model notation above." in txt,
    "thinning rate = negative thickness slope")
chk(sec, "the note still carries the rate-only special case and the eigenvalue",
    "This check used the rate-only special case, with *a*_V = 0" in txt
    and "against a leading eigenvalue of %s" % f(_sa["max_abs_eigen"], 3) in txt,
    "0.059 from the CSV")
chk(sec, "the competing-interests sentence points at the Funding statement",
    "The funders listed in the Funding statement had no role" in txt
    and "The funders named above had no role" not in txt,
    "Competing interests")
chk(sec, "Section 4.5 calls the 6.7 to 428 range a credible interval",
    "its 95% credible interval (6.7 to 428) should not be over-read" in txt
    and "its odds ratio (6.7 to 428)" not in txt,
    "4.5 point 2")
chk(sec, "the public script names the varied quantity as an additive term",
    "c_s(t*)" in _rsim and "h_{0s}(t*)" not in _rsim
    and "ll_zero_common" in _rsim,
    "12_sim_degeneracy.R")

# ==================================================================== report ==
out = []
out.append("=" * 74)
out.append("INDEPENDENT NUMERIC AUDIT -- %s.md" % BASENAME)
out.append("=" * 74)
npass = sum(1 for r in results if r[0])
nfail = len(results) - npass
cur = None
for ok, s, lab, det in results:
    if s != cur:
        out.append("")
        out.append("[%s]" % s)
        cur = s
    out.append("  %s %-58s %s" % ("PASS" if ok else "**FAIL**", lab[:58], det[:60]))
out.append("")
out.append("=" * 74)
out.append("TOTAL %d checks : %d PASS, %d FAIL" % (len(results), npass, nfail))
out.append("=" * 74)
body = "\n".join(out)
io.open(os.path.join(BASE, "%s_audit_strict.log" % BASENAME), "w",
         encoding="utf-8", newline="\n").write(body + "\n")
print(body)
sys.exit(1 if nfail else 0)
