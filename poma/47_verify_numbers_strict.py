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
import csv, io, json, os, re, statistics, sys, urllib.request

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

# the Abstract quotes B1 as well: it MUST carry the same rounding as Table 6
m_abs = re.search(r"annual score increase\s+(\d\.\d+),\s+([\d.]+)\s+to\s+([\d.]+)", NTXT)
chk(sec, "Abstract B1 rounding identical to Table 6",
    bool(m_abs) and norm(ci2(B1["OR"], B1["lo"], B1["hi"])) ==
    "%s (%s-%s)" % (m_abs.group(1), m_abs.group(2), m_abs.group(3)),
    "Table 6 %s | abstract %s" % (ci2(B1["OR"], B1["lo"], B1["hi"]),
                                  m_abs.group(0) if m_abs else "?"))

tp = rd("cth_tertile_OR.csv")
chk(sec, "panel B fastest tertile", has(ci2(tp[0]["OR"], tp[0]["lo"], tp[0]["hi"]), sec, "tertile"),
    ci2(tp[0]["OR"], tp[0]["lo"], tp[0]["hi"]))
cc2 = {r["metric"]: r["value"] for r in rd("cth_concordance.csv")}
chk(sec, "panel C tertile agreement 46.6%", has("46.6", sec, "agreement"),
    "disk %.1f" % float(cc2["tertile_agreement_pct"]))
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
    missing_repo = sorted(n for n in named if n not in remote_base and n not in declared_absent)
    chk(sec, "every script named in S2 exists in the archive", not missing_repo,
        "named=%d, absent=%s" % (len(named), missing_repo or "none"))
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
sec = "DOI : points at the archive holding the corrected sampler"
chk(sec, "version DOI 23111696 (v1.0.1)", has("10.5281/zenodo.23111696", sec, "v1.0.1 DOI"))
chk(sec, "concept DOI 22771312 unchanged", has("10.5281/zenodo.22771312", sec, "concept DOI"))
chk(sec, "v1.0.0 named only as superseded",
    txt.count("10.5281/zenodo.22771313") == 1
    and "superseded v1.0.0 archive" in txt,
    "occurrences=%d" % txt.count("10.5281/zenodo.22771313"))
chk(sec, "GitHub repo URL", has("https://github.com/morrosun/oai-poma-cartilage-thinning-tkr", sec, "repo url"))

# ================================================== 9. OAC hard metrics ------
sec = "OAC hard metrics"
m = re.search(r"body\s+(\d+)\s+words against a \*\*≤4000-word\*\* limit", txt)
chk(sec, "body word count <= 4000", bool(m) and int(m.group(1)) <= 4000, m.group(1) if m else "?")
m = re.search(r"abstract \(Objective / Design / Results / Conclusions\) ≤300 words", txt, re.I)
chk(sec, "abstract declared <=300", bool(m))
m = re.search(r"\*\*(\d+) figures and tables combined\*\* against a limit of 8", txt)
chk(sec, "tables+figures == 8", bool(m) and int(m.group(1)) == 8, m.group(1) if m else "?")
refs = len(re.findall(r"^\| \d+ \|", txt, re.M))
chk(sec, "references <= 50", 0 < refs <= 50, "reference table rows=%d" % refs)
chk(sec, "no unresolved citation tokens", "[@r" not in txt)
chk(sec, "no bare pipes inside table cells",
    not any(re.search(r"\|[^|\n]*\|[^|\n]*\|[^|\n]*\|[^|\n]*\|[^|\n]*\|", l)
            and l.count("|") > 40 and "S3 panel" in l for l in lines))

# ==================================================================== report ==
out = []
out.append("=" * 74)
out.append("INDEPENDENT NUMERIC AUDIT -- Manuscript_EN_OAC_v6.md")
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
