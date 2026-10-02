# ============================================================================
# 46 : export the fitted-object read-outs that the CSV summaries do not carry.
#
# WHY
#   `jm_multimetric.csv` stores the point estimate for the within-pair
#   adjusted joint model (JM_Aadj_ORSD / JM_Aadj_OR01) but NOT their credible
#   intervals, because 10_jm_multimetric.R never wrote those two columns.  The
#   manuscript quotes them (Table 2, panel B), and 45_verify_v5_numbers.R
#   recomputes them from the fitted object -- so the numbers are correct, but
#   the CSV trail is incomplete and an automated audit cannot close the loop.
#
#   This script writes the missing intervals (and a few neighbouring read-outs
#   that exist only inside the RDS) to `audit_values.csv`, so the manuscript
#   can be checked against a stored artefact rather than against a log.
#
#   -> audit_values.csv
# ============================================================================
BASE <- "D:/BaiduSyncdisk/OAI/Analysis/poma_pilot"

q3 <- function(x) sprintf("%.3f (%.3f-%.3f)",
                          median(x), quantile(x, .025), quantile(x, .975))

f <- readRDS(file.path(BASE, "jmfit_cMFTC_ThCtAB_aMe.rds"))

rows <- list()
add <- function(key, value) rows[[length(rows) + 1L]] <<- data.frame(key = key, value = value)

for (nm in c("A", "Aadj")) {
  d <- f[[nm]]$draws
  add(sprintf("%s_ORSD", nm),  q3(exp(-d$aS * d$SDslo)))
  add(sprintf("%s_OR01",  nm),  q3(exp(-0.1 * d$aS)))
  add(sprintf("%s_ORSDlevel", nm), q3(exp(-d$aV * d$SDval)))
  add(sprintf("%s_aS_median", nm), sprintf("%.4f", median(d$aS)))
  add(sprintf("%s_aV_median", nm), sprintf("%.4f", median(d$aV)))
}

## gradient of the rate association, from the same posterior
if (!is.null(f$A$draws$g_med_lat)) {
  d <- f$A$draws
  add("cMFTC_vs_cLF_ratio",   q3(exp(d$g_med_lat)))
  add("cMFTC_vs_cLFTC_ratio", q3(exp(d$g_med_lat2)))
}

out <- do.call(rbind, rows)
write.csv(out, file.path(BASE, "audit_values.csv"), row.names = FALSE)

cat("wrote audit_values.csv with", nrow(out), "rows\n")
print(out, row.names = FALSE)
