# Validate released Online Resource files (schema + key ranges)
library(data.table)

or1 <- fread("supplementary/OnlineResource1_full_scenario_results.csv")
or2 <- fread("supplementary/OnlineResource2_target_population_shift.csv")
or3 <- fread("supplementary/OnlineResource3_sensitivity_results.csv")

range2 <- function(x) sprintf("%.3f–%.3f", min(x, na.rm=TRUE), max(x, na.rm=TRUE))

cat("=== Online Resource 1 (Primary) ===\n")
cat("Rows:", nrow(or1), "Cols:", ncol(or1), "\n")
if ("bias" %in% names(or1)) cat("Bias:", range2(or1$bias), "\n")
if ("cover" %in% names(or1)) cat("Coverage:", range2(or1$cover), "\n")
if ("fpr_ci" %in% names(or1)) cat("CI excl null:", range2(or1$fpr_ci), "\n")

cat("\n=== Online Resource 2 (Target population shift) ===\n")
cat("Rows:", nrow(or2), "Cols:", ncol(or2), "\n")
if ("frac_confirm_median" %in% names(or2)) cat("Confirmation fraction:", range2(or2$frac_confirm_median), "\n")
if ("confirm_rate_diff_median" %in% names(or2)) cat("Confirm rate diff:", range2(or2$confirm_rate_diff_median), "\n")
if ("max_abs_smd_pop_median" %in% names(or2)) cat("Max abs SMD (pop):", range2(or2$max_abs_smd_pop_median), "\n")
if ("n_confirm_median" %in% names(or2)) cat("Confirmed N:", range2(or2$n_confirm_median), "\n")

cat("\n=== Online Resource 3 (Sensitivity) ===\n")
cat("Rows:", nrow(or3), "Cols:", ncol(or3), "\n")
if ("bias" %in% names(or3)) cat("Bias:", range2(or3$bias), "\n")
if ("cover" %in% names(or3)) cat("Coverage:", range2(or3$cover), "\n")
if ("fpr_ci" %in% names(or3)) cat("CI excl null:", range2(or3$fpr_ci), "\n")

cat("\nDone. Compare these ranges to Tables 2–3 in the manuscript.\n")
