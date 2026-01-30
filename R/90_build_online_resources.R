# Build Online Resource CSVs from raw output files (perf + balance + target shift)
# Raw outputs are expected in OUT_DIR (default: output/raw)

if (!dir.exists("output/raw")) dir.create("output/raw", recursive = TRUE)
if (!dir.exists("supplementary")) dir.create("supplementary", recursive = TRUE)

library(data.table)

out_dir  <- Sys.getenv("OUT_DIR", file.path("output","raw"))
supp_dir <- "supplementary"

B <- as.integer(Sys.getenv("B","10"))
N <- as.integer(Sys.getenv("N","100"))

# Primary
perf1 <- fread(file.path(out_dir, sprintf("perf_ZERO_ext_B%d_N%d.csv", B, N)))
bal1  <- fread(file.path(out_dir, sprintf("balance_summary_ZERO_ext_B%d_N%d.csv", B, N)))

or1 <- merge(
  perf1,
  bal1[, .(scenario, tau2nd, p2nd_delta, spec,
           ess_over_n_median, w_max_median, w_max_over_med_median,
           treat_rate_median, max_smd_median)],
  by = c("scenario","tau2nd","p2nd_delta","spec"),
  all.x = TRUE
)

fwrite(or1, file.path(supp_dir, "OnlineResource1_full_scenario_results.csv"))

# Target population shift
shift_sum <- fread(file.path(out_dir, sprintf("target_shift_summary_ZERO_ext_B%d_N%d.csv", B, N)))
fwrite(shift_sum, file.path(supp_dir, "OnlineResource2_target_population_shift.csv"))

# Sensitivity
perf3 <- fread(file.path(out_dir, sprintf("perf_ZERO_ext_refitPS_trunc_B%d_N%d.csv", B, N)))
bal3  <- fread(file.path(out_dir, sprintf("balance_summary_ZERO_ext_refitPS_trunc_B%d_N%d.csv", B, N)))

or3 <- merge(
  perf3,
  bal3[, .(scenario, tau2nd, p2nd_delta, spec,
           ess_over_n_median, w_max_median, w_med_median, w_max_over_med_median,
           pct_smd_le_0.1_median, max_smd_median)],
  by = c("scenario","tau2nd","p2nd_delta","spec"),
  all.x = TRUE
)

fwrite(or3, file.path(supp_dir, "OnlineResource3_sensitivity_results.csv"))

cat("Built Online Resource CSVs in ./supplementary\n")
