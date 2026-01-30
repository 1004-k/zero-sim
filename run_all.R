# Run all analyses (primary + sensitivity + non-null) and then rebuild/validate Online Resources
# Full run:
#   B=500 N=50000 N_CORES=3 Rscript R/run_all.R

source("R/01_primary_null_simulation.R")
source("R/02_sensitivity_analysis.R")
source("R/03_nonnull_extension.R")
source("R/90_build_online_resources.R")
source("R/99_validate_against_manuscript_tables.R")
