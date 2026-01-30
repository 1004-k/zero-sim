# When Time-Zero Confirmation Changes the Target Population (Simulation code)

This repository contains code and outputs for the simulation study in the manuscript:

**“When Time-Zero Confirmation Changes the Target Population: Operating Characteristics of ATE and Overlap-Weighted Cox Models Under a Causal Null.”**

## Contents

- `R/` – scripts to run simulations and rebuild/validate Online Resource files
- `supplementary/` – released aggregated Online Resource CSVs (1–3)(ignored by git).
- `manuscript/` – submitted manuscript (ignored by git).
- `docs/original_scripts/` – archived development scripts (ignored by git).

## Quick validation

```r
source("R/99_validate_against_manuscript_tables.R")
```

## Re-running the simulation

Full run:

```bash
B=500 N=50000 N_CORES=3 Rscript R/run_all.R
```

Raw outputs go to `output/raw/`, logs to `output/logs/`.

## Notes

- The manuscript focuses on operating characteristics under a known causal null. 
- Additional exploratory code (e.g., non-null scenarios) is provided for completeness but is not used in the main analyses.

