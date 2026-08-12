# -----------------------------------------------------------------------------
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# run_all.R -- reproduce the whole caffeine PBPK pipeline end to end.
#
#   source("run_all.R")
#
# Runs from the project root. Each numbered script sources _setup.R and writes
# its outputs to tables/ and figures/. The mrgsolve model compiles C++ on first
# build (needs a working compiler toolchain -- see R/00_install_packages.R).
#
# 03b_partition.R needs data/tissue_composition_ws.csv (already committed) and,
# to REGENERATE that file, the ospsuite route documented in R/03b_partition.R;
# the committed CSV lets the pipeline run without ospsuite.
# -----------------------------------------------------------------------------

steps <- c(
  "R/01_physiology.R",       # generic human physiology + balance checks
  "R/02_drug_params.R",      # caffeine compound parameters (provenance-tagged)
  "R/03b_partition.R",       # tissue:plasma Kp from PK-Sim composition (Poulin-Theil)
  "R/03_build_model.R",      # assemble + compile the mrgsolve PBPK model
  "R/04_simulate.R",         # simulate 200/250 mg oral, overlay observed
  "R/05_calibrate_validate.R", # calibrate CYP1A2 rate + ka; validate (AAFE, 2-fold)
  "R/06_sensitivity.R"       # local sensitivity (Cmax/AUC/half-life)
)

t0 <- Sys.time()
for (s in steps) {
  message("\n===== ", s, " =====")
  source(s)
}
message(sprintf("\nAll steps complete in %.1f s.", as.numeric(Sys.time() - t0, units = "secs")))
