# Run the unit tests.
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
#
#   source("run_tests.R")
#
# Tests cover the provenance data files (physiology balance, drug-param
# plausibility, Kp -> Vss), and the mrgsolve model itself (compiles, conserves
# mass with elimination off, reproduces the validated fit within 2-fold). They
# do NOT run the optim calibration, so they finish quickly. The mrgsolve model
# compiles C++ on first use, so the first run takes a little longer.

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Package 'testthat' is not installed. install.packages('testthat')")
}

Sys.setenv(CAFFEINE_ROOT = normalizePath("."))

testthat::test_dir(
  "tests/testthat",
  reporter = testthat::default_reporter(),
  stop_on_failure = TRUE
)
