# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# Shared test fixtures. Locate the project root (set by run_tests.R, else walk up
# from the test dir), load the data files, and build the mrgsolve model once.

ROOT <- Sys.getenv("CAFFEINE_ROOT", unset = NA)
if (is.na(ROOT)) ROOT <- normalizePath(file.path("..", ".."))
setwd(ROOT)

suppressPackageStartupMessages({
  library(mrgsolve); library(readr); library(dplyr)
})

DIR_DAT <- file.path(ROOT, "data")
DIR_DIG <- file.path(DIR_DAT, "digitized")
DIR_MOD <- file.path(ROOT, "models")

phys <- readr::read_csv(file.path(DIR_DAT, "physiology", "human_physiology.csv"),
                        show_col_types = FALSE)
cp   <- readr::read_csv(file.path(DIR_DAT, "caffeine_params.csv"), show_col_types = FALSE)
kp   <- readr::read_csv(file.path(DIR_DAT, "partition_coefficients.csv"), show_col_types = FALSE)

BODY_WEIGHT_KG <- 70; CARDIAC_OUTPUT_LH <- 312

# Build the model once and cache it for the model tests (compiles C++ on first use)
.build_model <- local({
  m <- NULL
  function() {
    if (is.null(m)) m <<- mrgsolve::mread("caffeine_pbpk", project = DIR_MOD)
    m
  }
})
