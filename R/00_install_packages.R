# -----------------------------------------------------------------------------
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# 00_install_packages.R -- one-time environment setup and toolchain check.
#
#   source("R/00_install_packages.R")
#
# mrgsolve compiles its models with C++, so a working toolchain is required:
# Rtools on Windows, `xcode-select --install` on macOS, build-essential on Linux.
# -----------------------------------------------------------------------------

pkgs <- c("mrgsolve", "dplyr", "tidyr", "readr", "ggplot2")

to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) {
  message("Installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install, repos = "https://cloud.r-project.org")
}
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) stop("Failed to install: ", paste(missing, collapse = ", "))

# --- Verify the compile-and-simulate pipeline --------------------------------
# A tiny one-compartment IV bolus model, written with syntax valid across
# mrgsolve versions:
#   - concentration is derived in $TABLE (not assigned inside $CAPTURE, which is
#     a name list, not an assignment block),
#   - `capture CP = ...` both computes and captures CP for the output.
suppressPackageStartupMessages(library(mrgsolve))
message("Compiling a tiny test model to verify the C++ toolchain...")

code <- "
$PARAM CL = 1, V = 20
$CMT CENT
$ODE dxdt_CENT = -(CL/V) * CENT;
$TABLE capture CP = CENT/V;
"

test <- mrgsolve::mcode("test_onecmt", code)

out <- mrgsolve::mrgsim(
  test,
  events = mrgsolve::ev(amt = 100, cmt = "CENT"),
  end = 24, delta = 1
)
od <- as.data.frame(out)
stopifnot(nrow(od) > 0,
          "CP" %in% names(od),
          all(is.finite(od$CP)),
          max(od$CP) > min(od$CP))   # decays over time

message("mrgsolve compiled and simulated a test model successfully.")
message("\nEnvironment ready. Next: source(\"R/01_physiology.R\") once the ",
        "physiology table is populated.")
