# -----------------------------------------------------------------------------
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# 02_drug_params.R -- load and sanity-check the caffeine compound parameters.
#
# Source of truth: data/caffeine_params.csv, extracted from the OSP "caffeine
# template model" (literature/osp_models/Example_Caffeine/Caffeine.pksim5) by
# tools/extract_caffeine_pksim5.py. Reference: Britz et al. 2019 (ref 15 = the
# OSP Example_Caffeine repository). Every value carries a status + source, as in
# the warfarin project. Nothing is typed from memory.
#
# status values:
#   verified          - value + unit confirmed from the model XML
#   confirm-units     - value read, but PK-Sim stores base units; the display
#                       conversion needs an ospsuite-R cross-check
#   confirm-alternative - parameter has >1 stored alternative; active one TBD
#   todo              - present in the model but not yet extracted
# -----------------------------------------------------------------------------

source("R/_setup.R")

path <- file.path(DIR_DAT, "caffeine_params.csv")
if (!file.exists(path)) {
  stop("Missing ", path, "\nRun tools/extract_caffeine_pksim5.py first ",
       "(needs literature/osp_models/Example_Caffeine).")
}
cp <- readr::read_csv(path, show_col_types = FALSE)

# --- Report provenance status ------------------------------------------------
tab <- table(cp$status)
message("Caffeine parameters loaded: ", nrow(cp), " rows -> ",
        paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = ", "))

needs <- cp$parameter[cp$status %in% c("confirm-units", "confirm-alternative", "todo")]
if (length(needs)) {
  message("Still to confirm (ospsuite-R / OSP evaluation report): ",
          paste(needs, collapse = ", "))
}

# --- Physical-plausibility guards on the values we will actually use ----------
val <- function(p) suppressWarnings(as.numeric(cp$value[cp$parameter == p]))
stopifnot(
  abs(val("molecular_weight") - 194.19) < 1,        # caffeine MW
  val("fraction_unbound_plasma") > 0 && val("fraction_unbound_plasma") <= 1,
  is.finite(val("lipophilicity_logMA")),
  val("CYP1A2_Km") > 0,
  val("CYP1A2_Vmax_invitro") > 0
)
message("Plausibility checks passed (MW, fu in (0,1], logP finite, Km>0, Vmax>0).")

save_tab(cp, "caffeine_params_used")
message("02_drug_params.R complete.")
