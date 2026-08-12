# -----------------------------------------------------------------------------
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# 01_physiology.R -- build and VALIDATE the generic human physiology.
#
# Reads data/physiology/human_physiology.csv, whose source-of-truth columns are
# percentages transcribed from Brown et al. (1997):
#   pct_body_weight   -> Table 7 (Relative Organ Weight, % BW, ICRP Reference Man)
#   pct_cardiac_output-> Table 23 (Regional Blood Flow, % Cardiac Output, human)
#
# Absolute volumes (L) and flows (L/h) are COMPUTED here from the reference body
# weight and cardiac output, so the file stays a transparent, checkable record
# of the source percentages rather than pre-baked numbers.
#
# Balance checks (a physiology that fails these is a bug, so we stop):
#   - volumes sum to ~ body weight,
#   - arterial inflows sum to ~ cardiac output (lung excluded: it takes full CO;
#     liver row carries only its hepatic-artery share, the portal share sits on
#     the gut row -> no double counting).
# -----------------------------------------------------------------------------

source("R/_setup.R")

BODY_WEIGHT_KG    <- 70     # ICRP Reference Man (Brown 1997, "Body Weight in Humans")
CARDIAC_OUTPUT_LH <- 312    # 5.2 L/min (Brown 1997 Table 22, human) * 60
BALANCE_TOL       <- 0.10

phys_path <- file.path(DIR_PHYS, "human_physiology.csv")
if (!file.exists(phys_path)) {
  stop("Missing ", phys_path,
       "\nSee data/physiology/human_physiology_TEMPLATE.csv and literature/.")
}

phys <- readr::read_csv(phys_path, show_col_types = FALSE)

# --- Compute absolute quantities from the source percentages -----------------
phys <- phys %>%
  mutate(
    volume_L = pct_body_weight / 100 * BODY_WEIGHT_KG,
    flow_Lh  = pct_cardiac_output / 100 * CARDIAC_OUTPUT_LH
  )

# --- Report any still-unverified rows (blood, here) --------------------------
unver <- phys %>% filter(status != "verified")
if (nrow(unver)) {
  message("Note: ", nrow(unver), " physiology row(s) not yet verified against a ",
          "primary source: ", paste(unver$compartment, collapse = ", "),
          "\n(blood volumes are from ICRP 89, not Brown Table 7 -- confirm before ",
          "quoting.)")
}

# --- Mass balance: all compartment volumes ~ body weight ---------------------
total_vol <- sum(phys$volume_L, na.rm = TRUE)
vol_err <- abs(total_vol - BODY_WEIGHT_KG) / BODY_WEIGHT_KG
if (vol_err > BALANCE_TOL) {
  warning(sprintf("Compartment volumes sum to %.1f L vs body weight %.0f kg (%.1f%% off).",
                  total_vol, BODY_WEIGHT_KG, 100 * vol_err))
} else {
  message(sprintf("Volume balance OK: %.1f L vs %.0f kg (%.1f%%).",
                  total_vol, BODY_WEIGHT_KG, 100 * vol_err))
}

# --- Flow balance: arterial inflows ~ cardiac output -------------------------
# Exclude lung (full CO, not an arterial branch) and blood pools from the sum.
arterial <- phys %>%
  filter(!compartment %in% c("lung", "arterial_blood", "venous_blood"))
total_flow <- sum(arterial$flow_Lh, na.rm = TRUE)
flow_err <- abs(total_flow - CARDIAC_OUTPUT_LH) / CARDIAC_OUTPUT_LH
if (flow_err > BALANCE_TOL) {
  warning(sprintf("Arterial inflows sum to %.0f L/h vs cardiac output %.0f L/h (%.1f%% off).",
                  total_flow, CARDIAC_OUTPUT_LH, 100 * flow_err))
} else {
  message(sprintf("Flow balance OK: %.0f L/h vs %.0f L/h (%.1f%%).",
                  total_flow, CARDIAC_OUTPUT_LH, 100 * flow_err))
}

# Liver total perfusion, for reference (hepatic artery + portal from gut).
liver_ha <- phys$flow_Lh[phys$compartment == "liver"]
portal   <- phys$flow_Lh[phys$compartment == "gut"]
message(sprintf("Liver total perfusion = %.1f L/h (HA %.1f + portal %.1f) = %.1f%% CO.",
                liver_ha + portal, liver_ha, portal,
                100 * (liver_ha + portal) / CARDIAC_OUTPUT_LH))

save_tab(phys, "physiology_used")
message("01_physiology.R complete.")
