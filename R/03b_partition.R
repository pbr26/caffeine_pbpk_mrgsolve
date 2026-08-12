# -----------------------------------------------------------------------------
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# 03b_partition.R -- compute caffeine tissue:plasma partition coefficients (Kp)
# from PK-Sim's own Willmann-Schmitt tissue composition, using the Poulin &
# Theil (2002) homogeneous-tissue model.
#
# WHY this is faithful to Britz et al. 2019:
#   Britz's caffeine model uses PK-Sim's Schmitt/Willmann partitioning. We could
#   not have PK-Sim recompute Kp on macOS (ospsuite::loadSimulation is broken in
#   the current dev build, and PK-Sim desktop is Windows-only). Instead we took
#   PK-Sim's exact Willmann-Schmitt (-WS) TISSUE COMPOSITION -- the same numbers
#   PK-Sim feeds its partition calculation -- and applied the standard, published
#   Poulin-Theil equation. Composition provenance (data/tissue_composition_ws.csv):
#     * Fat & Muscle  -> ospsuite::createIndividual (European ICRP 2002)
#     * all other organs -> PK-Sim's bundled Aciclovir.pkml (ospsuite extdata),
#       parsed directly (composition is drug-independent physiology).
#   For a neutral, hydrophilic drug like caffeine (pKa 0.8 base, logP ~ -0.07)
#   Poulin-Theil and Schmitt converge; the cross-check below is Vss ~ 0.5 L/kg.
#
# Output: data/partition_coefficients.csv  (consumed by 03_build_model.R)
# -----------------------------------------------------------------------------

source("R/_setup.R")

# --- Drug inputs (from the verified compound table) --------------------------
cp    <- readr::read_csv(file.path(DIR_DAT, "caffeine_params.csv"), show_col_types = FALSE)
gv    <- function(p) suppressWarnings(as.numeric(cp$value[cp$parameter == p]))
logP  <- gv("lipophilicity_logMA")            # ~ -0.07
fu_p  <- gv("fraction_unbound_plasma")        # 0.70
P     <- 10 ^ logP

# --- Plasma reference composition (Poulin & Theil 2002, human) ---------------
Vw_p <- 0.945; Vnl_p <- 0.0023; Vph_p <- 0.0013

# --- Tissue composition (PK-Sim Willmann-Schmitt) ----------------------------
comp <- readr::read_csv(file.path(DIR_DAT, "tissue_composition_ws.csv"),
                        show_col_types = FALSE)

# --- Poulin & Theil (2002) partition, neutral drug ---------------------------
# Kpu = [P*Vnl + (0.3P+0.7)*Vph + Vw] / [same for plasma];  Kp = Kpu * fu_p
# (fu_tissue ~ 1 for low-binding caffeine -> only plasma binding corrects).
kpu <- function(Vw, Vnl, Vph) {
  (P * Vnl + (0.3 * P + 0.7) * Vph + Vw) /
  (P * Vnl_p + (0.3 * P + 0.7) * Vph_p + Vw_p)
}
comp$Kp <- round(kpu(comp$Vw, comp$Vnl, comp$Vph) * fu_p, 4)

# rest_of_body: mean of the well-perfused soft tissues
soft <- comp$Kp[comp$compartment %in% c("br","he","ki","li","lu","mu","sk","gu")]
kp <- dplyr::bind_rows(
  comp[, c("compartment","pksim_organ","Kp","status","source")],
  tibble::tibble(compartment = "re", pksim_organ = "(soft-tissue mean)",
                 Kp = round(mean(soft), 4), status = "derived",
                 source = "mean of well-perfused soft-tissue Kp (rest-of-body closure)")
)
readr::write_csv(kp, file.path(DIR_DAT, "partition_coefficients.csv"))
message("Wrote data/partition_coefficients.csv (", nrow(kp), " compartments).")

# --- Cross-check: predicted steady-state volume of distribution --------------
BODY <- 70
pct  <- c(ad=21.42,bo=14.29,br=2.0,he=0.47,ki=0.44,li=2.57,lu=0.76,
          mu=40.0,sk=3.71,th=0.03,gu=2.13,re=3.20)
V    <- pct / 100 * BODY
Kpv  <- setNames(kp$Kp, kp$compartment)[names(V)]
Vpl  <- (1.8 + 3.9) / 100 * BODY * (1 - 0.47)     # plasma volume (blood x (1-Hct))
Vss  <- Vpl + sum(V * Kpv)
message(sprintf("Cross-check: predicted Vss = %.1f L (%.2f L/kg); caffeine target ~0.5 L/kg.",
                Vss, Vss / BODY))
if (abs(Vss / BODY - 0.5) > 0.15)
  warning("Predicted Vd is >0.15 L/kg from caffeine's ~0.5 L/kg -- inspect Kp.")
message("03b_partition.R complete. Re-run 03_build_model.R then 04_simulate.R.")
