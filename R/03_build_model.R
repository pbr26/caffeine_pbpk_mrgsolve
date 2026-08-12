# -----------------------------------------------------------------------------
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# 03_build_model.R -- assemble the whole-body PBPK parameter set from the
# VERIFIED data files and compile the mrgsolve model.
#
# Inputs (all provenance-tagged upstream):
#   data/physiology/human_physiology.csv  -> volumes (L) & flows (L/h)   [01]
#   data/caffeine_params.csv              -> drug disposition + CYP1A2    [02]
# Output:
#   a compiled mrgsolve model object `mod`, and tables/pbpk_param_set.csv
#   recording every value pushed into the model, with its source.
#
# Design rule (same as warfarin): the .cpp holds only neutral defaults; every
# number the model actually runs on is injected here from a sourced file, so the
# provenance chain stays unbroken.
# -----------------------------------------------------------------------------

source("R/_setup.R")
library(mrgsolve)

BODY_WEIGHT_KG    <- 70
CARDIAC_OUTPUT_LH <- 312

# --- Physiology: absolute volumes & flows from the source percentages --------
phys <- readr::read_csv(file.path(DIR_PHYS, "human_physiology.csv"),
                        show_col_types = FALSE)
phys$volume_L <- phys$pct_body_weight    / 100 * BODY_WEIGHT_KG
phys$flow_Lh  <- phys$pct_cardiac_output / 100 * CARDIAC_OUTPUT_LH
V <- setNames(phys$volume_L, phys$compartment)
Q <- setNames(phys$flow_Lh,  phys$compartment)

# --- Drug parameters ---------------------------------------------------------
cp  <- readr::read_csv(file.path(DIR_DAT, "caffeine_params.csv"),
                       show_col_types = FALSE)
gv  <- function(p) suppressWarnings(as.numeric(cp$value[cp$parameter == p]))
MW  <- gv("molecular_weight")            # g/mol
fup <- gv("fraction_unbound_plasma")     # -

# --- CYP1A2 whole-liver Vmax scaling (in vitro -> in vivo) --------------------
# Vmax_invitro [pmol/min/mg microsomal protein] x MPPGL [mg/g liver]
#   x liver mass [g] -> pmol/min; convert to mg/h with MW.
# NOTE: Vmax_invitro is still flagged "confirm-units" in caffeine_params.csv,
# and MPPGL below needs its own primary-source confirmation -> the whole
# metabolism term is PROVISIONAL until both are locked (see tables/ output).
MPPGL_mg_per_g   <- 40      # Barter et al. 2007 (mean ~39.8); CONFIRM from PDF
LIVER_DENSITY    <- 1.0     # g/mL ~ liver ~ 1
Vmax_invitro     <- gv("CYP1A2_Vmax_invitro")   # pmol/min/mg
Km_umol_L        <- gv("CYP1A2_Km")             # umol/L (unbound)

liver_mass_g <- V[["liver"]] * 1000 * LIVER_DENSITY
Vmax_pmol_min <- Vmax_invitro * MPPGL_mg_per_g * liver_mass_g   # pmol/min whole liver
Vmax_mg_h     <- Vmax_pmol_min * 1e-6 * (MW / 1000) * 60        # pmol/min -> mg/h
Km_mg_L       <- Km_umol_L * MW / 1000                          # umol/L -> mg/L

# --- Renal clearance ---------------------------------------------------------
# caffeine_params renal_clearance is currently NA (todo); keep 0 until sourced.
CLr <- ifelse(is.na(gv("renal_clearance")), 0, gv("renal_clearance"))

# --- Partition coefficients --------------------------------------------------
# Preferred: PK-Sim's own Schmitt-method Kp, extracted by 03b_partition_ospsuite.R
# into data/partition_coefficients.csv (faithful to Britz 2019). If that file is
# absent, fall back to unity Kp (a defensible first pass for hydrophilic caffeine)
# so the pipeline still runs -- but the fit will over-predict Vd ~1.8x.
Kp_tissues <- c("ad","bo","br","he","ki","li","lu","mu","sk","th","gu","re")
kp_path <- file.path(DIR_DAT, "partition_coefficients.csv")
if (file.exists(kp_path)) {
  kpf <- readr::read_csv(kp_path, show_col_types = FALSE)
  Kp  <- setNames(kpf$Kp[match(Kp_tissues, kpf$compartment)], paste0("Kp_", Kp_tissues))
  if (any(!is.finite(Kp))) stop("NA Kp in ", kp_path, " -- fix 03b organ mapping.")
  KP_SOURCE <- "PK-Sim Schmitt via ospsuite (data/partition_coefficients.csv)"
  message("Using PK-Sim Schmitt partition coefficients from ", kp_path, ".")
} else {
  Kp <- setNames(rep(1.0, length(Kp_tissues)), paste0("Kp_", Kp_tissues))
  KP_SOURCE <- "FIRST-PASS unity (run 03b_partition_ospsuite.R for Schmitt Kp)"
  message("No partition_coefficients.csv -> unity Kp fallback (Vd will be ~1.8x high).")
}

# --- Blood:plasma ------------------------------------------------------------
BP <- 1.0   # ASSUMED ~1.0 for hydrophilic caffeine; CONFIRM from OSP B:P value

# --- Absorption --------------------------------------------------------------
ka   <- 1.5   # 1/h, first-order oral absorption; PROVISIONAL (fit/confirm in step 4)
Fbio <- 1.0   # caffeine oral bioavailability ~ complete

# --- Map physiology -> model $PARAM names ------------------------------------
pars <- c(
  V_ad = V[["adipose"]], V_bo = V[["bone"]],  V_br = V[["brain"]],
  V_he = V[["heart"]],   V_ki = V[["kidney"]],V_li = V[["liver"]],
  V_lu = V[["lung"]],    V_mu = V[["muscle"]],V_sk = V[["skin"]],
  V_th = V[["thyroid"]], V_gu = V[["gut"]],   V_re = V[["rest_of_body"]],
  V_ar = V[["arterial_blood"]], V_ve = V[["venous_blood"]],

  Q_ad = Q[["adipose"]], Q_bo = Q[["bone"]],  Q_br = Q[["brain"]],
  Q_he = Q[["heart"]],   Q_ki = Q[["kidney"]],Q_ha = Q[["liver"]],
  Q_mu = Q[["muscle"]],  Q_sk = Q[["skin"]],  Q_th = Q[["thyroid"]],
  Q_gu = Q[["gut"]],     Q_re = Q[["rest_of_body"]], Q_co = CARDIAC_OUTPUT_LH,

  BP = BP, fup = fup, ka = ka, Fbio = Fbio,
  Vmax = Vmax_mg_h, Km = Km_mg_L, CLr = CLr
)
pars <- c(pars, Kp)

# --- Compile the model and apply the parameter set ---------------------------
mod <- mread("caffeine_pbpk", project = DIR_MOD) %>% param(pars)

# --- Balance guards (a broken flow/volume map is a bug, so stop) -------------
arterial_Q <- sum(Q[setdiff(names(Q), c("lung","arterial_blood","venous_blood"))],
                  na.rm = TRUE)
stopifnot(abs(arterial_Q - CARDIAC_OUTPUT_LH) / CARDIAC_OUTPUT_LH < 0.10)
stopifnot(abs(sum(V) - BODY_WEIGHT_KG) / BODY_WEIGHT_KG < 0.10)

# --- Record exactly what went into the model, with provenance ----------------
pset <- tibble::tibble(
  model_param = names(pars),
  value       = as.numeric(pars),
  source      = dplyr::case_when(
    grepl("^V_", model_param) ~ "physiology (Brown 1997 Table 7)",
    grepl("^Q_", model_param) ~ "physiology (Brown 1997 Table 23)",
    grepl("^Kp_", model_param) ~ KP_SOURCE,
    model_param == "fup"  ~ "caffeine_params.csv (Lelo 1986) - verified",
    model_param == "Vmax" ~ "CYP1A2 in vitro scaled x MPPGL - PROVISIONAL (confirm-units)",
    model_param == "Km"   ~ "caffeine_params.csv CYP1A2_Km - confirm-units",
    model_param == "BP"   ~ "ASSUMED ~1.0 - confirm from OSP",
    model_param == "ka"   ~ "PROVISIONAL - fit/confirm in simulation step",
    model_param == "CLr"  ~ "0 until renal_clearance sourced (todo)",
    TRUE ~ "assembled in 03_build_model.R"
  )
)
save_tab(pset, "pbpk_param_set")

message(sprintf("Model compiled. Whole-liver Vmax = %.1f mg/h, Km = %.3f mg/L, CLr = %.2f L/h.",
                Vmax_mg_h, Km_mg_L, CLr))
message("03_build_model.R complete. Parameter provenance -> tables/pbpk_param_set.csv")
