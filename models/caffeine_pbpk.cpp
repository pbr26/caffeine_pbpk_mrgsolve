// ---------------------------------------------------------------------------
// Author: Pramod B R  |  caffeine_pbpk_mrgsolve
// caffeine_pbpk.cpp -- whole-body, perfusion-limited (well-stirred) PBPK model
// for caffeine, in mrgsolve.
//
// 14 tissue/blood compartments + oral depot + mass-balance sinks. Every
// physiological and drug parameter is supplied from R (see R/03_build_model.R)
// out of the VERIFIED data files; nothing numeric is hard-coded here except
// neutral defaults, so the C++ stays a pure structural description.
//
// Structure (perfusion-limited, blood-flow driven):
//   depot --ka--> gut lumen absorbed into GUT tissue
//   arterial blood -> every tissue -> venous blood -> lung -> arterial blood
//   GUT drains via the portal vein into LIVER (not straight to venous)
//   LIVER also gets the hepatic-artery share; CYP1A2 Michaelis-Menten removes drug
//   KIDNEY carries renal clearance of unbound drug
//
// Concentration bookkeeping (mg/L == ug/mL; amounts in mg):
//   C_tissue,total   = A_i / V_i
//   C_blood,leaving  = BP * (A_i/V_i) / Kp_i        (Kp = tissue:plasma)
//   C_unbound,tissue = fup * (A_i/V_i) / Kp_i       (plasma-unbound-equivalent)
// Clinical samples are venous PLASMA, so the reported output CP = C_ven / BP.
// ---------------------------------------------------------------------------

[ PROB ]
Caffeine whole-body PBPK (perfusion-limited). Concentrations in mg/L (= ug/mL).

[ PARAM ] @annotated
// --- Tissue volumes (L) ---------------------------------------------------
V_ad :  18.0 : adipose volume (L)
V_bo :  10.0 : bone volume (L)
V_br :   1.4 : brain volume (L)
V_he :   0.33: heart volume (L)
V_ki :   0.31: kidney volume (L)
V_li :   1.8 : liver volume (L)
V_lu :   0.53: lung volume (L)
V_mu :  28.0 : muscle volume (L)
V_sk :   2.6 : skin volume (L)
V_th :   0.02: thyroid volume (L)
V_gu :   1.5 : gut (splanchnic) volume (L)
V_re :   2.2 : rest-of-body volume (L)
V_ar :   1.26: arterial blood volume (L)
V_ve :   2.73: venous blood volume (L)

// --- Blood flows (L/h) ----------------------------------------------------
Q_ad :  16.2 : adipose flow (L/h)
Q_bo :  13.1 : bone flow (L/h)
Q_br :  35.6 : brain flow (L/h)
Q_he :  12.5 : heart flow (L/h)
Q_ki :  54.6 : kidney flow (L/h)
Q_ha :  14.4 : liver hepatic-artery flow (L/h)
Q_mu :  59.6 : muscle flow (L/h)
Q_sk :  18.1 : skin flow (L/h)
Q_th :   5.0 : thyroid flow (L/h)
Q_gu :  56.5 : gut/portal flow (L/h)
Q_re :  26.5 : rest-of-body flow (L/h)
Q_co : 312.0 : cardiac output (L/h)

// --- Partition coefficients (tissue:plasma, unitless) ---------------------
Kp_ad : 1.0 : adipose:plasma
Kp_bo : 1.0 : bone:plasma
Kp_br : 1.0 : brain:plasma
Kp_he : 1.0 : heart:plasma
Kp_ki : 1.0 : kidney:plasma
Kp_li : 1.0 : liver:plasma
Kp_lu : 1.0 : lung:plasma
Kp_mu : 1.0 : muscle:plasma
Kp_sk : 1.0 : skin:plasma
Kp_th : 1.0 : thyroid:plasma
Kp_gu : 1.0 : gut:plasma
Kp_re : 1.0 : rest-of-body:plasma

// --- Drug / disposition ---------------------------------------------------
BP    : 1.04 : blood:plasma ratio (unitless)
fup   : 0.70 : fraction unbound in plasma
ka    : 1.5  : first-order oral absorption rate (1/h)
Fbio  : 1.0  : oral bioavailability fraction
Vmax  : 0.0  : whole-liver max metabolic rate (mg/h)  [set in R from CYP1A2]
Km    : 2.85 : Michaelis constant, unbound (mg/L)      [14.7 umol/L * MW/1000]
CLr   : 0.0  : renal clearance of unbound drug (L/h)

[ CMT ] @annotated
DEPOT : oral gut-lumen depot (mg)
AD    : adipose (mg)
BO    : bone (mg)
BR    : brain (mg)
HE    : heart (mg)
KI    : kidney (mg)
LI    : liver (mg)
LU    : lung (mg)
MU    : muscle (mg)
SK    : skin (mg)
TH    : thyroid (mg)
GU    : gut (mg)
RE    : rest of body (mg)
AR    : arterial blood (mg)
VE    : venous blood (mg)
MET   : cumulative amount metabolized (mg)   [mass-balance sink]
URINE : cumulative amount excreted renally (mg) [mass-balance sink]

[ ODE ]
// Blood concentrations entering/leaving each compartment (mg/L)
double CA = AR / V_ar;                 // arterial blood conc
double CV = VE / V_ve;                 // venous blood conc

// Total systemic flow = the exact sum of the arterial branch flows. Using this
// (rather than the separate Q_co parameter) for the lung/arterial/venous pools
// guarantees the circulation is closed and mass is conserved even if the branch
// flows and Q_co disagree by rounding (they sum to ~Q_co by construction).
double Qsum = Q_ad + Q_bo + Q_br + Q_he + Q_ki + Q_ha
            + Q_mu + Q_sk + Q_th + Q_gu + Q_re;

// Blood conc leaving each tissue = BP * plasma-equiv leaving = BP*(A/V)/Kp
double out_ad = BP * (AD / V_ad) / Kp_ad;
double out_bo = BP * (BO / V_bo) / Kp_bo;
double out_br = BP * (BR / V_br) / Kp_br;
double out_he = BP * (HE / V_he) / Kp_he;
double out_ki = BP * (KI / V_ki) / Kp_ki;
double out_li = BP * (LI / V_li) / Kp_li;
double out_lu = BP * (LU / V_lu) / Kp_lu;
double out_mu = BP * (MU / V_mu) / Kp_mu;
double out_sk = BP * (SK / V_sk) / Kp_sk;
double out_th = BP * (TH / V_th) / Kp_th;
double out_gu = BP * (GU / V_gu) / Kp_gu;
double out_re = BP * (RE / V_re) / Kp_re;

// Hepatic CYP1A2 metabolism on unbound liver concentration
double Cu_li = fup * (LI / V_li) / Kp_li;
double Rmet  = Vmax * Cu_li / (Km + Cu_li);

// Renal clearance of unbound drug from kidney tissue
double Rren  = CLr * fup * (KI / V_ki) / Kp_ki;

// Oral absorption
double Rabs  = ka * DEPOT * Fbio;

// --- Compartment ODEs -----------------------------------------------------
dxdt_DEPOT = -ka * DEPOT;

dxdt_AD = Q_ad * (CA - out_ad);
dxdt_BO = Q_bo * (CA - out_bo);
dxdt_BR = Q_br * (CA - out_br);
dxdt_HE = Q_he * (CA - out_he);
dxdt_KI = Q_ki * (CA - out_ki) - Rren;
dxdt_MU = Q_mu * (CA - out_mu);
dxdt_SK = Q_sk * (CA - out_sk);
dxdt_TH = Q_th * (CA - out_th);
dxdt_RE = Q_re * (CA - out_re);

// Gut: arterial in, portal out (to liver), plus oral absorption
dxdt_GU = Q_gu * (CA - out_gu) + Rabs;

// Liver: hepatic artery + portal inflow, venous out, minus metabolism
dxdt_LI = Q_ha * CA + Q_gu * out_gu - (Q_ha + Q_gu) * out_li - Rmet;

// Lung: full systemic flow, venous in -> arterial out
dxdt_LU = Qsum * CV - Qsum * out_lu;

// Arterial blood: lung output distributed to tissues
dxdt_AR = Qsum * out_lu - Qsum * CA;

// Venous blood: collects every tissue outflow (gut goes via liver)
dxdt_VE = Q_ad*out_ad + Q_bo*out_bo + Q_br*out_br + Q_he*out_he
        + Q_ki*out_ki + Q_mu*out_mu + Q_sk*out_sk + Q_th*out_th
        + Q_re*out_re + (Q_ha + Q_gu)*out_li - Qsum * CV;

// Mass-balance sinks
dxdt_MET   = Rmet;
dxdt_URINE = Rren;

[ TABLE ]
// Reported prediction = venous PLASMA concentration (mg/L = ug/mL)
capture CP     = (VE / V_ve) / BP;   // venous plasma
capture CP_art = (AR / V_ar) / BP;   // arterial plasma
capture Cblood = VE / V_ve;          // venous whole blood
