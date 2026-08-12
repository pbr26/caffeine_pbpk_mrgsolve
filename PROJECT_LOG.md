# Project log — caffeine_pbpk_mrgsolve

Author: Pramod B R. A running record of decisions, dead-ends, and fixes, so the
reasoning behind the final model is auditable.

## Scope

Reproduce a published whole-body PBPK model of caffeine in `mrgsolve` (a
simulation engine, not an estimator), validated against digitized clinical
plasma profiles. Sibling to the warfarin PK/PD reproduction; same provenance
discipline (every value carries a status + source; nothing typed from memory).

## Timeline of decisions

1. **Drug + reference.** Chose caffeine (open PBPK reference model, abundant
   public plasma data). Reference: Britz et al. 2019 (CPT:PSP) + the OSP
   caffeine and Fluvoxamine-Caffeine-DDI projects.

2. **Physiology.** Transcribed organ % body weight (Brown 1997 Table 7) and
   % cardiac output (Table 23) into `data/physiology/`; absolute volumes/flows
   are computed, not pre-baked, so the source percentages stay checkable. Liver
   row carries hepatic-artery flow only; gut carries the portal share (no double
   count). Volume and flow balance enforced.

3. **Drug parameters.** Extracted from the OSP caffeine `.pksim5` (SQLite →
   streamed-ZIP XML). PK-Sim stores values in internal BASE units, which caused
   an early MW of 1.94e-7 (kg/µmol) — flagged and converted rather than shipped
   wrong. Several kinetic values kept a `confirm-units` flag.

4. **Observed data.** Decoded the two caffeine-alone profiles (Culm-Merdek 2005
   250 mg; Jeppesen 1996 200 mg) from base64 .NET float32 arrays in the DDI
   project; converted from base units (kg/L, min) to µg/mL and time-after-dose.

5. **Model.** Perfusion-limited whole-body model in `mrgsolve`; all numbers
   injected from data files, only neutral defaults in the `.cpp`.

6. **Distribution — the hard part.** First pass used unity Kp → Vd ~66 L, ~1.8×
   too high, so predictions ran ~2-fold low. Chose the faithful route: PK-Sim's
   Schmitt method. `ospsuite` installed on macOS (via `.NET 8` + rSharp;
   `DOTNET_ROOT` had to be set by hand). `createIndividual` worked, but
   `loadSimulation` is broken in the dev build and PK-Sim desktop is Windows-only.
   Recovered by (a) reading Fat/Muscle composition from the created individual
   and (b) parsing all other organs' Willmann-Schmitt composition straight out
   of `ospsuite`'s bundled `Aciclovir.pkml` (composition is drug-independent).
   Computed Kp with Poulin-Theil → **Vss 0.50 L/kg**, essentially exact.

7. **Clearance / calibration.** In-vitro-scaled CYP1A2 Vmax gave a ~1.8 h
   half-life (clearance ~2.7× too high). Britz's paper states they *fit* the
   CYP1A2 catalytic rate, so we did the same: jointly calibrated a Vmax scale and
   `ka` to both profiles. Result: AAFE ~1.3, 100% within 2-fold, t½ ≈ 3.8 h.

8. **Sensitivity.** Individual branch flows can't be perturbed alone (they must
   sum to cardiac output — doing so blew the metrics to 1e37); replaced with a
   balanced global-perfusion probe. Final analysis confirms low-extraction
   behavior (AUC flow-insensitive; driven by Vmax/fu/Km).

## Provenance status (final)

- physiology: 12/14 verified (blood volumes from ICRP, flagged)
- drug params: 9 verified / 3 confirm-units / 1 confirm-alternative / 1 todo
- partition coefficients: 11 verified (PK-Sim composition) + 1 derived (rest)
