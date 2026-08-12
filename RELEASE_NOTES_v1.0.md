# Release v1.0 — Whole-body caffeine PBPK in mrgsolve

Author: Pramod B R

First complete, validated release of the caffeine whole-body PBPK reproduction.

## What's in it

A 15-compartment, perfusion-limited (well-stirred) whole-body PBPK model of
caffeine in `mrgsolve`, driven end-to-end by a provenance-tagged parameter
pipeline and validated against two published clinical profiles.

- **Structure** — arterial/venous blood, lung, 12 tissues, oral depot; gut
  drains through the portal vein into liver; hepatic CYP1A2 Michaelis-Menten
  metabolism; renal clearance term; mass-balance sinks.
- **Physiology** — tissue volumes and blood flows from Brown et al. (1997),
  with volume and cardiac-output balance checks.
- **Distribution** — tissue:plasma partition coefficients from PK-Sim's own
  Willmann-Schmitt tissue composition (extracted via `ospsuite`), computed with
  the Poulin-Theil method. Predicted **Vss = 35 L ≈ 0.50 L/kg**, matching
  caffeine's established value.
- **Clearance** — CYP1A2 catalytic rate and oral `ka` calibrated to the observed
  profiles (as Britz et al. 2019 did), terminal half-life ≈ 3.8 h.
- **Validation** — **AAFE 1.25 / 1.38** for the 250 mg / 200 mg arms; **100% of
  points within 2-fold**; Cmax 5.12 vs 5.18 µg/mL observed (250 mg).
- **Sensitivity** — local one-at-a-time analysis confirms caffeine is a
  low-extraction drug (AUC insensitive to perfusion; governed by Vmax/fu/Km).

## Provenance highlights

Nothing was typed from memory. Caffeine parameters were extracted from the OSP
caffeine model, observed profiles decoded from the Fluvoxamine-Caffeine-DDI
project, and tissue composition pulled from PK-Sim (running in-process via
`ospsuite`) — each with a reproducible extractor in `tools/`.

## Known limitations / next

- The 200 mg (Jeppesen) arm runs slightly high in the absorption phase; a
  first-order depot is the simplest absorption model and could be refined.
- Terminal half-life (3.8 h) sits at the low end of caffeine's 4–6 h range.
- `ospsuite::loadSimulation` is broken in the current dev build on macOS ARM;
  Kp were computed from extracted composition rather than read from a PK-Sim
  simulation. Documented in `R/03b_partition.R`.
