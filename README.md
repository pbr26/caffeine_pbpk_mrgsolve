# Whole-body PBPK model of caffeine in `mrgsolve`

**Author: Pramod B R**

An independent reproduction of a published whole-body physiologically based
pharmacokinetic (PBPK) model of caffeine, implemented in
[`mrgsolve`](https://mrgsolve.org), and validated by overlaying the simulated
plasma concentration–time profile on digitized clinical data.

> **Status: v1.0 — complete and validated.** AAFE ≈ 1.3 with 100% of points
> within 2-fold across two clinical studies; predicted Vss ≈ 0.50 L/kg matches
> caffeine's established value. Same provenance discipline as the sibling
> [warfarin PK/PD reproduction](../warfarin_pkpd_reproduction): every parameter
> carries a status and a source, and nothing is typed from memory.

![Calibrated fit](figures/calibrated_fit.png)

---

## What "reproduction" means here

`mrgsolve` is an ODE **simulation** engine, not an estimation engine. A PBPK
reproduction therefore means: implement the published model *structure*, enter
the *parameters* from cited sources, *simulate* a dose, and *validate* the
predicted curve against observed data.

1. **Structure** — perfusion-limited (well-stirred) whole body: arterial and
   venous blood, lung, 12 tissues, an oral depot, and mass-balance sinks. The
   gut drains through the portal vein into the liver; the liver also receives
   its hepatic-artery share and clears caffeine by CYP1A2 Michaelis-Menten
   metabolism; the kidney carries a renal-clearance term.
2. **Parameters** — generic human physiology (Brown 1997), caffeine physchem and
   CYP1A2 kinetics (OSP caffeine model), and tissue:plasma partition coefficients
   computed from PK-Sim's own Willmann-Schmitt tissue composition.
3. **Simulate** — 200 mg and 250 mg single oral doses.
4. **Validate** — overlay on two digitized profiles; report AAFE, % within
   2-fold, Cmax, and terminal half-life.

## Headline results

| Metric | 250 mg (Culm-Merdek 2005) | 200 mg (Jeppesen 1996) |
|---|---|---|
| AAFE | 1.25 | 1.38 |
| % within 2-fold | 100% | 100% |
| Cmax obs → pred (µg/mL) | 5.18 → 5.12 | 4.27 → 3.98 |

Predicted **Vss = 35 L (0.50 L/kg)**; terminal **half-life ≈ 3.8 h**; caffeine
confirmed as a **low-extraction drug** (AUC insensitive to blood flow).

## Pipeline

```
run_all.R                        # reproduce everything end to end
R/
├── _setup.R                     # shared paths + helpers (fold_error, aafe)
├── 00_install_packages.R        # packages + C++ toolchain check
├── 01_physiology.R              # generic human physiology + balance checks
├── 02_drug_params.R             # caffeine compound parameters (provenance-tagged)
├── 03b_partition.R              # tissue:plasma Kp (PK-Sim composition, Poulin-Theil)
├── 03_build_model.R             # assemble + compile the mrgsolve model
├── 04_simulate.R                # simulate 200/250 mg, overlay observed
├── 05_calibrate_validate.R      # calibrate CYP1A2 rate + ka; validate
└── 06_sensitivity.R             # local sensitivity (Cmax / AUC / half-life)
models/caffeine_pbpk.cpp         # the mrgsolve model (structure only)
data/
├── physiology/human_physiology.csv       # Brown 1997 (source %s)
├── caffeine_params.csv                    # drug params, each flagged + sourced
├── tissue_composition_ws.csv              # PK-Sim Willmann-Schmitt composition
├── partition_coefficients.csv             # computed Kp (Vss cross-checked)
└── digitized/                             # observed plasma profiles
tools/                           # reproducible extractors (pksim5, observed, composition)
tests/testthat/                  # physiology balance, params, Kp→Vss, model mass balance
```

Run it:

```r
source("run_all.R")     # full pipeline (compiles C++ on first build)
source("run_tests.R")   # unit tests
quarto::quarto_render("analysis.qmd")   # the report -> reports/
```

## How the partition coefficients were obtained (the interesting part)

The faithful route is PK-Sim's Schmitt method, which Britz's caffeine model uses.
PK-Sim desktop is Windows-only and `ospsuite::loadSimulation` is broken in the
current macOS dev build, so instead of having PK-Sim recompute Kp, this project
extracts **PK-Sim's own Willmann-Schmitt tissue composition** — Fat and Muscle
from `ospsuite::createIndividual` (which runs PK-Sim in-process), and the other
organs parsed directly from `ospsuite`'s bundled `Aciclovir.pkml` (tissue
composition is drug-independent physiology). The published Poulin-Theil equation
is then applied to those fractions. For a neutral, hydrophilic drug like caffeine
(pKa 0.8 base, logP ≈ −0.07) Poulin-Theil and Schmitt converge, and the result is
cross-checked by Vss ≈ 0.5 L/kg. Every step is reproducible via `tools/` and
`R/03b_partition.R`.

## Provenance

- **Physiology:** Brown et al. (1997), *Toxicol Ind Health* 13:407–484; ICRP 89.
- **Drug parameters & tissue composition:** OSP caffeine model / PK-Sim
  (Britz et al. 2019, DOI 10.1002/psp4.12397), extracted via `ospsuite`.
- **Observed profiles:** Culm-Merdek et al. (2005); Jeppesen et al. (1996), via
  the OSP Fluvoxamine-Caffeine-DDI project.

Third-party OSP models and paper PDFs are **not** committed (see
`literature/DOWNLOADS.md`); only provenance-tagged extracted values are included.

## License

MIT (see `LICENSE`). Reference data belong to their original sources; cite them.
