# Observed caffeine plasma profiles

Digitized reference data used to validate the mrgsolve PBPK model (step 5).
Regenerate with `python3 tools/extract_observed_pksim5.py` — nothing here is
typed by hand.

| File | Study | Dose / route | Points | Cmax (ug/mL) |
|---|---|---|---|---|
| `culm_merdek_2005_250mg.csv` | Culm-Merdek et al. 2005 | 250 mg single oral | 12 | 5.18 |
| `jeppesen_1996_200mg.csv` | Jeppesen et al. 1996 (EM) | 200 mg single oral | 11 | ~4.3 |

Both are the **caffeine-alone** arms (no fluvoxamine), extracted from the OSP
`Fluvoxamine-Caffeine-DDI.pksim5` project (Britz et al. 2019). Columns:

- `time_h` — time after dose (h). Stored in PK-Sim base units (min); converted `(t - t0)/60`.
- `conc_ug_mL` — plasma caffeine (ug/mL). Stored as float32 in base units kg/L; converted `x 1e6`.
- `sd_ng_mL` — reported SD (ng/mL); `0` where the study reported none.
- `study`, `source` — provenance carried on every row.

The PK-Sim base-unit convention (mass concentration = kg/L) was confirmed by the
plausibility of the decoded curves and is what lets us treat the compound
`solubility` (0.0216 kg/L -> 21.6 mg/mL) as verified rather than confirm-units.
