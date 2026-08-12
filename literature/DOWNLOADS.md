# Downloaded resources

Compiled 2026-08-06. The candidate list is in `caffeine_pbpk_literature_search.xlsx`.

## Downloaded automatically (GitHub, public) -> `osp_models/`

These are third-party Open Systems Pharmacology repositories, cloned (not
committed to this repo — they're gitignored and reproducible with the commands
below). They carry their own OSP licence.

| Folder | What it contains | Use here |
|---|---|---|
| `Fluvoxamine-Caffeine-DDI/` | `*.pksim5` (SQLite) — **the caffeine PBPK model + digitized observed plasma data** | **Primary source.** Caffeine compound parameters and validation profiles both live here. |
| `Theophylline-Model/` | `*.pksim5` — whole-body model of theophylline (caffeine's CYP1A2-victim twin) | Structural analog / cross-check |
| `Fluvoxamine-Model/` | `*.json` (readable) + `Evaluation/workflow.R` | Example of the readable OSP JSON export format |

**Caffeine content confirmed inside `Fluvoxamine-Caffeine-DDI.pksim5`:**
- caffeine compound building block (partition coefficients, fu, B:P, clearance)
- observed plasma data: `CulmMerdek2005_250mgSD` and `Jeppesen1996_200mgSD`
  (caffeine alone + DDI arms)

Re-clone if the folder is missing:

```bash
cd literature/osp_models
git clone --depth 1 https://github.com/Open-Systems-Pharmacology/Fluvoxamine-Caffeine-DDI.git
git clone --depth 1 https://github.com/Open-Systems-Pharmacology/Theophylline-Model.git
git clone --depth 1 https://github.com/Open-Systems-Pharmacology/Fluvoxamine-Model.git
```

## NOT downloadable here — please fetch manually (paywall / reCAPTCHA / PDF binary)

I cannot pull journal or PMC PDFs (they're paywalled, behind reCAPTCHA, or the
fetcher only returns rendered text, not the PDF file). Download these yourself
and drop the PDFs into `literature/` (they're gitignored):

| ID | Reference | Link |
|---|---|---|
| C1 | Britz et al. 2019, CPT:PSP (open access) | https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.12397 |
| C4 | CYP1A2/2C19 DDGI network, Pharmaceutics 2020 (OA) | https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7764797/ |
| C5 | Ginsberg et al. 2004 (paywall) | https://www.tandfonline.com/doi/abs/10.1080/15287390490273550 |
| C6 | Simplified PBPK rat/human 2023 (OA) | https://pmc.ncbi.nlm.nih.gov/articles/PMC10014059/ |
| C7 | Baron et al., mrgsolve tutorial | https://www.researchgate.net/publication/336814947 |
| C8 | Jones & Rowland-Yeo 2013, CPT:PSP | https://ascpt.onlinelibrary.wiley.com/doi/10.1038/psp.2013.41 |

## Extraction status

- Caffeine compound parameters -> `data/caffeine_params.csv` (via
  `tools/extract_caffeine_pksim5.py`). 9 verified / 3 confirm-units / 1
  confirm-alternative / 1 todo.
- Observed plasma profiles -> `data/digitized/*.csv` (via
  `tools/extract_observed_pksim5.py`): Culm-Merdek 2005 (250 mg) and
  Jeppesen 1996 (200 mg), both caffeine-alone arms.

Next: build the mrgsolve whole-body PBPK model (`models/*.cpp`, step 3).
