# Literature

## Reproduction target

**Open Systems Pharmacology CYP1A2 / CYP2C19 drug–drug–gene interaction PBPK
network** (caffeine as the CYP1A2 substrate). Open access; the whole-body PBPK
model and all parameters are published on GitHub, which is why it was chosen —
the drug-specific parameters can be transcribed from a public source rather than
recalled.

- Article (open access): PMC7764797 —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7764797/
- The associated qualified PBPK models are hosted by Open Systems Pharmacology
  on GitHub (search "OSP PBPK model library caffeine" / the paper's data
  availability statement for the exact repository).

> To obtain: open the PMC link, download the PDF into this folder as
> `caffeine_pbpk.pdf` (gitignored — copyrighted), and note the GitHub repo URL
> below once confirmed. The parameter tables + the GitHub model files are what
> this project reads.

GitHub model repo (fill in once confirmed): `________`

## Supporting / methods references

- Brown RP, et al. (1997) *Physiological parameter values for physiologically
  based pharmacokinetic models.* Toxicol Ind Health 13:407–484. — generic human
  tissue volumes and blood flows.
- ICRP Publication 89 (2002) — reference organ masses and blood flows.
- Poulin P, Theil FP (2002) *Prediction of pharmacokinetics prior to in vivo
  studies. II. Generic PBPK model.* J Pharm Sci 91:1358–1370. — tissue:plasma
  partition-coefficient prediction.
- Rodgers T, Rowland M (2006) — mechanistic tissue-composition Kp prediction.
- Baron K, et al. — *QSP and PBPK Modeling with mrgsolve: A Hands-on Tutorial.*
  CPT:PSP. — the mrgsolve PBPK implementation pattern this project follows.

## Provenance policy (same as the warfarin project)

Every drug-specific parameter carries a `status` and a `source`. Values are
`unverified` until checked against the paper or the GitHub model, then set to
`verified`. Generic physiology is cited to Brown 1997 / ICRP 89. No paper PDF is
committed (copyright); DOIs and links live here.
