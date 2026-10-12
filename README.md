# Comparative genome-scale metabolic modelling of PHB production from sugarcane molasses

Code and derived data for the manuscript:

> **Comparative genome-scale metabolic modelling identifies optimal microbial hosts for PHB production from sugarcane molasses.** Muhammad Fauzi Arif. *(under review).*

This repository reproduces the in-silico comparison of seven microbial hosts
(*Saccharomyces cerevisiae*, *Yarrowia lipolytica*, *Rhodotorula toruloides*,
*Pseudomonas putida*, *Escherichia coli*, *Cupriavidus necator*, *Corynebacterium glutamicum*)
for poly(3-hydroxybutyrate) (PHB) and PHBV production, using flux balance analysis
with an identical heterologous *phaCAB* pathway inserted into every model.

## 1. Software environment

The analysis was performed with:

| Software | Version |
|---|---|
| Python | 3.11.15 |
| COBRApy | 0.32.1 |
| GLPK (solver) | 5.0 |
| NumPy | 2.4.4 |
| SciPy | 1.17.1 |
| Matplotlib | 3.10.9 |
| RDKit | (for the chemical-structure panels in the pathway figure only) |

Install the dependencies into a fresh environment:

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```

GLPK is used through COBRApy's bundled `swiglpk`; no separate solver licence is required.

## 2. Metabolic models (obtain separately)

The seven genome-scale models and the three enzyme-constrained models are **not
redistributed here** because they are maintained, and licensed, by their original
authors. Download each model from the source listed below and place it in the
`models/` directory under the exact filename given, which is the filename the
scripts expect.

| Host (strain) | Model | Expected filename | Source / accession |
|---|---|---|---|
| *S. cerevisiae* | Yeast8 v9.1.1 | `models/yeast-GEM.xml` | SysBioChalmers/yeast-GEM |
| *Y. lipolytica* | W29-GEM | `models/iYali.xml` | SysBioChalmers/Yarrowia_lipolytica_W29-GEM (Guo et al., 2022) |
| *R. toruloides* | rhto-GEM v1.3.0 | `models/rhto.xml` | SysBioChalmers/rhto-GEM |
| *P. putida* KT2440 | iJN1463 | `models/iJN1463.xml` | BiGG Models |
| *E. coli* | iJO1366 | `models/iJO1366.xml` | BiGG Models; BioModels MODEL1108160000 |
| *C. necator* H16 | RehMBEL1391 | `models/cnecator.xml` | Park et al. (2011), supplementary material |
| *C. glutamicum* ATCC 13032 | iCW773 | `models/iCW773.xml` | BioModels MODEL2110010001 |
| *S. cerevisiae* (ec) | ecYeast8 | `models/ecYeastGEM_batch.xml` | SysBioChalmers/ecModels |
| *Y. lipolytica* (ec) | eciYali | `models/eciYali.xml` | SysBioChalmers/ecModels |
| *E. coli* (ec) | eciML1515 | `models/eciML1515.xml` | SysBioChalmers/ecModels |

See `models/README.md` for direct repository notes.

## 3. How to reproduce

Run from the repository root (so that `import cmp6` and the relative `models/`
paths resolve). `cmp6.py` is the shared core: it loads each model, inserts the
identical cytosolic *phaCAB* pathway, and defines the culture media.

Pipeline order (each step writes a JSON result and the corresponding figure to
the working directory):

```bash
# 1. core yields, substrate scope, cofactor swap  -> cmp6_fast.json
python cmp6_fast.py
# 2. growth vs production envelope                 -> cmp6_env.json
python cmp6_env.py
# 3. PHBV copolymer capacity                       -> cmp6_phbv.json
python cmp6_phbv.py
# 4. FSEOF engineering targets                     -> cmp6_fseof.json
python cmp6_fseof.py
# 5. CCE + ranking-robustness (Spearman)          -> cmp6_sens.json
python cmp6_sens.py
# 6. assemble the multi-metric table              -> cmp6_metrics.json
python cmp6_metrics.py
# 7. fingerprint, grouping, PCA/clustering, matrix -> cmp6_decision.json
python cmp6_synth.py
# enzyme-constrained check (3 hosts, standalone)   -> ec_combined.json
python ec_combined.py
# Methods schematics (need RDKit + the struct_*.png assets)
python fig_pathway_sci.py
python fig_fseof_sci.py
```

Steps 6 and 7 read the JSON produced by the earlier steps. To regenerate the
synthesis figures directly from the deposited snapshot without rerunning the
model steps, copy the provided results first:

```bash
cp results/*.json .
python cmp6_metrics.py   # optional; results/cmp6_metrics.json is already provided
python cmp6_synth.py
```

## 4. Repository contents

```
cmp6.py            core: model loading, phaCAB insertion, media definitions
cmp6_fast.py       theoretical yield by sugar, substrate scope, phaB cofactor swap
cmp6_env.py        growth vs PHB production envelopes
cmp6_phbv.py       PHBV copolymer capacity (C5 branch + propionate co-feed)
cmp6_fseof.py      FSEOF: acetyl-CoA / NADPH / phaA module slopes
cmp6_sens.py       carbon conversion efficiency + Spearman ranking robustness
cmp6_metrics.py    assembles the per-host multi-metric table
cmp6_synth.py      fingerprint heatmap, functional grouping, PCA + Ward clustering, decision matrix
ec_combined.py     enzyme-constrained analysis (ecYeast8, eciYali, eciML1515)
fig_pathway_sci.py Methods Figure 1 (phaCAB + PHBV pathway)
fig_fseof_sci.py   Methods Figure 2 (FSEOF schematic)
struct_*.png       RDKit chemical structures used by fig_pathway_sci.py
results/*.json     deposited snapshot of the derived data behind the manuscript tables/figures
models/            place the downloaded GEMs here (see Section 2)
```

> **Note on script naming.** The `cmp6` prefix is the pipeline's version tag (the sixth iteration of the comparison script), not a host count. Every script analyses the full set of **seven hosts** defined by `HOSTS` in `cmp6.py`.

### Mapping of results to the manuscript

| Result file | Manuscript element |
|---|---|
| `cmp6_fast.json` | Fig. 3 / Table 2 (yield by sugar), Fig. 4 (substrate scope), Fig. 6 (cofactor) |
| `cmp6_sens.json` | Fig. 5 (CCE and ranking robustness) |
| `cmp6_env.json` | Fig. 7 (growth-production envelopes) |
| `cmp6_fseof.json` | Fig. 8 / Table 3 (FSEOF slopes) |
| `cmp6_phbv.json` | Fig. 9 (PHBV capacity) |
| `ec_combined.json` | Fig. 11 / Table 4 (enzyme-constrained analysis) |
| `cmp6_metrics.json` | Fig. 12, Table 5 (multi-metric fingerprint) |
| `cmp6_decision.json` | Fig. 13–15, Table 6 (grouping, clustering, decision matrix) |

## 5. Licence

The code in this repository is released under the MIT Licence (see `LICENSE`).
The third-party genome-scale models are governed by the licences of their
respective source repositories and are not covered by this licence.

## 6. How to cite

If you use this code or data, please cite the manuscript above. A machine-readable
citation is provided in `CITATION.cff`.

## 7. Contact

Muhammad Fauzi Arif — ORCID: https://orcid.org/0000-0001-9964-1294
Repository: https://github.com/mfauziarif94/phb-molasses-gems
