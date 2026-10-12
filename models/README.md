# Metabolic models

The genome-scale models are not included in this repository. Download each one
from its source and save it here under the exact filename the scripts expect.

## Stoichiometric models (used by the `cmp6_*` scripts)

| Save as | Model | Where to get it |
|---|---|---|
| `yeast-GEM.xml` | Yeast8 v9.1.1 (*S. cerevisiae*) | github.com/SysBioChalmers/yeast-GEM |
| `iYali.xml` | W29-GEM (*Y. lipolytica*) | github.com/SysBioChalmers/Yarrowia_lipolytica_W29-GEM |
| `rhto.xml` | rhto-GEM v1.3.0 (*R. toruloides*) | github.com/SysBioChalmers/rhto-GEM |
| `iJN1463.xml` | iJN1463 (*P. putida* KT2440) | bigg.ucsd.edu |
| `iJO1366.xml` | iJO1366 (*E. coli*) | bigg.ucsd.edu ; BioModels MODEL1108160000 |
| `cnecator.xml` | RehMBEL1391 (*C. necator* H16) | Park et al. (2011), BMC Syst Biol 5:101, supplementary material |
| `iCW773.xml` | iCW773 (*C. glutamicum* ATCC 13032) | BioModels MODEL2110010001 |

## Enzyme-constrained models (used by `ec_combined.py`)

| Save as | Model | Where to get it |
|---|---|---|
| `ecYeastGEM_batch.xml` | ecYeast8 | github.com/SysBioChalmers/ecModels |
| `eciYali.xml` | eciYali | github.com/SysBioChalmers/ecModels |
| `eciML1515.xml` | eciML1515 | github.com/SysBioChalmers/ecModels |

Notes:
- Filenames may differ in the source repositories; rename the downloaded file to
  the "Save as" name above.
- The models are redistributed under the licences of their source repositories,
  not under this repository's MIT licence.
