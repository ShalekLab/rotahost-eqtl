# RotaHost sc-eQTL + sc-pQTL Pipeline

WDL workflow for the RotaHost single-cell eQTL and protein QTL analysis.

Implements the sc-eQTLGen WG3 approach ([Kaptijn et al. 2026](https://doi.org/10.64898/2026.01.20.700519)) for cell-type-specific genetic regulatory mapping.

## Analyses

| Mode | Script | Description |
|------|--------|-------------|
| `global` | run_QTL_analysis.py | Baseline cis-eQTL, all donors |
| `longitudinal` | run_QTL_analysis.py | Longitudinal donor-visit pairs, timepoint covariate |
| `iga_interaction` | run_interaction_QTL_analysis.py | SNP × IgA seroconversion interaction |
| `temporal_igay` | run_interaction_QTL_analysis.py | SNP × timepoint in seroconverters |
| `temporal_igan` | run_interaction_QTL_analysis.py | SNP × timepoint in non-seroconverters |
| `adt` | run_QTL_analysis.py | Protein QTL using TotalVI-denoised CITE-seq ADT |

## Docker Image

`ghcr.io/shaleklab/rotahost-eqtl:latest`

## Terra / Dockstore

Import this workflow into Terra from Dockstore:
`https://dockstore.org/workflows/github.com/ShalekLab/rotahost-eqtl/rotahost_eqtl`

## Quick Start (Toy Test)

```
cell_types: ["B"]
chromosomes: ["22"]
n_permutations: 10
analysis_mode: "global"
```
Expected: ~35-40 min, ~$0.03

## Full Run

```
cell_types: ["B","CD4_T","CD8_T","NK","Mono","DC","gdT"]
chromosomes: ["1"..."22"]
n_permutations: 1000
analysis_mode: "global"  # run separately for each of 6 modes
```
Expected: ~2-3 hours per mode, ~$7-10 per mode on preemptible VMs

## Citation

If you use this pipeline, please cite:
- Kaptijn et al. 2026 (sc-eQTLGen consortium)
- Cuomo et al. 2021 (limix_qtl)
- Guzman, Triana et al. (RotaHost, in preparation)
