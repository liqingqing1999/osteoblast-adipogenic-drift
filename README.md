# Adipogenic Drift of the Osteoblast Lineage — Single-Cell Landscape Analysis

**Single-cell evidence for a continuous adipogenic drift of osteoblast-lineage cells in osteoporosis, validated against a healthy mouse reference.**

This repository contains the complete, reproducible R analysis pipeline for the study:

> *Single-Cell Landscape Reveals a Continuous Adipogenic Drift of the Osteoblast Lineage in Osteoporosis: Cross-Species Validation Against a Healthy Mouse Reference*
> (target journal: *Journal of Bone and Mineral Research*)

## Summary of the study

Using single-cell RNA-seq of human femoral heads stratified by lumbar spine T-score (osteoporosis / osteopenia / normal-range) plus a healthy mouse metaphysis reference, we show that:

1. **45.8%** of human bone-lineage cells exhibit a hybrid (osteogenic⁺ / adipogenic⁺) intermediate state — the largest single category, indicating that adipogenic drift is a *continuous* process rather than a discrete binary shift;
2. The cell-level **DriftIndex** (= UCell adipogenic-program score − UCell osteogenic-program score) decreases continuously along the osteoblast differentiation trajectory (Spearman ρ = −0.321, p = 3.31 × 10⁻⁹³);
3. Drift **peaks in the osteopenic stage** (56.0% hybrid) rather than in established osteoporosis, suggesting drift is an early driver of bone-quality deterioration;
4. Healthy mouse bone-lineage cells show only **9.6%** hybrid cells and a negative DriftIndex, providing a cross-species healthy baseline.

## Datasets (all public, from GEO)

| Dataset | Content | Use |
|---|---|---|
| **GSE169396** | Human femoral-head scRNA-seq, 4 donors (10x v3), stratified by lumbar T-score: S1 = OP (61y F, T = −3.0), S2 = osteopenia (45y F, T = −1.3), S4 = normal-range (31y M, T = +0.6); S3 (66y M, T = NA) retained for age comparisons only | Main human analysis (steps 1–3) |
| **GSE317069** | Healthy adult C57BL/6 mouse metaphysis scRNA-seq, 5 replicates (GSM6552950–54, `Metaphysis_1–5`) | Cross-species healthy reference (step 6) |

Download (example for the mouse metaphysis samples):

```bash
# Human data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169396/suppl/GSE169396_RAW.tar
# Mouse metaphysis data (5 samples × 3 files)
for i in 50 51 52 53 54; do
  base=https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6552nnn/GSM65529${i}/suppl
  rep=$((i - 49))
  wget $base/GSM65529${i}_Metaphysis_${rep}_barcodes.tsv.gz
  wget $base/GSM65529${i}_Metaphysis_${rep}_genes.tsv.gz
  wget $base/GSM65529${i}_Metaphysis_${rep}_matrix.mtx.gz
done
```

## Directory layout required to run

```
Lu-2026-mesenchymal-drift/
├── data/
│   ├── raw/
│   │   ├── GSM5201883_S1_barcodes.tsv.gz      # GSE169396 RAW.tar extracted,
│   │   ├── ...                                 # GSM<id>_<S1..S4>_*.gz files
│   │   └── GSE317069_metaphysis/
│   │       ├── GSM6552950_Metaphysis_1_*.gz    # mouse 10x files
│   │       └── ...
│   └── processed/                              # created by the scripts
├── scripts/                                    # this repository
├── results/
│   ├── figures/
│   └── tables/
```

The scripts handle all preprocessing automatically (per-sample standard-10x directories, QC, integration). Set `base_dir` at the top of each script to your local project path.

## Dependencies

- R ≥ 4.4 (tested on R 4.6.1)
- [Seurat](https://satijalab.org/seurat/) v5 (`Seurat`, `sctransform`, `glmGamPoi`)
- [harmony](https://github.com/immunogenomics/harmony)
- [UCell](https://github.com/carmonalab/UCell) (Bioconductor)
- [slingshot](https://github.com/kstreet13/slingshot) (Bioconductor)
- `dplyr`, `tidyr`, `ggplot2`, `patchwork`

## Pipeline (run in order)

Scripts are numbered by **execution order** (01–05); comments in each script header map it back to the corresponding step of the parent study protocol (protocol steps 1, 2, 2b/sensitivity, 3, 6).

| Script | Analysis step | Input | Key output |
|---|---|---|---|
| `01_qc_cluster.R` | QC + SCTransform + Harmony integration + clustering + **bone-lineage extraction** | GSE169396 raw 10x (4 samples) | `seurat_gse169396_all.rds`, `seurat_gse169396_bone.rds`; QC/UMAP/DotPlot figures |
| `02_drift_analysis.R` | **Core drift analysis**: UCell scoring, cell-level DriftIndex, hybrid-state classification, disease-gradient + age-dimension trends | `seurat_gse169396_bone.rds` | `seurat_gse169396_drift.rds`; `02_drift_summary.csv`; drift figures |
| `03_sensitivity.R` | Threshold-sweep sensitivity of the hybrid-state classification (informative thresholds ≥ 0) | `seurat_gse169396_drift.rds` | `03_sensitivity_summary.csv`, `03_sensitivity_trend.pdf` |
| `04_trajectory.R` | slingshot pseudotime; DriftIndex–pseudotime Spearman correlation; decile trend; lineage-gene dynamics | `seurat_gse169396_drift.rds` | `seurat_gse169396_traj.rds`; `04_trajectory_summary.csv`; trajectory figures |
| `05_cross_species.R` | Healthy mouse metaphysis processing + same scoring; cross-species DriftIndex & hybrid-proportion comparison | GSE317069 mouse 10x (5 samples) + `seurat_gse169396_drift.rds` | `seurat_mouse_metaphysis_bone.rds`; `05_cross_species_summary.csv`; cross-species figures |

```bash
Rscript scripts/01_qc_cluster.R
Rscript scripts/02_drift_analysis.R
Rscript scripts/03_sensitivity.R
Rscript scripts/04_trajectory.R
Rscript scripts/05_cross_species.R
```

## Sample classification (human, GSE169396)

Per the original study's Supplementary Table 1 (Aging 2021, PMID 34111027), using lumbar spine T-scores and WHO criteria:

| Sample | Age/Sex | Lumbar T | Hip T | Group |
|---|---|---|---|---|
| S1 | 61 y / F | −3.0 | −1.9 | Osteoporosis |
| S2 | 45 y / F | −1.3 | −1.2 | Osteopenia |
| S3 | 66 y / M | NA | NA | Age-dimension only (excluded from disease strata) |
| S4 | 31 y / M | +0.6 | −1.1 | Normal-range (lumbar) |

## Key results tables (pre-computed, in `results/tables/`)

- `02_drift_summary.csv` — per-sample DriftIndex (mean ± SD), hybrid/osteo/adipo proportions, MD & SASP scores
- `03_sensitivity_summary.csv` — hybrid % across 7 thresholds
- `04_trajectory_summary.csv` — mean DriftIndex/AD/OS/MD per pseudotime decile
- `05_cross_species_summary.csv` — human vs healthy-mouse comparison

## Methodological notes

- **DriftIndex** (cell level) = `UCell(AD_drift)` − `UCell(OS_identity)`; hybrid = both UCell scores > 0 (rank-based baseline).
- UCell scores are rank-normalized; **cross-species** comparisons therefore rely on hybrid proportions and relative trends, not absolute scores.
- Cells within a donor are non-independent; per-sample summaries (mean ± SD) are reported alongside cell-level statistics.

## License & citation

Data: GSE169396 (Aging 2021, PMID 34111027); GSE317069 (Nat Genet 2026, PMID 42432248). Please cite the original data papers and this repository if you reuse the code.

**License: All rights reserved.** This repository is released for academic review and reproducibility; reuse requires written permission from the corresponding author (see `LICENSE`).

---

*Analysis performed 2026-09. Contact: [your email / GitHub profile]*
