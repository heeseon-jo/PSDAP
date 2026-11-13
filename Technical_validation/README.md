# Technical validation

R scripts that reproduce the technical validation of perturbagen similarity and gene-level agreement reported in the PSDAP Data Descriptor.

## Input files

| File | Description |
|---|---|
| `Perturbagen_similarity/Similarity_*.csv` | Distributed perturbagen similarity files (GO-BP scores) |
| `Differential_expression_signatures/` | Distributed DE signature files (`DEsignature_SIG*.csv`) |
| `signature_info.csv` | Distributed signature metadata |
| Shared raw-data pairs (csv) | Signature pairs generated from the same raw data (TCM Atlas vs KORE-Map and GSE85871), provided as a supplementary table. Save as csv with columns `cvc1` and `cvc2` (Signature_id of each signature). |

Set the paths in the `Settings` block at the top of each script.

## Scripts

Run in this order. Each script reads the outputs of the previous ones from `output/`.

| Script | Analysis | Output |
|---|---|---|
| `01_top100_lists.R` | Top-100 most similar profiles for each query profile | `top100_{up,down}.rds` (intermediate) |
| `02_recovery_profile.R` | Same-perturbagen recovery at the profile level (all-source and cross-source, top 10/30/50/100), random expectation and random iterations | `recovery_profile_level.csv`, `recovery_profile_level_random_iterations.csv`, `recovery_{up,down}.rds` (intermediate) |
| `03_recovery_perturbagen.R` | Same-perturbagen recovery at the perturbagen level (one query profile per perturbagen) | `recovery_perturbagen_level.csv`, `recovery_perturbagen_level_iterations.csv` |
| `04_recovery_source_pairs.R` | Cross-source recovery by pair of data sources, profile and perturbagen levels | `recovery_by_source_pair.csv` |
| `05_condition_matching.R` | Cross-source recovery by matched experimental conditions (model, dose, duration), pooled by experimental model, and composition of conditions | `recovery_by_matched_conditions.csv` |
| `06_gene_level.R` | Gene-level agreement (Spearman and Pearson correlations) of same-perturbagen signatures across data sources, overall and by matched conditions | `gene_level_agreement.csv`, `gene_level_pairs.rds` (intermediate) |

`06_gene_level.R` reads all DE signature files and takes several hours. Files are cached in `output/gene_cache/`, so an interrupted run continues from the cache.

## Settings

| Setting | Value |
|---|---|
| Top N | 10, 30, 50, 100 (conditions: 100) |
| Random iterations / subsampling / bootstrap | 10,000 |
| Empirical P | (b + 1) / (B + 1) |
| Dose matched | same unit family, no more than twofold difference after conversion to a common unit |
| Duration matched | identical after conversion to hours |
| Background candidates (gene level) | up to 500 signatures per data source |
| Random seeds | given at the top of each script |

Results based on random iterations may differ slightly in the last digits from the reported values, depending on the R version and the random number sequence.

## Requirements

R (>= 4.3) with `data.table` and `matrixStats`.
