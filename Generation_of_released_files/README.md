# Generation of released files

Generation of the released metadata files and the per-signature result files from
the internal PSDAP database. The database is not distributed, so these scripts
document how the released files were generated and cannot be run without it.

| Folder | Description |
|---|---|
| [`Metadata`](Metadata) | Builds the metadata files `signature_info.csv`, `model_info.csv`, `perturbagen_info.csv` and `source_info.csv`. |
| [`Differential_expression_signature`](Differential_expression_signature) | Writes the differential expression signature files, one file per signature. |
| [`Enrichment_analysis`](Enrichment_analysis) | Writes the pathway and disease enrichment files (KS test and hypergeometric test), one file per signature. |
| [`Perturbagen_similarity`](Perturbagen_similarity) | Writes the perturbagen similarity tables. |

## Usage

The scripts of the `Metadata` folder are run first, because the scripts of the
other folders read `output/signature_info.csv` created by
`Build_signature_info.R`. The database credentials are read from the environment
variables `DB_HOST`, `DB_PORT`, `DB_NAME`, `DB_USER` and `DB_PW`. All files are
written to `output/`. See the README of each folder for the scripts.

## Requirements

R (v4.4.1). The packages are listed in the README of each folder.
