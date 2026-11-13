# Metadata

Generation of the released metadata files `signature_info.csv`, `model_info.csv`,
`perturbagen_info.csv` and `source_info.csv` from the internal PSDAP database,
which is not distributed. Only comparisons marked for open release are included.

| Script | Description |
|---|---|
| `Build_signature_info.R` | Builds `signature_info.csv` (one row per signature) with the perturbagen, model, dose, duration, extraction, control condition and source of each case versus control comparison. |
| `Build_model_info.R` | Builds `model_info.csv` (one row per model used in `signature_info.csv`) with the standard identifier (Cellosaurus or UBERON), species, disease, tissue and synonyms. |
| `Build_perturbagen_info.R` | Builds `perturbagen_info.csv` (one row per standardized perturbagen) with the standard identifier, its type and synonyms. |
| `Build_source_info.R` | Builds `source_info.csv` (one row per data source) with the category, assay, device, DE method and reference. |

## Usage

The database credentials are read from the environment variables `DB_HOST`,
`DB_PORT`, `DB_NAME`, `DB_USER` and `DB_PW` (for example set in `~/.Renviron`).
`Build_signature_info.R` must be run before `Build_model_info.R`, which reads
`output/signature_info.csv`.

```
Rscript Build_signature_info.R
Rscript Build_model_info.R
Rscript Build_perturbagen_info.R
Rscript Build_source_info.R
```

The files are written to `output/`.

> **Do not open `perturbagen_info.csv` in Excel and save it (Ctrl+S). The columns may be corrupted.**

## Requirements

R (v4.4.1).

- `Build_signature_info.R`, `Build_model_info.R`: data.table, DBI, RMariaDB.
- `Build_perturbagen_info.R`: data.table, DBI, RMariaDB, readr, dplyr.
- `Build_source_info.R`: DBI, RMariaDB, readr, dplyr.
