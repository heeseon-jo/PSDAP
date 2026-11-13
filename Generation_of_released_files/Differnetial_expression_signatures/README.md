# Differential expression signature

Generation of the released differential expression signature files, one file per
signature, from the internal PSDAP database, which is not distributed. Genes are
annotated with gene symbols and Ensembl identifiers, and only signatures released
to the public are exported.

| Script | Description |
|---|---|
| `Export_designature_other_sources.R` | Writes the signature files of all sources except Connectivity Map with the log2 fold change, P value and adjusted P value of each gene. Genes with a missing value are removed, and rows are ordered by adjusted P value and absolute log2 fold change. |
| `Export_designature_cmap.R` | Writes the signature files of Connectivity Map with the moderated z-score (MODZ) and rank of each gene, ordered by rank. Existing files are skipped, and the result of each signature is saved in `export_log_cmap.csv`. |

## Usage

The database credentials are read from the environment variables `DB_HOST`,
`DB_PORT`, `DB_NAME`, `DB_USER` and `DB_PW` (for example set in `~/.Renviron`).
Both scripts read `output/signature_info.csv`, which is created by
`Build_signature_info.R` in `Metadata`.

```
Rscript Export_designature_other_sources.R
Rscript Export_designature_cmap.R
```

The files are written to
`output/Differential_expression_signatures/<Category>/DEsignature_<Signature_id>_<Source>.csv`.

## Requirements

R (v4.4.1). Both scripts: data.table, DBI, RMariaDB.
