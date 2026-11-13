# Perturbagen similarity

Generation of the released perturbagen similarity tables from the pairwise GO
semantic similarity scores stored in the internal PSDAP database, which is not
distributed. For each signature, the 100 most similar signatures per GO ontology
(BP, CC, MF) are kept, separately for the up and down direction. Only signatures
released to the public are included.

| Script | Description |
|---|---|
| `Export_perturbagen_similarity_non_CMAP.R` | Writes one similarity file per stringency level (log2 fold change and adjusted P value cutoffs) for all signatures except the Connectivity Map top 100 lists. |
| `Export_perturbagen_similarity_CMAP.R` | Writes the similarity file of the Connectivity Map signatures (top 100 gene lists), which also contains the stringency level of each signature of a pair. |

## Usage

The database credentials are read from the environment variables `DB_HOST`,
`DB_PORT`, `DB_NAME`, `DB_USER` and `DB_PW` (for example set in `~/.Renviron`).
Both scripts read `output/signature_info.csv`, which is created by
`Build_signature_info.R` in `Metadata`.

```
Rscript Export_perturbagen_similarity_non_CMAP.R
Rscript Export_perturbagen_similarity_CMAP.R
```

The files are written to `output/Perturbagen_similarity/` as
`Similarity_<Stringency_level>.csv` and `Similarity_Connectivity_map_top100.csv`.
A script stops if its output file already exists.

## Requirements

R (v4.4.1). Both scripts: DBI, RMariaDB, dplyr, data.table.
