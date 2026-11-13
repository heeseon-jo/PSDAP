# Enrichment analysis

Generation of the released pathway and disease enrichment files, one file per
signature, from the results stored in the internal PSDAP database, which is not
distributed. Only signatures released to the public are exported, and signatures
without any result have no file.

| Script | Description |
|---|---|
| `Export_pathway_ks_test.R` | Writes the pathway enrichment files of the KS test: gene sets with an adjusted P value <= 0.05 in the UP or DOWN direction. |
| `Export_disease_ks_test.R` | Writes the disease enrichment files of the KS test: gene sets with an adjusted P value <= 0.05 in the UP or DOWN direction. |
| `Export_pathway_hypergeometric_test.R` | Writes the pathway enrichment files of the hypergeometric test, with the direction and the stringency level of the input gene sets. |
| `Export_disease_hypergeometric_test.R` | Writes the disease enrichment files of the hypergeometric test, with the direction and the stringency level of the input gene sets. |

## Usage

The database credentials are read from the environment variables `DB_HOST`,
`DB_PORT`, `DB_NAME`, `DB_USER` and `DB_PW` (for example set in `~/.Renviron`).
All scripts read `output/signature_info.csv`, which is created by
`Build_signature_info.R` in `Metadata`.

```
Rscript Export_pathway_ks_test.R
Rscript Export_disease_ks_test.R
Rscript Export_pathway_hypergeometric_test.R
Rscript Export_disease_hypergeometric_test.R
```

The files are written to `output/<Pathway or Disease>_enrichment/<Kolmogorov_Smirnov_test or Hypergeometric_test>/<Category>/`

(`Category`: Bioactive_compound, Herbal_formula, Herbal_ingredient, Medicinal_herb)
and named `PathwayKS_`, `DiseaseKS_`, `PathwayHG_` or `DiseaseHG_` followed by
`<Signature_id>_<Source>.csv`. Signatures whose file already exists are skipped,
so an interrupted run can be resumed.

## Requirements

R (v4.4.1).

- `Export_pathway_ks_test.R`, `Export_disease_ks_test.R`: DBI, RMariaDB, dplyr, tidyr, pbapply.
- `Export_pathway_hypergeometric_test.R`, `Export_disease_hypergeometric_test.R`: DBI, RMariaDB, dplyr, pbapply.
