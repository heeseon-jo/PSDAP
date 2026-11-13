# Ortholog mapping

Conversion of rat and mouse gene symbols to human gene symbols. Rat symbols are
converted with mapping tables built from Ensembl BioMart, using
`convert_to_human()`. Mouse symbols are converted with the Ensembl annotation
package EnsDb.Mmusculus.v79 and the MGI file `HMD_HumanPhenotype.rpt`
(downloaded from MGI in 2020).

| Script | Description |
|---|---|
| `Build_ortholog_table_rat.R` | Builds the rat-to-human mapping tables (current and former rat symbols) from Ensembl BioMart. |
| `Convert_rat_to_human.R` | `convert_to_human()`: converts rat gene symbols to human gene symbols using the mapping tables, with NCBI HomoloGene as a last resort. |
| `Build_ortholog_table_mouse.R` | Builds the mouse gene table with the human ortholog (symbol and Entrez ID) of each mouse gene, using EnsDb.Mmusculus.v79 and `HMD_HumanPhenotype.rpt`. |

## Usage
### Build_ortholog_table_rat.R

```
Rscript Build_ortholog_table_rat.R
```

The tables are saved to `gene_mapping/MAPPING_orth_tbl_human_rat_final.Rdata`. The
script needs internet access to Ensembl BioMart, and its result depends on the
Ensembl release at the time of the query.


### Convert_rat_to_human.R

```r
load("gene_mapping/MAPPING_orth_tbl_human_rat_final.Rdata")
source("Convert_rat_to_human.R")

# rat_symbols: character vector of rat gene symbols, already cleaned of
#              dataset-specific prefixes or suffixes
res <- convert_to_human(
  symbols       = rat_symbols,
  mapping       = map_use_all_human_rat,
  mapping_alias = map_use_all_human_rat_alias
)
```

`res` has one row per input symbol: `orig_symbol`, `human_symbol` (NA if
unmapped) and `mapping_class` (`homologene` for symbols resolved by HomoloGene).


### Build_ortholog_table_mouse.R

```r
source("Build_ortholog_table_mouse.R")   # creates the mouse gene table ginfo
```

`HMD_HumanPhenotype.rpt` is read from `gene_mapping/` (change `gene_dir` in the
script as needed). `ginfo` has the columns `Ensembl`, `Entrez`, `Symbol`,
`GENEBIOTYPE`, `Human_Symbol` and `Human_Entrez` (`-` if missing).

Mouse gene symbols are converted to human gene symbols by matching them to ginfo. If a mouse symbol matches more than one row of ginfo, the first matching row is used.

## Requirements

R (v4.4.1).

- `Build_ortholog_table_rat.R`: biomaRt, dplyr.
- `Convert_rat_to_human.R`: dplyr, homologene.
- `Build_ortholog_table_mouse.R`: EnsDb.Mmusculus.v79.
