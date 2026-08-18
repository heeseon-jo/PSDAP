# PSDAP

Code accompanying the data descriptor "Perturbagen-induced transcriptomic signatures with functional, disease, and similarity annotations."

## Contents

- `hypergeometric_test.R` — Hypergeometric testing of up- and down-regulated gene sets against a gene set collection (pathway or disease). The background gene universe is defined per data source as the genes measured in that resource.

## Usage

```r
source("hypergeometric_test.R")

res <- run_hypergeometric(
  up_gene_list = up,
  dn_gene_list = dn,
  gene_sets    = gene_sets,
  universe     = universe,
  out_dir      = "results/"
)
```

## Requirements

R (v4.4.1); base packages only (stats, parallel).

## Data

The dataset is available on Zenodo: [DOI].

