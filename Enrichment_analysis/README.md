# Enrichment analysis

Enrichment of each signature against gene set collections (e.g. pathway or
disease gene sets). Two approaches are used: a KS test on the ranked gene list of
a signature, and a hypergeometric test on its up- and down-regulated gene sets.

| Script | Description |
|---|---|
| `KS_test.R` | KS-test enrichment of a differential expression signature. Genes are scored by `sign(log2FC) * -log10(P)`, and each gene set is tested against all other genes of the signature, separately for UP and DOWN. Gene sets with fewer than 5 genes in the signature are not tested. |
| `Hypergeometric_test.R` | Hypergeometric test of the up- and down-regulated gene sets of each signature against the gene set collection. Gene sets are restricted to the background gene universe (genes measured in each data source) before testing; p values are adjusted by the Benjamini-Hochberg method. |

## Usage

### KS test

```r
source("KS_test.R")

gs  <- make_gs(nested)                    
sc  <- read_score("signature.xlsx")       # columns: Symbol, Log2FoldChange, P_value
raw <- run_ks(sc, gs)                     
res <- make_processed(raw, padj_cutoff = 0.05)  
```

`make_processed()` applies the Benjamini-Hochberg adjustment separately for UP
and DOWN, so a gene set significant in both directions appears in both.

### Hypergeometric test

```r
source("Hypergeometric_test.R")

run_hypergeometric(
  up_gene_list = cupL,    
  dn_gene_list = cdnL,
  gene_sets    = gene_sets,
  universe     = universe, 
  out_dir      = "hypergeometric_results/source_name"
)
```

`cupL` and `cdnL` are the directional gene sets created in
`Signature_standardization`. One `.rds` file (list with `up` and `dn` result
tables) is written per signature to `out_dir`; existing files are skipped, so an
interrupted run can be resumed.

## Requirements

R (v4.4.1).

- `KS_test.R`: readxl, data.table, parallel.
- `Hypergeometric_test.R`: base R only (stats, parallel).
