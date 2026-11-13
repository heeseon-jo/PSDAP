# Signature standardization

Standardization of each signature into directional (up / down) gene sets. Ranked
gene lists are not built here; they are created when the signatures are read for
the KS test (`read_score()` in `KS_enrichment.R`), where genes are scored by
`sign(log2FC) * -log10(P)` and ordered by this score.

| Script | Description |
|---|---|
| `Signature_standardization.R` | Builds the up- and down-regulated gene sets of each signature. For signatures with differential expression results, six stringency levels are used (2 padj cutoffs: 0.05, 0.01 x 3 absolute log2 fold change cutoffs: 0, 1, 2; a cutoff of 0 means no fold change cutoff). For CMap, the top and bottom 100 genes of the Level 5 moderated z-score are used. |

## Usage

The script has two independent sections. Run them separately, because both
create `cupL` and `cdnL`.

**Differential expression signatures (section 2)**

- Input
  - `degl_all`: named list of per-signature differential expression tables
    (names = signature IDs; columns `Gene`, `log2FoldChange` and `padj`)
- Output
  - `compare_up_dn_list_log2fc.Rdata` containing `cupL` (up) and `cdnL` (down),
    named lists of unique gene symbols. Names are
    `<signature ID>@<log2FC cutoff>@<padj cutoff>`, the format used by
    `run_hypergeometric()` in `hypergeometric_test.R`.

**CMap signatures (section 3)**

- Input
  - `m`: Level 5 moderated z-score matrix (genes x signatures)
  - `ginfo`: gene annotation (`gene_id`, `gene_symbol`)
- Output
  - `compare_up_dn_list_cmap.Rdata` containing `cupL` (top 100 genes) and `cdnL`
    (bottom 100 genes), named by signature.

## Requirements

R (v4.4.1); base R only.
