# Differential expression analysis

Differential expression (DE) analysis of the individual experiments. The DE
method used for each source is listed in `source_info.csv` (DESeq2 for RNA-seq
sources, limma for microarray sources). Sources that provide precomputed
results (TCM Atlas, DrugMatrix) are not re-analysed here.

| Script | Description |
|---|---|
| `RNA_seq_DEG_analysis.R` | DE analysis of RNA-seq count data with DESeq2.
| `Microarray_DEG_analysis.R` | DE analysis of microarray data with limma.

## Output

TODO: columns of the result tables. The downstream steps expect, for each
signature, a data frame with `Gene`, `log2FoldChange` and `padj`
(see `C_Signature_standardization`).


## Usage
### RNA-seq
```r
# case_ids    : unique_ID of the perturbed samples
# ctrl_ids    : unique_ID of the control samples
# meta        : sample metadata with the columns unique_ID and plate_rep
# count_mat_t : raw read count matrix, genes x samples (colnames = unique_ID)

source("RNA_seq_DEG_analysis.R")
head(res)
```
The result is the data frame `res`, overwritten for each treatment. Run the
script once per treatment, for example in a loop that sets `case_ids` and
`ctrl_ids`.


### Microarray

```r
source("Microarray_DEG_analysis.R")

case_cels <- parse_cels(case_cel_string)   # e.g. "GSM001.CEL, GSM002.CEL"
ctrl_cels <- parse_cels(ctrl_cel_string)

out <- run_limma(expr_mat, case_cels, ctrl_cels)
out$reason    # "OK" or the reason why the treatment was skipped
out$result    # limma results of all genes
```
`expr_mat` is the normalized expression matrix (genes x samples) whose column
names are the CEL file names; names are compared in upper case.

## Requirements

R (v4.4.1). 

- RNA-seq: DESeq2, edgeR.
- Microarray: limma.
