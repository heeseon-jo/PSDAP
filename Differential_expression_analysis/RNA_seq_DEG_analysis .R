# ==============================================================================
# RNA_seq_DEG_analysis.R
#
# Purpose
#   Per-treatment differential expression analysis of RNA-seq count data with
#   DESeq2 (treated vs. matched control)
#
# Model
#   design               : ~ plate_rep + group (design may differ between datasets)
#   contrast             : group "2" vs "1" ("2" for perturbed sample, "1" for control sample)
#   multiple testing     : Benjamini-Hochberg (pAdjustMethod = "fdr")
#
# Dependencies: R 4.4.1; DESeq2, edgeR
# ==============================================================================



# ==============================================================================
# Differential expression analysis of one treatment (perturbed samples) vs its control samples.
# Input
#   case_ids    : unique_ID of the perturbed samples.
#   ctrl_ids    : unique_ID of the control samples.
#   meta        : sample metadata.
#   count_mat_t : raw read count matrix, genes x samples.
# Output
#   res         : DESeq2 results for all tested genes, including BH-adjusted p-values
# ==============================================================================



library(DESeq2)
library(edgeR)


# Factorize group
an <- rbind(
  data.frame(ID=ctrl_ids, fac=1, stringsAsFactors=FALSE),
  data.frame(ID=case_ids, fac=2, stringsAsFactors=FALSE)
)
an$plate_rep <- meta$plate_rep[match(an$ID, meta$unique_ID)]

group     <- factor(an$fac, levels=c(1,2))
plate_rep <- factor(an$plate_rep)



# Extract count matrix
xm <- count_mat_t[, an$ID, drop=FALSE]
storage.mode(xm) <- "integer"

keep <- filterByExpr(xm, group=group)
xm_use <- xm[keep, , drop=FALSE]



# DESeq2
## Create coldata
coldata <- data.frame(
  plate_rep = plate_rep,
  group     = group,
  row.names = colnames(xm_use)
)



## Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = xm_use,
  colData   = coldata,
  design    = ~ plate_rep + group
)
dds <- dds[rowSums(counts(dds)) > 0, ]

if (nrow(dds) < 10) {
  stop(sprintf("too few genes: %d", nrow(dds)))
}



## size factor: first for ratio. If unavailable, poscount
sf_method <- "ratio"
dds <- tryCatch({
  estimateSizeFactors(dds, type="ratio")
}, error = function(e) {
  message("sparse counts: using poscounts size factors")
  sf_method <<- "poscounts"
  estimateSizeFactors(dds, type="poscounts")
})



## Run DEseq2
ddsMF <- DESeq(dds, quiet=TRUE, fitType="local")
resMF <- results(
  ddsMF,
  contrast             = c("group","2","1"),
  cooksCutoff          = FALSE,
  pAdjustMethod        = "fdr",
  independentFiltering = FALSE
)



## DESeq2 results
res <- as.data.frame(resMF)
res <- data.frame(
  gene = rownames(res),
  res,
  keep = TRUE,
  sizeFactor_method = sf_method,
  check.names = FALSE
)


## Prevent zero-value of adjusted P-value.
res$logFDR <- ifelse(is.na(res$padj)|res$padj<=0, NA, -log10(res$padj))
res <- na.omit(res)