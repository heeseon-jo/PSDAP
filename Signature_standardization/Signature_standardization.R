# =============================================================================
# Signature_standardization.R
#
# Directional (up / down) gene sets from differential expression results.
#
# For every signature, upregulated and downregulated genes are defined
# separately at six stringency levels: adjusted P value <= 0.05 or <= 0.01,
# each combined with no fold change cutoff (|log2FC| >= 0) or an absolute
# log2 fold change cutoff of >= 1 or >= 2. Genes without a gene symbol are
# excluded. Deterministic; no random seed required.
#
# Input
#   degl_all : named list of data frames (names = signature IDs), one per
#              signature, with the columns Gene, log2FoldChange and padj.
#              Source-specific preprocessing (reading each resource, renaming
#              columns, mapping to human gene symbols) must be done beforehand.
#
# Output
#   cupL, cdnL : named lists of unique gene symbols (upregulated and
#                downregulated genes) with identical names
#                "<signature ID>@<log2FC cutoff>@<padj cutoff>"; the format
#                expected by run_hypergeometric() (Hypergeometric_test.R).
#
# Dependencies: base R only.
# =============================================================================


# 1. Gene ranked list
# Ranked gene lists are not built in this script. For signatures with
# differential expression statistics, genes are scored by
# sign(log2FC) * -log10(P) in read_score() (KS_test.R) and ordered by
# this score.



# 2. Directional gene sets at six stringency levels
#    up   : log2FC >=  l2fc_cut & padj <= padj_cut
#    down : log2FC <= -l2fc_cut & padj <= padj_cut
#    (l2fc_cut = 0 means no fold change cutoff)
cut_comb <- expand.grid(l2fc = c(0,1,2), padj = c(0.05, 0.01))
cupL <- list()
cdnL <- list()
for(k in seq_len(nrow(cut_comb))){
  l2fc_cut <- cut_comb$l2fc[k]
  padj_cut <- cut_comb$padj[k]
  
  for(sn in names(degl_all)){
    df  <- degl_all[[sn]]
    key <- paste0(sn, "@", l2fc_cut, "@", padj_cut)
    
    cupL[[key]] <- unique(df$Gene[df$log2FoldChange >=  l2fc_cut & df$padj <= padj_cut & !is.na(df$Gene)])
    cdnL[[key]] <- unique(df$Gene[df$log2FoldChange <= -l2fc_cut & df$padj <= padj_cut & !is.na(df$Gene)])
  }
  cat(sprintf("[%d/%d] Done\n", k, nrow(cut_comb)))
}

save(cupL, cdnL, file ="compare_up_dn_list_log2fc.Rdata")



# 3. CMap directional gene sets (top / bottom 100)
# Input:
#   m     : Level 5 moderated z-score matrix (genes x signatures)
#   ginfo : Gene annotation (gene_id, gene_symbol)

topn <- 100

cupL <- vector("list", ncol(m))
cdnL <- vector("list", ncol(m))

for (i in seq_len(ncol(m))) {
  v <- m[, i]
  v <- v[order(v, decreasing = TRUE)]
  
  names(v) <- ginfo$gene_symbol[match(names(v), ginfo$gene_id)]
  
  cupL[[i]] <- names(v)[v >= v[topn]]
  cdnL[[i]] <- names(v)[v <= v[length(v) - topn + 1]]
}

names(cupL) <- names(cdnL) <- colnames(m)

# Save
save(cupL, cdnL, file = "compare_up_dn_list_cmap.Rdata")