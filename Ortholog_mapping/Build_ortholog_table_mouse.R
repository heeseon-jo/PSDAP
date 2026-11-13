# =============================================================================
# Build_ortholog_table_mouse.R
#
# Build the mouse gene table `ginfo` with human orthologs.
#
# Mouse genes (Ensembl gene ID, Entrez ID, gene symbol and gene biotype) are
# taken from the Ensembl annotation package EnsDb.Mmusculus.v79. Each mouse gene
# is linked by its gene symbol to its human ortholog (human symbol and Entrez ID)
# in the MGI file HMD_HumanPhenotype.rpt, which was downloaded from MGI in 2020.
#
# Mouse only. If a mouse symbol appears more than once in the MGI file, the
# first entry is used. Mouse genes without an entry in the MGI file, and missing
# values in general, are set to "-".
#
# Input (in `gene_dir`):
#   HMD_HumanPhenotype.rpt   MGI human-mouse ortholog report
#
# Result:
#   ginfo   data frame with the columns Ensembl, Entrez, Symbol, GENEBIOTYPE,
#           Human_Symbol and Human_Entrez
#
# Dependencies: EnsDb.Mmusculus.v79
# =============================================================================

library(EnsDb.Mmusculus.v79)

gene_dir <- "gene_mapping"   # change as needed


# 1. Mouse genes from EnsDb
edb <- EnsDb.Mmusculus.v79
gid <- keys(edb, keytype = "GENEID")
gdb <- ensembldb::select(edb, keys = gid,
                         columns = c("ENTREZID", "SYMBOL", "GENEBIOTYPE"),
                         keytype = "GENEID")
colnames(gdb) <- c("Ensembl", "Entrez", "Symbol", "GENEBIOTYPE")


# 2. MGI human-mouse ortholog file
ort <- read.csv(file = file.path(gene_dir, "HMD_HumanPhenotype.rpt"),
                header = FALSE, sep = "\t", check.names = FALSE)
ort <- ort[, -c(4, 8)]
colnames(ort) <- c("Human_Symbol", "Human_Entrez", "HomoloGene_ID",
                   "Mouse_Symbol", "MGI_AccessionID", "Phenotype ID")


# 3. Add the human ortholog to each mouse gene
gdb <- cbind(gdb, ort[match(gdb$Symbol, ort$Mouse_Symbol), 1:2])
ginfo <- gdb
ginfo[is.na(ginfo)] <- "-"

# 4. Save
save(ginfo, file = file.path(gene_dir, "MAPPING_orth_tbl_human_mouse_final.Rdata"))