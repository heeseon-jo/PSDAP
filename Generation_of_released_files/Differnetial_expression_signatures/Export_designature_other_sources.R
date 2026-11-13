# =============================================================================
# export_designature_other_sources.R
#
# Writes the released differential-expression signature files for all sources
# except Connectivity Map, one file per signature:
#   Differential_expression_signatures/<Category>/DEsignature_<Signature_id>_<Source>.csv
#
# Gene-level differential-expression results (log2 fold change, P value,
# adjusted P value) are read from the internal PSDAP database, which is not
# distributed, and genes are annotated with gene symbols and Ensembl
# identifiers.
#
# Output columns (rows ordered by adjusted P value, then by decreasing
# absolute log2 fold change):
#   Symbol            gene symbol
#   Ensembl           Ensembl gene identifier
#   Log2FoldChange    log2 fold change
#   P_value           P value
#   Adjusted_p_value  adjusted P value
#   Symbol_species    gene symbol in the species of the experiment
#
# Genes with a missing Symbol, Log2FoldChange, P_value or Adjusted_p_value are
# removed.
#
# Database credentials are read from environment variables (DB_HOST, DB_PORT,
# DB_NAME, DB_USER, DB_PW), e.g. set in ~/.Renviron.
#
# Dependencies: data.table, DBI, RMariaDB.
# =============================================================================

library(data.table)
library(DBI)

out_dir  <- "output/Differential_expression_signatures"
key_cols <- c("Symbol", "Log2FoldChange", "P_value", "Adjusted_p_value")
clean_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

con <- dbConnect(RMariaDB::MariaDB(),
                 host = Sys.getenv("DB_HOST"), port = as.integer(Sys.getenv("DB_PORT")),
                 dbname = Sys.getenv("DB_NAME"),
                 user = Sys.getenv("DB_USER"), password = Sys.getenv("DB_PW"))

# 1) Signatures other than Connectivity Map -------------------------------------
signature_info <- fread("output/signature_info.csv", data.table = FALSE)
target <- signature_info[signature_info$Source != "Connectivity Map",
                         c("Signature_id", "Source", "Category")]
target$cvc_id <- as.integer(sub("^SIG", "", target$Signature_id))
target <- target[order(target$cvc_id), ]

gene_info <- dbGetQuery(con, "SELECT * FROM ptn_ginfo")

for (cat in unique(target$Category)) {
  dir.create(file.path(out_dir, clean_name(cat)), recursive = TRUE, showWarnings = FALSE)
}

# 2) One signature -------------------------------------------------------------
export_one <- function(cvc_id, sig_id, source, category) {
  path <- file.path(out_dir, clean_name(category),
                    sprintf("DEsignature_%s_%s.csv", sig_id, clean_name(source)))
  
  deg <- dbGetQuery(con, sprintf("SELECT * FROM deg WHERE deg_cvc_id = %d", cvc_id))
  if (nrow(deg) == 0) stop("no genes")
  
  idx <- match(deg$deg_gene_id, gene_info$gene_id)
  out <- data.frame(Symbol           = gene_info$Symbol[idx],
                    Ensembl          = gene_info$Ensembl[idx],
                    Log2FoldChange   = deg$log2FoldChange,
                    P_value          = deg$pvalue,
                    Adjusted_p_value = deg$padj,
                    Symbol_species   = deg$Symbol)
  
  out <- out[complete.cases(out[, key_cols]), ]          # remove genes with missing values
  out <- out[order(out$Adjusted_p_value, -abs(out$Log2FoldChange)), ]
  
  fwrite(out, path, na = "", bom = TRUE)
}

# 3) All signatures --------------------------------------------------------------
for (i in seq_len(nrow(target))) {
  export_one(target$cvc_id[i], target$Signature_id[i], target$Source[i], target$Category[i])
}
dbDisconnect(con)