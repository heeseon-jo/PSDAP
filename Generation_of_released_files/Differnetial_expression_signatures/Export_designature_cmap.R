# =============================================================================
# export_designature_cmap.R
#
# Writes the released differential-expression signature files for
# Connectivity Map, one file per signature:
#   Differential_expression_signatures/<Category>/DEsignature_<Signature_id>_Connectivity_Map.csv
#
# Gene-level moderated z-scores (MODZ) and their within-signature ranks are
# read from the internal PSDAP database, which is not distributed, and genes
# are annotated with gene symbols and Ensembl identifiers.
#
# Output columns (rows ordered by Gene_rank):
#   Symbol     gene symbol
#   Ensembl    Ensembl gene identifier
#   MODZ       moderated z-score
#   Gene_rank  rank of the gene within the signature
#
# Database credentials are read from environment variables (DB_HOST, DB_PORT,
# DB_NAME, DB_USER, DB_PW), e.g. set in ~/.Renviron.
#
# Dependencies: data.table, DBI, RMariaDB.
# =============================================================================

library(data.table)
library(DBI)

out_dir <- "output/Differential_expression_signatures"
clean_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

con <- dbConnect(RMariaDB::MariaDB(),
                 host = Sys.getenv("DB_HOST"), port = as.integer(Sys.getenv("DB_PORT")),
                 dbname = Sys.getenv("DB_NAME"),
                 user = Sys.getenv("DB_USER"), password = Sys.getenv("DB_PW"))

# 1) Connectivity Map signatures -------------------------------------------------
signature_info <- fread("output/signature_info.csv", data.table = FALSE)
target <- signature_info[signature_info$Source == "Connectivity Map",
                         c("Signature_id", "Source", "Category")]
target$cvc_id <- as.integer(sub("^SIG", "", target$Signature_id))
target <- target[order(target$cvc_id), ]

cmap_ids <- dbGetQuery(con, "SELECT DISTINCT cvc_id FROM cmap_mode_z")$cvc_id
stopifnot(all(target$cvc_id %in% cmap_ids))            # every signature has scores

gene_info <- dbGetQuery(con, "SELECT * FROM ptn_ginfo")

for (cat in unique(target$Category)) {
  dir.create(file.path(out_dir, clean_name(cat)), recursive = TRUE, showWarnings = FALSE)
}

# 2) One signature -------------------------------------------------------------
export_one <- function(cvc_id, sig_id, source, category) {
  path <- file.path(out_dir, clean_name(category),
                    sprintf("DEsignature_%s_%s.csv", sig_id, clean_name(source)))
  if (file.exists(path)) return("SKIP")
  
  cz <- dbGetQuery(con, sprintf("SELECT * FROM cmap_mode_z WHERE cvc_id = %d", cvc_id))
  if (nrow(cz) == 0) stop("no scores")
  
  idx <- match(cz$gene_id, gene_info$gene_id)
  out <- data.frame(Symbol    = gene_info$Symbol[idx],
                    Ensembl   = gene_info$Ensembl[idx],
                    MODZ      = cz$mode_z,
                    Gene_rank = cz$mode_z_rank)
  out <- out[order(out$Gene_rank), ]
  
  n_na <- colSums(is.na(out[, c("Symbol", "MODZ", "Gene_rank")]))
  if (any(n_na > 0)) stop("NA: ", paste(names(n_na)[n_na > 0], n_na[n_na > 0],
                                        sep = "=", collapse = ","))
  
  fwrite(out, path, na = "", bom = TRUE)
  "OK"
}

# 3) All signatures --------------------------------------------------------------
log <- target[, c("Signature_id", "Source", "Category")]
log$status <- ""
log$reason <- ""
for (i in seq_len(nrow(target))) {
  r <- tryCatch(export_one(target$cvc_id[i], target$Signature_id[i],
                           target$Source[i], target$Category[i]),
                error = function(e) conditionMessage(e))
  log$status[i] <- if (r %in% c("OK", "SKIP")) r else "ERROR"
  log$reason[i] <- if (log$status[i] == "ERROR") r else ""
}
dbDisconnect(con)

print(table(log$status))
fwrite(log, file.path(out_dir, "export_log_cmap.csv"), bom = TRUE)