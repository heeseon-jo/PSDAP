# =============================================================================
# Export_disease_hypergeometric_test.R
#
# Writes the released disease enrichment (hypergeometric test) results, one
# file per signature:
#   Disease_enrichment/Hypergeometric_test/<Category>/DiseaseHG_<Signature_id>_<Source>.csv
#
# Hypergeometric test results (P value, adjusted P value, overlap
# count) are read from the internal PSDAP database, which is not distributed,
# and annotated with disease gene set names, databases and stringency levels.
# Only signatures released to the public are exported; signatures without any
# enrichment result in the database have no file.
#
# Output columns:
#   Gene_set_name     disease gene set name
#   Database          database of the gene set
#   P_value           P value of the hypergeometric test
#   Adjusted_P_value  adjusted P value
#   Direction         UP or DOWN
#   Count             number of overlapping genes
#   Stringency_level  cutoffs used to define the input genes
#                     ("Log2FoldChange_<|FC cutoff|>_AdjustedPvalue_<FDR cutoff>")
#
# Requires output/signature_info.csv (build_signature_info.R).
#
# Database credentials are read from environment variables (DB_HOST, DB_PORT,
# DB_NAME, DB_USER, DB_PW), e.g. set in ~/.Renviron.
#
# Dependencies: DBI, RMariaDB, dplyr, pbapply.
# =============================================================================

library(DBI)
library(dplyr)
library(pbapply)



con <- dbConnect(RMariaDB::MariaDB(),
                host = Sys.getenv("DB_HOST"), port = as.integer(Sys.getenv("DB_PORT")),       
                dbname = Sys.getenv("DB_NAME"),
                user = Sys.getenv("DB_USER"), password = Sys.getenv("DB_PW"))



# 1) Signature IDs -------------------------------------------------------------
# Signature_id (e.g. SIG123) -> numeric cvc_id used in the database
sig_info <- read.csv("output/signature_info.csv")                                                  
sig <- sig_info[,c("Signature_id"), drop = FALSE]
sig$cvc_id <- gsub("SIG", "",sig$Signature_id) %>% as.integer()



# 2) Signatures to export ------------------------------------------------------
odn <-  dbGetQuery(con, "select * from option_db_name;")
pert <- dbGetQuery(con, "select * from option_perturbagen;")
cvc <- dbGetQuery(con, "select * from c_vs_c;")
opt <- dbGetQuery(con, "select * from option_pert_type;")

# Keep signatures released to the public
cvc <- cvc[is.na(cvc$exclude_from_open) | cvc$exclude_from_open != 1, ]
cvc <- cvc[,c("cvc_id", "case_id","pertur_id", "exclude_from_open")]

# Attach source name, perturbagen category and Signature_id
cvc$odn_id   <- pert$odn_id[match(cvc$pertur_id, pert$pertur_id)]
cvc$db_name  <- odn$db_name[match(cvc$odn_id, odn$odn_id)]
cvc$opt_id   <- odn$opt_id[match(cvc$odn_id, odn$odn_id)]
cvc$pert_type <- opt$pert_type[match(cvc$opt_id, opt$opt_id)]
cvc$sig_id <- sig$Signature_id[match(cvc$cvc_id, sig$cvc_id)]

cvc <- cvc[,c("cvc_id", "db_name", "pert_type", "sig_id")]



# 3) Disease gene set annotation -----------------------------------------------
# disease_id -> gene set name and database
dsub <- dbGetQuery(con, "select * from disease_subcat;") 
dgan <- dbGetQuery(con, "select * from disease_gannot;")

d <- left_join(dgan, dsub, by = "dsubcat_id")
d <- d[,c("disease_id", "description", "gs_name")]



# 4) Stringency levels ---------------------------------------------------------
# FF_id -> direction (UP/DOWN) and cutoff label (|log2FC| and adjusted P-value)
fdc <- dbGetQuery(con, "select * from fdrfc;")

fdc$FDR_cut_rev <- gsub(pattern = "\\.",replacement = "_", x = fdc$FDR_cut)
fdc$Stringency_level <- paste0("Log2FoldChange_", abs(fdc$FC_cut),"_AdjustedPvalue_",fdc$FDR_cut_rev) 
ff <- fdc[,c("FF_id", "up_down","Stringency_level")]








# 5) Export hypergeometric test results per category ------------------------------
base_dir <- "output/Disease_enrichment/Hypergeometric_test"
opt_name <- c("Bioactive compound", "Herbal formula", "Herbal ingredient", "Medicinal herb")
chunk_size <- 1000

for (i in opt_name) {
  # Output directory for this category
  target_dir <- file.path(base_dir, gsub(" ", "_", i))
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("target file path : %s\n", target_dir))
  
  # All signatures of this category
  cvc_id_all <- cvc$cvc_id[cvc$pert_type == i]
  n_total <- length(cvc_id_all)
  cat(sprintf("total cvc_id count : %d\n", n_total))
  
  # Skip signatures whose file already exists
  sig_ids_all <- cvc$sig_id[match(cvc_id_all, cvc$cvc_id)]
  sours_all   <- cvc$db_name[match(cvc_id_all, cvc$cvc_id)]
  expected_files_all <- sprintf("DiseaseHG_%s_%s.csv", sig_ids_all, sours_all)
  
  existing_files <- list.files(target_dir)
  already_done <- expected_files_all %in% existing_files
  
  cvc_id_todo <- cvc_id_all[!already_done]
  cat(sprintf("already done : %d, remaining : %d\n", sum(already_done), length(cvc_id_todo)))
  
  # Split the remaining signatures into chunks
  chunks <- split(cvc_id_todo, ceiling(seq_along(cvc_id_todo) / chunk_size))
  if (length(cvc_id_todo) == 0) chunks <- list()
  cat(sprintf("total chunks : %d\n\n", length(chunks)))
  
  for (ci in seq_along(chunks)) {
    cvc_id_chunk <- chunks[[ci]]
    
    # Upload the chunk's cvc_ids to a temporary table
    dbExecute(con, "DROP TEMPORARY TABLE IF EXISTS tmp_cvc_ids")
    dbExecute(con, "CREATE TEMPORARY TABLE tmp_cvc_ids (cvc_id BIGINT PRIMARY KEY)")
    dbWriteTable(con, "tmp_cvc_ids", data.frame(cvc_id = cvc_id_chunk), append = TRUE, row.names = FALSE)
    
    # Fetch HG results of the chunk
    ora_chunk <- dbGetQuery(con, "
      SELECT o.cvc_id, o.FF_id, o.disease_id, o.pvalue, o.padj, o.count
      FROM ora_disease o
      INNER JOIN tmp_cvc_ids t ON o.cvc_id = t.cvc_id
    ")
    dbExecute(con, "DROP TEMPORARY TABLE tmp_cvc_ids")
    
    
    # Add gene set name, database, direction and stringency level
    ora_chunk$gene_set <- d$gs_name[match(ora_chunk$disease_id, d$disease_id)]
    ora_chunk$db <- d$description[match(ora_chunk$disease_id, d$disease_id)]
    ora_chunk$direc <- toupper(ff$up_down[match(ora_chunk$FF_id, ff$FF_id)])
    ora_chunk$strin <- ff$Stringency_level[match(ora_chunk$FF_id, ff$FF_id)]
    
    # Select and rename output columns (cvc_id is kept to split by signature)
    ora_df <- ora_chunk[, c("cvc_id", "gene_set", "db", "pvalue", "padj", "direc", "count", "strin")]
    colnames(ora_df) <- c("cvc_id", "Gene_set_name", "Database", "P_value", "Adjusted_P_value", "Direction", "Count", "Stringency_level")
    
    # Write one CSV per signature
    ora_split <- split(ora_df, ora_df$cvc_id)
    
    invisible(pblapply(names(ora_split), function(id) {
      sour <- cvc$db_name[match(as.numeric(id), cvc$cvc_id)] #Source name
      sig_id <- cvc$sig_id[match(as.numeric(id), cvc$cvc_id)] #SIG ID
      df_out <- ora_split[[id]][, -1]
      out_file <- file.path(target_dir, sprintf("DiseaseHG_%s_%s.csv", sig_id, sour))
      write.csv(df_out, out_file, row.names = FALSE)
    }))
    
    cat(sprintf("  chunk %d/%d : %d rows saved\n", ci, length(chunks), nrow(ora_chunk)))
    
    rm(ora_chunk, ora_df, ora_split)
    gc()
  }
}

dbDisconnect(con)