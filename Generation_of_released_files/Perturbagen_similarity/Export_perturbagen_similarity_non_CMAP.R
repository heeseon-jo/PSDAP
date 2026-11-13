# =============================================================================
# Export_perturbagen_similarity_non_CMAP.R
#
# Writes the released perturbagen similarity tables for signatures defined by
# log2 fold-change and adjusted P-value cutoffs (all stringency levels except
# the Connectivity Map top-100 lists), one file per stringency level:
#   Perturbagen_similarity/Similarity_<Stringency_level>.csv
#
# Pairwise GO semantic similarity scores between signatures (tables go_up and
# go_dn, for up- and down-regulated genes) are read from the internal PSDAP
# database, which is not distributed. For each signature of a stringency level,
# the 100 most similar signatures per GO ontology (BP, CC, MF) are kept,
# separately for the up and down direction. Only signatures released to the
# public are included. The Connectivity Map signatures are handled by
# Export_perturbagen_similarity_CMAP.R.
#
# Output columns (one row per signature pair):
#   Signature_id_1, Signature_id_2          signature identifiers
#   BP_term_count_1, BP_term_count_2        number of GO terms used (BP)
#   MF_term_count_1, MF_term_count_2        number of GO terms used (MF)
#   CC_term_count_1, CC_term_count_2        number of GO terms used (CC)
#   BP_score, MF_score, CC_score            GO semantic similarity score
#   Direction                               Up or Down
#
# The script stops if an output file already exists.
#
# Requires output/signature_info.csv (build_signature_info.R).
#
# Database credentials are read from environment variables (DB_HOST, DB_PORT,
# DB_NAME, DB_USER, DB_PW), e.g. set in ~/.Renviron.
#
# Dependencies: DBI, RMariaDB, dplyr, data.table.
# =============================================================================

library(DBI)
library(RMariaDB)
library(dplyr)
library(data.table)



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
cvc <- dbGetQuery(con, "select * from c_vs_c;")

# Keep signatures released to the public
cvc <- cvc[is.na(cvc$exclude_from_open) | cvc$exclude_from_open != 1, ]
cvc <- cvc[,c("cvc_id", "case_id","pertur_id", "exclude_from_open")]

# Attach Signature_id
cvc$sig_id <- sig$Signature_id[match(cvc$cvc_id, sig$cvc_id)]

cvc <- cvc[,c("cvc_id", "sig_id")]



# 3) GO term counts per signature ----------------------------------------------
# gtc_id identifies a signature at one stringency level in the similarity tables
gcount <- dbGetQuery(con, "select * from goterm_count;") 

gc <- gcount[,c("gtc_id", "cvc_id", "FF_id", "BPO_used_go_n", "CCO_used_go_n", "MFO_used_go_n")]
gc <- gc[gc$cvc_id %in% cvc$cvc_id,]

# Temporary table of the gtc_ids of released signatures
dbExecute(con, "DROP TEMPORARY TABLE IF EXISTS tmp_open_gtc_ids")
dbExecute(con, "CREATE TEMPORARY TABLE tmp_open_gtc_ids (gtc_id BIGINT PRIMARY KEY)")
dbWriteTable(con, "tmp_open_gtc_ids", data.frame(gtc_id = unique(gc$gtc_id)), append = TRUE, row.names = FALSE)



# 4) Stringency levels ---------------------------------------------------------
# FF_id -> direction (up/down) and cutoff label (|log2FC| and adjusted P-value)
fdc <- dbGetQuery(con, "select * from fdrfc;")

fdc$FDR_cut_rev <- gsub(pattern = "\\.",replacement = "_", x = fdc$FDR_cut)
fdc$Stringency_level <- paste0("Log2FoldChange_", abs(fdc$FC_cut),"_AdjustedPvalue_",fdc$FDR_cut_rev) 
ff <- fdc[,c("FF_id", "up_down","Stringency_level")]



# 5) Stringency levels other than Connectivity Map (FF_id 1-12) ----------------
fdcl <- ff[!grepl("_NA_", ff$Stringency_level),] %>% split(., .$Stringency_level)
ffl <- lapply(fdcl, function(x) x[, c("FF_id", "up_down")])

chunk_size <- 500

for (j in names(ffl)) {
  
  i <- ffl[[j]]
  ni <- names(ffl[j])
  
  cat(sprintf("Stringency_level : %s\n", ni))
  
  # Output directory
  base_dir <- "output/Perturbagen_similarity/"
  target_dir <- file.path(base_dir)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("target file path : %s\n", target_dir))
  
  # gtc_ids of the up (larger FF_id) and down (smaller FF_id) signatures
  up_gtc <- gc$gtc_id[gc$FF_id == max(i$FF_id)]
  dn_gtc <- gc$gtc_id[gc$FF_id == min(i$FF_id)]
  cat(sprintf("up_gtc count : %d\ndn_gtc_count : %d\n\n", length(up_gtc), length(dn_gtc)))
  
  if (length(up_gtc) == 0 || length(dn_gtc) == 0) {
    stop(sprintf("'%s': up_gtc or dn_gtc is empty.", ni))
  }
  
  # Stop if the output file already exists
  expected_files_all <- sprintf("Similarity_%s.csv", ni)
  existing_files <- list.files(target_dir)
  already_done <- expected_files_all %in% existing_files
  
  cat(sprintf("'%s' exist : %s ",expected_files_all , already_done))
  if(already_done){
    stop(sprintf("File already exists: %s", expected_files_all))
  }
  
  
  # 6) Most similar signatures, up direction (go_up) -----------------------------
  up_chunks <- split(up_gtc, ceiling(seq_along(up_gtc) / chunk_size))
  if (length(up_gtc) == 0) up_chunks <- list()
  cat(sprintf("total chunks : %d\n\n", length(up_chunks)))

  up_result_list <- list()
  for (ci in seq_along(up_chunks)) {
    gtc_up_chunk <- up_chunks[[ci]]
    
    # Upload the chunk's gtc_ids to a temporary table
    dbExecute(con, "DROP TEMPORARY TABLE IF EXISTS tmp_gtc_ids")
    dbExecute(con, "CREATE TEMPORARY TABLE tmp_gtc_ids (gtc_id BIGINT PRIMARY KEY)")
    dbWriteTable(con, "tmp_gtc_ids", data.frame(gtc_id = gtc_up_chunk), append = TRUE, row.names = FALSE)
    
    cat(sprintf("chunk %d / %d - querying...\n", ci, length(up_chunks)))
    
    # Pairs in which the chunk's signature is sam1 or sam2, with a released partner
    part1 <- dbGetQuery(con, "
      SELECT gu.sam1_id, gu.sam2_id, gu.GO_BP, gu.GO_CC, gu.GO_MF,
             t.gtc_id AS matched_gtc
      FROM go_up gu
      JOIN tmp_gtc_ids t ON gu.sam1_id = t.gtc_id
      JOIN tmp_open_gtc_ids o ON gu.sam2_id = o.gtc_id
    ")
    part2 <- dbGetQuery(con, "
      SELECT gu.sam1_id, gu.sam2_id, gu.GO_BP, gu.GO_CC, gu.GO_MF,
             t.gtc_id AS matched_gtc
      FROM go_up gu
      JOIN tmp_gtc_ids t ON gu.sam2_id = t.gtc_id
      JOIN tmp_open_gtc_ids o ON gu.sam1_id = o.gtc_id
    ")
    combined <- rbind(part1, part2)
    setDT(combined)
    
    # Top 100 per ontology (BP, CC, MF) for each matched signature
    chunk_result <- combined[, {
      bp <- .SD[!is.na(GO_BP)][frank(-GO_BP, ties.method = "min") <= 100]
      cc <- .SD[!is.na(GO_CC)][frank(-GO_CC, ties.method = "min") <= 100]
      mf <- .SD[!is.na(GO_MF)][frank(-GO_MF, ties.method = "min") <= 100]
      rbind(bp, cc, mf)
    }, by = matched_gtc]
    
    up_result_list[[ci]] <- chunk_result
  }
    up_result_df <- do.call(rbind, up_result_list)
    up_result_df <- unique(up_result_df)
    
    # Add direction, GO term counts and Signature_id of both signatures
    up_result_df$direc <- "Up"
    up_result_df$bp_c1 <- gc$BPO_used_go_n[match(up_result_df$sam1_id, gc$gtc_id)]
    up_result_df$cc_c1 <- gc$CCO_used_go_n[match(up_result_df$sam1_id, gc$gtc_id)]
    up_result_df$mf_c1 <- gc$MFO_used_go_n[match(up_result_df$sam1_id, gc$gtc_id)]
    up_result_df$bp_c2 <- gc$BPO_used_go_n[match(up_result_df$sam2_id, gc$gtc_id)]
    up_result_df$cc_c2 <- gc$CCO_used_go_n[match(up_result_df$sam2_id, gc$gtc_id)]
    up_result_df$mf_c2 <- gc$MFO_used_go_n[match(up_result_df$sam2_id, gc$gtc_id)]

    # gtc_id -> cvc_id -> Signature_id
    cvc_matched1 <- gc$cvc_id[match(up_result_df$sam1_id, gc$gtc_id)]
    up_result_df$sam1_sig <- sig$Signature_id[match(cvc_matched1, sig$cvc_id)]
    
    cvc_matched2 <- gc$cvc_id[match(up_result_df$sam2_id, gc$gtc_id)]
    up_result_df$sam2_sig <- sig$Signature_id[match(cvc_matched2, sig$cvc_id)]
    
    
    
    
    
    # 7) Most similar signatures, down direction (go_dn) ---------------------------
    dn_chunks <- split(dn_gtc, ceiling(seq_along(dn_gtc) / chunk_size))
    if (length(dn_gtc) == 0) dn_chunks <- list()
    cat(sprintf("total chunks : %d\n\n", length(dn_chunks)))
    
    dn_result_list <- list()
    for (ci in seq_along(dn_chunks)) {
      gtc_dn_chunk <- dn_chunks[[ci]]
      
      # Upload the chunk's gtc_ids to a temporary table
      dbExecute(con, "DROP TEMPORARY TABLE IF EXISTS tmp_gtc_ids")
      dbExecute(con, "CREATE TEMPORARY TABLE tmp_gtc_ids (gtc_id BIGINT PRIMARY KEY)")
      dbWriteTable(con, "tmp_gtc_ids", data.frame(gtc_id = gtc_dn_chunk), append = TRUE, row.names = FALSE)
      
      cat(sprintf("chunk %d / %d - querying...\n", ci, length(dn_chunks)))
      
      # Pairs in which the chunk's signature is sam1 or sam2, with a released partner
    part1 <- dbGetQuery(con, "
      SELECT gd.sam1_id, gd.sam2_id, gd.GO_BP, gd.GO_CC, gd.GO_MF,
             t.gtc_id AS matched_gtc
      FROM go_dn gd
      JOIN tmp_gtc_ids t ON gd.sam1_id = t.gtc_id
      JOIN tmp_open_gtc_ids o ON gd.sam2_id = o.gtc_id
    ")
    part2 <- dbGetQuery(con, "
      SELECT gd.sam1_id, gd.sam2_id, gd.GO_BP, gd.GO_CC, gd.GO_MF,
             t.gtc_id AS matched_gtc
      FROM go_dn gd
      JOIN tmp_gtc_ids t ON gd.sam2_id = t.gtc_id
      JOIN tmp_open_gtc_ids o ON gd.sam1_id = o.gtc_id
    ")
      
      combined <- rbind(part1, part2)
      setDT(combined)
      
      # Top 100 per ontology (BP, CC, MF) for each matched signature
      chunk_result_dn <- combined[, {
        bp <- .SD[!is.na(GO_BP)][frank(-GO_BP, ties.method = "min") <= 100]
        cc <- .SD[!is.na(GO_CC)][frank(-GO_CC, ties.method = "min") <= 100]
        mf <- .SD[!is.na(GO_MF)][frank(-GO_MF, ties.method = "min") <= 100]
        rbind(bp, cc, mf)
      }, by = matched_gtc]
      
      dn_result_list[[ci]] <- chunk_result_dn
    }
    dn_result_df <- do.call(rbind, dn_result_list)
    dn_result_df <- unique(dn_result_df)
    
    # Add direction, GO term counts and Signature_id of both signatures
    dn_result_df$direc <- "Down"
    dn_result_df$bp_c1 <- gc$BPO_used_go_n[match(dn_result_df$sam1_id, gc$gtc_id)]
    dn_result_df$cc_c1 <- gc$CCO_used_go_n[match(dn_result_df$sam1_id, gc$gtc_id)]
    dn_result_df$mf_c1 <- gc$MFO_used_go_n[match(dn_result_df$sam1_id, gc$gtc_id)]
    dn_result_df$bp_c2 <- gc$BPO_used_go_n[match(dn_result_df$sam2_id, gc$gtc_id)]
    dn_result_df$cc_c2 <- gc$CCO_used_go_n[match(dn_result_df$sam2_id, gc$gtc_id)]
    dn_result_df$mf_c2 <- gc$MFO_used_go_n[match(dn_result_df$sam2_id, gc$gtc_id)]

    # gtc_id -> cvc_id -> Signature_id
    cvc_matched1 <- gc$cvc_id[match(dn_result_df$sam1_id, gc$gtc_id)]
    dn_result_df$sam1_sig <- sig$Signature_id[match(cvc_matched1, sig$cvc_id)]
    
    cvc_matched2 <- gc$cvc_id[match(dn_result_df$sam2_id, gc$gtc_id)]
    dn_result_df$sam2_sig <- sig$Signature_id[match(cvc_matched2, sig$cvc_id)]

    
    # 8) Combine up and down results ----------------------------------------------
    result_df <- rbind(up_result_df, dn_result_df)
    result_df <- unique(result_df)
    
 
    if (!identical(result_df, rbind(up_result_df, dn_result_df))) {
      stop("result_df differs from rbind(up_result_df, dn_result_df).")
    }
    if (any(is.na(result_df$sam1_sig)) | any(is.na(result_df$sam2_sig))) {
      stop("result_df has rows with a missing Signature_id.")
    }
    
    # Select and rename output columns
    fdf <- result_df[, c(
      "sam1_sig",
      "bp_c1",
      "mf_c1",
      "cc_c1",
      "sam2_sig",
      "bp_c2",
      "mf_c2",
      "cc_c2",
      "GO_BP",
      "GO_MF",
      "GO_CC",
      "direc"
    )] %>% as.data.frame()
    
    colnames(fdf) <- c(
      "Signature_id_1",
      "BP_term_count_1",
      "MF_term_count_1",
      "CC_term_count_1",
      "Signature_id_2",
      "BP_term_count_2",
      "MF_term_count_2",
      "CC_term_count_2",
      "BP_score",
      "MF_score",
      "CC_score",
      "Direction"
    )
    
    # Remove duplicated rows
    n_before <- nrow(fdf)
    fdf <- unique(fdf)
    n_after <- nrow(fdf)
    cat(sprintf("rows before unique: %d, after: %d (removed: %d)\n", n_before, n_after, n_before - n_after))
    



    # 9) Check for missing values and write the file -------------------------------
    # All columns except the scores must have no missing values
    check_cols <- setdiff(colnames(fdf), grep("_score", colnames(fdf), value = TRUE))
    
    na_summary <- sapply(fdf[, check_cols], function(x) sum(is.na(x)))
    print(na_summary)
    
    if (any(na_summary > 0)) {
      stop(sprintf("Missing values in columns : %s",
                   paste(names(na_summary[na_summary > 0]), collapse = ", ")))
    }
    
    
    # write csv file
    out_path <- file.path(target_dir, sprintf("Similarity_%s.csv", ni))
    fwrite(fdf, out_path)
    
    cat(sprintf("saved: %s (%d rows)\n\n", out_path, nrow(fdf)))
}
dbDisconnect(con)