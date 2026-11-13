# =============================================================================
# check_metadata_integrity.R
#
# Automated integrity check of the four PSDAP metadata tables
# (signature_info, perturbagen_info, model_info, source_info).
#
# Checks, each reported as the number of problem records:
#   1. Table format: unnamed columns and fully duplicated rows.
#   2. Primary IDs (Signature_id, Pertur_id, Model_id, Source_id): none empty,
#      none duplicated within their own table.
#   3. Referential integrity: every Pertur_id, Model_id and Source_id used in
#      signature_info is defined in perturbagen_info, model_info and
#      source_info, respectively.
#   4. Consistency: the Source name stored in signature_info matches the name
#      given for its Source_id in source_info.
#   5. Unused records (note only): perturbagens, models and sources that no
#      signature refers to.
# The check passes when all "error" rows report 0 problems. Deterministic.
#
# Input : the four CSV files in `meta_dir` (UTF-8, with or without BOM).
# Output: data frame (check, problems, type), printed with a PASS / FAIL line;
#         returned invisibly.
#
# Usage : Rscript check_metadata_integrity.R
#         (the metadata directory is set in the last line of this script)
#
# Dependencies: base R only.
# =============================================================================


# metadata file of each table and its primary ID column
meta_files <- c(signature   = "signature_info.csv",
                perturbagen = "perturbagen_info.csv",
                model       = "model_info.csv",
                source      = "source_info.csv")

meta_keys  <- c(signature   = "Signature_id",
                perturbagen = "Pertur_id",
                model       = "Model_id",
                source      = "Source_id")

# ID columns of signature_info that point to another table
meta_fk    <- c(perturbagen = "Pertur_id",
                model       = "Model_id",
                source      = "Source_id")


#' Flag missing values (NA or empty string).
#'
#' @param x Character vector.
#' @return Logical vector, TRUE where the value is missing.
is_missing <- function(x) {
  is.na(x) | trimws(x) == ""
}


#' Check the integrity of the metadata tables.
#'
#' @param meta_dir Directory containing signature_info.csv, perturbagen_info.csv,
#'   model_info.csv and source_info.csv.
#' @return Data frame with one row per check: `check` (description), `problems`
#'   (number of problem records) and `type` ("error" or "note"). Returned
#'   invisibly; the table and an overall PASS / FAIL line are printed.
check_metadata_integrity <- function(meta_dir) {
  
  tabs <- lapply(meta_files, function(f) {
    read.csv(file.path(meta_dir, f), colClasses = "character",
             fileEncoding = "UTF-8-BOM", check.names = FALSE)
  })
  names(tabs) <- names(meta_files)
  
  res <- list()
  add <- function(check, n, type = "error") {
    res[[length(res) + 1]] <<- data.frame(check = check, problems = n,
                                          type = type, stringsAsFactors = FALSE)
  }
  
  # 1-2. format and primary IDs
  for (nm in names(tabs)) {
    d   <- tabs[[nm]]
    key <- meta_keys[[nm]]
    stopifnot(key %in% names(d))
    add(paste0(nm, ": unnamed columns"),    sum(is.na(names(d)) | names(d) == ""))
    add(paste0(nm, ": duplicated rows"),    sum(duplicated(d)))
    add(paste0(nm, ": empty ", key),        sum(is_missing(d[[key]])))
    add(paste0(nm, ": duplicated ", key),   sum(duplicated(d[[key]])))
  }
  
  # 3 and 5. signature_info -> parent tables
  sig <- tabs$signature
  for (parent in names(meta_fk)) {
    col <- meta_fk[[parent]]
    ids <- tabs[[parent]][[col]]
    add(paste0("signature: ", col, " not defined in ", parent, "_info"),
        sum(!(sig[[col]] %in% ids)))
    add(paste0(parent, ": records not used by any signature"),
        sum(!(ids %in% sig[[col]])), type = "note")
  }
  
  # 4. Source name in signature_info vs source_info
  src_name <- tabs$source$Source[match(sig$Source_id, tabs$source$Source_id)]
  add("signature: Source differs from source_info",
      sum(is.na(src_name) | is.na(sig$Source) | sig$Source != src_name))
  
  res <- do.call(rbind, res)
  print(res, row.names = FALSE)
  
  ok <- all(res$problems[res$type == "error"] == 0)
  cat(sprintf("\nSignatures: %d | Perturbagens: %d | Models: %d | Sources: %d\n",
              nrow(tabs$signature), nrow(tabs$perturbagen),
              nrow(tabs$model), nrow(tabs$source)))
  cat(if (ok) "Integrity check: PASS\n" else "Integrity check: FAIL\n")
  
  invisible(res)
}



check_metadata_integrity("output/")
