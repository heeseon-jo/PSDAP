# =============================================================================
# hypergeometric_test.R
#
# Hypergeometric testing of up- and down-regulated gene sets against a gene
# set collection (e.g. pathway or disease gene sets).
#
# The background gene universe is supplied by the caller and is defined per
# data source as the genes measured in that resource; gene sets are restricted
# to this universe before testing. Deterministic; no random seed required.
#
# Dependencies: base R only (stats, parallel).
# =============================================================================


#' Test one gene set against a gene set collection.
#'
#' @param query      Character vector of query genes.
#' @param gene_sets  Named list, already restricted to `universe`.
#' @param universe   Character vector of background genes.
#' @return data.frame with one row per gene set, including BH-adjusted p-values.
hypergeometric_test <- function(query, gene_sets, universe) {
  
  set_sizes <- lengths(gene_sets)
  overlaps  <- vapply(gene_sets,
                      function(g) length(intersect(query, g)),
                      numeric(1L))
  
  p_value <- phyper(q = overlaps - 1L,
                    m = set_sizes,
                    n = length(universe) - set_sizes,
                    k = length(query),
                    lower.tail = FALSE)
  
  data.frame(
    pathway    = names(gene_sets),
    p.value    = p_value,
    gene_ratio = ifelse(set_sizes > 0L, overlaps / set_sizes, 0),
    count      = overlaps,
    padj_BH    = p.adjust(p_value, method = "BH"),
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
}


#' Run the test over many signatures from one data source.
#'
#' Gene sets are restricted to the universe once, then each signature is tested
#' in the up and down directions. One .rds file is written per signature;
#' existing files are skipped so an interrupted run can be resumed.
#'
#' @param up_gene_list,dn_gene_list  Named lists of character vectors, with
#'   matching names identifying each signature.
#' @param gene_sets  Unrestricted gene set collection (named list).
#' @param universe   Genes measured in this data source.
#' @param out_dir    Output directory for per-signature .rds files.
#' @param n_cores    Number of worker processes.
run_hypergeometric <- function(up_gene_list, dn_gene_list,
                               gene_sets, universe, out_dir,
                               n_cores = max(1L, parallel::detectCores() - 1L)) {
  
  stopifnot(identical(names(up_gene_list), names(dn_gene_list)))
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(out_dir, "run.log")
  
  universe  <- unique(as.character(universe))
  gene_sets <- lapply(gene_sets, function(g) {
    intersect(unique(as.character(g)), universe)
  })
  
  cat(sprintf("[%s] start | signatures: %d | gene sets: %d | universe: %d\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              length(up_gene_list), length(gene_sets), length(universe)),
      file = log_file, append = TRUE)
  
  run_one <- function(sig_id) {
    out_file <- file.path(out_dir, paste0(sig_id, ".rds"))
    if (file.exists(out_file)) return(out_file)
    
    res <- list(
      up = hypergeometric_test(up_gene_list[[sig_id]], gene_sets, universe),
      dn = hypergeometric_test(dn_gene_list[[sig_id]], gene_sets, universe)
    )
    saveRDS(res, file = out_file)
    
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sig_id),
        file = log_file, append = TRUE)
    out_file
  }
  
  ids <- names(up_gene_list)
  
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      varlist = c("up_gene_list", "dn_gene_list", "gene_sets", "universe",
                  "out_dir", "log_file", "hypergeometric_test"),
      envir   = environment()
    )
    files <- unlist(parallel::parLapply(cl, ids, run_one))
  } else {
    files <- unlist(lapply(ids, run_one))
  }
  
  cat(sprintf("[%s] done | %d files\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), length(files)),
      file = log_file, append = TRUE)
  
  invisible(files)
}