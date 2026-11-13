# =============================================================================
# KS_enrichment.R
#
# KS-test gene-set enrichment of a differential expression signature against
# a gene set collection (e.g. pathway or disease gene sets).
# Genes are scored by sign(log2FC) * -log10(P) and each gene set is tested
# against all other genes in the signature, separately for UP and DOWN.
# Gene sets with fewer than MIN_GENE genes in the signature are not tested.
# Deterministic; no random seed required.
#
# Dependencies: R 4.4.1; readxl, parallel, data.table
# =============================================================================


library(parallel)
library(readxl)
library(data.table)

MIN_GENE <- 5




#' Flatten a nested gene-set list (database > gene set > genes).
#' Database and gene-set names are kept as separate vectors
#' @param nested Named list of databases. Each element is a named list of
#'        gene-symbol vectors.
#' @return List with `sets` (gene vectors), `db` and `name`.
make_gs <- function(nested) {
  stopifnot(is.list(nested), !is.null(names(nested)),
            all(vapply(nested, is.list, logical(1))))
  list(
    sets = unlist(nested, recursive = FALSE, use.names = FALSE),
    db   = rep(names(nested), lengths(nested)),
    name = unlist(lapply(nested, names), use.names = FALSE)
  )
}



#' Read a DE signature and convert it to a named score vector.
#' @param path Path to an .xlsx file (Symbol, Log2FoldChange, P_value).
#' @return Named numeric vector (names = gene symbols), one entry per symbol.
read_score <- function(path) {
  df <- as.data.frame(read_excel(path, sheet = 1))
  
  need <- c("Symbol", "Log2FoldChange", "P_value")
  if (!all(need %in% names(df))) {
    stop(sprintf("No such column [%s]: %s", basename(path),
                 paste(setdiff(need, names(df)), collapse = ", ")))
  }
  
  sym <- df$Symbol; lfc <- df$Log2FoldChange; pv <- df$P_value
  ok  <- !is.na(sym) & nzchar(sym) & is.finite(lfc) & is.finite(pv)
  sym <- sym[ok]; lfc <- lfc[ok]; pv <- pv[ok]
  
  s <- sign(lfc) * -log10(pmax(pv, .Machine$double.xmin))
  
  o   <- order(-abs(s)) # duplicated symbols: keep larger |score|
  sym <- sym[o]; s <- s[o]
  dup <- duplicated(sym)
  setNames(s[!dup], sym[!dup])
}



#' KS test of every gene set against the background, both directions.
#' @param sc Named score vector from read_score().
#' @param g  Gene-set collection from make_gs().
#' @return Data frame, one row per gene set: n_genes, n_bg, D and P-value
#'   for UP and DOWN (NA if fewer than MIN_GENE genes are in the signature).
run_ks <- function(sc, g) {
  
  N <- length(sc)
  
  res <- lapply(seq_along(g$sets), function(k) {
    hit   <- names(sc) %in% g$sets[[k]]
    in_s  <- sc[ hit]
    out_s <- sc[!hit]
    
    if (length(in_s) < MIN_GENE || length(out_s) < 1) {
      return(c(n_genes = length(in_s),
               D_up = NA_real_, pval_up = NA_real_,
               D_down = NA_real_, pval_down = NA_real_))
    }
    r_up <- ks.test(in_s, out_s, alternative = "less")
    r_dn <- ks.test(in_s, out_s, alternative = "greater")
    
    c(n_genes = length(in_s),
      D_up    = as.numeric(r_up$statistic), pval_up   = r_up$p.value,
      D_down  = as.numeric(r_dn$statistic), pval_down = r_dn$p.value)
  })
  
  m  <- do.call(rbind, res)
  df <- data.frame(
    Gene_set_name = g$name,
    Database      = g$db,
    n_genes       = as.integer(m[, "n_genes"]),
    n_bg          = N,
    D_up          = m[, "D_up"],
    pval_up       = m[, "pval_up"],
    D_down        = m[, "D_down"],
    pval_down     = m[, "pval_down"],
    stringsAsFactors = FALSE
  )
  rownames(df) <- NULL
  df
}



## Empty result with the processed-table columns
EMPTY_PRC <- data.frame(
  Gene_set_name = character(), Database = character(),
  D_statistic = numeric(), P_value = numeric(),
  Adjusted_P_value = numeric(), Direction = character(),
  stringsAsFactors = FALSE
)



#' Convert raw KS results to the processed table.
#' BH is applied per direction; each direction is judged independently
#' (a gene set significant in both directions appears as UP and DOWN rows).
#' @param df Output of run_ks().
#' @param padj_cutoff Adjusted P-value cutoff.
#' @return Data frame of significant gene sets (empty if none).
make_processed <- function(df, padj_cutoff = 0.05) {
  
  padj_up   <- p.adjust(df$pval_up,   method = "BH")
  padj_down <- p.adjust(df$pval_down, method = "BH")
  
  ok_up <- !is.na(padj_up)   & padj_up   <= padj_cutoff
  ok_dn <- !is.na(padj_down) & padj_down <= padj_cutoff
  if (!any(ok_up) && !any(ok_dn)) return(EMPTY_PRC)
  
  mk <- function(idx, D, P, PA, dir) {
    if (!any(idx)) return(EMPTY_PRC)
    data.frame(
      Gene_set_name    = df$Gene_set_name[idx],
      Database         = df$Database[idx],
      D_statistic      = D[idx],
      P_value          = P[idx],
      Adjusted_P_value = PA[idx],
      Direction        = dir,
      stringsAsFactors = FALSE
    )
  }
  
  out <- rbind(
    mk(ok_up, df$D_up,   df$pval_up,   padj_up,   "UP"),
    mk(ok_dn, df$D_down, df$pval_down, padj_down, "DOWN")
  )
  
  out <- out[order(out$Adjusted_P_value, -out$D_statistic), ]
  rownames(out) <- NULL
  out
}