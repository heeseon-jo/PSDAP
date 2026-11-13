# =============================================================================
# convert_to_human.R
#
# Convert rat gene symbols to human gene symbols.
# Each symbol is resolved in three steps; a symbol is only passed to the next
# step if the previous one found no human symbol:
#   1. direct match to the current rat symbol   (mapping)
#   2. match to a former rat symbol (alias)     (mapping_alias)
#   3. NCBI HomoloGene                          (homologene)
#
# When a rat symbol has several human orthologs in steps 1-2, one is chosen by,
# in order: one-to-one orthology first, then higher orthology confidence, then
# higher % identity (rat to human). Ties are broken by row order of the mapping
# table, so use the same mapping tables to reproduce the result.
# HomoloGene (step 3) has not been updated since 2014 and is used only as a
# last resort; the first returned human symbol is taken.
#
# Rat only. Dataset-specific preprocessing (e.g. removal of platform-specific
# suffixes from symbols) must be done before calling this function.
# The mapping tables are produced by build_ortholog_table.R
# (`map_use_all_human_rat`, `map_use_all_human_rat_alias`).
#
# Dependencies: R 4.4.1; dplyr, homologene
# =============================================================================

library(dplyr)
library(homologene)


#' Convert rat gene symbols to human gene symbols.
#'
#' @param symbols Character vector of rat gene symbols, already cleaned of
#'   dataset-specific suffixes or prefixes.
#' @param mapping `map_use_all_human_rat` from build_ortholog_table.R
#'   (rat_symbol, human_symbol, orthology_type, orthology_confidence,
#'   perc_id_rat_to_human, mapping_class, ...).
#' @param mapping_alias `map_use_all_human_rat_alias` from
#'   build_ortholog_table.R (same columns, with `alias` in place of
#'   `rat_symbol`).
#' @param in_tax NCBI taxonomy ID of the input species, passed to HomoloGene
#'   in step 3. Default 10116 (rat); the mapping tables are rat-specific.
#' @return Data frame with one row per input symbol: `orig_symbol`,
#'   `human_symbol` (NA if unmapped) and `mapping_class`. `mapping_class` is
#'   the class given in the mapping tables, or "homologene" for symbols
#'   resolved in step 3.
convert_to_human <- function(symbols, mapping, mapping_alias, in_tax = 10116) {
  
  sym_col <- if (in_tax == 10116) "rat_symbol" else "mouse_symbol"
  
  result <- data.frame(
    orig_symbol   = symbols,
    human_symbol  = NA_character_,
    mapping_class = NA_character_,
    stringsAsFactors = FALSE
  )
  
  m1 <- mapping %>%
    dplyr::filter(.data[[sym_col]] %in% symbols) %>%
    dplyr::mutate(priority = ifelse(orthology_type == "ortholog_one2one", 1L, 0L)) %>%
    dplyr::group_by(.data[[sym_col]]) %>%
    dplyr::arrange(desc(priority), desc(orthology_confidence),
                   desc(perc_id_rat_to_human), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(orig_symbol = !!sym_col, human_symbol, mapping_class)
  
  result <- result %>%
    dplyr::left_join(m1, by = "orig_symbol") %>%
    dplyr::mutate(
      human_symbol  = dplyr::coalesce(human_symbol.y,  human_symbol.x),
      mapping_class = dplyr::coalesce(mapping_class.y, mapping_class.x)
    ) %>%
    dplyr::select(orig_symbol, human_symbol, mapping_class)
  
  still_na <- is.na(result$human_symbol)
  if (any(still_na)) {
    alias_syms <- result$orig_symbol[still_na]
    m2 <- mapping_alias %>%
      dplyr::filter(alias %in% alias_syms) %>%
      dplyr::mutate(priority = ifelse(orthology_type == "ortholog_one2one", 1L, 0L)) %>%
      dplyr::group_by(alias) %>%
      dplyr::arrange(desc(priority), desc(orthology_confidence),
                     desc(perc_id_rat_to_human), .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(orig_symbol = alias, human_symbol, mapping_class)
    
    result <- result %>%
      dplyr::left_join(m2, by = "orig_symbol") %>%
      dplyr::mutate(
        human_symbol  = dplyr::coalesce(human_symbol.y,  human_symbol.x),
        mapping_class = dplyr::coalesce(mapping_class.y, mapping_class.x)
      ) %>%
      dplyr::select(orig_symbol, human_symbol, mapping_class)
  }
  
  still_na <- is.na(result$human_symbol)
  if (any(still_na)) {
    hom_syms <- result$orig_symbol[still_na]
    hom_syms <- hom_syms[
      !grepl("^RGD\\d+|^LOC\\d+|^MGC\\d+|^AABR\\d+|^Gm\\d+", hom_syms) &
        !is.na(hom_syms) & nchar(hom_syms) > 2
    ]
    if (length(hom_syms) > 0) {
      hom_res <- tryCatch(
        homologene::homologene(genes = hom_syms, inTax = in_tax, outTax = 9606),
        error = function(e) NULL
      )
      if (!is.null(hom_res) && nrow(hom_res) > 0) {
        m3 <- hom_res %>%
          dplyr::select(orig_symbol = !!as.character(in_tax), human_symbol = "9606") %>%
          dplyr::filter(!is.na(human_symbol), human_symbol != "") %>%
          dplyr::distinct() %>%
          dplyr::group_by(orig_symbol) %>%
          dplyr::slice(1) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(mapping_class = "homologene")
        
        result <- result %>%
          dplyr::left_join(m3, by = "orig_symbol") %>%
          dplyr::mutate(
            human_symbol  = dplyr::coalesce(human_symbol.y,  human_symbol.x),
            mapping_class = dplyr::coalesce(mapping_class.y, mapping_class.x)
          ) %>%
          dplyr::select(orig_symbol, human_symbol, mapping_class)
      }
    }
  }
  
  return(result)
}