# =============================================================================
# Build_ortholog_table.R
#
# Build the rat-to-human gene symbol mapping tables from Ensembl BioMart.
#
# Three BioMart queries are made against the rat gene dataset
# (rnorvegicus_gene_ensembl):
#   1. rat-human orthologs (orthology type, confidence, % identity)
#   2. rat-side annotation (Entrez ID, description)
#   3. rat symbol aliases (former symbols used in older resources)
# They are merged, each rat symbol is classified by its mapping pattern
# (`mapping_class`), and two tables are derived:
#   map_use_all_human_rat       : current rat symbol -> human symbol(s)
#   map_use_all_human_rat_alias : former rat symbol (alias) -> current rat
#                                 symbol -> human symbol(s)
# The alias table covers outdated rat symbols that appear in older resources
# (e.g. microarray annotations) but are no longer the current Ensembl symbol.
# Both tables are the input of convert_to_human() (convert_to_human.R).
#
# Rat only. All orthologs are kept here (one-to-one, one-to-many, many-to-many);
# choosing a single human symbol per rat symbol is done in convert_to_human().
# Queries are made against the live Ensembl BioMart, so the result depends on
# the Ensembl release at the time of the query (queried on 2026.04.09).
#
# Output (written to `out_dir`):
#   MAPPING_orth_tbl_human_rat_final.Rdata    both mapping tables
#
# Dependencies: R 4.4.1 ; biomaRt, dplyr
# =============================================================================

library(biomaRt)
library(dplyr)

############################################################
## 1. Output directory
############################################################
out_dir <- "gene_mapping"   # change as needed
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


############################################################
## 2. Connect to Ensembl
############################################################
mirrors <- c("www.ensembl.org", "useast.ensembl.org", "asia.ensembl.org", "uswest.ensembl.org")
#' Connect to the rat Ensembl mart, trying several mirrors in turn.
#' @param mirrors Character vector of Ensembl hosts, tried in order.
#' @param dataset BioMart dataset name.
#' @param max_try Number of attempts per mirror (waits 5 * attempt seconds
#'   between attempts).
#' @return A Mart object. Stops if every mirror fails.
connect_mart <- function(mirrors, dataset = "rnorvegicus_gene_ensembl", max_try = 3) {
  for (mirror in mirrors) {
    for (attempt in seq_len(max_try)) {
      message(sprintf("[%s] Attempt %d/%d ...", mirror, attempt, max_try))
      mart <- tryCatch(
        useEnsembl(biomart = "genes", dataset = dataset, host = mirror),
        error = function(e) {
          message("  -> Error: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(mart)) {
        message(sprintf("  -> Connected: %s", mirror))
        return(mart)
      }
      Sys.sleep(5 * attempt)
    }
  }
  stop("All mirrors failed. Please try again later.")
}
rat_mart <- connect_mart(mirrors)


############################################################
## 3. getBM wrapper (with retries)
############################################################

#' getBM() with retries.
#' @param attributes Character vector of BioMart attributes.
#' @param mart Mart object from connect_mart().
#' @param max_try Maximum number of attempts.
#' @param wait Seconds to wait between attempts.
#' @return Data frame returned by getBM(). Stops if all attempts fail.
safe_getBM <- function(attributes, mart, max_try = 3, wait = 10) {
  for (attempt in seq_len(max_try)) {
    result <- tryCatch(
      getBM(attributes = attributes, mart = mart),
      error = function(e) {
        message(sprintf("getBM attempt %d failed: %s", attempt, conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(result)) return(result)
    Sys.sleep(wait)
  }
  stop("getBM failed after ", max_try, " attempts.")
}



############################################################
## 4. Query 1: rat-human ortholog table
############################################################
orth_core <- safe_getBM(
  attributes = c(
    "external_gene_name",
    "ensembl_gene_id",
    "hsapiens_homolog_associated_gene_name",
    "hsapiens_homolog_ensembl_gene",
    "hsapiens_homolog_orthology_type",
    "hsapiens_homolog_perc_id",
    "hsapiens_homolog_perc_id_r1",
    "hsapiens_homolog_orthology_confidence"
  ),
  mart = rat_mart
)
colnames(orth_core) <- c(
  "rat_symbol",
  "rat_ensembl_gene_id",
  "human_symbol",
  "human_ensembl_gene_id",
  "orthology_type",
  "perc_id_rat_to_human",
  "perc_id_human_to_rat",
  "orthology_confidence"
)



############################################################
## 5. Query 2: rat-side annotation only
############################################################
rat_annot <- safe_getBM(
  attributes = c(
    "ensembl_gene_id",
    "entrezgene_id",
    "description"
  ),
  mart = rat_mart
)

colnames(rat_annot) <- c(
  "rat_ensembl_gene_id",
  "rat_entrez_id",
  "rat_description"
)

############################################################
## 6. Query 3: rat alias table (for former symbols)
############################################################
rat_alias_raw <- safe_getBM(
  attributes = c(
    "external_gene_name", 
    "external_synonym"     
  ),
  mart = rat_mart
)
colnames(rat_alias_raw) <- c("current_rat_symbol", "alias")

rat_alias <- rat_alias_raw %>%
  dplyr::filter(
    !is.na(alias), alias != "",
    !is.na(current_rat_symbol), current_rat_symbol != ""
  ) %>%
  dplyr::distinct() %>%
  dplyr::filter(alias != current_rat_symbol)



############################################################
## 7. Merge
############################################################
orth_tbl <- orth_core %>%
  left_join(rat_annot, by = "rat_ensembl_gene_id") %>%
  distinct() %>%
  mutate(
    rat_symbol            = trimws(as.character(rat_symbol)),
    rat_ensembl_gene_id   = trimws(as.character(rat_ensembl_gene_id)),
    rat_entrez_id         = suppressWarnings(as.numeric(rat_entrez_id)),
    rat_description       = as.character(rat_description),
    human_symbol          = trimws(as.character(human_symbol)),
    human_ensembl_gene_id = trimws(as.character(human_ensembl_gene_id)),
    orthology_type        = trimws(as.character(orthology_type)),
    perc_id_rat_to_human  = suppressWarnings(as.numeric(perc_id_rat_to_human)),
    perc_id_human_to_rat  = suppressWarnings(as.numeric(perc_id_human_to_rat)),
    orthology_confidence  = suppressWarnings(as.numeric(orthology_confidence))
  ) %>%
  filter(!is.na(rat_symbol), rat_symbol != "")



############################################################
## 8. Summarise mapping status per rat_symbol
##    mapping_class: unmapped / one2one_unique / one2one_plus_others /
##                   non_one2one_but_unique / one_to_many
############################################################
map_summary <- orth_tbl %>%
  group_by(rat_symbol) %>%
  summarise(
    n_rows          = n(),
    n_human_symbol  = n_distinct(human_symbol[!is.na(human_symbol) & human_symbol != ""]),
    n_human_ensembl = n_distinct(human_ensembl_gene_id[!is.na(human_ensembl_gene_id) & human_ensembl_gene_id != ""]),
    has_one2one      = any(orthology_type == "ortholog_one2one", na.rm = TRUE),
    has_human_symbol = any(!is.na(human_symbol) & human_symbol != ""),
    human_symbols_all = paste(sort(unique(human_symbol[!is.na(human_symbol) & human_symbol != ""])), collapse = ";"),
    human_ensembl_all = paste(sort(unique(human_ensembl_gene_id[!is.na(human_ensembl_gene_id) & human_ensembl_gene_id != ""])), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    mapping_class = case_when(
      !has_human_symbol                        ~ "unmapped",
      has_one2one & n_human_symbol == 1        ~ "one2one_unique",
      has_one2one & n_human_symbol > 1         ~ "one2one_plus_others",
      !has_one2one & n_human_symbol == 1       ~ "non_one2one_but_unique",
      n_human_symbol > 1                       ~ "one_to_many",
      TRUE                                     ~ "other"
    )
  )

orth_tbl_annot <- orth_tbl %>%
  left_join(map_summary, by = "rat_symbol")



############################################################
## 9. Derived mapping tables
############################################################

# (A) all mappings preserved
map_use_all <- orth_tbl_annot %>%
  dplyr::filter(!is.na(human_symbol), human_symbol != "") %>%
  dplyr::select(
    rat_symbol, human_symbol, human_ensembl_gene_id,
    orthology_type, orthology_confidence,
    perc_id_rat_to_human, perc_id_human_to_rat,
    mapping_class, n_human_symbol, human_symbols_all
  ) %>%
  dplyr::distinct()

# (B) former rat symbol (alias) -> current_rat_symbol -> human symbol(s)
map_use_all_with_alias <- rat_alias %>%
  dplyr::inner_join(map_use_all, by = c("current_rat_symbol" = "rat_symbol")) %>%
  dplyr::select(
    alias,                  
    current_rat_symbol,    
    human_symbol, human_ensembl_gene_id,
    orthology_type, orthology_confidence,
    perc_id_rat_to_human, perc_id_human_to_rat,
    mapping_class, n_human_symbol, human_symbols_all
  ) %>%
  dplyr::distinct()


############################################################
## 10. Save
############################################################
map_use_all_human_rat         <- map_use_all
map_use_all_human_rat_alias   <- map_use_all_with_alias

save(map_use_all_human_rat, map_use_all_human_rat_alias,
     file = file.path(out_dir, "MAPPING_orth_tbl_human_rat_final.Rdata"))
