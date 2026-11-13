# =============================================================================
# standardize_model_info.R
#
# Reference-identifier lookup for experimental model names.
#   In vitro models -> Cellosaurus accession
#   In vivo models  -> UBERON term
#
# Not run by this script (results are in the released metadata):
#   - manual selection of the parental cell line when several Cellosaurus
#     entries matched a name (derivatives such as Cas9, GFP or resistant
#     sublines were not chosen)
#   - manual pre-classification of in vivo names (cell types vs tissues) and
#     name corrections before lookup
#   - manual review of names without an identifier (identifier assigned, or
#     the reason recorded, e.g. patient-derived or primary specimen)
#
# Cellosaurus and UBERON (EBI OLS) web services are queried; results depend on
# the query date. Model names are read from the internal PSDAP database, which
# is not distributed. Credentials are read from environment variables
# (DB_HOST, DB_PORT, DB_NAME, DB_USER, DB_PW).
#
# Output: output/model_identifier_mapping.csv
#
# Dependencies: data.table, DBI, RMariaDB, httr, jsonlite.
# =============================================================================

library(data.table)
library(DBI)
library(httr)
library(jsonlite)

out_file <- "output/model_identifier_mapping.csv"

# Input: model names (internal database) --------------------------------------------
con <- dbConnect(RMariaDB::MariaDB(),
                 host = Sys.getenv("DB_HOST"), port = as.integer(Sys.getenv("DB_PORT")),
                 dbname = Sys.getenv("DB_NAME"),
                 user = Sys.getenv("DB_USER"), password = Sys.getenv("DB_PW"))
invitro_names <- unique(dbGetQuery(con, "SELECT name FROM option_cell_line WHERE omt_id IN (1, 3, 5)")$name)
invivo_names  <- unique(dbGetQuery(con, "SELECT name FROM option_cell_line WHERE omt_id IN (2, 4, 6)")$name)
dbDisconnect(con)

# STEP 1. In vitro: name -> Cellosaurus accession ---------------------------------------
# A name is mapped only when it equals the Cellosaurus identifier, or else one of
# its synonyms, after normalization (lower case; "-", "_", "/", ".", spaces removed).
normalize_cell_name <- function(x) tolower(gsub("[-_/. ]", "", x))

get_cvcl_by_name <- function(name) {
  url <- sprintf("https://api.cellosaurus.org/search/cell-line?q=%s&fields=id,ac,sy&format=json&rows=1000",
                 URLencode(name, reserved = TRUE))
  res <- tryCatch(GET(url), error = function(e) NULL)
  if (is.null(res) || status_code(res) != 200) return(c(ac = NA_character_, cvcl_id = NA_character_))
  cells <- fromJSON(content(res, as = "text", encoding = "UTF-8"))$Cellosaurus$`cell-line-list`
  if (is.null(cells) || length(cells) == 0) return(c(ac = NA_character_, cvcl_id = NA_character_))
  ac_list        <- sapply(cells$`accession-list`, function(x) x$value[x$type == "primary"][1])
  name_list      <- sapply(cells$`name-list`,      function(x) x$value[x$type == "identifier"][1])
  all_names_list <- lapply(cells$`name-list`,      function(x) x$value)
  # 1. identifier
  exact <- which(normalize_cell_name(name_list) == normalize_cell_name(name))
  # 2. synonym
  if (length(exact) == 0)
    exact <- which(sapply(all_names_list, function(nms) any(normalize_cell_name(nms) == normalize_cell_name(name))))
  if (length(exact) == 0) return(c(ac = NA_character_, cvcl_id = NA_character_))
  idx <- exact[1]
  c(ac = ac_list[idx], cvcl_id = name_list[idx])
}

get_cvcl_details <- function(ac) {
  url <- sprintf("https://api.cellosaurus.org/cell-line/%s?fields=id,ac,ca,di,derived-from-site,ox,sy&format=json", ac)
  empty <- c(ca = NA_character_, disease = NA_character_, tissue = NA_character_,
             species = NA_character_, synonyms = NA_character_)
  res <- tryCatch(GET(url), error = function(e) NULL)
  if (is.null(res) || status_code(res) != 200) return(empty)
  cell <- tryCatch(fromJSON(content(res, as = "text", encoding = "UTF-8"))$Cellosaurus$`cell-line-list`,
                   error = function(e) NULL)
  if (is.null(cell)) return(empty)
  ca <- tryCatch(cell$category[1], error = function(e) NA_character_)
  disease <- tryCatch({                                      # NCIt label preferred
    di <- cell$`disease-list`[[1]]
    if (!is.null(di) && "label" %in% names(di)) {
      ncit <- di[di$database == "NCIt", ]
      if (nrow(ncit) > 0) ncit$label[1] else di$label[1]
    } else NA_character_
  }, error = function(e) NA_character_)
  tissue <- tryCatch({
    site <- cell$`derived-from-site-list`[[1]]
    if (!is.null(site) && "site" %in% names(site)) site$site$value[1] else NA_character_
  }, error = function(e) NA_character_)
  species <- tryCatch({
    ox <- cell$`species-list`[[1]]
    if (!is.null(ox) && "label" %in% names(ox)) ox$label[1] else NA_character_
  }, error = function(e) NA_character_)
  synonyms <- tryCatch({
    sy <- cell$`name-list`[[1]]
    syns <- sy$value[sy$type == "synonym"]
    if (length(syns) > 0) paste(syns, collapse = "; ") else NA_character_
  }, error = function(e) NA_character_)
  c(ca = ca, disease = disease, tissue = tissue, species = species, synonyms = synonyms)
}

cvcl_map <- rbindlist(lapply(invitro_names, function(nm) {
  Sys.sleep(0.3)
  info <- get_cvcl_by_name(nm)
  data.table(name = nm, cvcl_ac = info[["ac"]], cvcl_id = info[["cvcl_id"]])
}))
details <- rbindlist(lapply(cvcl_map[!is.na(cvcl_ac), unique(cvcl_ac)], function(ac) {
  Sys.sleep(0.3)
  as.data.table(as.list(c(cvcl_ac = ac, get_cvcl_details(ac))))
}))
cvcl_map <- merge(cvcl_map, details, by = "cvcl_ac", all.x = TRUE, sort = FALSE)

# STEP 2. In vitro: merge names by accession ----------------------------------------
# Standard name: Cellosaurus identifier; names without an accession keep their own name.
cvcl_map[, model_key := fifelse(!is.na(cvcl_ac), cvcl_ac, name)]
invitro <- cvcl_map[, .(
  model_system  = "in vitro",
  standard_id   = cvcl_ac[!is.na(cvcl_ac)][1],
  standard_name = if (any(!is.na(cvcl_id))) cvcl_id[!is.na(cvcl_id)][1] else name[1],
  ca = ca[1], disease = disease[1], tissue = tissue[1], species = species[1], synonyms = synonyms[1],
  names = paste(sort(unique(name)), collapse = "; ")
), by = model_key]

# STEP 3. In vivo: name -> UBERON term ------------------------------------------------
# A name is mapped when it equals a UBERON label, or else a UBERON synonym
# (case-insensitive).
get_uberon <- function(term) {
  url <- sprintf("https://www.ebi.ac.uk/ols/api/search?q=%s&ontology=uberon&fieldList=iri,label,short_form,synonym,description&rows=20",
                 URLencode(tolower(term), reserved = TRUE))
  empty <- c(uberon_id = NA_character_, uberon_label = NA_character_, synonyms = NA_character_)
  res <- tryCatch(GET(url), error = function(e) NULL)
  if (is.null(res) || status_code(res) != 200) return(empty)
  docs <- tryCatch(fromJSON(content(res, as = "text", encoding = "UTF-8"))$response$docs, error = function(e) NULL)
  if (is.null(docs) || !is.data.frame(docs) || nrow(docs) == 0) return(empty)
  exact <- which(tolower(docs$label) == tolower(term))
  if (length(exact) == 0)
    exact <- which(sapply(docs$synonym, function(s) length(s) > 0 && any(tolower(unlist(s)) == tolower(term))))
  if (length(exact) == 0) return(empty)
  idx <- exact[1]
  c(uberon_id = docs$short_form[idx], uberon_label = docs$label[idx],
    synonyms = tryCatch(paste(unlist(docs$synonym[[idx]]), collapse = "; "), error = function(e) NA_character_))
}

uberon_map <- rbindlist(lapply(invivo_names, function(nm) {
  Sys.sleep(0.3)
  info <- get_uberon(nm)
  data.table(name = nm, uberon_id = info[["uberon_id"]], uberon_label = info[["uberon_label"]],
             synonyms = info[["synonyms"]])
}))

# STEP 4. In vivo: merge names by UBERON term -----------------------------------------
# Standard name: UBERON label; names without a term keep their own name.
uberon_map[, model_key := fifelse(!is.na(uberon_id), uberon_id, name)]
invivo <- uberon_map[, .(
  model_system  = "in vivo",
  standard_id   = uberon_id[!is.na(uberon_id)][1],
  standard_name = if (any(!is.na(uberon_label))) uberon_label[!is.na(uberon_label)][1] else name[1],
  ca = if (any(!is.na(uberon_id))) "tissue" else NA_character_,
  disease = NA_character_,
  tissue  = if (any(!is.na(uberon_label))) uberon_label[!is.na(uberon_label)][1] else NA_character_,
  species = NA_character_, synonyms = synonyms[!is.na(synonyms)][1],
  names = paste(sort(unique(name)), collapse = "; ")
), by = model_key]

# Output --------------------------------------------------------------------------------
models <- rbind(invitro, invivo)
dir.create("output", showWarnings = FALSE)
fwrite(models, out_file, na = "", bom = TRUE)