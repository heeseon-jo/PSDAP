# =============================================================================
# standardize_perturbagen_info.R
#
# InChIKey lookup for perturbagen names (bioactive compounds, herbal ingredients).
#
# Not run by this script (results are in the released metadata):
#   - further matching of unmatched herbal ingredients: PubChem synonyms with
#     name variants, PubChem PUG-REST and ChEMBL web searches
#   - manual review of names without an InChIKey
#   - Latin names of medicinal herbs and names of herbal formulas
#
# Inputs (third-party files, not included in this repository; place in data/):
#   PubChem      https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/
#                CID-Synonym-filtered, CID-InChI-Key, CID-Title.gz
#   CMap (GEO)   GSE92742_Broad_LINCS_pert_info.txt, GSE70138_Broad_LINCS_pert_info.txt
#   CIGS         supplementary table 41592_2025_2781_MOESM3_ESM.xlsx
#   GtoPdb       https://www.guidetopharmacology.org/DATA/ligands.csv
# Perturbagen names are read from the internal PSDAP database, which is not
# distributed. Credentials are read from environment variables (DB_HOST,
# DB_PORT, DB_NAME, DB_USER, DB_PW).
#
# Output: output/perturbagen_inchikey_mapping.csv
#   name, category, CID, InChIKey, match_source, perturbagen_key, standard_name
#
# Dependencies: data.table, DBI, RMariaDB, readxl.
# =============================================================================

library(data.table)
library(DBI)

pubchem_dir <- "data/pubchem"
cmap_files  <- c("data/GSE92742_Broad_LINCS_pert_info.txt",
                 "data/GSE70138_Broad_LINCS_pert_info.txt")
cigs_file   <- "data/41592_2025_2781_MOESM3_ESM.xlsx"
gtopdb_file <- "data/ligands.csv"
out_file    <- "output/perturbagen_inchikey_mapping.csv"

lower  <- function(x) tolower(trimws(x))
squish <- function(x) trimws(gsub("\\s+", " ", gsub("[\r\n\t]+", " ", as.character(x))))

# Input: perturbagen names (internal database) ---------------------------------------
con <- dbConnect(RMariaDB::MariaDB(),
                 host = Sys.getenv("DB_HOST"), port = as.integer(Sys.getenv("DB_PORT")),
                 dbname = Sys.getenv("DB_NAME"),
                 user = Sys.getenv("DB_USER"), password = Sys.getenv("DB_PW"))
pert <- setDT(dbGetQuery(con, "
  SELECT op.pertur_id, op.name, op.odn_id, odn.db_name, odn.opt_id, op.exclude_from_open
  FROM option_perturbagen op JOIN option_db_name odn ON op.odn_id = odn.odn_id"))
dbDisconnect(con)

# Input: PubChem files -------------------------------------------------------------------
synonym <- fread(file.path(pubchem_dir, "CID-Synonym-filtered"), header = FALSE, sep = "\t",
                 col.names = c("CID", "synonym"), colClasses = "character", quote = "")
n_col    <- ncol(fread(file.path(pubchem_dir, "CID-InChI-Key"), header = FALSE, sep = "\t", nrows = 1))
inchikey <- fread(file.path(pubchem_dir, "CID-InChI-Key"), header = FALSE, sep = "\t",
                  select = c(1, n_col), col.names = c("CID", "InChIKey"),
                  colClasses = "character", quote = "")
inchikey <- inchikey[CID != "" & InChIKey != ""]

# CID-Title.gz is filtered with awk to the queried CIDs or names
title_by_cid <- function(cids) {
  f <- tempfile(); writeLines(unique(cids), f)
  fread(cmd = sprintf("gzip -dc %s | awk -F '\\t' 'NR==FNR {a[$1]=1; next} ($1 in a) {print $1 \"\\t\" $2}' %s -",
                      shQuote(file.path(pubchem_dir, "CID-Title.gz")), shQuote(f)),
        sep = "\t", header = FALSE, col.names = c("CID", "title"), colClasses = "character", quote = "")
}
title_by_name <- function(keys) {
  f <- tempfile(); writeLines(unique(keys), f)
  fread(cmd = sprintf("gzip -dc %s | awk -F '\\t' 'NR==FNR {a[$1]=1; next} {k=tolower($2); gsub(/^[ \\t]+|[ \\t]+$/, \"\", k); if (k in a) print k \"\\t\" $1}' %s -",
                      shQuote(file.path(pubchem_dir, "CID-Title.gz")), shQuote(f)),
        sep = "\t", header = FALSE, col.names = c("name_key", "CID"), colClasses = "character", quote = "")
}

# STEP 1. Bioactive compounds: name -> InChIKey ------------------------------------
# Resources are tried in order; a name keeps the first InChIKey found.
cp <- unique(pert[opt_id == 1, .(name, odn_id, db_name)])
cp_names <- unique(cp$name)

# 1-1. PubChem synonym (first CID per synonym)
synonym[, synonym_lower := lower(synonym)]
hit <- synonym[data.table(name = cp_names, synonym_lower = lower(cp_names)),
               on = "synonym_lower", .(name, CID), mult = "first", nomatch = 0]
m_pubchem <- inchikey[hit, on = "CID", .(name, CID, InChIKey), nomatch = 0]

# 1-2. CMap compound information from GEO (CMap names only): pert_iname, then pert_id
cmap <- unique(rbindlist(lapply(cmap_files, fread), fill = TRUE)[
  pert_type == "trt_cp" & !inchi_key %in% c("", "-666"), .(pert_id, pert_iname, inchi_key)])
cmap_names <- unique(cp[odn_id == 15, name])                 # Connectivity Map
m_cmap <- rbind(
  data.table(name = cmap_names)[, InChIKey := cmap$inchi_key[match(lower(name), lower(cmap$pert_iname))]],
  data.table(name = cmap_names)[!lower(name) %in% lower(cmap$pert_iname)][
    , InChIKey := cmap$inchi_key[match(name, cmap$pert_id)]]
)[!is.na(InChIKey)]

# 1-3. CIGS CAS number -> PubChem synonym (exact) -> InChIKey
cigs <- as.data.table(readxl::read_excel(cigs_file, sheet = 3, skip = 1))
cas  <- data.table(name = cigs[["Compound name"]], synonym = trimws(cigs[["CAS Number"]]))
hit  <- synonym[cas, on = "synonym", .(name, CID), mult = "first", nomatch = 0]
m_cigs <- inchikey[hit, on = "CID", .(name, CID, InChIKey), nomatch = 0]

# 1-1 > 1-2 > 1-3
cp_map <- data.table(name = cp_names)
cp_map[m_pubchem, on = "name", `:=`(CID = i.CID, InChIKey = i.InChIKey, match_source = "pubchem_synonym")]
cp_map[is.na(InChIKey), c("InChIKey", "match_source") :=
         .(m_cmap$InChIKey[match(name, m_cmap$name)], "cmap")]
cp_map[is.na(InChIKey), c("CID", "InChIKey", "match_source") :=
         .(m_cigs$CID[match(name, m_cigs$name)], m_cigs$InChIKey[match(name, m_cigs$name)], "cigs_cas")]

# 1-4. GtoPdb for names still unmatched: ligand name, then ligand synonym
gtopdb <- as.data.table(read.csv(gtopdb_file, skip = 1, stringsAsFactors = FALSE))[InChIKey != ""]
gtopdb_syn <- gtopdb[Synonyms != ""][, .(InChIKey, synonym = lower(unlist(strsplit(Synonyms, "\\|")))),
                                     by = .(row = seq_len(sum(gtopdb$Synonyms != "")))]   # keeps file order
todo <- cp_map[is.na(InChIKey), name]
todo <- todo[!duplicated(lower(todo))]                       # one name per case-insensitive key
m_gtopdb <- unique(rbind(
  data.table(name = todo, InChIKey = gtopdb$InChIKey[match(lower(todo), lower(gtopdb$Name))],
             source = "gtopdb_name"),
  data.table(name = todo, InChIKey = gtopdb_syn$InChIKey[match(lower(todo), gtopdb_syn$synonym)],
             source = "gtopdb_synonym")
)[!is.na(InChIKey)], by = "name")
cp_map[is.na(InChIKey), c("InChIKey", "match_source") :=
         .(m_gtopdb$InChIKey[match(name, m_gtopdb$name)], m_gtopdb$source[match(name, m_gtopdb$name)])]
cp_map[is.na(InChIKey), match_source := NA_character_]

# STEP 2. Bioactive compounds: InChIKey -> standard name ---------------------------
# PubChem title of the InChIKey; otherwise the most frequent source name
rows   <- merge(pert[opt_id == 1, .(name)], cp_map[!is.na(InChIKey), .(name, InChIKey)], by = "name")
local  <- rows[, .(local_representative = names(sort(table(squish(name)), decreasing = TRUE))[1]), by = InChIKey]
ik_cid <- inchikey[InChIKey %in% local$InChIKey]                       # PubChem file order
ti     <- title_by_cid(ik_cid$CID)
titles <- ti[ik_cid, on = "CID"][
  , .(pubchem_cid_primary = CID[1], pubchem_title = title[!is.na(title) & title != ""][1]), by = InChIKey]
std <- merge(local, titles, by = "InChIKey", all.x = TRUE)
std[, standard_name := fifelse(!is.na(pubchem_title), pubchem_title, local_representative)]
cp_map <- merge(cp_map, std[, .(InChIKey, standard_name, local_representative)], by = "InChIKey", all.x = TRUE)

# STEP 3. Herbal ingredients: name -> InChIKey --------------------------------------
# Candidates from STEP 1-2 and PubChem are pooled; per name the highest priority is kept.
ig_names <- unique(pert[opt_id == 4 & (is.na(exclude_from_open) | exclude_from_open != 1), name])
q <- data.table(name = ig_names, name_key = tolower(squish(ig_names)))

cand <- rbind(
  cp_map[!is.na(InChIKey), .(name_key = tolower(squish(name)), CID, InChIKey,
                             source = "compound_final", priority = 1)],
  std[, .(name_key = tolower(squish(standard_name)), CID = pubchem_cid_primary, InChIKey,
          source = "standard_name", priority = 2)],
  std[, .(name_key = tolower(squish(local_representative)), CID = pubchem_cid_primary, InChIKey,
          source = "standard_local_representative", priority = 3)]
)[name_key %in% q$name_key]

synonym[, name_key := tolower(squish(synonym))]
syn_hit <- synonym[q, on = "name_key", .(name_key, CID), nomatch = 0, allow.cartesian = TRUE]
cand <- rbind(cand,
              inchikey[syn_hit, on = "CID", .(name_key, CID, InChIKey, source = "pubchem_synonym_local", priority = 10),
                       nomatch = 0, allow.cartesian = TRUE],
              inchikey[title_by_name(q$name_key), on = "CID", .(name_key, CID, InChIKey, source = "pubchem_title_local", priority = 11),
                       nomatch = 0, allow.cartesian = TRUE])
cand <- unique(cand[!is.na(InChIKey) & InChIKey != ""])
setorder(cand, name_key, priority, source, CID, InChIKey)
best <- cand[, .SD[1], by = name_key]

ig_map <- merge(q, best[, .(name_key, CID, InChIKey, match_source = source)], by = "name_key", all.x = TRUE)

# STEP 4. Merge names by InChIKey: one perturbagen per InChIKey -----------------------
# Standard name: STEP 2 for InChIKeys found among bioactive compounds; otherwise
# the first herbal-ingredient name (by source and perturbagen ID). Names without
# an InChIKey are kept as separate perturbagens (reviewed manually, see header).
ig_rows  <- pert[opt_id == 4 & (is.na(exclude_from_open) | exclude_from_open != 1)][order(odn_id, pertur_id)]
ig_first <- merge(ig_rows[, .(name, ord = .I)], ig_map[!is.na(InChIKey), .(name, InChIKey)], by = "name")[
  order(ord), .(ig_name = name[1]), by = InChIKey]

mapping <- rbind(
  cp_map[, .(name, category = "Bioactive compound", CID, InChIKey, match_source)],
  ig_map[, .(name, category = "Herbal ingredient", CID, InChIKey, match_source)])
mapping <- merge(mapping, std[, .(InChIKey, standard_name)], by = "InChIKey", all.x = TRUE, sort = FALSE)
mapping[is.na(standard_name) & !is.na(InChIKey),
        standard_name := ig_first$ig_name[match(InChIKey, ig_first$InChIKey)]]
mapping[is.na(InChIKey), standard_name := name]
mapping[, perturbagen_key := fifelse(!is.na(InChIKey), InChIKey, paste0("NO_INCHIKEY:", name))]

# Output ------------------------------------------------------------------------------
dir.create("output", showWarnings = FALSE)
fwrite(mapping[, .(name, category, CID, InChIKey, match_source, perturbagen_key, standard_name)],
       out_file, na = "", bom = TRUE)