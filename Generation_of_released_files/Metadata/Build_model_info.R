# =============================================================================
# build_model_info.R
#
# Builds the released metadata file model_info.csv, with one row per
# experimental model (Model_id) used in signature_info.csv.
#
# Model annotations come from an internal curation table of the PSDAP database
# (option_cell_line_standard), which maps cell line and tissue names to
# Cellosaurus or UBERON identifiers with disease, tissue and synonym
# information. The table itself is not distributed; model_info.csv contains
# its curated annotations for every model used in PSDAP. Species is taken
# from signature_info.csv.
#
# Database credentials are read from environment variables (DB_HOST, DB_PORT,
# DB_NAME, DB_USER, DB_PW), e.g. set in ~/.Renviron.
#
# Dependencies: data.table, DBI, RMariaDB.
# =============================================================================

library(data.table)
library(DBI)

# 1) Models used in PSDAP --------------------------------------------------------
signature_info <- fread("output/signature_info.csv", data.table = FALSE)
species <- unique(signature_info[, c("Model_id", "Experimental_system_species")])
stopifnot(!anyDuplicated(species$Model_id))            # one species per model

# 2) Curation table from the PSDAP database --------------------------------------
con <- dbConnect(RMariaDB::MariaDB(),
                 host = Sys.getenv("DB_HOST"), port = as.integer(Sys.getenv("DB_PORT")),
                 dbname = Sys.getenv("DB_NAME"),
                 user = Sys.getenv("DB_USER"), password = Sys.getenv("DB_PW"))
curation <- dbGetQuery(con, "SELECT * FROM option_cell_line_standard")
dbDisconnect(con)

stopifnot(!anyDuplicated(curation$ocls_id),
          all(species$Model_id %in% curation$ocls_id))  # every model is curated
m <- curation[curation$ocls_id %in% species$Model_id, ]

# 3) Columns of model_info.csv -------------------------------------------------
model_info <- data.frame(
  Model_id       = m$ocls_id,
  Model_name     = m$standard_name2,
  Model_standard = m$samekey,
  Model_standard_type = ifelse(grepl("CVCL", m$samekey), "Cellosaurus",
                               ifelse(grepl("UBERON", m$samekey), "UBERON", NA)),
  Experimental_system_species =
    species$Experimental_system_species[match(m$ocls_id, species$Model_id)],
  Disease  = ifelse(m$disease == "", NA, m$disease),
  Tissue   = m$tissue,
  Synonyms = m$synonyms
)

# 4) Save --------------------------------------------------------------------------
dir.create("output", showWarnings = FALSE)
fwrite(model_info, "output/model_info.csv", na = "", bom = TRUE)