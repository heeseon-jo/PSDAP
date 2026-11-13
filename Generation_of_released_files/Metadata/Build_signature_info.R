# =============================================================================
# build_signature_info.R
#
# Builds the released metadata file signature_info.csv, with one row per
# differential-expression signature (case versus control comparison).
#
# Signatures and their sample annotations are read from the internal PSDAP
# database, which is not distributed. Only comparisons marked for open release
# are included.
#
# Output columns:
#   Signature_id                 "SIG" + zero-padded comparison identifier
#   Perturbagen, Pertur_id       perturbagen of the case sample
#   Model_type, Model, Model_id  experimental model of the case sample
#   Experimental_system_species  species of the model
#   Dose, Duration, Extraction   condition of the case sample ("-" as NA)
#   Control_condition            perturbagen of the control sample (empty as NA)
#   Source, Category, Source_id  data source of the case perturbagen
#
# Database credentials are read from environment variables (DB_HOST, DB_PORT,
# DB_NAME, DB_USER, DB_PW), e.g. set in ~/.Renviron.
#
# Dependencies: data.table, DBI, RMariaDB.
# =============================================================================

library(data.table)
library(DBI)

con <- dbConnect(RMariaDB::MariaDB(),
                 host = Sys.getenv("DB_HOST"), port = as.integer(Sys.getenv("DB_PORT")),
                 dbname = Sys.getenv("DB_NAME"),
                 user = Sys.getenv("DB_USER"), password = Sys.getenv("DB_PW"))
read_table <- function(name) dbGetQuery(con, paste("SELECT * FROM", name))

# 1) Signatures released openly ------------------------------------------------
cvc <- dbGetQuery(con, "SELECT * FROM c_vs_c
                         WHERE exclude_from_open IS NULL OR exclude_from_open <> 1")
w <- nchar(max(cvc$cvc_id))
cvc$Signature_id <- sprintf("SIG%0*d", w, cvc$cvc_id)

# 2) Case and control sample annotations ---------------------------------------
expc <- read_table("exp_condition2")
case <- expc[match(cvc$case_id,    expc$exp_id), ]
ctrl <- expc[match(cvc$control_id, expc$exp_id), ]

# 3) Option tables (curated vocabularies) --------------------------------------
pert   <- read_table("option_perturbagen")
cell   <- read_table("option_cell_line")
cls    <- read_table("option_cell_line_standard")
mtype  <- read_table("option_model_type")
dose   <- read_table("option_dose")
otime  <- read_table("option_time")
extr   <- read_table("option_extract")
dbname <- read_table("option_db_name")
ptype  <- read_table("option_pert_type")
dbDisconnect(con)

lookup   <- function(key, table, key_col, value_col) table[[value_col]][match(key, table[[key_col]])]
dash_na  <- function(x) { x[x == "-"] <- NA; x }
blank_na <- function(x) { x[x == ""]  <- NA; x }

# 4) Columns of signature_info.csv ---------------------------------------------
ocls_id   <- lookup(case$cell_id,   cell,   "cell_id",  "ocls_id")
source_id <- lookup(case$pertur_id, pert,   "pertur_id", "odn_id")

signature_info <- data.frame(
  Signature_id = cvc$Signature_id,
  Perturbagen  = lookup(case$pertur_id, pert, "pertur_id", "name"),
  Model_type   = lookup(ocls_id, cls, "ocls_id", "ca"),
  Model        = lookup(case$cell_id, cell, "cell_id", "name_origin"),
  Dose         = dash_na(lookup(case$dose_id,    dose, "dose_id",    "name")),
  Duration     = dash_na(lookup(case$ot_id,      otime, "ot_id",     "time")),
  Extraction   = dash_na(lookup(case$extract_id, extr, "extract_id", "name")),
  Control_condition = blank_na(lookup(ctrl$pertur_id, pert, "pertur_id", "name")),
  Source       = lookup(source_id, dbname, "odn_id", "db_name"),
  Category     = lookup(lookup(source_id, dbname, "odn_id", "opt_id"), ptype, "opt_id", "pert_type"),
  Experimental_system_species = lookup(lookup(case$cell_id, cell, "cell_id", "omt_id"),
                                       mtype, "omt_id", "model_type"),
  Pertur_id    = lookup(case$pertur_id, pert, "pertur_id", "ops_id"),
  Model_id     = ocls_id,
  Source_id    = source_id
)

stopifnot(!anyDuplicated(signature_info$Signature_id),
          !anyNA(signature_info[, c("Perturbagen", "Pertur_id", "Source", "Category", "Source_id")]))

# 5) Save ----------------------------------------------------------------------
dir.create("output", showWarnings = FALSE)
fwrite(signature_info, "output/signature_info.csv", na = "", bom = TRUE)