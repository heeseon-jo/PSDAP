# =============================================================================
# build_source_info.R
#
# Builds the released metadata file source_info.csv, with one row per data
# source (Source_id) of the PSDAP database.
#
# Source annotations are read from the internal PSDAP database, which is not
# distributed.
#
# Output columns:
#   Source_id, Source     identifier and name of the data source
#   Categories            perturbagen category of the source
#   Assay, Assay_id       assay of the source and its identifier
#   Device, Device_id     measurement device and its identifier
#   DE_method             differential expression method
#   Reference             reference of the source
#
# Database credentials are read from environment variables (DB_HOST, DB_PORT,
# DB_NAME, DB_USER, DB_PW), e.g. set in ~/.Renviron.
#
# Dependencies: DBI, RMariaDB, dplyr, readr.
# =============================================================================

library(DBI)
library(readr)


con <- dbConnect(RMariaDB::MariaDB(),
                 host = Sys.getenv("DB_HOST"), port = as.integer(Sys.getenv("DB_PORT")),
                 dbname = Sys.getenv("DB_NAME"),
                 user = Sys.getenv("DB_USER"), password = Sys.getenv("DB_PW"))



#data load
od <- dbGetQuery(con, "select * from option_db_name;")
opt <- dbGetQuery(con, "select * from option_pert_type;")
assay <- dbGetQuery(con, "select * from assay_identifier;")
de <- dbGetQuery(con, "select * from de_method;")
dbDisconnect(con)


# 1) Data sources ------------------------------------------------------------------
od <- od[,c("odn_id", "opt_id", "db_name", "ai_id")]

# Category (perturbagen type)
od$cat <- opt$pert_type[match(od$opt_id, opt$opt_id)]

# Assay, device and reference, matched by assay identifier (ai_id).
od$assay <- assay$assay[match(od$ai_id, assay$ai_id)]
od$assay_id <- assay$assay_id[match(od$ai_id, assay$ai_id)]
od$device <- assay$device[match(od$ai_id, assay$ai_id)]
od$device_id <- assay$device_id[match(od$ai_id, assay$ai_id)]
od$refer <- assay$reference[match(od$ai_id, assay$ai_id)]

# DE method of each assay (ai_id -> de_id -> de_meth).
tmp_od <- unique(od[, c("db_name", "ai_id")])
tmp_assay <- unique(assay[,c("ai_id", "de_id")])
tmp <- dplyr::left_join(tmp_od, tmp_assay, by = "ai_id")
tmp_de <- dplyr::left_join(tmp, de, by = "de_id")
de_clear <- tmp_de[,c("ai_id", "de_meth")]

od$de <- de_clear$de_meth[match(od$ai_id, de_clear$ai_id)]



# 2) Output columns and save -------------------------------------------------------
source_df <- od[,c("odn_id","db_name","cat", "assay","assay_id","device", "device_id", "de", "refer")]
colnames(source_df) <- c("Source_id", "Source", "Categories", "Assay","Assay_id","Device", "Device_id", "DE_method", "Reference")

dir.create("output", showWarnings = FALSE)
write_excel_csv(source_df, "output/source_info.csv", na = "")
