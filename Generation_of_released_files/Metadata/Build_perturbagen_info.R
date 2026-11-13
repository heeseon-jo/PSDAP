# =============================================================================
# build_perturbagen_info.R
#
# Builds the released metadata file perturbagen_info.csv, with one row per
# standardized perturbagen used in the openly released signatures.
#
# Perturbagens and their standard identifiers are read from the internal PSDAP
# database, which is not distributed. Only perturbagens of comparisons marked
# for open release are included.
#
# Output columns:
#   Pertur_id, Pertur_name       identifier and name of the standardized perturbagen
#   Perturbagen_standard         standard identifier (InChIKey, formula name or
#                                Latin name; empty if none)
#   Perturbagen_standard_type    type of the standard identifier
#   Synonyms                     synonyms of the perturbagen
#
# Database credentials are read from environment variables (DB_HOST, DB_PORT,
# DB_NAME, DB_USER, DB_PW), e.g. set in ~/.Renviron.
#
# Dependencies: DBI, RMariaDB, dplyr, readr.
# =============================================================================

library(DBI)
library(readr)
library(data.table)



con <- dbConnect(RMariaDB::MariaDB(),
                 host = Sys.getenv("DB_HOST"), port = as.integer(Sys.getenv("DB_PORT")),
                 dbname = Sys.getenv("DB_NAME"),
                 user = Sys.getenv("DB_USER"), password = Sys.getenv("DB_PW"))



# data load
c_vs_c <- dbGetQuery(con, "select * from c_vs_c;")
exp_condition2 <- dbGetQuery(con, "select * from exp_condition2;")
op <- dbGetQuery(con, "select * from option_perturbagen;")
ops <- dbGetQuery(con, "select * from option_perturbagen_standard;")
opts <- dbGetQuery(con, "select * from perturbagen_standard_type;")
dbDisconnect(con)


# 1) Signatures released openly ------------------------------------------------
c_vs_c_open <- c_vs_c[is.na(c_vs_c$exclude_from_open) | c_vs_c$exclude_from_open != 1, ]
case_id_list <- unique(c_vs_c_open$case_id)



# 2) Case samples of the open signatures -----------------------------------------
# Match case samples (exp_id) and keep those not excluded from open release.
exp_matched <- exp_condition2[exp_condition2$exp_id %in% case_id_list, ] # case_id & exp_id match
exp_matched_open <- exp_matched[is.na(exp_matched$exclude_from_open) | exp_matched$exclude_from_open != 1, ]
pertur_id_list <- unique(exp_matched_open$pertur_id) # select unique perturbagen id



# 3) Perturbagens of those samples -----------------------------------------------
# Match perturbagens (pertur_id), keep open ones, and get their standard ID (ops_id).
op_matched <- op[op$pertur_id %in% pertur_id_list, ]
op_matched_open <- op_matched[is.na(op_matched$exclude_from_open) | op_matched$exclude_from_open != 1, ]
ops_id_list <- unique(op_matched_open$ops_id[!is.na(op_matched_open$ops_id)])
 


# 4) Standardized perturbagens ---------------------------------------------------
# Match standard perturbagens (ops_id) and keep those not excluded from open release.
ops_matched <- ops[ops$ops_id %in% ops_id_list, ]
ops_matched_open <- ops_matched[is.na(ops_matched$exclude_from_open) | ops_matched$exclude_from_open != 1, ]



# 5) Standard type ---------------------------------------------------------------
df <- ops_matched_open[, c("ops_id", "name", "samekey","pst_id" ,"synonyms_all", "opt_id_standard", "exclude_from_open")]

# Map opt_id_standard to the standard type id (pst_id).
df$pst_id <- dplyr::case_when(
  df$opt_id_standard %in% c("1", "1;4", "4") ~ 1,
  df$opt_id_standard == "2" ~ 2,
  df$opt_id_standard == "3" ~ 3,
  TRUE ~ NA_real_
)

# No standard identifier (samekey NA or empty) -> no standard type.
df$pst_id[is.na(df$samekey) | df$samekey == ""] <- NA_real_

# Add the standard type name.
df$standard_type <- opts$pertur_stand_type[match(df$pst_id, opts$pst_id)]



# 6) Output columns and save -----------------------------------------------------
# Select and rename the columns of perturbagen_info.csv.
fdf <- df[,c("ops_id","name", "samekey", "standard_type", "synonyms_all")]
colnames(fdf) <- c("Pertur_id", "Pertur_name", "Perturbagen_standard", "Perturbagen_standard_type", "Synonyms")

dir.create("output", showWarnings = FALSE)
write_excel_csv(fdf, "output/perturbagen_info.csv", na = "")