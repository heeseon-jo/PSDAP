# =============================================================================
# 01. Top-100 most similar profiles for each query profile (GO-BP similarity)
#
#   Profile : (Signature_id, direction, threshold). Connectivity Map signatures
#             have the threshold "CMap".
#   Query   : the signature on the side of the file threshold
#             (logFC files: non-CMap side; CMap file: CMap side). Both columns
#             of each file are checked.
#   Partner threshold : CMap file -> Stringency columns; logFC files -> "CMap"
#             if the partner is a CMap signature, otherwise the file threshold.
#   Ranking : pairs with a GO-BP score; score rounded to three decimals; ordered
#             by (-score, partner Signature_id, partner threshold); top 100.
#
#   Input   : Perturbagen_similarity/Similarity_*.csv, signature_info.csv
#   Output  : <OUT>/top100_{up,down}.rds  (intermediate, read by 02-05)
# =============================================================================
rm(list = ls()); gc()
suppressPackageStartupMessages(library(data.table))

## ---- Settings: edit the paths ----------------------------------------------
SIM_DIR  <- "PSDAP/Perturbagen_similarity"
SIG_FILE <- "PSDAP/signature_info.csv"
OUT      <- "output"
## ---------------------------------------------------------------------------
TOPN <- 100L
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

sig  <- fread(SIG_FILE, select = c("Signature_id", "Pertur_id", "Source"), colClasses = "character")
cmap <- sig[Source == "Connectivity Map", Signature_id]
fs   <- list.files(SIM_DIR, "^Similarity_.*\\.csv$", full.names = TRUE)
thr_of <- function(f) if (grepl("Connectivity_map", f)) "CMap" else
  sub("^Similarity_Log2FoldChange_([0-9]+)_AdjustedPvalue_0_([0-9]+)\\.csv$", "LFC\\1_P0.\\2", basename(f))
lab <- function(v) fifelse(grepl("Connectivity|CMap|_NA", v), "CMap",
                           sub("^Log2FoldChange_([0-9]+)_AdjustedPvalue_0_([0-9]+)$", "LFC\\1_P0.\\2", v))
LABS <- c("CMap", outer(c("LFC0", "LFC1", "LFC2"), c("P0.05", "P0.01"), paste, sep = "_"))

# ---- 1. All rows containing a query, from both columns of each file ---------
log <- list()
L <- rbindlist(lapply(fs, function(f) {
  th <- thr_of(f); h <- names(fread(f, nrows = 0))
  tc <- if (th == "CMap") c(grep("^(Stringency|Threshold).*_1$", h, value = TRUE)[1], grep("^(Stringency|Threshold).*_2$", h, value = TRUE)[1])
  if (th == "CMap" && anyNA(tc)) stop("No threshold columns in the CMap file. Columns: ", paste(h, collapse = ", "))
  x <- fread(f, select = c("Signature_id_1", "Signature_id_2", "BP_score", "Direction", tc))
  setnames(x, c("s1", "s2", "bp", "dir", if (length(tc)) c("h1", "h2")))
  if (length(tc)) x[, `:=`(h1 = lab(h1), h2 = lab(h2))] else
    x[, `:=`(h1 = fifelse(s1 %in% cmap, "CMap", th), h2 = fifelse(s2 %in% cmap, "CMap", th))]
  log[[th]] <<- data.table(file = basename(f), thr = th, rows = nrow(x), bp_na = x[, sum(is.na(bp))],
                           bad_label = x[, sum(!h1 %in% LABS | !h2 %in% LABS)])
  x <- x[!is.na(bp)][, `:=`(bp = round(bp, 3), dir = tolower(dir))]
  rbind(x[h1 == th, .(q_sig = s1, t_sig = s2, t_thr = h2, bp, dir)],
        x[h2 == th, .(q_sig = s2, t_sig = s1, t_thr = h1, bp, dir)])[, q_thr := th]
}))
cat("\n[Rows per file, missing GO-BP scores, unrecognized threshold labels]\n"); print(rbindlist(log))

# ---- 2. Profile table --------------------------------------------------------
prof <- unique(L[, .(Signature_id = q_sig, dir, thr = q_thr)])
prof <- merge(prof, sig, by = "Signature_id", all.x = TRUE)
prof[, `:=`(cmap = thr == "CMap", pid = .I)]
cat(sprintf("\n[Profiles] %s\n", format(nrow(prof), big.mark = ",")))
L[prof, on = .(q_sig = Signature_id, dir, q_thr = thr), q := i.pid]

# ---- 3. Top-100 list per direction ------------------------------------------
for (D in c("up", "down")) {
  T <- L[dir == D, .(q, t_sig, t_thr, bp)]
  setorder(T, q, -bp, t_sig, t_thr); T[, rnk := rowid(q)]
  TOP <- T[rnk <= TOPN]; rm(T); gc()
  saveRDS(list(TOP = TOP, prof = prof[dir == D], sig = sig), file.path(OUT, sprintf("top100_%s.rds", D)))
  cat(sprintf("[%s] queries %s, saved %s\n", D, format(uniqueN(TOP$q), big.mark = ","), file.path(OUT, sprintf("top100_%s.rds", D))))
}
