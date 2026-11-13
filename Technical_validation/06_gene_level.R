# =============================================================================
# 06. Gene-level agreement of same-perturbagen signatures across data sources
#
#   Input values : DE signature files (all genes, matched by Symbol; duplicated
#                  Symbols averaged). MODZ for Connectivity Map, Log2FoldChange
#                  for other sources.
#   Pairs        : all signature pairs of the same perturbagen from different
#                  data sources. Pairs generated from the same raw data
#                  (SHARED_FILE) are excluded.
#   Genes        : for each correlation, genes measured in both signatures.
#   Measures     : Spearman correlation (all pairs);
#                  Pearson correlation (pairs in which both signatures report
#                  log2 fold changes).
#   Background   : for a pair (a from source s, b from source t), the mean
#                  correlation of a with signatures of other perturbagens from
#                  t and of b with those from s, averaged. Candidates: POOL_N
#                  signatures selected at random once per source (all if fewer).
#   Difference   : same-perturbagen minus background correlation.
#   Conditions   : model, dose and duration matched as in 05.
#   Summary      : pair level = mean over pairs;
#                  perturbagen level = B bootstrap iterations, each resampling
#                  perturbagens with replacement and sampling one pair per
#                  perturbagen. Mean difference, 2.5th and 97.5th percentiles,
#                  Empirical P = (b + 1) / (B + 1), b = iterations with a
#                  difference of zero or less.
#                  Rows: all pairs, eight combinations of matched conditions
#                  (and not comparable), pooled by experimental model.
#
#   Cache   : <OUT>/gene_cache/ (signature files are read once; an interrupted
#             run continues from the cache)
#   Input   : DE_DIR/**/DEsignature_SIG*.csv, SIG_FILE, SHARED_FILE
#   Output  : <OUT>/gene_level_agreement.csv
#             <OUT>/gene_level_pairs.rds  (intermediate, correlations of every pair)
#   Note    : reading all signature files takes several hours.
# =============================================================================
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(data.table); library(matrixStats) })
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), sprintf(...), "\n")

## ---- Settings: edit the paths ----------------------------------------------
DE_DIR      <- "PSDAP/Differential_expression_signatures"
SIG_FILE    <- "PSDAP/signature_info.csv"
SHARED_FILE <- "shared_raw_signature_pairs.csv"     # columns cvc1, cvc2
OUT         <- "output"
## ---------------------------------------------------------------------------
B         <- 10000L
SEED      <- 20261019
DOSE_FOLD <- 2
CHUNK     <- 1000L                                   # signatures per cache file
BATCH     <- 2000L                                   # pairs per block in correlation calculation (memory)
POOL_N    <- 500L                                    # background candidates per source
CDIR      <- file.path(OUT, "gene_cache"); dir.create(CDIR, showWarnings = FALSE, recursive = TRUE)
fmt <- function(x) format(x, big.mark = ",")
hm  <- function(s) if (s < 3600) sprintf("%.0f min", s / 60) else sprintf("%.1f h", s / 3600)

# ---- 1. Signatures and conditions (as in 05) -----------------------------------
S <- fread(SIG_FILE, encoding = "UTF-8", colClasses = "character", select = c("Signature_id", "Source", "Pertur_id", "Model_id", "Dose", "Duration"))
setnames(S, c("sig", "source", "drug", "model", "dose_raw", "time_raw"))
du <- c(nm = 1e-3, um = 1, mm = 1e3, "mg/kg" = 1, "ug/kg" = 1e-3, mg_per_kg = 1, "ug/ml" = 1, "ng/ml" = 1e-3, "mg/ml" = 1e3)
df <- c(nm = "M", um = "M", mm = "M", "mg/kg" = "kg", "ug/kg" = "kg", mg_per_kg = "kg", "ug/ml" = "ml", "ng/ml" = "ml", "mg/ml" = "ml")
d0 <- tolower(trimws(sub("\\s+(low|middle|high)$", "", S$dose_raw, ignore.case = TRUE)))
dun <- trimws(sub("^[0-9.eE+-]+\\s*", "", d0))
S[, `:=`(dose_val = suppressWarnings(as.numeric(sub("^([0-9.eE+-]+)\\s*.*$", "\\1", d0))) * du[dun], dose_fam = df[dun])]
tm <- tolower(trimws(S$time_raw)); tv <- suppressWarnings(as.numeric(sub("^([0-9.]+)\\s*(hr|h|hours?|days?)$", "\\1", tm)))
S[, time_h := fifelse(grepl("days?$", tm), tv * 24, tv)]
S[model %in% c("", "NA"), model := NA]
S <- S[!is.na(drug) & drug != ""]

IDX <- file.path(OUT, "de_index.rds")
if (file.exists(IDX)) f <- readRDS(IDX) else { f <- list.files(DE_DIR, pattern = "^DEsignature_SIG[0-9]+", recursive = TRUE, full.names = TRUE); saveRDS(f, IDX) }
S <- merge(S, data.table(sig = regmatches(basename(f), regexpr("SIG[0-9]+", basename(f))), file = f), by = "sig")
setorder(S, source, drug, sig)
msg("signatures with a DE file and a perturbagen: %s, sources: %d", fmt(nrow(S)), uniqueN(S$source))

# ---- 2. Pairs: same perturbagen, different sources ------------------------------
sh <- fread(SHARED_FILE, colClasses = "character")
EX <- unique(rbind(sh[, .(a = cvc1, b = cvc2)], sh[, .(a = cvc2, b = cvc1)]))
md <- S[, uniqueN(source), by = drug][V1 >= 2, drug]
X  <- S[drug %in% md, .(sig, source, drug, model, dose_val, dose_fam, time_h)]
PR <- merge(X, X, by = "drug", allow.cartesian = TRUE, suffixes = c("_a", "_b"))[source_a < source_b]
n0 <- nrow(PR); PR <- PR[!EX, on = .(sig_a = a, sig_b = b)]
PR[, `:=`(m_model = fifelse(is.na(model_a) | is.na(model_b), NA, model_a == model_b),
          m_dose  = fifelse(is.na(dose_val_a) | is.na(dose_val_b) | dose_fam_a != dose_fam_b, NA, abs(log2(dose_val_a / dose_val_b)) <= log2(DOSE_FOLD) + 1e-9),
          m_time  = fifelse(is.na(time_h_a) | is.na(time_h_b), NA, abs(time_h_a - time_h_b) < 1e-9))]
PR[, combo := fifelse(is.na(m_model) | is.na(m_dose) | is.na(m_time), "not comparable",
                      sprintf("model %s | dose %s | time %s", fifelse(m_model, "same", "diff"), fifelse(m_dose, "same", "diff"), fifelse(m_time, "same", "diff")))]
PR <- PR[, .(a = sig_a, b = sig_b, drug, s = source_a, t = source_b, combo)]
msg("cross-source pairs %s (perturbagens %s), excluded pairs from the same raw data %s", fmt(nrow(PR)), fmt(uniqueN(PR$drug)), fmt(n0 - nrow(PR)))
SPS <- unique(PR[, .(s, t)])

# ---- 3. Read signature files -> cache (by source and perturbagen) ---------------
rd <- function(ff) {
  x  <- suppressWarnings(fread(ff, showProgress = FALSE, select = c("Symbol", "MODZ", "Log2FoldChange")))
  st <- intersect(c("MODZ", "Log2FoldChange"), names(x))[1]
  x  <- x[!is.na(Symbol) & Symbol != "" & is.finite(get(st))][, .(v = mean(get(st))), by = Symbol]
  list(sym = x$Symbol, v = x$v, stat = st)
}
MAN <- file.path(CDIR, "manifest.rds")
man <- if (file.exists(MAN)) readRDS(MAN) else data.table(sig = character(), source = character(), drug = character(), stat = character(), chunk = character())
todo <- S[!sig %in% man$sig]
if (nrow(todo)) {
  msg("reading %s files (cached %s)", fmt(nrow(todo)), fmt(nrow(man))); t0 <- Sys.time()
  k <- nrow(man)
  for (src in unique(todo$source)) {
    x <- todo[source == src]; brk <- c(which(diff(match(x$drug, unique(x$drug))) != 0), nrow(x))   # perturbagen boundaries
    st <- 1L
    while (st <= nrow(x)) {
      en <- brk[brk >= min(st + CHUNK - 1L, nrow(x))][1]
      y  <- x[st:en]; L <- lapply(y$file, function(ff) tryCatch(rd(ff), error = function(e) NULL)); names(L) <- y$sig
      ok <- !vapply(L, is.null, TRUE)
      cf <- sprintf("%s_%s.rds", gsub("[^A-Za-z0-9]", "", src), y$sig[1]); saveRDS(L[ok], file.path(CDIR, cf))
      man <- rbind(man, data.table(sig = y$sig[ok], source = src, drug = y$drug[ok], stat = vapply(L[ok], `[[`, "", "stat"), chunk = cf)); saveRDS(man, MAN)
      k <- k + length(y$sig); el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      msg("  %s / %s  (%s elapsed)", fmt(k), fmt(nrow(S)), hm(el))
      st <- en + 1L
    }
  }
}
man <- man[sig %in% S$sig]

# ---- 4. Gene index over all cached Symbols ---------------------------------------
GF <- file.path(CDIR, "gene_universe.rds")
GU <- if (file.exists(GF)) readRDS(GF) else NULL
if (is.null(GU) || !identical(GU$key, nrow(man))) {
  GN <- character()
  for (cf in unique(man$chunk)) GN <- unique(c(GN, unlist(lapply(readRDS(file.path(CDIR, cf)), `[[`, "sym"), use.names = FALSE)))
  GU <- list(GN = GN, key = nrow(man)); saveRDS(GU, GF)
}
GN <- GU$GN
ld <- function(sigs) {                               # signatures -> list(g = sorted gene index, v, k = gene-set key)
  out <- list()
  for (cf in unique(man[.(sigs), on = "sig", chunk])) {
    L <- readRDS(file.path(CDIR, cf)); L <- L[names(L) %in% sigs]
    out <- c(out, lapply(L, function(e) { g <- match(e$sym, GN); o <- order(g); g <- g[o]
    list(g = g, v = e$v[o], k = sprintf("%d_%.0f_%.0f", length(g), sum(as.numeric(g)), sum(as.numeric(g)^2 %% 1e9))) }))
  }
  out
}

# ---- 5. Correlations over genes measured in both signatures ----------------------
#   Signatures with the same gene set are processed together as a matrix
stdz <- function(M) { M <- M - rowMeans(M); M / (rowSds(M) * sqrt(ncol(M) - 1)) }   # row dot product = correlation
zmat <- function(Lg, G) {                            # signatures with the same gene set, genes G -> standardized ranks and values
  idx <- match(G, Lg[[1]]$g); M <- matrix(unlist(lapply(Lg, function(e) e$v[idx]), use.names = FALSE), length(Lg), byrow = TRUE)
  R <- stdz(rowRanks(M, ties.method = "average")); P <- stdz(M); rownames(R) <- rownames(P) <- names(Lg); list(r = R, p = P)
}
pdot <- function(A, ia, Bm, ib) { out <- numeric(length(ia)); for (st in seq(1, length(ia), BATCH)) { j <- st:min(st + BATCH - 1, length(ia)); out[j] <- rowSums(A[ia[j], , drop = FALSE] * Bm[ib[j], , drop = FALSE]) }; out }
keyv <- function(L) vapply(L, `[[`, "", "k")

set.seed(SEED)                                       # background candidates, fixed per source
POOL <- lapply(split(man$sig, man$source), function(x) if (length(x) > POOL_N) sample(x, POOL_N) else x)
POOL <- lapply(POOL, function(x) { L <- ld(x); list(L = L, drug = man[.(names(L)), on = "sig", drug], k = keyv(L)) })

pairs1 <- function(pr, LA, LB) {                     # pair correlations, grouped by gene-set combination
  pr[, `:=`(ka = keyv(LA)[a], kb = keyv(LB)[b])]
  pr[, c("n_gene", "rho", "r") := { G <- intersect(LA[[a[1]]]$g, LB[[b[1]]]$g)
  if (length(G) < 3) list(length(G), NA_real_, NA_real_) else {
    ua <- unique(a); ub <- unique(b); ZA <- zmat(LA[ua], G); ZB <- zmat(LB[ub], G); ia <- match(a, ua); ib <- match(b, ub)
    list(length(G), pdot(ZA$r, ia, ZB$r, ib), pdot(ZA$p, ia, ZB$p, ib)) } }, by = .(ka, kb)]
  pr[, c("ka", "kb") := NULL][]
}
base1 <- function(L, pool) {                         # mean correlation of each signature with candidates of other perturbagens
  d <- man[.(names(L)), on = "sig", drug]; k <- keyv(L)
  sr <- sp <- nn <- setNames(numeric(length(L)), names(L)); memo <- new.env()
  for (ks in unique(k)) {
    i <- which(k == ks)
    for (kp in unique(pool$k)) {
      j <- which(pool$k == kp); G <- intersect(L[[i[1]]]$g, pool$L[[j[1]]]$g); if (length(G) < 3) next
      mk <- paste(kp, length(G), sum(as.numeric(G)))
      ZP <- memo[[mk]]; if (is.null(ZP)) { ZP <- zmat(pool$L[j], G); memo[[mk]] <- ZP }
      ZS <- zmat(L[i], G); m <- outer(d[i], pool$drug[j], "!=")
      ok <- is.finite(ZP$r[, 1]); m[, !ok] <- FALSE; ZPr <- ZP$r; ZPr[!ok, ] <- 0; ZPp <- ZP$p; ZPp[!ok, ] <- 0
      sr[i] <- sr[i] + rowSums((ZS$r %*% t(ZPr)) * m); sp[i] <- sp[i] + rowSums((ZS$p %*% t(ZPp)) * m); nn[i] <- nn[i] + rowSums(m)
    }
  }
  list(r = sr / nn, p = sp / nn)
}

RES <- list(); t0 <- Sys.time()
setorder(SPS, s, t)
for (i in seq_len(nrow(SPS))) {
  s <- SPS$s[i]; t <- SPS$t[i]; pr <- PR[SPS[i], on = .(s, t)][a %in% man$sig & b %in% man$sig]; if (!nrow(pr)) next
  sm <- if (uniqueN(pr$a) <= uniqueN(pr$b)) "a" else "b"; lg <- setdiff(c("a", "b"), sm)   # smaller side loaded once, larger side by cache file
  Lsm <- ld(unique(pr[[sm]])); bsm <- base1(Lsm, POOL[[if (sm == "a") t else s]])
  out <- list()
  for (cf in unique(man[sig %in% pr[[lg]], chunk])) {
    ss <- intersect(man[chunk == cf, sig], pr[[lg]]); Llg <- ld(ss)
    p1 <- pr[get(lg) %in% names(Llg) & get(sm) %in% names(Lsm)]; if (!nrow(p1)) next
    p1 <- if (sm == "a") pairs1(p1, Lsm, Llg) else pairs1(p1, Llg, Lsm)
    blg <- base1(Llg, POOL[[if (lg == "a") t else s]])
    ba <- if (sm == "a") bsm else blg; bb <- if (sm == "a") blg else bsm
    p1[, `:=`(rho_base = (ba$r[a] + bb$r[b]) / 2, r_base = (ba$p[a] + bb$p[b]) / 2)]; out[[cf]] <- p1
  }
  pr <- rbindlist(out)
  if (!nrow(pr)) next
  pr[, lfc := man[.(a), on = "sig", stat] == "Log2FoldChange" & man[.(b), on = "sig", stat] == "Log2FoldChange"]
  RES[[i]] <- pr; gc()
  msg("  [%d/%d] %s - %s : pairs %s  (%s elapsed)", i, nrow(SPS), s, t, fmt(nrow(pr)), hm(as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
P <- rbindlist(RES)[is.finite(rho) & is.finite(rho_base)]
saveRDS(list(pairs = P, pool = lapply(POOL, function(p) names(p$L)), B = B, seed = SEED, pool_n = POOL_N), file.path(OUT, "gene_level_pairs.rds"))

# ---- 6. Summary: pair level and perturbagen level -------------------------------
CMB <- c("model same | dose same | time same", "model same | dose same | time diff", "model same | dose diff | time same", "model same | dose diff | time diff",
         "model diff | dose same | time same", "model diff | dose same | time diff", "model diff | dose diff | time same", "model diff | dose diff | time diff",
         "not comparable")
MCL <- c("Same model", "Different model", "Not comparable")
P[, `:=`(combo = factor(combo, CMB),
         model = factor(fifelse(combo == "not comparable", MCL[3], fifelse(grepl("^model same", combo), MCL[1], MCL[2])), MCL))]
summ <- function(z, v, bv) {
  x <- z[[v]]; y <- z[[bv]]; d <- x - y; ok <- is.finite(d); z <- z[ok]; x <- x[ok]; y <- y[ok]; d <- d[ok]
  if (length(d) < 2) return(NULL)
  o <- order(z$drug); g <- z[o, .(st = .I[1], sz = .N), by = drug]; dd <- d[o]
  D  <- matrix(sample.int(nrow(g), B * nrow(g), replace = TRUE), B)                # resample perturbagens with replacement
  ix <- matrix(g$st[D] + floor(runif(B * nrow(g)) * g$sz[D]), B)                  # one pair per sampled perturbagen
  ds <- rowMeans(matrix(dd[ix], B))
  data.table(n_pair = length(d), n_drug = nrow(g), same = mean(x), base = mean(y), diff = mean(d),
             diff_drug = mean(ds), diff_lo = unname(quantile(ds, .025)), diff_hi = unname(quantile(ds, .975)), p_drug = (sum(ds <= 0) + 1) / (B + 1))
}
both <- function(by = NULL) {                        # Spearman (all pairs), Pearson (log2 fold-change pairs)
  f <- function(z, v, bv) if (is.null(by)) summ(z, v, bv) else z[, summ(.SD, v, bv), by = by]
  rbind(f(P, "rho", "rho_base")[, measure := "Spearman"], f(P[lfc == TRUE], "r", "r_base")[, measure := "Pearson"], fill = TRUE)
}
set.seed(SEED)
A  <- both()[, `:=`(Comparison = "All cross-source pairs", Model = "Any", Dose = "Any", Duration = "Any", row = 0L)]
Bc <- both("combo")
Cm <- both("model")[model != "Not comparable"]

yn <- function(cb, w) fifelse(grepl(paste(w, "same"), cb), "Yes", "No")
Bc[, combo := as.character(combo)]
Bc[, `:=`(Comparison = fifelse(combo == "not comparable", "Not comparable", "Combination"),
          Model = fifelse(combo == "not comparable", NA_character_, yn(combo, "model")),
          Dose = fifelse(combo == "not comparable", NA_character_, yn(combo, "dose")),
          Duration = fifelse(combo == "not comparable", NA_character_, yn(combo, "time")),
          row = fifelse(combo == "not comparable", 11L, match(combo, CMB)))]
Cm[, `:=`(Comparison = fifelse(model == "Same model", "Model matched (pooled)", "Model unmatched (pooled)"),
          Model = fifelse(model == "Same model", "Yes", "No"), Dose = "Any", Duration = "Any",
          row = fifelse(model == "Same model", 9L, 10L))]
d <- rbindlist(list(A, Bc, Cm), fill = TRUE)
d[, m_o := match(measure, c("Spearman", "Pearson"))]; setorder(d, m_o, row)
tab <- d[, .(Measure = measure, Comparison, Model_matched = Model, Dose_matched = Dose, Duration_matched = Duration,
             Perturbagens = n_drug, Pairs = n_pair,
             Pair_same_perturbagen = same, Pair_background = base, Pair_difference = diff,
             Perturbagen_difference = diff_drug, Perturbagen_lower_95 = diff_lo, Perturbagen_upper_95 = diff_hi, Perturbagen_P = p_drug)]
fwrite(tab, file.path(OUT, "gene_level_agreement.csv"))
print(tab[, .(Measure, Comparison, Model_matched, Dose_matched, Duration_matched, Perturbagens, Pair_difference, Perturbagen_difference, Perturbagen_P)])
msg("saved: %s", file.path(OUT, "gene_level_agreement.csv"))
