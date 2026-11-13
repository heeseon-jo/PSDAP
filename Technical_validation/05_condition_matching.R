# =============================================================================
# 05. Cross-source recovery by matched experimental conditions
#
#   Comparisons : each pair of a query profile and a comparable profile of the
#                 same perturbagen from a different data source. Queries and
#                 candidate numbers M are as in 02.
#   Conditions  : Model  = same Model_id
#                 Dose   = same unit family and differing by no more than twofold
#                          after conversion to a common unit
#                 Duration = identical after conversion to hours
#                 -> eight combinations of matched / unmatched conditions.
#                 Comparisons with a missing or non-comparable value for any
#                 condition are "not comparable" and excluded from recovery.
#   Recovery    : for each (query, group of comparisons), k = number of partner
#                 profiles in the group, hit = the best-ranked one within the
#                 top 100, expected = 1 - C(M - k, 100) / C(M, 100).
#                 Profile level: all queries. Perturbagen level: one query
#                 profile per perturbagen in each of B iterations; fold
#                 enrichment = mean observed / mean expected recovery over
#                 iterations, 95% range = 2.5th and 97.5th percentiles of the
#                 per-iteration fold.
#     A. Eight combinations of matched conditions
#     B. Pooled by experimental model (matched vs unmatched). Empirical P at
#        the perturbagen level = (b + 1) / (B + 1), b = number of iterations in
#        which recovery drawn with the expected probabilities is at least the
#        observed recovery.
#   Composition : cross-source profile pairs of the same perturbagen by
#                 combination, and (perturbagen, source pair) units assigned to
#                 their best-matched combination; percentages include not comparable.
#   Signature pairs generated from the same raw data are excluded as in 02.
#
#   Input   : <OUT>/top100_{up,down}.rds (01), <OUT>/recovery_{up,down}.rds (02), SIG_FILE, SHARED_FILE
#   Output  : <OUT>/recovery_by_matched_conditions.csv
# =============================================================================
rm(list = ls()); gc()
suppressPackageStartupMessages(library(data.table))
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), sprintf(...), "\n")

## ---- Settings: edit the paths ----------------------------------------------
SIG_FILE    <- "PSDAP/signature_info.csv"
SHARED_FILE <- "shared_raw_signature_pairs.csv"     # columns cvc1, cvc2
OUT         <- "output"
## ---------------------------------------------------------------------------
TOPN      <- 100L
B         <- 10000L
SEED_A    <- 20261016                                # A. eight combinations
SEED_B    <- 20261019                                # B. pooled by model
DOSE_FOLD <- 2
pany <- function(k, M, n) ifelse(k <= 0, 0, 1 - exp(lchoose(pmax(M - k, 0), n) - lchoose(M, n)))
sub_index <- function(drug, B) {                     # B x (number of perturbagens) row indices, one query per perturbagen
  g <- data.table(drug)[, .(st = .I[1], sz = .N), by = drug]
  matrix(g$st, B, nrow(g), byrow = TRUE) + floor(matrix(runif(B * nrow(g)), B) * rep(g$sz, each = B))
}
sh <- fread(SHARED_FILE, colClasses = "character")
EX <- unique(rbind(sh[, .(a = cvc1, b = cvc2)], sh[, .(a = cvc2, b = cvc1)]))
CMB <- CJ(model = c(TRUE, FALSE), dose = c(TRUE, FALSE), time = c(TRUE, FALSE), sorted = FALSE)
CMB[, combo := sprintf("model %s | dose %s | time %s", fifelse(model, "same", "diff"), fifelse(dose, "same", "diff"), fifelse(time, "same", "diff"))]
lv <- CMB$combo

# ---- Conditions: dose in a common unit within each unit family, duration in hours
S <- fread(SIG_FILE, encoding = "UTF-8", colClasses = "character", select = c("Signature_id", "Model_id", "Dose", "Duration"))
setnames(S, c("Signature_id", "model", "dose_raw", "time_raw"))
du <- c(nm = 1e-3, um = 1, mm = 1e3, "mg/kg" = 1, "ug/kg" = 1e-3, mg_per_kg = 1, "ug/ml" = 1, "ng/ml" = 1e-3, "mg/ml" = 1e3)
df <- c(nm = "M", um = "M", mm = "M", "mg/kg" = "kg", "ug/kg" = "kg", mg_per_kg = "kg", "ug/ml" = "ml", "ng/ml" = "ml", "mg/ml" = "ml")
d0 <- tolower(trimws(sub("\\s+(low|middle|high)$", "", S$dose_raw, ignore.case = TRUE)))
dn <- suppressWarnings(as.numeric(sub("^([0-9.eE+-]+)\\s*.*$", "\\1", d0)))
dun <- trimws(sub("^[0-9.eE+-]+\\s*", "", d0))
S[, `:=`(dose_val = dn * du[dun], dose_fam = df[dun])]
tm <- tolower(trimws(S$time_raw)); tv <- suppressWarnings(as.numeric(sub("^([0-9.]+)\\s*(hr|h|hours?|days?)$", "\\1", tm)))
S[, time_h := fifelse(grepl("days?$", tm), tv * 24, tv)]
S[model %in% c("", "NA"), model := NA]

A <- list(); G <- list(); PR <- list(); PRD <- list()
for (D in c("up", "down")) {
  r <- readRDS(file.path(OUT, sprintf("top100_%s.rds", D)))
  P <- r$prof[!is.na(Pertur_id), .(pid, sig = Signature_id, drug = Pertur_id, source = Source, thr, cmap)]
  P <- merge(P, S, by.x = "sig", by.y = "Signature_id", all.x = TRUE); setkey(P, pid)

  # ---- 1. Queries and M (as in 02) -------------------------------------------
  qt <- readRDS(file.path(OUT, sprintf("recovery_%s.rds", D)))$qtab
  Q  <- P[pid %in% qt$pid]
  nthr <- P[, .N, by = thr]; ncm <- P[cmap == TRUE, .N]
  Q[cmap == FALSE, M := nthr[.(Q[cmap == FALSE, thr]), on = "thr", N] + ncm - 1L]
  Q[cmap == TRUE,  M := nrow(P) - 1L]

  # ---- 2. Cross-source comparisons, matched conditions, ranks ----------------
  X <- merge(Q[, .(q = pid, q_sig = sig, drug, s_src = source, q_thr = thr, q_cmap = cmap, q_model = model, q_dv = dose_val, q_df = dose_fam, q_t = time_h, M)],
             P[, .(t = pid, drug, t_sig = sig, t_src = source, t_thr = thr, t_cmap = cmap, model, dose_val, dose_fam, time_h)],
             by = "drug", allow.cartesian = TRUE)[t_src != s_src & (q_cmap | t_cmap | q_thr == t_thr)]
  X <- X[!EX, on = .(q_sig = a, t_sig = b)]
  X[, `:=`(m_model = fifelse(is.na(q_model) | is.na(model), NA, q_model == model),
           m_dose  = fifelse(is.na(q_dv) | is.na(dose_val) | q_df != dose_fam, NA, abs(log2(q_dv / dose_val)) <= log2(DOSE_FOLD) + 1e-9),
           m_time  = fifelse(is.na(q_t) | is.na(time_h), NA, abs(q_t - time_h) < 1e-9))]
  X[CMB, on = .(m_model = model, m_dose = dose, m_time = time), combo := i.combo]
  X[r$TOP, on = .(q, t_sig, t_thr), rnk := i.rnk]

  # ---- 3. Composition ----------------------------------------------------------
  PR[[D]] <- X[, .(n_pair = .N), by = .(combo)][, `:=`(direction = D, share = n_pair / sum(n_pair))]
  PRD[[D]] <- X[s_src < t_src, .(best = if (all(is.na(combo))) "not comparable" else lv[min(match(combo, lv), na.rm = TRUE)]),
                by = .(drug, s_src, t_src)][, .(n_unit = .N), by = .(combo = best)][, `:=`(direction = D, share = n_unit / sum(n_unit))]

  # ---- 4. A. Eight combinations ------------------------------------------------
  Y <- X[!is.na(combo), .(k = .N, first = suppressWarnings(as.numeric(min(rnk, na.rm = TRUE)))), by = .(q, drug, M, combo)]
  Y[is.infinite(first), first := NA]
  Y[, `:=`(h = 1 * (!is.na(first) & first <= TOPN), p = pany(k, M, TOPN))]
  setorder(Y, combo, drug); set.seed(SEED_A)
  for (cb in unique(Y$combo)) {
    y <- Y[combo == cb]
    ix <- sub_index(y$drug, B); oh <- rowMeans(matrix(y$h[ix], B)); ep <- rowMeans(matrix(y$p[ix], B)); fi <- oh / ep
    A[[paste(D, cb)]] <- data.table(direction = D, combo = cb, n_query = nrow(y), n_drug = uniqueN(y$drug),
                                    obs = mean(y$h), exp = mean(y$p), fold_profile = mean(y$h) / mean(y$p),
                                    fold = mean(oh) / mean(ep), fold_lo = unname(quantile(fi, .025)), fold_hi = unname(quantile(fi, .975)), p = NA_real_)
  }

  # ---- 5. B. Pooled by experimental model --------------------------------------
  Z <- X[!is.na(combo)][, mcls := fifelse(m_model, "Model same", "Model different")][
    , .(k = .N, first = suppressWarnings(as.numeric(min(rnk, na.rm = TRUE)))), by = .(q, drug, M, mcls)]
  Z[is.infinite(first), first := NA]
  setorder(Z, mcls, drug); set.seed(SEED_B)
  for (m in unique(Z$mcls)) {
    z <- Z[mcls == m]
    h <- 1 * (!is.na(z$first) & z$first <= TOPN); p <- pany(z$k, z$M, TOPN)
    ix <- sub_index(z$drug, B); oh <- rowMeans(matrix(h[ix], B)); ep <- rowMeans(matrix(p[ix], B))
    dr <- rowMeans(matrix(runif(length(ix)) < p[ix], B))                           # recovery drawn with the expected probabilities
    fi <- oh / ep
    G[[paste(D, m)]] <- data.table(direction = D, combo = m, n_query = nrow(z), n_drug = uniqueN(z$drug),
                                   obs = mean(h), exp = mean(p), fold_profile = mean(h) / mean(p),
                                   fold = mean(oh) / mean(ep), fold_lo = unname(quantile(fi, .025, na.rm = TRUE)), fold_hi = unname(quantile(fi, .975, na.rm = TRUE)),
                                   p = (sum(dr >= oh) + 1) / (B + 1))
  }
  msg("[%s] queries %s | cross-source comparisons %s (comparable for all conditions %s)", D, format(nrow(Q), big.mark = ","),
      format(nrow(X), big.mark = ","), format(X[!is.na(combo), .N], big.mark = ","))
  rm(r, X, Y, Z); gc()
}

# ---- 6. Output -----------------------------------------------------------------
comp <- merge(rbindlist(PR)[, .(direction, combo = fifelse(is.na(combo), "not comparable", combo), n_pair, pct_pair = 100 * share)],
              rbindlist(PRD)[, .(direction, combo, n_unit, pct_unit = 100 * share)], by = c("direction", "combo"), all = TRUE)
pool <- comp[combo != "not comparable", .(n_pair = sum(n_pair), pct_pair = sum(pct_pair), n_unit = sum(n_unit), pct_unit = sum(pct_unit)),
             by = .(direction, combo = fifelse(grepl("^model same", combo), "Model same", "Model different"))]
yn <- function(cb, w) fifelse(grepl(paste(w, "same"), cb), "Yes", "No")
a <- merge(rbindlist(A), comp, by = c("direction", "combo"), all.x = TRUE)[, `:=`(Comparison = "Combination", Model = yn(combo, "model"), Dose = yn(combo, "dose"), Duration = yn(combo, "time"), row = match(combo, lv))]
g <- merge(rbindlist(G), pool, by = c("direction", "combo"), all.x = TRUE)[, `:=`(Comparison = fifelse(combo == "Model same", "Model matched (pooled)", "Model unmatched (pooled)"),
                                                                                  Model = fifelse(combo == "Model same", "Yes", "No"), Dose = "Any", Duration = "Any",
                                                                                  row = fifelse(combo == "Model same", 9L, 10L))]
nc <- comp[combo == "not comparable"][, `:=`(Comparison = "Not comparable", Model = NA_character_, Dose = NA_character_, Duration = NA_character_, row = 11L)]
d <- rbindlist(list(a, g, nc), fill = TRUE)
d[, dir_o := match(direction, c("up", "down"))]; setorder(d, dir_o, row)
tab <- d[, .(Direction = direction, Comparison, Model_matched = Model, Dose_matched = Dose, Duration_matched = Duration,
             Perturbagens = n_drug, Query_profiles = n_query,
             Profile_observed = obs, Profile_expected = exp, Profile_fold = fold_profile,
             Perturbagen_fold = fold, Perturbagen_lower_95 = fold_lo, Perturbagen_upper_95 = fold_hi, Perturbagen_P = p,
             Profile_pairs = n_pair, Profile_pairs_pct = pct_pair, Perturbagen_source_pairs = n_unit, Perturbagen_source_pairs_pct = pct_unit)]
fwrite(tab, file.path(OUT, "recovery_by_matched_conditions.csv"))
print(tab[, .(Direction, Comparison, Model_matched, Dose_matched, Duration_matched, Perturbagens, Perturbagen_fold, Perturbagen_P, Perturbagen_source_pairs_pct)])
msg("saved: %s", file.path(OUT, "recovery_by_matched_conditions.csv"))
