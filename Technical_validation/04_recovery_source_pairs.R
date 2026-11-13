# =============================================================================
# 04. Cross-source recovery by pair of data sources
#
#   Splits the cross-source recovery of 02 by (query source s -> target source t).
#   Queries, top-N lists, hit definition and random expectation are as in 02/03.
#   For each pair (s -> t):
#     queries  = queries from s with at least one comparable profile of the same
#                perturbagen from t (k_t > 0)
#     hit      = a profile of the same perturbagen from t among the top N
#     expected = mean over queries of 1 - C(M - k_t, N) / C(M, N); fold = observed / expected
#   Profile level     : all queries. Empirical P = (b + 1) / (B + 1), b = number of
#                       random iterations (a hit drawn for each query with
#                       probability p) with at least as many hits as observed.
#   Perturbagen level : one query profile sampled per perturbagen in each of B
#                       iterations, as in 03.
#   Signature pairs generated from the same raw data are excluded from k_t and
#   from observed hits, as in 02.
#
#   Input   : <OUT>/top100_{up,down}.rds (01), <OUT>/recovery_{up,down}.rds (02), SHARED_FILE
#   Output  : <OUT>/recovery_by_source_pair.csv
# =============================================================================
rm(list = ls()); gc()
suppressPackageStartupMessages(library(data.table))
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), sprintf(...), "\n")

## ---- Settings: edit the paths ----------------------------------------------
SHARED_FILE <- "shared_raw_signature_pairs.csv"     # columns cvc1, cvc2
OUT         <- "output"
## ---------------------------------------------------------------------------
TOPS  <- c(10L, 30L, 50L, 100L); k4 <- length(TOPS)
B     <- 10000L
SEED  <- 20261016
CHUNK <- 250L                                     # iterations per block in the profile-level simulation (memory)
pany <- function(k, M, n) ifelse(k <= 0, 0, 1 - exp(lchoose(M - k, n) - lchoose(M, n)))

sh <- fread(SHARED_FILE, colClasses = "character")
EX <- unique(rbind(sh[, .(a = cvc1, b = cvc2)], sh[, .(a = cvc2, b = cvc1)]))
A <- list(); BB <- list()
for (D in c("up", "down")) {
  r <- readRDS(file.path(OUT, sprintf("top100_%s.rds", D)))
  P <- r$prof[!is.na(Pertur_id), .(pid, sig = Signature_id, drug = Pertur_id, source = Source, thr, cmap)]; setkey(P, pid)

  # ---- 0. Excluded (query, partner) profile pairs (as in 02) ----------------
  XP <- merge(P[, .(q = pid, a = sig, dq = drug, sq = source, hq = thr, cq = cmap)], EX, by = "a", allow.cartesian = TRUE)
  XP <- merge(XP, P[, .(t = pid, b = sig, dt = drug, st = source, ht = thr, ct = cmap)], by = "b", allow.cartesian = TRUE)
  XP <- unique(XP[dq == dt & sq != st & (cq | ct | hq == ht), .(q, t, t_src = st)])

  # ---- 1. Queries (as in 02) -------------------------------------------------
  qt  <- readRDS(file.path(OUT, sprintf("recovery_%s.rds", D)))$qtab
  Q   <- P[pid %in% qt$pid]
  nthr <- P[, .N, by = thr]; ncm <- P[cmap == TRUE, .N]
  Q[cmap == FALSE, M := nthr[.(Q[cmap == FALSE, thr]), on = "thr", N] + ncm - 1L]
  Q[cmap == TRUE,  M := nrow(P) - 1L]

  # ---- 2. Query x target source, number of same-perturbagen profiles k_t -----
  ds <- P[, .(n = .N), by = .(drug, thr, source)]
  QL <- Q[cmap == FALSE, .(pid, drug, thr)]; QC <- Q[cmap == TRUE, .(pid, drug)]
  K <- rbind(ds[QL, on = .(drug, thr), nomatch = 0L, allow.cartesian = TRUE, .(pid = i.pid, t_src = x.source, n = x.n)],                                  # same threshold
             ds[thr == "CMap"][QL, on = "drug", nomatch = 0L, allow.cartesian = TRUE, .(pid = i.pid, t_src = x.source, n = x.n)],  # + CMap
             ds[QC, on = "drug", nomatch = 0L, allow.cartesian = TRUE, .(pid = i.pid, t_src = x.source, n = x.n)])               # CMap query: all
  K <- K[, .(k = sum(n)), by = .(pid, t_src)]
  K[XP[, .(nx = .N), by = .(pid = q, t_src)], on = .(pid, t_src), k := k - i.nx]
  K <- K[k > 0]
  C <- merge(K, Q[, .(pid, drug, s_src = source, M)], by = "pid")[t_src != s_src]

  # ---- 3. Observed: first rank of a same-perturbagen profile from t ---------
  T1 <- r$TOP[q %in% Q$pid]
  T1[, `:=`(dq = P[.(q), drug], sq = P[.(q), source])]
  T1[r$sig, on = .(t_sig = Signature_id), `:=`(dt = i.Pertur_id, st = i.Source)]
  T1[P, on = .(t_sig = sig, t_thr = thr), t := i.pid]
  T1[, ex := FALSE][XP, on = .(q, t), ex := TRUE]
  fr <- T1[!ex & !is.na(dt) & dt == dq & st != sq, .(first = min(rnk)), by = .(pid = q, t_src = st)]
  C[fr, on = .(pid, t_src), first := i.first]
  for (N in TOPS) {
    set(C, j = paste0("h", N), value = !is.na(C$first) & C$first <= N)
    set(C, j = paste0("p", N), value = pany(C$k, C$M, N))
  }

  # ---- 4. Profile level --------------------------------------------------------
  A[[D]] <- rbindlist(lapply(TOPS, function(N) {
    h <- paste0("h", N); p <- paste0("p", N)
    C[, .(n_drug = uniqueN(drug), n_query = .N, obs = mean(get(h)), exp = mean(get(p))), by = .(s_src, t_src)][, `:=`(direction = D, topN = N)]
  }))
  C[, cell := .GRP, by = .(s_src, t_src)]
  cid <- C[, .(cell = cell[1]), by = .(s_src, t_src)]
  Hm <- as.matrix(C[, paste0("h", TOPS), with = FALSE]) * 1; Pm <- as.matrix(C[, paste0("p", TOPS), with = FALSE])
  nh <- rowsum(Hm, C$cell)                                                          # observed hits per pair x top N
  bc <- matrix(0, nrow(nh), k4)
  set.seed(SEED + 1L)
  for (b0 in seq(1L, B, by = CHUNK)) {
    nb <- min(CHUNK, B - b0 + 1L)
    for (j in seq_len(k4)) {
      U <- matrix(runif(nrow(C) * nb), nrow(C)) < Pm[, j]
      S <- rowsum(U * 1, C$cell)
      bc[, j] <- bc[, j] + rowSums(S >= nh[, j])
    }
  }
  PE <- rbindlist(lapply(seq_len(k4), function(j) data.table(cell = as.integer(rownames(nh)), topN = TOPS[j], p = (bc[, j] + 1) / (B + 1))))
  PE <- merge(PE, cid, by = "cell")[, cell := NULL]
  A[[D]][PE, on = .(s_src, t_src, topN), p := i.p]
  rm(U, S); gc()

  # ---- 5. Perturbagen level ----------------------------------------------------
  cells <- C[, .(n_drug = uniqueN(drug)), by = .(s_src, t_src)]
  set.seed(SEED)
  for (i in seq_len(nrow(cells))) {
    x  <- C[s_src == cells$s_src[i] & t_src == cells$t_src[i]]; setorder(x, drug)
    g  <- x[, .(st = .I[1], sz = .N), by = drug]; ng <- nrow(g)
    Hm <- as.matrix(x[, paste0("h", TOPS), with = FALSE]) * 1; Pm <- as.matrix(x[, paste0("p", TOPS), with = FALSE])
    it <- matrix(0, B, 3 * k4)
    for (b in seq_len(B)) {
      ix <- g$st + floor(runif(ng) * g$sz)
      pp <- Pm[ix, , drop = FALSE]
      it[b, ] <- c(colMeans(Hm[ix, , drop = FALSE]), colMeans(pp), colMeans(matrix(runif(length(pp)), ng) < pp))
    }
    ob <- it[, 1:k4, drop = FALSE]; ex <- it[, k4 + 1:k4, drop = FALSE]; dr <- it[, 2 * k4 + 1:k4, drop = FALSE]
    BB[[paste(D, i)]] <- data.table(direction = D, s_src = cells$s_src[i], t_src = cells$t_src[i], topN = TOPS,
                                    obs = colMeans(ob), obs_lo = apply(ob, 2, quantile, .025), obs_hi = apply(ob, 2, quantile, .975),
                                    exp = colMeans(ex), fold = colMeans(ob) / colMeans(ex), p = (colSums(dr >= ob) + 1) / (B + 1))
  }
  msg("[%s] queries %s | source pairs %d", D, format(nrow(Q), big.mark = ","), nrow(cells))
  rm(r, T1); gc()
}

# ---- 6. Output -----------------------------------------------------------------
a <- rbindlist(A); b <- rbindlist(BB)
d <- merge(a, b, by = c("direction", "topN", "s_src", "t_src"), suffixes = c("_pf", "_pt"), all.x = TRUE)
d[, dir_o := match(direction, c("up", "down"))]; setorder(d, dir_o, -n_drug, s_src, t_src, topN)
tab <- d[, .(Direction = direction, Query_source = s_src, Target_source = t_src, Shared_perturbagens = n_drug, Query_profiles = n_query, Top_N = topN,
             Profile_observed = obs_pf, Profile_expected = exp_pf, Profile_fold = obs_pf / exp_pf, Profile_P = p_pf,
             Perturbagen_observed = obs_pt, Perturbagen_lower_95 = obs_lo, Perturbagen_upper_95 = obs_hi,
             Perturbagen_expected = exp_pt, Perturbagen_fold = fold, Perturbagen_P = p_pt)]
fwrite(tab, file.path(OUT, "recovery_by_source_pair.csv"))
msg("saved: %s (%d rows)", file.path(OUT, "recovery_by_source_pair.csv"), nrow(tab))
