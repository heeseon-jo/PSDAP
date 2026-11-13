# =============================================================================
# 02. Same-perturbagen recovery at the profile level
#
#   Query   : profile with at least one comparable profile of the same
#             perturbagen in another data source (comparable = same threshold,
#             or either profile from CMap, i.e. pairs for which similarity is
#             calculated).
#   Hit     : a profile of the same perturbagen among the top N (10/30/50/100)
#             most similar profiles. All-source: any data source;
#             cross-source: other data sources only.
#   Random expectation : probability that N profiles drawn at random without
#             replacement from the M comparable candidates include at least one
#             of the k same-perturbagen profiles, p = 1 - C(M - k, N) / C(M, N).
#             M: logFC query = same threshold + CMap; CMap query = all profiles;
#             the query itself is excluded.
#             Expected recovery = mean p over queries; fold = observed / expected.
#   Random iterations and P : in each of N_ITER iterations, a hit is drawn for
#             each query with probability p and the hits are summed.
#             Empirical P = (b + 1) / (N_ITER + 1), b = number of iterations with
#             at least as many hits as observed.
#   Signature pairs generated from the same raw data (SHARED_FILE) are not
#             counted as same-perturbagen profiles of each other (observed hits
#             and k); their ranks in the top-N lists are kept.
#
#   Input   : <OUT>/top100_{up,down}.rds (01), SHARED_FILE
#   Output  : <OUT>/recovery_profile_level.csv
#             <OUT>/recovery_profile_level_random_iterations.csv
#             <OUT>/recovery_{up,down}.rds  (intermediate, per-query table read by 03-05)
# =============================================================================
rm(list = ls()); gc()
suppressPackageStartupMessages(library(data.table))
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), sprintf(...), "\n")

## ---- Settings: edit the paths ----------------------------------------------
SHARED_FILE <- "shared_raw_signature_pairs.csv"     # columns cvc1, cvc2
OUT         <- "output"
## ---------------------------------------------------------------------------
TOPS   <- c(10L, 30L, 50L, 100L); k4 <- length(TOPS)
N_ITER <- 10000
SEED   <- 20260709
pany <- function(k, M, n) ifelse(k <= 0, 0, 1 - exp(lchoose(M - k, n) - lchoose(M, n)))   # 1 - C(M-k, n) / C(M, n)
sh <- fread(SHARED_FILE, colClasses = "character")
EX <- unique(rbind(sh[, .(a = cvc1, b = cvc2)], sh[, .(a = cvc2, b = cvc1)]))              # both orientations

RES <- list(); ITER <- list()
for (DIR in c("up", "down")) {
  r <- readRDS(file.path(OUT, sprintf("top100_%s.rds", DIR)))
  P <- r$prof[!is.na(Pertur_id), .(pid, sig = Signature_id, drug = Pertur_id, source = Source, thr, cmap)]; setkey(P, pid)

  # ---- 0. Excluded (query, partner) profile pairs: same raw data, same perturbagen, different source, comparable
  XP <- merge(P[, .(q = pid, a = sig, dq = drug, sq = source, hq = thr, cq = cmap)], EX, by = "a", allow.cartesian = TRUE)
  XP <- merge(XP, P[, .(t = pid, b = sig, dt = drug, st = source, ht = thr, ct = cmap)], by = "b", allow.cartesian = TRUE)
  XP <- unique(XP[dq == dt & sq != st & (cq | ct | hq == ht), .(q, t)])
  nx <- XP[, .(nx = .N), by = .(pid = q)]

  # ---- 1. Queries --------------------------------------------------------
  cnt <- P[, .N, by = .(drug, thr_t = thr, cmap_t = cmap, src_t = source)]
  E   <- merge(P, cnt, by = "drug", allow.cartesian = TRUE)[cmap | cmap_t | thr == thr_t]
  K   <- E[, .(k = sum(N[src_t != source])), by = pid]
  K[nx, on = "pid", k := k - i.nx]
  Q   <- P[K[k > 0], on = "pid"]; nq <- nrow(Q)
  msg("[%s] profiles %s | queries %s | excluded profile pairs %s", DIR, format(nrow(P), big.mark = ","),
      format(nq, big.mark = ","), format(nrow(XP), big.mark = ","))

  # ---- 2. Observed hits ----------------------------------------------------
  T1 <- r$TOP[q %in% Q$pid]
  T1[, `:=`(dq = P[.(q), drug], sq = P[.(q), source])]
  T1[r$sig, on = .(t_sig = Signature_id), `:=`(dt = i.Pertur_id, st = i.Source)]
  T1[P, on = .(t_sig = sig, t_thr = thr), t := i.pid]
  T1[, ex := FALSE][XP, on = .(q, t), ex := TRUE]                                    # rank kept, not counted as a hit
  for (N in TOPS) {
    h <- T1[rnk <= N & !ex, .(a = any(dt == dq, na.rm = TRUE), c = any(dt == dq & st != sq, na.rm = TRUE)), by = q]
    Q[h, on = .(pid = q), c(paste0("ha", N), paste0("hc", N)) := .(i.a, i.c)]
  }
  for (v in grep("^h[ac][0-9]+$", names(Q), value = TRUE)) set(Q, which(is.na(Q[[v]])), v, FALSE)

  # ---- 3. Random expectation: M, k and p for each query --------------------
  nthr <- P[, .N, by = thr]; ncm <- P[cmap == TRUE, .N]
  Q[cmap == FALSE, M := nthr[.(Q[cmap == FALSE, thr]), on = "thr", N] + ncm - 1L]
  Q[cmap == TRUE,  M := nrow(P) - 1L]
  look <- function(tab, keys, on) { v <- tab[keys, on = on, N]; v[is.na(v)] <- 0L; v }
  dthr <- P[, .N, by = .(drug, thr)]; dsrc <- P[, .N, by = .(drug, thr, source)]
  nd_all <- P[, .N, by = drug]; nd_cm <- P[cmap == TRUE, .N, by = drug]
  ns_all <- P[, .N, by = .(drug, source)]; ns_cm <- P[cmap == TRUE, .N, by = .(drug, source)]
  Q[, `:=`(d_h = look(dthr, .(drug = drug, thr = thr), c("drug", "thr")), d_cm = look(nd_cm, .(drug = drug), "drug"), d_all = look(nd_all, .(drug = drug), "drug"),
           s_h = look(dsrc, .(drug = drug, thr = thr, source = source), c("drug", "thr", "source")),
           s_cm = look(ns_cm, .(drug = drug, source = source), c("drug", "source")), s_all = look(ns_all, .(drug = drug, source = source), c("drug", "source")))]
  Q[, n_same := fifelse(cmap, d_all, d_h + d_cm)][, n_src := fifelse(cmap, s_all, s_h + s_cm)]
  Q[nx, on = "pid", n_ex := i.nx][is.na(n_ex), n_ex := 0L]
  Q[, `:=`(k_any = n_same - 1L - n_ex, k_cross = n_same - n_src - n_ex)]
  for (N in TOPS) { set(Q, j = paste0("pa", N), value = pany(Q$k_any, Q$M, N)); set(Q, j = paste0("pc", N), value = pany(Q$k_cross, Q$M, N)) }

  # ---- 4. Random iterations ----------------------------------------------
  set.seed(SEED)
  pa <- as.matrix(Q[, paste0("pa", TOPS), with = FALSE]); pc <- as.matrix(Q[, paste0("pc", TOPS), with = FALSE])
  NL <- matrix(0, N_ITER, 2 * k4)
  for (i in seq_len(N_ITER)) NL[i, ] <- c(colSums(matrix(runif(length(pa)), nq) < pa), colSums(matrix(runif(length(pc)), nq) < pc))

  # ---- 5. Summary ------------------------------------------------------------
  res <- rbindlist(lapply(seq_len(k4), function(k) {
    N <- TOPS[k]; a <- NL[, k]; c_ <- NL[, k4 + k]; oa <- sum(Q[[paste0("ha", N)]]); oc <- sum(Q[[paste0("hc", N)]])
    ea <- sum(pa[, k]); ec <- sum(pc[, k])
    data.table(direction = DIR, topN = N, n_query = nq,
               any = oa / nq, null_any = ea / nq, fold_any = oa / ea, p_any = (sum(a >= oa) + 1) / (N_ITER + 1),
               cross = oc / nq, null_cross = ec / nq, fold_cross = oc / ec, p_cross = (sum(c_ >= oc) + 1) / (N_ITER + 1))
  }))
  RES[[DIR]] <- res
  ITER[[DIR]] <- rbindlist(lapply(seq_len(k4), function(k) data.table(direction = DIR, topN = TOPS[k], iteration = seq_len(N_ITER),
                                                                     all_source = NL[, k] / nq, cross_source = NL[, k4 + k] / nq)))
  qtab <- Q[, c("pid", "sig", "drug", "source", "thr", "cmap", "M", "k_any", "k_cross", paste0(rep(c("ha", "hc", "pa", "pc"), each = k4), TOPS)), with = FALSE]
  saveRDS(list(res = res, n_query = nq, tops = TOPS, qtab = qtab, seed = SEED, n_iter = N_ITER),
          file.path(OUT, sprintf("recovery_%s.rds", DIR)))
  rm(r, T1, E); gc()
}

# ---- 6. Output -----------------------------------------------------------------
res <- rbindlist(RES)
tab <- res[, .(Direction = direction, Top_N = topN, Query_profiles = n_query,
               All_source_observed = any, All_source_expected = null_any, All_source_fold = fold_any, All_source_P = p_any,
               Cross_source_observed = cross, Cross_source_expected = null_cross, Cross_source_fold = fold_cross, Cross_source_P = p_cross)]
fwrite(tab, file.path(OUT, "recovery_profile_level.csv"))
fwrite(rbindlist(ITER), file.path(OUT, "recovery_profile_level_random_iterations.csv"))   # recovery rate in each random iteration
print(tab)
msg("saved: %s", file.path(OUT, "recovery_profile_level.csv"))
