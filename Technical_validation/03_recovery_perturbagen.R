# =============================================================================
# 03. Same-perturbagen recovery at the perturbagen level
#
#   Uses the per-query table from 02 (queries, observed hits, random-expectation
#   probabilities p, exclusion of signature pairs generated from the same raw
#   data). Only the unit of aggregation changes.
#   In each of B iterations, one query profile is sampled per perturbagen and
#   observed recovery, expected recovery and fold enrichment are calculated.
#   Observed recovery is the mean over iterations, with the 2.5th and 97.5th
#   percentiles as the 95% range.
#   Random iterations and P : in each iteration, a hit is drawn for each sampled
#   query with probability p. Empirical P = (b + 1) / (B + 1), b = number of
#   iterations in which this random recovery is at least the observed recovery.
#
#   Input   : <OUT>/recovery_{up,down}.rds (02)
#   Output  : <OUT>/recovery_perturbagen_level.csv
#             <OUT>/recovery_perturbagen_level_iterations.csv
# =============================================================================
rm(list = ls()); gc()
suppressPackageStartupMessages(library(data.table))
msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S]"), sprintf(...), "\n")

## ---- Settings: edit the path -----------------------------------------------
OUT  <- "output"
## ---------------------------------------------------------------------------
B    <- 10000L
SEED <- 20261015

SUB <- list(); IT <- list()
for (D in c("up", "down")) {
  x <- readRDS(file.path(OUT, sprintf("recovery_%s.rds", D)))
  Q <- copy(x$qtab); TOPS <- x$tops

  set.seed(SEED)
  gi <- split(seq_len(nrow(Q)), Q$drug); ng <- length(gi)
  pick1 <- function(ix) if (length(ix) == 1L) ix else ix[sample.int(length(ix), 1L)]
  it <- rbindlist(lapply(seq_len(B), function(b) {
    ix <- vapply(gi, pick1, 0L)
    rbindlist(lapply(TOPS, function(N) {
      pa <- Q[[paste0("pa", N)]][ix]; pc <- Q[[paste0("pc", N)]][ix]
      data.table(b = b, topN = N, any = mean(Q[[paste0("ha", N)]][ix]), null_any = mean(pa), draw_any = mean(runif(ng) < pa),
                 cross = mean(Q[[paste0("hc", N)]][ix]), null_cross = mean(pc), draw_cross = mean(runif(ng) < pc))
    }))
  }))
  it[, direction := D]
  IT[[D]] <- it
  SUB[[D]] <- it[, .(n_unit = ng, any = mean(any), any_lo = quantile(any, .025), any_hi = quantile(any, .975),
                     null_any = mean(null_any), fold_any = mean(any / null_any), p_any = (sum(draw_any >= any) + 1) / (B + 1),
                     cross = mean(cross), cross_lo = quantile(cross, .025), cross_hi = quantile(cross, .975),
                     null_cross = mean(null_cross), fold_cross = mean(cross / null_cross), p_cross = (sum(draw_cross >= cross) + 1) / (B + 1)),
                 by = .(direction, topN)]
  msg("[%s] queries %s, perturbagens %s", D, format(nrow(Q), big.mark = ","), format(ng, big.mark = ","))
}

sub <- rbindlist(SUB)
tab <- sub[, .(Direction = direction, Top_N = topN, Perturbagens = n_unit,
               All_source_observed = any, All_source_lower_95 = any_lo, All_source_upper_95 = any_hi,
               All_source_expected = null_any, All_source_fold = fold_any, All_source_P = p_any,
               Cross_source_observed = cross, Cross_source_lower_95 = cross_lo, Cross_source_upper_95 = cross_hi,
               Cross_source_expected = null_cross, Cross_source_fold = fold_cross, Cross_source_P = p_cross)]
fwrite(tab, file.path(OUT, "recovery_perturbagen_level.csv"))
fwrite(rbindlist(IT)[, .(direction, topN, iteration = b, all_source_observed = any, all_source_random = draw_any,
                         cross_source_observed = cross, cross_source_random = draw_cross)],
       file.path(OUT, "recovery_perturbagen_level_iterations.csv"))   # observed and random recovery in each iteration
print(tab)
msg("saved: %s", file.path(OUT, "recovery_perturbagen_level.csv"))
