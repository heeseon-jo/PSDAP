# ==============================================================================
# Microarray_DEG_analysis.R
#
# Purpose
#   Per-treatment differential expression analysis of microarray expression data with
#   limma (treated vs. matched control)
#
# Model
#   comparison       : case = perturbed samples with several conditions (e.g., perturbagen, dose, time)
#                      control = control samples sharing the conditions of the case samples.
#   input expression : Expression matrix after normalization of raw expression data
#   design           : ~ group (levels "ctrl", "case"); coefficient "groupcase" = case vs. ctrl
#   multiple testing : Benjamini-Hochberg (`adj.P.Val`)
#   minimum samples  : MIN_CASE = 2 perturbed and MIN_CTRL = 2 control samples
#
# Dependencies: R 4.4.1; limma
# ==============================================================================

library(limma)


#==============================================================================
# Set minimum count of perturbed, control samples.
MIN_CASE  <- 2
MIN_CTRL  <- 2


#==============================================================================
#' Split a comma-separated string of CEL file names.
#'
#' @param cel_str  Comma-separated CEL file names.
#' @return Character vector of upper-case CEL file names.
parse_cels <- function(cel_str) {
  toupper(trimws(unlist(strsplit(as.character(cel_str), ","))))
}


#=============================================================================
#' Differential expression of one treatment (perturbed vs. control samples).
#'
#' @param expr_mat    Normalized expression matrix, genes x samples
#'                     (colnames = CEL names).
#' @param case_cels   CEL names of the perturbed samples (see `parse_cels`).
#' @param ctrl_cels   CEL names of the control samples (see `parse_cels`).
#' @return List with
#'          `result`               data.frame of limma results for all genes, including BH-adjusted p-values; NULL if skipped
#'          `reason`               "OK" or the skip reason
#'          `n_case` and `n_ctrl`  samples used.
run_limma <- function(expr_mat, case_cels, ctrl_cels) {
  expr_cols <- toupper(colnames(expr_mat))
  case_cols <- intersect(case_cels, expr_cols)
  ctrl_cols <- intersect(ctrl_cels, expr_cols)
  
  if (length(case_cols) < MIN_CASE || length(ctrl_cols) < MIN_CTRL) {
    return(list(
      result = NULL,
      reason = sprintf("Too few samples (case=%d, ctrl=%d)",
                       length(case_cols), length(ctrl_cols)),
      n_case = length(case_cols),
      n_ctrl = length(ctrl_cols)
    ))
  }
  
  sub_mat <- expr_mat[, c(ctrl_cols, case_cols), drop = FALSE]
  group   <- factor(c(rep("ctrl", length(ctrl_cols)),
                      rep("case", length(case_cols))),
                    levels = c("ctrl", "case"))
  design  <- model.matrix(~ group)
  fit     <- lmFit(sub_mat, design)
  fit     <- eBayes(fit)
  res     <- topTable(fit, coef = "groupcase", number = Inf, sort.by = "none")
  res     <- data.frame(gene = rownames(res), res, stringsAsFactors = FALSE)
  rownames(res) <- NULL
  
  list(result = res, reason = "OK",
       n_case = length(case_cols),
       n_ctrl = length(ctrl_cols))
}

