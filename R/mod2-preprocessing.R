#' @title Module 2: Quality Control and Preprocessing
#' @description Functions for quality control filtering, normalization, and
#'   variable gene selection of spatial transcriptomic data.
#' @name preprocessing
NULL

#' Run Quality Control on Spatial Transcriptomic Data
#'
#' Performs quality control filtering on a Seurat object using configurable
#' thresholds. Filters cells based on transcript count, gene count, and
#' negative probe metrics. Optionally filters FOVs with too few cells.
#'
#' This function generalizes the QC steps from the original CellDynamicST
#' pipeline (generalized_clustering.R lines 40-55), replacing hardcoded
#' percentile cutoffs with configurable parameters from the CdstConfig object.
#'
#' @param seurat_obj A Seurat object.
#' @param config A \code{\link{CdstConfig}} object with QC thresholds.
#' @param verbose Logical. Print QC statistics. Default: TRUE.
#'
#' @return A filtered Seurat object with QC statistics added to metadata.
#'
#' @details
#' The QC pipeline applies the following filters in order:
#' \enumerate{
#'   \item Remove cells below the \code{min_count_percentile} for nCount_RNA.
#'   \item Remove cells below the \code{min_feature_percentile} for nFeature_RNA.
#'   \item Remove cells with fewer than \code{min_transcripts_per_cell} transcripts.
#'   \item Remove cells with negative probe count above \code{max_neg_probe_per_cell}
#'     (if the column exists).
#'   \item Flag FOVs with fewer than \code{min_cells_per_fov} cells.
#' }
#'
#' QC summary statistics are stored in \code{obj@misc$cdst_qc_summary}.
#'
#' @examples
#' \dontrun{
#' config <- cdst_load_config("experiment_config.yaml")
#' obj <- cdst_load_seurat("data.rds", config)
#' obj <- cdst_run_qc(obj, config)
#' }
#'
#' @export
cdst_run_qc <- function(seurat_obj, config, verbose = TRUE) {
  checkmate::assert_class(config, "CdstConfig")
  validate_seurat(seurat_obj)
  start_time <- Sys.time()

  qc <- config@qc
  n_before <- ncol(seurat_obj)

  if (verbose) {
    cli::cli_inform("Starting QC filtering on {.val {n_before}} cells...")
  }

  # ---- Filter by nCount_RNA percentile ----
  counts <- seurat_obj$nCount_RNA
  count_cutoff <- stats::quantile(counts, probs = qc$min_count_percentile,
                                   na.rm = TRUE)
  cells_pass_count <- counts > count_cutoff

  # ---- Filter by nFeature_RNA percentile ----
  features <- seurat_obj$nFeature_RNA
  feature_cutoff <- stats::quantile(features, probs = qc$min_feature_percentile,
                                     na.rm = TRUE)
  cells_pass_feature <- features > feature_cutoff

  # ---- Filter by absolute transcript minimum ----
  cells_pass_min <- counts >= qc$min_transcripts_per_cell

  # ---- Filter by negative probe count (if available) ----
  neg_col <- intersect(
    c("nCount_negprobes", "nCount_NegProbe", "neg_probe_count"),
    colnames(seurat_obj@meta.data)
  )
  if (length(neg_col) > 0) {
    neg_counts <- seurat_obj@meta.data[[neg_col[1]]]
    cells_pass_neg <- neg_counts <= qc$max_neg_probe_per_cell
  } else {
    cells_pass_neg <- rep(TRUE, ncol(seurat_obj))
  }

  # ---- Combine filters ----
  cells_keep <- cells_pass_count & cells_pass_feature &
                cells_pass_min & cells_pass_neg

  seurat_obj <- seurat_obj[, cells_keep]
  n_after <- ncol(seurat_obj)

  # ---- Flag low-cell FOVs ----
  fov_flagged <- character(0)
  if ("fov_id" %in% colnames(seurat_obj@meta.data)) {
    fov_counts <- table(seurat_obj$fov_id)
    low_fovs <- names(fov_counts[fov_counts < qc$min_cells_per_fov])
    if (length(low_fovs) > 0) {
      seurat_obj$cdst_low_fov <- seurat_obj$fov_id %in% low_fovs
      fov_flagged <- low_fovs
      if (verbose) {
        cli::cli_warn("{length(low_fovs)} FOV(s) flagged with < {qc$min_cells_per_fov} cells.")
      }
    }
  }

  # ---- Store QC summary ----
  qc_summary <- list(
    cells_before       = n_before,
    cells_after        = n_after,
    cells_removed      = n_before - n_after,
    pct_removed        = round((n_before - n_after) / n_before * 100, 2),
    count_cutoff       = unname(count_cutoff),
    feature_cutoff     = unname(feature_cutoff),
    fovs_flagged       = fov_flagged,
    parameters         = qc,
    timestamp          = Sys.time()
  )
  seurat_obj@misc$cdst_qc_summary <- qc_summary

  if (verbose) {
    cli::cli_inform(c(
      "v" = "QC complete.",
      "i" = "Cells: {.val {n_before}} -> {.val {n_after}} ({qc_summary$pct_removed}% removed)",
      "i" = "nCount cutoff (p{qc$min_count_percentile * 100}): {.val {round(count_cutoff)}}",
      "i" = "nFeature cutoff (p{qc$min_feature_percentile * 100}): {.val {round(feature_cutoff)}}"
    ))
  }

  seurat_obj
}


#' Find Variable Genes
#'
#' Identifies highly variable genes using Seurat's FindVariableFeatures.
#' This is a thin wrapper that records parameters for provenance tracking.
#'
#' @param seurat_obj A Seurat object (post-QC).
#' @param n_features Integer. Number of variable features to select. Default: 2000.
#' @param selection_method Character. Method for variable feature selection.
#'   Default: "vst".
#' @param verbose Logical. Default: TRUE.
#'
#' @return The Seurat object with variable features identified.
#'
#' @export
cdst_find_variable_genes <- function(seurat_obj,
                                     n_features = 2000L,
                                     selection_method = "vst",
                                     verbose = TRUE) {
  validate_seurat(seurat_obj)

  if (verbose) cli::cli_inform("Finding {.val {n_features}} variable genes...")

  seurat_obj <- Seurat::FindVariableFeatures(
    seurat_obj,
    selection.method = selection_method,
    nfeatures = n_features,
    verbose = FALSE
  )

  if (verbose) {
    n_var <- length(Seurat::VariableFeatures(seurat_obj))
    cli::cli_inform(c("v" = "Identified {.val {n_var}} variable features."))
  }

  seurat_obj
}


#' Normalize and Scale Data
#'
#' Performs log-normalization and scaling of the expression data. Also computes
#' log2CPM normalization used by the neurotransmitter classification module.
#'
#' This function wraps Seurat's NormalizeData and ScaleData, adding the log2CPM
#' assay required by downstream CellDynamicST modules.
#'
#' @param seurat_obj A Seurat object (post-QC).
#' @param normalization_method Character. Normalization method for Seurat.
#'   Default: "LogNormalize".
#' @param scale_factor Numeric. Scale factor for normalization. Default: 10000.
#' @param scale_all Logical. Scale all genes (TRUE) or only variable features
#'   (FALSE). Default: TRUE.
#' @param verbose Logical. Default: TRUE.
#'
#' @return The Seurat object with normalized and scaled data, plus a log2CPM
#'   assay.
#'
#' @export
cdst_normalize <- function(seurat_obj,
                           normalization_method = "LogNormalize",
                           scale_factor = 10000,
                           scale_all = TRUE,
                           verbose = TRUE) {
  validate_seurat(seurat_obj)

  if (verbose) cli::cli_inform("Normalizing data...")

  # Standard Seurat normalization
  seurat_obj <- Seurat::NormalizeData(
    seurat_obj,
    normalization.method = normalization_method,
    scale.factor = scale_factor,
    verbose = FALSE
  )

  # Scale data
  if (scale_all) {
    features <- rownames(seurat_obj)
  } else {
    features <- Seurat::VariableFeatures(seurat_obj)
    if (length(features) == 0) {
      cli::cli_warn("No variable features found. Scaling all genes instead.")
      features <- rownames(seurat_obj)
    }
  }

  if (verbose) cli::cli_inform("Scaling {.val {length(features)}} features...")

  seurat_obj <- Seurat::ScaleData(
    seurat_obj,
    features = features,
    verbose = FALSE
  )

  # Compute log2CPM if not already present
  if (!"log2CPM" %in% Seurat::Assays(seurat_obj)) {
    if (verbose) cli::cli_inform("Computing log2CPM normalization...")
    seurat_obj <- .compute_log2cpm(seurat_obj)
  }

  if (verbose) {
    cli::cli_inform(c("v" = "Normalization and scaling complete."))
  }

  seurat_obj
}
