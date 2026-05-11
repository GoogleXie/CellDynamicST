#' @title Shared Utility Functions
#' @description Internal helper functions used across multiple CellDynamicST
#'   modules.
#' @name utils
#' @keywords internal
NULL

#' Create a timestamped output directory
#'
#' @param base_dir Character. Base output directory.
#' @param module_name Character. Module identifier for subdirectory.
#' @param create Logical. Create the directory if it doesn't exist. Default: TRUE.
#'
#' @return Character. Path to the output directory.
#' @keywords internal
cdst_output_dir <- function(base_dir, module_name, create = TRUE) {
  dir_path <- file.path(base_dir, module_name)
  if (create && !dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  dir_path
}


#' Log a provenance entry
#'
#' Records the execution of a pipeline module with timestamp, parameters,
#' and package version for FAIR4RS compliance.
#'
#' @param module_name Character. Name of the module.
#' @param params List. Parameters used.
#' @param start_time POSIXct. When the module started.
#'
#' @return A list with provenance information.
#' @keywords internal
cdst_provenance_entry <- function(module_name, params = list(),
                                  start_time = Sys.time()) {
  list(
    module       = module_name,
    started_at   = format(start_time, "%Y-%m-%d %H:%M:%S %Z"),
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    duration_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
    parameters   = params,
    package_version = as.character(utils::packageVersion("CellDynamicST")),
    r_version    = paste0(R.version$major, ".", R.version$minor),
    platform     = R.version$platform
  )
}


#' Safely subset a Seurat object by group
#'
#' @param seurat_obj A Seurat object.
#' @param group_name Character. The group name to subset.
#' @param group_column Character. Metadata column containing group labels.
#'
#' @return A subsetted Seurat object.
#' @keywords internal
subset_by_group <- function(seurat_obj, group_name,
                            group_column = "experiment_group") {
  cells <- which(seurat_obj@meta.data[[group_column]] == group_name)
  if (length(cells) == 0) {
    cli::cli_abort("No cells found for group {.val {group_name}} in column {.val {group_column}}.")
  }
  seurat_obj[, cells]
}


#' Convert log2CPM expression to binary presence/absence
#'
#' Used by neurotransmitter and glial classification to determine whether
#' a gene is "expressed" in a cell based on a log2CPM threshold.
#'
#' @param expr_matrix A matrix or dgCMatrix of expression values (log-normalized).
#' @param genes Character vector of gene names to evaluate.
#' @param threshold Numeric. Expression threshold (log2CPM scale). Default: 3.
#'
#' @return A logical matrix (genes x cells) indicating expression above threshold.
#' @keywords internal
expression_above_threshold <- function(expr_matrix, genes, threshold = 3) {
  # Filter to available genes
  available <- intersect(genes, rownames(expr_matrix))
  if (length(available) == 0) {
    cli::cli_warn("None of the specified genes found in expression matrix: {.val {genes}}")
    return(matrix(FALSE, nrow = length(genes), ncol = ncol(expr_matrix),
                  dimnames = list(genes, colnames(expr_matrix))))
  }

  result <- as.matrix(expr_matrix[available, , drop = FALSE]) > threshold

  # Add rows of FALSE for missing genes
  missing <- setdiff(genes, available)
  if (length(missing) > 0) {
    missing_mat <- matrix(FALSE, nrow = length(missing), ncol = ncol(expr_matrix),
                          dimnames = list(missing, colnames(expr_matrix)))
    result <- rbind(result, missing_mat)
  }

  result[genes, , drop = FALSE]
}


#' Compute cluster-level summary statistics
#'
#' @param seurat_obj A Seurat object.
#' @param cluster_column Character. Metadata column with cluster assignments.
#' @param group_column Character. Metadata column with group labels.
#'
#' @return A data.frame with cluster-level statistics.
#' @keywords internal
cluster_summary <- function(seurat_obj, cluster_column = "cdst_cluster",
                            group_column = "experiment_group") {
  meta <- seurat_obj@meta.data
  if (!cluster_column %in% colnames(meta)) {
    cli::cli_abort("Cluster column {.val {cluster_column}} not found in metadata.")
  }

  stats <- meta |>
    dplyr::group_by(.data[[cluster_column]]) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_groups = dplyr::n_distinct(.data[[group_column]]),
      mean_nCount = mean(.data[["nCount_RNA"]], na.rm = TRUE),
      mean_nFeature = mean(.data[["nFeature_RNA"]], na.rm = TRUE),
      .groups = "drop"
    )

  as.data.frame(stats)
}


#' Safe file path construction
#'
#' Constructs a file path with optional date prefix and ensures the parent
#' directory exists.
#'
#' @param dir Character. Directory path.
#' @param filename Character. File name (without date prefix).
#' @param ext Character. File extension (without dot).
#' @param date_prefix Logical. Add date prefix. Default: TRUE.
#'
#' @return Character. Full file path.
#' @keywords internal
safe_path <- function(dir, filename, ext, date_prefix = TRUE) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (date_prefix) {
    filename <- paste0(format(Sys.Date(), "%Y%m%d"), "_", filename)
  }

  file.path(dir, paste0(filename, ".", ext))
}
