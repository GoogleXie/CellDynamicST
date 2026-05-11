#' @title Input Validation Utilities
#' @description Internal functions for validating Seurat objects, metadata
#'   columns, and analysis prerequisites before computation.
#' @name validation
#' @keywords internal
NULL

#' Validate that a Seurat object has the required structure
#'
#' Checks that the Seurat object contains the expected assay, data slots,
#' and metadata columns. Throws informative errors if any requirement is
#' not met.
#'
#' @param seurat_obj A Seurat object.
#' @param required_metadata Character vector of required metadata column names.
#' @param required_assay Character. Required assay name. Default: "RNA".
#' @param require_normalized Logical. Check for normalized data slot. Default: FALSE.
#' @param caller Character. Name of the calling function for error messages.
#'
#' @return Invisibly returns TRUE if all checks pass.
#' @keywords internal
validate_seurat <- function(seurat_obj,
                            required_metadata = character(0),
                            required_assay = "RNA",
                            require_normalized = FALSE,
                            caller = "CellDynamicST") {
  # Check class
  if (!inherits(seurat_obj, "Seurat")) {
    cli::cli_abort(c(
      "Expected a Seurat object, got {.cls {class(seurat_obj)}}.",
      "i" = "Load your data with {.fn Seurat::CreateSeuratObject} or {.fn cdst_load_seurat}."
    ), call = rlang::caller_env())
  }

  # Check assay exists
  available_assays <- Seurat::Assays(seurat_obj)
  if (!required_assay %in% available_assays) {
    cli::cli_abort(c(
      "Required assay {.val {required_assay}} not found in Seurat object.",
      "i" = "Available assays: {.val {available_assays}}"
    ), call = rlang::caller_env())
  }

  # Check normalized data
  if (require_normalized) {
    assay_obj <- seurat_obj[[required_assay]]
    if (all(dim(Seurat::GetAssayData(seurat_obj, assay = required_assay,
                                      layer = "data")) == 0)) {
      cli::cli_abort(c(
        "Normalized data not found in assay {.val {required_assay}}.",
        "i" = "Run {.fn cdst_normalize} or {.fn Seurat::NormalizeData} first."
      ), call = rlang::caller_env())
    }
  }

  # Check required metadata columns
  if (length(required_metadata) > 0) {
    meta_cols <- colnames(seurat_obj@meta.data)
    missing <- setdiff(required_metadata, meta_cols)
    if (length(missing) > 0) {
      cli::cli_abort(c(
        "Required metadata column{?s} missing: {.val {missing}}",
        "i" = "Available columns: {.val {head(meta_cols, 20)}}",
        "i" = "Ensure your Seurat object has the required metadata. See {.fn cdst_load_seurat}."
      ), call = rlang::caller_env())
    }
  }

  invisible(TRUE)
}


#' Validate that experimental groups exist in metadata
#'
#' @param seurat_obj A Seurat object.
#' @param config A CdstConfig object.
#' @param group_column Character. Name of the metadata column containing group labels.
#'   Default: "experiment_group".
#'
#' @return Invisibly returns TRUE.
#' @keywords internal
validate_groups <- function(seurat_obj, config,
                            group_column = "experiment_group") {
  validate_seurat(seurat_obj, required_metadata = group_column)

  available <- unique(seurat_obj@meta.data[[group_column]])
  expected <- vapply(config@groups, function(g) g$name, character(1))
  missing <- setdiff(expected, available)

  if (length(missing) > 0) {
    cli::cli_abort(c(
      "Experimental group{?s} not found in column {.val {group_column}}: {.val {missing}}",
      "i" = "Available groups in data: {.val {available}}",
      "i" = "Check your configuration file or metadata column name."
    ), call = rlang::caller_env())
  }

  invisible(TRUE)
}


#' Validate comparison definitions against available groups
#'
#' @param config A CdstConfig object.
#'
#' @return Invisibly returns TRUE.
#' @keywords internal
validate_comparisons <- function(config) {
  group_names <- vapply(config@groups, function(g) g$name, character(1))

  for (comp in config@comparisons) {
    if (!comp$group1 %in% group_names) {
      cli::cli_abort(c(
        "Comparison {.val {comp$name}}: group1 {.val {comp$group1}} not in defined groups.",
        "i" = "Available groups: {.val {group_names}}"
      ), call = rlang::caller_env())
    }
    if (!comp$group2 %in% group_names) {
      cli::cli_abort(c(
        "Comparison {.val {comp$name}}: group2 {.val {comp$group2}} not in defined groups.",
        "i" = "Available groups: {.val {group_names}}"
      ), call = rlang::caller_env())
    }
  }

  invisible(TRUE)
}


#' Check if a package is available, with informative error
#'
#' @param pkg Character. Package name.
#' @param reason Character. Why the package is needed.
#'
#' @return Invisibly returns TRUE.
#' @keywords internal
require_package <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- paste0("Package {.pkg ", pkg, "} is required")
    if (!is.null(reason)) msg <- paste0(msg, " for ", reason)
    msg <- paste0(msg, " but is not installed.")

    cli::cli_abort(c(
      msg,
      "i" = "Install with: {.code install.packages(\"{pkg}\")}"
    ), call = rlang::caller_env())
  }
  invisible(TRUE)
}
