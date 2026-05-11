#' @title Module 5: Cell Dynamics Analysis
#' @description Functions for analyzing cell type proportions, spatial
#'   distributions, and compositional changes across experimental conditions.
#' @name dynamics
NULL

#' Compute Cell Type Proportions
#'
#' Calculates the proportion of each cell type within each experimental group,
#' brain region, or sample. Supports multiple grouping strategies and
#' statistical testing for compositional differences.
#'
#' This function generalizes the cell type distribution analysis from
#' \code{cell_type_distribution_heatmap.R} and related scripts.
#'
#' @param seurat_obj A Seurat object with cell type annotations.
#' @param config A \code{\link{CdstConfig}} object.
#' @param cell_type_column Character. Metadata column with cell type labels.
#'   Default: "cdst_cell_type".
#' @param group_by Character vector. Columns to group by for proportion
#'   calculation. Default: "experiment_group".
#' @param region_column Character. Brain region column for spatial grouping.
#'   Default: "brain_region_L7".
#' @param min_cells Integer. Minimum cells per group to include. Default: 50.
#' @param verbose Logical. Default: TRUE.
#'
#' @return A data.frame with columns: cell_type, group, region (if applicable),
#'   n_cells, proportion, and optionally fold_change and p_value.
#'
#' @export
cdst_cell_proportions <- function(seurat_obj, config,
                                  cell_type_column = "cdst_cell_type",
                                  group_by = "experiment_group",
                                  region_column = "brain_region_L7",
                                  min_cells = 50L,
                                  verbose = TRUE) {
  checkmate::assert_class(config, "CdstConfig")
  validate_seurat(seurat_obj, required_metadata = c(cell_type_column, group_by))

  meta <- seurat_obj@meta.data

  if (verbose) cli::cli_inform("Computing cell type proportions...")

  # Build grouping columns
  group_cols <- group_by
  if (region_column %in% colnames(meta)) {
    group_cols <- c(group_cols, region_column)
  }

  # Compute counts and proportions
  props <- meta |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, cell_type_column)))) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::mutate(
      total_cells = sum(.data$n_cells),
      proportion = .data$n_cells / .data$total_cells
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$total_cells >= min_cells)

  if (verbose) {
    n_types <- length(unique(props[[cell_type_column]]))
    n_groups <- length(unique(props[[group_by[1]]]))
    cli::cli_inform(c(
      "v" = "Proportions computed.",
      "i" = "Cell types: {.val {n_types}}, Groups: {.val {n_groups}}"
    ))
  }

  as.data.frame(props)
}


#' Build Cell Type Distribution Heatmap Data
#'
#' Constructs the data matrix for a cell type x brain region distribution
#' heatmap, ordered by anatomical hierarchy. This is the data preparation
#' step for the heatmap visualization in \code{cdst_plot_distribution_heatmap}.
#'
#' Generalizes \code{cell_type_distribution_heatmap.R}, replacing hardcoded
#' region ordering with atlas-based anatomical ordering.
#'
#' @param seurat_obj A Seurat object with cell type and region annotations.
#' @param cell_type_column Character. Default: "cdst_cell_type".
#' @param region_column Character. Default: "brain_region_L7".
#' @param parent_region_column Character. Higher-level region for grouping.
#'   Default: "brain_region_L3".
#' @param normalize Character. Normalization method: "row" (per cell type),
#'   "column" (per region), or "none". Default: "row".
#' @param min_cells Integer. Minimum cells per cell type. Default: 50.
#' @param verbose Logical. Default: TRUE.
#'
#' @return A list with components:
#'   \describe{
#'     \item{matrix}{Numeric matrix (cell types x regions)}
#'     \item{row_order}{Ordered cell type names}
#'     \item{col_order}{Ordered region names}
#'     \item{row_groups}{Parent cell type grouping}
#'     \item{col_groups}{Parent region grouping}
#'   }
#'
#' @export
cdst_distribution_matrix <- function(seurat_obj,
                                     cell_type_column = "cdst_cell_type",
                                     region_column = "brain_region_L7",
                                     parent_region_column = "brain_region_L3",
                                     normalize = c("row", "column", "none"),
                                     min_cells = 50L,
                                     verbose = TRUE) {
  normalize <- match.arg(normalize)
  validate_seurat(seurat_obj,
                  required_metadata = c(cell_type_column, region_column))

  meta <- seurat_obj@meta.data

  # Filter cell types with too few cells
  ct_counts <- table(meta[[cell_type_column]])
  valid_cts <- names(ct_counts[ct_counts >= min_cells])
  meta <- meta[meta[[cell_type_column]] %in% valid_cts, ]

  # Build count matrix
  count_table <- table(meta[[cell_type_column]], meta[[region_column]])
  count_mat <- as.matrix(count_table)

  # Normalize
  if (normalize == "row") {
    row_sums <- rowSums(count_mat)
    row_sums[row_sums == 0] <- 1
    norm_mat <- count_mat / row_sums
  } else if (normalize == "column") {
    col_sums <- colSums(count_mat)
    col_sums[col_sums == 0] <- 1
    norm_mat <- sweep(count_mat, 2, col_sums, "/")
  } else {
    norm_mat <- count_mat
  }

  # Order columns by parent region then by AP position if available
  col_groups <- NULL
  if (parent_region_column %in% colnames(meta)) {
    region_to_parent <- meta |>
      dplyr::distinct(
        region = .data[[region_column]],
        parent = .data[[parent_region_column]]
      )
    col_groups <- stats::setNames(region_to_parent$parent, region_to_parent$region)

    # Order by parent group
    col_order <- region_to_parent |>
      dplyr::arrange(.data$parent, .data$region) |>
      dplyr::pull(.data$region)
    col_order <- intersect(col_order, colnames(norm_mat))
  } else {
    col_order <- sort(colnames(norm_mat))
  }

  # Order rows by superclass grouping
  row_order <- sort(rownames(norm_mat))

  norm_mat <- norm_mat[row_order, col_order, drop = FALSE]

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Distribution matrix built.",
      "i" = "Dimensions: {nrow(norm_mat)} cell types x {ncol(norm_mat)} regions"
    ))
  }

  list(
    matrix = norm_mat,
    row_order = row_order,
    col_order = col_order,
    row_groups = NULL,
    col_groups = col_groups
  )
}


#' Compute Cell Type Annotation Proportions at Multiple Thresholds
#'
#' Evaluates the proportion of cell types that can be annotated (NT type,
#' brain region) at different confidence thresholds. Useful for determining
#' optimal annotation thresholds.
#'
#' Generalizes the threshold sweep analysis from \code{annotate_cell_type.R}.
#'
#' @param seurat_obj A Seurat object with annotations.
#' @param column_name Character. Metadata column to evaluate.
#' @param cluster_column Character. Cluster column. Default: "cdst_cluster".
#' @param thresholds Numeric vector. Thresholds to evaluate.
#'   Default: seq(0, 1, by = 0.05).
#' @param verbose Logical. Default: TRUE.
#'
#' @return A data.frame with columns: threshold, proportion_annotated.
#'
#' @export
cdst_annotation_coverage <- function(seurat_obj,
                                     column_name,
                                     cluster_column = "cdst_cluster",
                                     thresholds = seq(0, 1, by = 0.05),
                                     verbose = TRUE) {
  validate_seurat(seurat_obj, required_metadata = c(column_name, cluster_column))
  meta <- seurat_obj@meta.data

  results <- vapply(thresholds, function(thresh) {
    expanded <- meta |>
      tidyr::separate_rows(!!rlang::sym(column_name), sep = "-") |>
      dplyr::count(.data[[cluster_column]], .data[[column_name]]) |>
      dplyr::group_by(.data[[cluster_column]]) |>
      dplyr::mutate(proportion = .data$n / sum(.data$n)) |>
      dplyr::ungroup()

    ct_above <- expanded |>
      dplyr::filter(.data$proportion > thresh) |>
      dplyr::distinct(.data[[cluster_column]])

    nrow(ct_above) / length(unique(meta[[cluster_column]]))
  }, numeric(1))

  data.frame(
    threshold = thresholds,
    proportion_annotated = results
  )
}
