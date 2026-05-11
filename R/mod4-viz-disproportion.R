# =============================================================================
# CellDynamicST — Module 4b: Disproportion Score & Reduction Visualization
# Faithfully generalized from: gene_expression_reduction.py
# Produces:
#   (1) tSNE/UMAP with per-cell-type disproportion scores (coolwarm scale,
#       leader-line labels, symmetric limits centered at 0)
#   (2) Multi-panel tSNE/UMAP colored by different metadata columns
# =============================================================================

#' Calculate Disproportion Scores per Cell Type
#'
#' For a given comparison, computes a per-cell-type disproportion score based
#' on the log2 ratio of cell type proportions between two groups. Scores are
#' z-scored and optionally filtered by significance.
#' Faithfully follows gene_expression_reduction.py logic.
#'
#' @param data Seurat object or data.frame with group and cell type columns.
#' @param group1 Character string for the control group name.
#' @param group2 Character string for the treatment group name.
#' @param group_col Metadata column containing group labels (default
#'   "experiment_group").
#' @param celltype_col Metadata column for cell type annotations (default
#'   "cell_type_annotation").
#' @param min_cells Minimum cells per cell type per group (default 5).
#' @return data.frame with columns: cell_type, n_group1, n_group2,
#'   prop_group1, prop_group2, log2_ratio, zscore.
#' @export
cdst_disproportion_scores <- function(
    data, group1, group2,
    group_col = "experiment_group",
    celltype_col = "cell_type_annotation",
    min_cells = 5
) {
  meta <- if (inherits(data, "Seurat")) data@meta.data
         else if (is.data.frame(data)) data
         else stop("`data` must be a Seurat object or data.frame.", call. = FALSE)

  for (col in c(group_col, celltype_col))
    if (!col %in% colnames(meta))
      stop("Column '", col, "' not found.", call. = FALSE)

  meta_sub <- meta[meta[[group_col]] %in% c(group1, group2), ]
  g1 <- table(meta_sub[[celltype_col]][meta_sub[[group_col]] == group1])
  g2 <- table(meta_sub[[celltype_col]][meta_sub[[group_col]] == group2])
  all_types <- union(names(g1), names(g2))

  res <- data.frame(cell_type = all_types,
                    n_group1 = as.integer(g1[all_types]),
                    n_group2 = as.integer(g2[all_types]),
                    stringsAsFactors = FALSE)
  res$n_group1[is.na(res$n_group1)] <- 0L
  res$n_group2[is.na(res$n_group2)] <- 0L
  res <- res[res$n_group1 >= min_cells | res$n_group2 >= min_cells, ]
  if (nrow(res) == 0) { warning("No cell types pass min_cells."); return(data.frame()) }

  tot1 <- sum(res$n_group1); tot2 <- sum(res$n_group2)
  res$prop_group1 <- res$n_group1 / max(tot1, 1)
  res$prop_group2 <- res$n_group2 / max(tot2, 1)
  pseudo <- 1e-6
  res$log2_ratio <- log2((res$prop_group2 + pseudo) / (res$prop_group1 + pseudo))
  if (nrow(res) >= 3) {
    res$zscore <- (res$log2_ratio - mean(res$log2_ratio)) / stats::sd(res$log2_ratio)
  } else {
    res$zscore <- res$log2_ratio
  }
  res[order(-abs(res$zscore)), ]
}


#' Plot Disproportion Scores on tSNE/UMAP
#'
#' Faithfully reproduces the visualization from gene_expression_reduction.py:
#' - All cells plotted; scored cells colored by coolwarm diverging scale
#' - Unscored cells shown as light grey background
#' - Cell type labels with leader lines (ggrepel) for significant types
#' - Symmetric color limits centered at 0
#'
#' @param data Seurat object or data.frame with reduction coordinates.
#' @param scores data.frame from \code{cdst_disproportion_scores} with
#'   columns \code{cell_type} and \code{zscore}.
#' @param celltype_col Cell type column name (default "cell_type_annotation").
#' @param reduction_x X coordinate column (default "tSNE_1").
#' @param reduction_y Y coordinate column (default "tSNE_2").
#' @param reduction Seurat reduction name (default "tsne").
#' @param score_col Score column to map (default "zscore").
#' @param comparison_label Title label (default "").
#' @param point_size Point size (default 0.8).
#' @param label_size Label font size (default 3.5).
#' @param low_color Negative score color (default "#3B4CC0", blue).
#' @param mid_color Zero score color (default "grey85").
#' @param high_color Positive score color (default "#B40426", red).
#' @param bg_color Background cell color (default "grey90").
#' @param bg_alpha Background cell alpha (default 0.25).
#' @param output_file Path to save (NULL = return ggplot).
#' @param width Plot width inches (default 10).
#' @param height Plot height inches (default 8).
#' @return ggplot object.
#' @export
cdst_plot_disproportion <- function(
    data, scores,
    celltype_col = "cell_type_annotation",
    reduction_x = "tSNE_1", reduction_y = "tSNE_2",
    reduction = "tsne", score_col = "zscore",
    comparison_label = "",
    point_size = 0.8, label_size = 3.5,
    low_color = "#3B4CC0", mid_color = "grey85", high_color = "#B40426",
    bg_color = "grey90", bg_alpha = 0.25,
    output_file = NULL, width = 10, height = 8
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("ggrepel required.")

  # Extract coordinates
  if (inherits(data, "Seurat")) {
    red <- Seurat::Embeddings(data, reduction = reduction)
    meta <- data@meta.data
    plot_df <- data.frame(dim1 = red[, 1], dim2 = red[, 2],
                          cell_type = meta[[celltype_col]],
                          stringsAsFactors = FALSE)
  } else {
    plot_df <- data.frame(dim1 = data[[reduction_x]], dim2 = data[[reduction_y]],
                          cell_type = data[[celltype_col]],
                          stringsAsFactors = FALSE)
  }

  # Map scores to cells (matching Python: fillna(0))
  score_map <- stats::setNames(scores[[score_col]], scores$cell_type)
  plot_df$score <- score_map[plot_df$cell_type]
  plot_df$has_score <- !is.na(plot_df$score)
  plot_df$score[is.na(plot_df$score)] <- 0

  bg_df <- plot_df[!plot_df$has_score, ]
  fg_df <- plot_df[plot_df$has_score, ]

  # Symmetric limits

  vmax <- max(abs(fg_df$score), na.rm = TRUE)
  if (is.na(vmax) || vmax == 0) vmax <- 1

  # Label positions (centroids of scored cell types)
  sig_types <- scores$cell_type[abs(scores[[score_col]]) > 0]
  label_df <- do.call(rbind, lapply(sig_types, function(ct) {
    sub <- fg_df[fg_df$cell_type == ct, ]
    if (nrow(sub) == 0) return(NULL)
    data.frame(x = mean(sub$dim1), y = mean(sub$dim2), label = ct,
               stringsAsFactors = FALSE)
  }))

  p <- ggplot2::ggplot() +
    # Background cells
    ggplot2::geom_point(data = bg_df,
      ggplot2::aes(x = dim1, y = dim2),
      color = bg_color, size = point_size * 0.5, alpha = bg_alpha) +
    # Scored cells with coolwarm scale
    ggplot2::geom_point(data = fg_df,
      ggplot2::aes(x = dim1, y = dim2, color = score),
      size = point_size, alpha = 0.85) +
    ggplot2::scale_color_gradient2(
      low = low_color, mid = mid_color, high = high_color,
      midpoint = 0, limits = c(-vmax, vmax),
      name = "Disproportion\nScore (z)") +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      legend.title = ggplot2::element_text(face = "bold", size = 14),
      legend.text = ggplot2::element_text(size = 12),
      legend.position = "right") +
    ggplot2::labs(title = paste0("Disproportion Scores: ", comparison_label))

  # Leader-line labels
  if (!is.null(label_df) && nrow(label_df) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(x = x, y = y, label = label),
      size = label_size, fontface = "bold",
      segment.color = "grey40", segment.size = 0.4,
      segment.alpha = 0.7,
      max.overlaps = 30, force = 2, force_pull = 0.5,
      box.padding = 0.5, point.padding = 0.3,
      min.segment.length = 0)
  }

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height, dpi = 300)
    message("Disproportion plot saved to: ", output_file)
    return(invisible(p))
  }
  p
}


#' Plot Multi-Panel tSNE/UMAP Colored by Metadata
#'
#' Generates a multi-panel reduction plot where each panel colors cells by a
#' different metadata column (cell_type, genotype, treatment, region, etc.).
#' Matches the multi-panel tSNE style from manuscript Fig 1B-E.
#'
#' @param metadata data.frame with coordinate and metadata columns.
#' @param coord_cols Character vector of length 2 for x, y coordinates.
#' @param color_by Character vector of column names (one panel each).
#' @param color_maps Optional named list of named vectors: col -> (val -> color).
#' @param point_size Point size (default 0.4).
#' @param ncol Panel grid columns (default 2).
#' @param title Overall title (default NULL).
#' @param output_file Path to save (NULL = return ggplot).
#' @param width Plot width inches (default 14).
#' @param height Plot height inches (default 10).
#' @return ggplot object (patchwork composite or list).
#' @export
cdst_plot_reduction_panels <- function(
    metadata,
    coord_cols = c("tSNE_1", "tSNE_2"),
    color_by = c("cell_type_annotation", "genotype", "treatment", "brain_region"),
    color_maps = NULL,
    point_size = 0.4, ncol = 2, title = NULL,
    output_file = NULL, width = 14, height = 10
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")

  color_by <- intersect(color_by, colnames(metadata))
  if (length(color_by) == 0) stop("No valid color_by columns in metadata.")

  plots <- list()
  for (col in color_by) {
    df <- metadata[sample(nrow(metadata)), ]  # shuffle for overplotting
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data[[coord_cols[1]]], y = .data[[coord_cols[2]]],
      color = .data[[col]])) +
      ggplot2::geom_point(size = point_size, alpha = 0.6) +
      ggplot2::theme_void(base_size = 14) +
      ggplot2::theme(
        legend.title = ggplot2::element_text(face = "bold", size = 12),
        legend.text = ggplot2::element_text(size = 10),
        legend.key.size = ggplot2::unit(0.4, "cm"),
        plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0.5)) +
      ggplot2::labs(title = col, color = col)

    if (!is.null(color_maps) && col %in% names(color_maps)) {
      p <- p + ggplot2::scale_color_manual(values = color_maps[[col]])
    } else {
      n_vals <- length(unique(df[[col]]))
      if (n_vals <= 25) {
        p <- p + ggplot2::scale_color_manual(
          values = grDevices::hcl.colors(n_vals, palette = "Dynamic"))
      }
    }

    if (length(unique(df[[col]])) > 15) {
      p <- p + ggplot2::guides(color = ggplot2::guide_legend(
        ncol = 2, override.aes = list(size = 2, alpha = 1)))
    }
    plots[[col]] <- p
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(plots, ncol = ncol)
    if (!is.null(title))
      combined <- combined + patchwork::plot_annotation(title = title)
  } else {
    message("Install 'patchwork' for multi-panel layout. Returning list.")
    return(plots)
  }

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, combined, width = width, height = height, dpi = 300)
    message("Multi-panel reduction plot saved to: ", output_file)
    return(invisible(combined))
  }
  combined
}


#' Run Full Disproportion Analysis and Visualization
#'
#' Convenience wrapper that computes scores and generates plots for multiple
#' comparisons.
#'
#' @param data Seurat object or data.frame.
#' @param comparisons List of length-2 character vectors (group1, group2).
#' @param group_col Group column (default "experiment_group").
#' @param celltype_col Cell type column (default "cell_type_annotation").
#' @param reduction_x X coordinate column (default "tSNE_1").
#' @param reduction_y Y coordinate column (default "tSNE_2").
#' @param output_dir Output directory (default "disproportion_plots").
#' @param ... Additional arguments to \code{cdst_disproportion_scores}.
#' @return Named list of score data.frames.
#' @export
cdst_run_disproportion <- function(
    data, comparisons,
    group_col = "experiment_group",
    celltype_col = "cell_type_annotation",
    reduction_x = "tSNE_1", reduction_y = "tSNE_2",
    output_dir = "disproportion_plots", ...
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  results <- list()
  for (comp in comparisons) {
    label <- paste(comp[1], "vs", comp[2])
    message("\n=== Processing: ", label, " ===")
    scores <- cdst_disproportion_scores(data, comp[1], comp[2],
                                         group_col, celltype_col, ...)
    if (nrow(scores) == 0) { warning("No scores for ", label); next }
    safe <- gsub("[^A-Za-z0-9]", "_", label)
    cdst_plot_disproportion(data, scores, celltype_col, reduction_x, reduction_y,
                            comparison_label = label,
                            output_file = file.path(output_dir,
                                                     paste0("disproportion_", safe, ".png")))
    results[[label]] <- scores
  }
  results
}
