# =============================================================================
# Module 4e: Additional Visualization Functions
# =============================================================================
# Contains: enrichment lollipop, enrichment heatmap, standalone QC plot,
#   and standalone proportions plot (all accept plain data frames).
# =============================================================================

#' Bidirectional GO Enrichment Comparison Plot (2x2 Design)
#'
#' Creates a bidirectional dot-and-bar plot comparing GO enrichment between
#' two conditions, optimized for 2x2 factorial designs (e.g., genotype x
#' treatment). Each GO term has:
#'   - A background bar showing the difference in weighted scores between
#'     condition A and condition B (left = A stronger, right = B stronger).
#'   - Dots on the left (condition A) and right (condition B) showing the
#'     weighted score magnitude, with size = |log2FC| and color = -log10(qvalue).
#'
#' The weighted score for each term is: |−log10(qvalue) × mean(log2FC of genes)|.
#'
#' @param enrichment_A Data frame for condition A with columns:
#'   Description (GO term), qvalue, geneID (slash-separated gene list).
#'   This is the standard output from clusterProfiler enrichGO().
#' @param enrichment_B Data frame for condition B (same format as A).
#' @param log2fc_A Data frame with columns: gene, log2FC for condition A.
#' @param log2fc_B Data frame with columns: gene, log2FC for condition B.
#' @param label_A Label for condition A. Default "Condition A".
#' @param label_B Label for condition B. Default "Condition B".
#' @param top_n Number of top terms per condition to include. Default 15.
#' @param title Plot title.
#' @param output_file Optional file path to save the plot.
#' @param width Plot width. Default 18.
#' @param height Plot height. Default 10.
#' @return A ggplot object.
#' @export
cdst_plot_enrichment_lollipop <- function(
    enrichment_A,
    enrichment_B = NULL,
    log2fc_A = NULL,
    log2fc_B = NULL,
    label_A = "Condition A",
    label_B = "Condition B",
    top_n = 15,
    title = NULL,
    output_file = NULL,
    width = 18,
    height = 10
) {
  # ---- Helper: compute weighted score ----
  .calc_weighted <- function(go_df, l2fc_df) {
    if (is.null(go_df) || nrow(go_df) == 0) return(go_df)
    go_df$log2FC <- NA_real_
    go_df$weighted_score <- NA_real_
    for (i in seq_len(nrow(go_df))) {
      genes_str <- as.character(go_df$geneID[i])
      genes <- unlist(strsplit(genes_str, "/"))
      if (!is.null(l2fc_df) && nrow(l2fc_df) > 0) {
        matched <- l2fc_df$log2FC[l2fc_df$gene %in% genes]
        go_df$log2FC[i] <- if (length(matched) > 0) mean(matched, na.rm = TRUE) else 0
      } else {
        go_df$log2FC[i] <- 0
      }
      qv <- go_df$qvalue[i]
      if (is.na(qv) || qv <= 0) qv <- 1e-300
      go_df$weighted_score[i] <- abs(-log10(qv) * go_df$log2FC[i])
    }
    go_df
  }

  # ---- Single-condition fallback (simple lollipop) ----
  if (is.null(enrichment_B)) {
    df <- enrichment_A
    # Support both clusterProfiler format and simple format
    if (!"go_name" %in% colnames(df) && "Description" %in% colnames(df)) {
      df$go_name <- df$Description
    }
    if (!"p_value" %in% colnames(df) && "qvalue" %in% colnames(df)) {
      df$p_value <- df$qvalue
    }
    if (!"enrichment_score" %in% colnames(df)) {
      if ("geneID" %in% colnames(df) && !is.null(log2fc_A)) {
        df <- .calc_weighted(df, log2fc_A)
        df$enrichment_score <- df$weighted_score
      } else {
        df$enrichment_score <- -log10(df$p_value)
      }
    }
    if (!"n_genes" %in% colnames(df) && "Count" %in% colnames(df)) {
      df$n_genes <- df$Count
    } else if (!"n_genes" %in% colnames(df)) {
      df$n_genes <- 10
    }

    df$abs_score <- abs(df$enrichment_score)
    df <- df[order(-df$abs_score), ]
    df <- utils::head(df, top_n)
    df$go_name <- factor(df$go_name, levels = rev(df$go_name))

    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = enrichment_score, y = go_name, color = -log10(p_value)
    )) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = enrichment_score, y = go_name, yend = go_name),
        linewidth = 0.8
      ) +
      ggplot2::geom_point(ggplot2::aes(size = n_genes), alpha = 0.85) +
      ggplot2::scale_color_viridis_c(option = "plasma", name = "-log10(p)") +
      ggplot2::scale_size_continuous(range = c(3, 8), name = "Gene Count") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(
        title = title %||% "GO Enrichment",
        x = "Enrichment Score",
        y = ""
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = ggplot2::element_text(face = "bold", size = 15),
        axis.text.y = ggplot2::element_text(size = 14),
        axis.text.x = ggplot2::element_text(size = 13),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank()
      )

    if (!is.null(output_file)) {
      ggplot2::ggsave(output_file, p, width = width, height = height, dpi = 300)
      message("Saved enrichment lollipop to: ", output_file)
      return(invisible(p))
    }
    return(p)
  }

  # ---- Two-condition bidirectional comparison (2x2 design) ----
  df_A <- as.data.frame(enrichment_A)
  df_B <- as.data.frame(enrichment_B)

  # Calculate weighted scores
  df_A <- .calc_weighted(df_A, log2fc_A)
  df_B <- .calc_weighted(df_B, log2fc_B)

  # Rank within each condition
  df_A <- df_A[order(-df_A$weighted_score), ]
  df_A$rank_A <- seq_len(nrow(df_A))
  df_B <- df_B[order(-df_B$weighted_score), ]
  df_B$rank_B <- seq_len(nrow(df_B))

  # Full join on Description
  combined <- merge(df_A, df_B, by = "Description", all = TRUE,
                    suffixes = c("_A", "_B"))

  # Select top N from each side
  top_A <- utils::head(combined[order(combined$rank_A), ], top_n)
  top_B <- utils::head(combined[order(combined$rank_B), ], top_n)
  all_terms <- unique(c(top_A$Description, top_B$Description))

  combined <- combined[combined$Description %in% all_terms, ]

  # Fill NAs
  combined$weighted_score_A[is.na(combined$weighted_score_A)] <- 0
  combined$weighted_score_B[is.na(combined$weighted_score_B)] <- 0
  combined$log2FC_A[is.na(combined$log2FC_A)] <- 0
  combined$log2FC_B[is.na(combined$log2FC_B)] <- 0
  combined$qvalue_A[is.na(combined$qvalue_A)] <- 1
  combined$qvalue_B[is.na(combined$qvalue_B)] <- 1

  # Compute difference and combined score for ranking
  combined$diff <- combined$weighted_score_B - combined$weighted_score_A
  combined$combined_score <- pmax(combined$weighted_score_A,
                                   combined$weighted_score_B)
  combined <- combined[order(-combined$combined_score), ]
  combined$rank <- seq_len(nrow(combined))

  # Filter to terms present in both conditions for the bar
  bar_data <- combined[combined$weighted_score_A > 0 &
                         combined$weighted_score_B > 0, ]

  # Dot data for each side
  dot_A <- combined[combined$weighted_score_A > 0, ]
  dot_B <- combined[combined$weighted_score_B > 0, ]

  # Build the plot
  p <- ggplot2::ggplot(combined,
    ggplot2::aes(y = stats::reorder(Description, -rank))) +
    # Background bars showing difference
    ggplot2::geom_bar(
      data = bar_data,
      ggplot2::aes(x = diff, fill = diff > 0),
      stat = "identity", position = "identity",
      alpha = 0.4, color = "grey70", linewidth = 0.3
    ) +
    # Dots for condition A (left side, negative x)
    ggplot2::geom_point(
      data = dot_A,
      ggplot2::aes(
        x = -weighted_score_A,
        size = abs(log2FC_A),
        color = -log10(qvalue_A)
      )
    ) +
    # Dots for condition B (right side, positive x)
    ggplot2::geom_point(
      data = dot_B,
      ggplot2::aes(
        x = weighted_score_B,
        size = abs(log2FC_B),
        color = -log10(qvalue_B)
      )
    ) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#FFA500", "FALSE" = "#87CEEB"),
      labels = c("TRUE" = paste(label_A, "<", label_B),
                 "FALSE" = paste(label_A, ">", label_B)),
      name = "Difference"
    ) +
    ggplot2::scale_color_gradient(
      low = "darkblue", high = "darkorange",
      name = "-log10(qvalue)"
    ) +
    ggplot2::scale_size_continuous(
      range = c(2, 7), name = "|log2FC|"
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid",
                        color = "grey30", linewidth = 0.5) +
    ggplot2::labs(
      title = title %||% paste("GO Enrichment:", label_A, "vs", label_B),
      x = paste0("Weighted Score  (← ", label_A, "  |  ", label_B, " →)"),
      y = ""
    ) +
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold", size = 15),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 13),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = ggplot2::element_text(face = "bold", size = 13),
      legend.text = ggplot2::element_text(size = 12),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(
        color = "grey90", linewidth = 0.3
      )
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height, dpi = 300)
    message("Saved enrichment comparison plot to: ", output_file)
    return(invisible(p))
  }
  p
}


#' Enrichment Heatmap
#'
#' Creates a heatmap of GO enrichment scores across regions for a given
#' contrast, using pheatmap.
#'
#' @param enrichment_results Data frame with columns: go_name,
#'   enrichment_score, region, contrast.
#' @param contrast Which contrast to plot. Default NULL (all).
#' @param top_n Number of top GO terms to include. Default 20.
#' @param title Plot title.
#' @param output_file Optional file path to save the plot.
#' @param width Plot width. Default 12.
#' @param height Plot height. Default 10.
#' @return A pheatmap object.
#' @export
cdst_plot_enrichment_heatmap <- function(
    enrichment_results,
    contrast = NULL,
    top_n = 20,
    title = NULL,
    output_file = NULL,
    width = 12,
    height = 10
) {
  df <- enrichment_results
  if (!is.null(contrast) && "contrast" %in% colnames(df)) {
    df <- df[df$contrast == contrast, ]
  }
  if (nrow(df) == 0) {
    warning("No enrichment results after filtering.")
    return(invisible(NULL))
  }

  # Select top terms by mean absolute score
  term_scores <- stats::aggregate(
    enrichment_score ~ go_name, data = df,
    FUN = function(x) mean(abs(x))
  )
  term_scores <- term_scores[order(-term_scores$enrichment_score), ]
  top_terms <- head(term_scores$go_name, top_n)
  df <- df[df$go_name %in% top_terms, ]

  # Pivot to matrix
  mat <- stats::reshape(
    df[, c("go_name", "region", "enrichment_score")],
    idvar = "go_name", timevar = "region",
    direction = "wide"
  )
  rownames(mat) <- mat$go_name
  mat$go_name <- NULL
  colnames(mat) <- gsub("^enrichment_score\\.", "", colnames(mat))
  mat[is.na(mat)] <- 0
  mat <- as.matrix(mat)

  if (is.null(title)) title <- paste("Enrichment Heatmap:", contrast)

  colors <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)
  max_val <- max(abs(mat), na.rm = TRUE)
  if (max_val == 0) max_val <- 1

  hp <- pheatmap::pheatmap(
    mat, color = colors, main = title,
    border_color = "grey90",
    breaks = seq(-max_val, max_val, length.out = 101),
    cluster_rows = nrow(mat) >= 2,
    cluster_cols = ncol(mat) >= 2,
    fontsize_row = 11, fontsize_col = 12, fontsize = 14
  )

  if (!is.null(output_file)) {
    grDevices::png(output_file, width = width, height = height,
                    units = "in", res = 300)
    print(hp)
    grDevices::dev.off()
    message("Saved enrichment heatmap to: ", output_file)
  }

  hp
}


#' Standalone QC Plot (from data frame)
#'
#' Generates a multi-panel QC summary figure from a metadata data frame
#' (does not require Seurat).
#'
#' @param metadata Data frame with columns: nCount_RNA, nFeature_RNA,
#'   and the grouping column.
#' @param group_by Grouping column. Default "experiment_group".
#' @param output_file Optional file path to save.
#' @param width Plot width. Default 16.
#' @param height Plot height. Default 6.
#' @return A ggplot or patchwork object.
#' @export
cdst_plot_qc_standalone <- function(
    metadata,
    group_by = "experiment_group",
    output_file = NULL,
    width = 16,
    height = 6
) {
  meta <- metadata

  p1 <- ggplot2::ggplot(meta, ggplot2::aes_string(
    x = group_by, y = "nCount_RNA", fill = group_by
  )) +
    ggplot2::geom_violin(scale = "width", alpha = 0.7) +
    ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Transcripts per Cell", y = "nCount_RNA") +
    ggplot2::theme(legend.position = "none",
                    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p2 <- ggplot2::ggplot(meta, ggplot2::aes_string(
    x = group_by, y = "nFeature_RNA", fill = group_by
  )) +
    ggplot2::geom_violin(scale = "width", alpha = 0.7) +
    ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Genes per Cell", y = "nFeature_RNA") +
    ggplot2::theme(legend.position = "none",
                    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p3 <- ggplot2::ggplot(meta, ggplot2::aes_string(
    x = "nCount_RNA", y = "nFeature_RNA", color = group_by
  )) +
    ggplot2::geom_point(size = 0.1, alpha = 0.3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Genes vs Transcripts",
                  x = "nCount_RNA", y = "nFeature_RNA")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- p1 + p2 + p3 + patchwork::plot_layout(ncol = 3)
  } else {
    combined <- p1
  }

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, combined, width = width, height = height, dpi = 300)
    message("Saved QC plot to: ", output_file)
    return(invisible(combined))
  }
  combined
}


#' Standalone Proportions Bar Plot (from data frame)
#'
#' Creates a grouped or stacked bar plot showing cell type proportions
#' across experimental conditions from a plain data frame.
#'
#' @param metadata Data frame with cell type and group columns.
#' @param celltype_col Column with cell types. Default "cell_type_annotation".
#' @param group_col Column with groups. Default "experiment_group".
#' @param plot_type "stacked" or "grouped". Default "stacked".
#' @param colors Optional named color vector.
#' @param output_file Optional file path to save.
#' @param width Plot width. Default 12.
#' @param height Plot height. Default 8.
#' @return A ggplot object.
#' @export
cdst_plot_proportions_standalone <- function(
    metadata,
    celltype_col = "cell_type_annotation",
    group_col = "experiment_group",
    plot_type = "stacked",
    colors = NULL,
    output_file = NULL,
    width = 12,
    height = 8
) {
  counts <- as.data.frame(table(
    group = metadata[[group_col]], celltype = metadata[[celltype_col]]
  ))
  colnames(counts) <- c("group", "celltype", "count")
  totals <- stats::aggregate(count ~ group, data = counts, FUN = sum)
  counts <- merge(counts, totals, by = "group", suffixes = c("", "_total"))
  counts$proportion <- counts$count / counts$count_total
  position <- if (plot_type == "stacked") "stack" else "dodge"

  p <- ggplot2::ggplot(counts,
                        ggplot2::aes(x = group, y = proportion, fill = celltype)) +
    ggplot2::geom_bar(stat = "identity", position = position, width = 0.7) +
    ggplot2::labs(title = "Cell Type Proportions by Group",
                  x = "", y = "Proportion", fill = "Cell Type") +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      axis.title = ggplot2::element_text(face = "bold", size = 15),
      axis.text = ggplot2::element_text(size = 13),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 13),
      legend.title = ggplot2::element_text(face = "bold", size = 14),
      legend.text = ggplot2::element_text(size = 12)
    )

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  } else {
    n_types <- length(unique(counts$celltype))
    p <- p + ggplot2::scale_fill_manual(
      values = grDevices::hcl.colors(n_types, palette = "Set2")
    )
  }

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height, dpi = 300)
    message("Saved proportion plot to: ", output_file)
    return(invisible(p))
  }
  p
}


# =============================================================================
# Distribution Heatmap (from cdst_distribution_matrix output)
# =============================================================================

#' Plot Cell Type Distribution Heatmap
#'
#' Creates a publication-quality heatmap from the output of
#' \code{\link{cdst_distribution_matrix}}, showing cell type proportions
#' across brain regions. Rows are cell types, columns are regions, and
#' cells are colored by normalized proportion. Optional column grouping
#' annotations show parent brain region hierarchy.
#'
#' This is the R-native heatmap for the distribution matrix. For the
#' multi-dimensional Python-based heatmap with circle encoding, see
#' \code{\link{cdst_plot_multidim_heatmap}}.
#'
#' @param dist_data A list from \code{\link{cdst_distribution_matrix}} with
#'   components: matrix, row_order, col_order, row_groups, col_groups.
#' @param title Plot title. Default: "Cell Type Distribution".
#' @param color_palette Character vector of colors for the heatmap gradient.
#'   Default: blue-white-red.
#' @param show_values Logical. Show numeric values in cells. Default: FALSE.
#' @param fontsize_row Row label font size. Default: 8.
#' @param fontsize_col Column label font size. Default: 7.
#' @param output_file Optional file path to save the plot.
#' @param width Plot width in inches. Default: 14.
#' @param height Plot height in inches. Default: 10.
#'
#' @return A pheatmap object (invisibly if output_file is specified).
#'
#' @examples
#' \dontrun{
#' dist <- cdst_distribution_matrix(seurat_obj)
#' cdst_plot_distribution_heatmap(dist,
#'   output_file = "distribution_heatmap.png")
#' }
#'
#' @export
cdst_plot_distribution_heatmap <- function(
    dist_data,
    title = "Cell Type Distribution",
    color_palette = NULL,
    show_values = FALSE,
    fontsize_row = 11,
    fontsize_col = 10,
    output_file = NULL,
    width = 14,
    height = 10
) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Install with: ",
         "install.packages('pheatmap')", call. = FALSE)
  }

  # Accept either a list from cdst_distribution_matrix or a raw matrix
  if (is.list(dist_data) && "matrix" %in% names(dist_data)) {
    mat <- dist_data$matrix
    col_groups <- dist_data$col_groups
  } else if (is.matrix(dist_data)) {
    mat <- dist_data
    col_groups <- NULL
  } else {
    stop("dist_data must be a list from cdst_distribution_matrix() or a numeric matrix.",
         call. = FALSE)
  }

  if (nrow(mat) == 0 || ncol(mat) == 0) {
    warning("Empty matrix — nothing to plot.")
    return(invisible(NULL))
  }

  # Color palette
  if (is.null(color_palette)) {
    color_palette <- grDevices::colorRampPalette(
      c("#2166ac", "#67a9cf", "#d1e5f0", "#fddbc7", "#ef8a62", "#b2182b")
    )(100)
  }

  # Build column annotation if parent groups are available
  annotation_col <- NULL
  annotation_colors <- NULL
  if (!is.null(col_groups) && length(col_groups) > 0) {
    valid_cols <- intersect(colnames(mat), names(col_groups))
    if (length(valid_cols) > 0) {
      annotation_col <- data.frame(
        Parent = col_groups[valid_cols],
        row.names = valid_cols,
        stringsAsFactors = FALSE
      )
      # Generate colors for parent regions
      unique_parents <- unique(annotation_col$Parent)
      parent_colors <- grDevices::hcl.colors(length(unique_parents), palette = "Set2")
      annotation_colors <- list(
        Parent = stats::setNames(parent_colors, unique_parents)
      )
    }
  }

  # Display numbers?
  display_numbers <- if (show_values) round(mat, 2) else FALSE

  hp <- pheatmap::pheatmap(
    mat,
    color = color_palette,
    main = title,
    border_color = "grey90",
    display_numbers = display_numbers,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    cluster_rows = nrow(mat) >= 2,
    cluster_cols = ncol(mat) >= 2,
    angle_col = 45
  )

  if (!is.null(output_file)) {
    grDevices::png(output_file, width = width, height = height,
                    units = "in", res = 300)
    print(hp)
    grDevices::dev.off()
    message("Saved distribution heatmap to: ", output_file)
    return(invisible(hp))
  }

  hp
}


#' DEG Scatter Plot (MA-style)
#'
#' Creates a scatter plot of average expression vs log2 fold change,
#' highlighting significant DEGs. This is an MA-plot variant commonly
#' used alongside volcano plots.
#'
#' @param deg_results Data frame with gene, avg_log2FC, p_val_adj columns.
#'   Optionally includes avg_expression or pct_1/pct_2 columns.
#' @param contrast Character. Which contrast to plot. Default NULL (all).
#' @param region Character. Which region to filter. Default NULL (all).
#' @param logfc_threshold LogFC threshold for significance. Default 0.5.
#' @param p_threshold Adjusted p-value threshold. Default 0.05.
#' @param n_label Number of top genes to label. Default 10.
#' @param title Plot title.
#' @param output_file Optional file path to save.
#' @param width Plot width. Default 10.
#' @param height Plot height. Default 8.
#' @return A ggplot object.
#'
#' @export
cdst_plot_deg_scatter <- function(
    deg_results,
    contrast = NULL,
    region = NULL,
    logfc_threshold = 0.5,
    p_threshold = 0.05,
    n_label = 10,
    title = NULL,
    output_file = NULL,
    width = 10,
    height = 8
) {
  df <- deg_results
  if (!is.null(contrast) && "contrast" %in% colnames(df)) {
    df <- df[df$contrast == contrast, ]
  }
  if (!is.null(region) && "region" %in% colnames(df)) {
    df <- df[df$region == region, ]
  }
  if (!"gene" %in% colnames(df)) {
    if (!is.null(rownames(df))) df$gene <- rownames(df)
    else stop("DEG results must have a 'gene' column or row names.", call. = FALSE)
  }

  # Compute average expression proxy if not present
  if ("avg_expression" %in% colnames(df)) {
    df$avg_expr <- df$avg_expression
  } else if (all(c("pct.1", "pct.2") %in% colnames(df))) {
    df$avg_expr <- (df$pct.1 + df$pct.2) / 2
  } else {
    df$avg_expr <- abs(df$avg_log2FC) + runif(nrow(df), 0, 0.1)
  }

  df$category <- "ns"
  df$category[df$avg_log2FC > logfc_threshold & df$p_val_adj < p_threshold] <- "up"
  df$category[df$avg_log2FC < -logfc_threshold & df$p_val_adj < p_threshold] <- "down"

  colors <- c(up = "#D62728", down = "#1F77B4", ns = "#CCCCCC")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = avg_expr, y = avg_log2FC)) +
    ggplot2::geom_point(ggplot2::aes(color = category), size = 1, alpha = 0.6) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::geom_hline(yintercept = c(-logfc_threshold, logfc_threshold),
                        linetype = "dashed", color = "grey50", linewidth = 0.4) +
    ggplot2::labs(
      title = title %||% "DEG Scatter (MA Plot)",
      x = "Average Expression",
      y = expression(log[2]~"Fold Change"),
      color = "Category"
    ) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      axis.title = ggplot2::element_text(face = "bold", size = 15),
      axis.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(face = "bold", size = 14),
      legend.text = ggplot2::element_text(size = 12),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Label top genes
  sig_df <- df[df$category != "ns", ]
  sig_df <- sig_df[order(-abs(sig_df$avg_log2FC)), ]
  label_genes <- utils::head(sig_df$gene, n_label)
  if (length(label_genes) > 0) {
    label_df <- df[df$gene %in% label_genes, ]
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_df, ggplot2::aes(label = gene),
        size = 3, max.overlaps = 15
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = label_df, ggplot2::aes(label = gene),
        size = 3, check_overlap = TRUE
      )
    }
  }

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height, dpi = 300)
    message("Saved DEG scatter plot to: ", output_file)
    return(invisible(p))
  }
  p
}


#' WGCNA Module Network Plot (R-native)
#'
#' Creates a network visualization of WGCNA module hub genes using igraph.
#' This is the pure-R alternative to \code{\link{cdst_plot_wgcna_network}}
#' (which uses Python). Nodes are hub genes, colored by module, and edges
#' represent co-expression relationships.
#'
#' @param hub_genes Data frame from \code{\link{cdst_wgcna_hub_genes}} with
#'   columns: module, gene, kME.
#' @param wgcna_result List from \code{\link{cdst_wgcna_detect}} (used to
#'   compute adjacency for edges).
#' @param top_n Number of hub genes per module to include. Default 5.
#' @param min_adjacency Minimum adjacency weight to draw an edge. Default 0.1.
#' @param title Plot title.
#' @param output_file Optional file path to save.
#' @param width Plot width. Default 10.
#' @param height Plot height. Default 10.
#' @return Invisible NULL (plot is drawn to device).
#'
#' @export
cdst_plot_network <- function(
    hub_genes,
    wgcna_result = NULL,
    top_n = 5,
    min_adjacency = 0.1,
    title = "Hub Gene Network",
    output_file = NULL,
    width = 10,
    height = 10
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. Install with: ",
         "install.packages('igraph')", call. = FALSE)
  }

  # Select top genes per module
  modules <- unique(hub_genes$module)
  selected <- do.call(rbind, lapply(modules, function(m) {
    mod_df <- hub_genes[hub_genes$module == m, ]
    mod_df <- mod_df[order(-mod_df$kME), ]
    utils::head(mod_df, top_n)
  }))

  genes <- selected$gene
  n_genes <- length(genes)

  # Build adjacency matrix from expression data if available
  if (!is.null(wgcna_result) && !is.null(wgcna_result$datExpr)) {
    expr <- wgcna_result$datExpr
    common_genes <- intersect(genes, colnames(expr))
    if (length(common_genes) >= 2) {
      cor_mat <- stats::cor(expr[, common_genes], use = "pairwise.complete.obs")
      adj_mat <- abs(cor_mat)
      diag(adj_mat) <- 0
      adj_mat[adj_mat < min_adjacency] <- 0
    } else {
      adj_mat <- matrix(0, n_genes, n_genes)
    }
  } else {
    # Create a simple module-based adjacency
    adj_mat <- matrix(0, n_genes, n_genes,
                      dimnames = list(genes, genes))
    for (m in modules) {
      mod_genes <- selected$gene[selected$module == m]
      for (i in seq_along(mod_genes)) {
        for (j in seq_along(mod_genes)) {
          if (i != j) adj_mat[mod_genes[i], mod_genes[j]] <- 0.5
        }
      }
    }
  }

  # Build igraph
  g <- igraph::graph_from_adjacency_matrix(
    adj_mat, mode = "undirected", weighted = TRUE, diag = FALSE
  )

  # Node attributes
  node_modules <- stats::setNames(selected$module, selected$gene)
  node_kme <- stats::setNames(selected$kME, selected$gene)

  igraph::V(g)$color <- node_modules[igraph::V(g)$name]
  igraph::V(g)$size <- scales::rescale(
    node_kme[igraph::V(g)$name], to = c(5, 15)
  )
  igraph::V(g)$label.cex <- 0.9

  # Edge attributes
  if (igraph::ecount(g) > 0) {
    igraph::E(g)$width <- scales::rescale(igraph::E(g)$weight, to = c(0.3, 2))
    igraph::E(g)$color <- grDevices::adjustcolor("grey60", alpha.f = 0.5)
  }

  layout <- igraph::layout_with_fr(g)

  if (!is.null(output_file)) {
    grDevices::png(output_file, width = width, height = height,
                    units = "in", res = 300)
  }

  graphics::plot(g, layout = layout, main = title,
       vertex.label = igraph::V(g)$name,
       vertex.frame.color = "grey30",
       edge.curved = 0.1)

  # Legend
  unique_modules <- unique(selected$module)
  graphics::legend("bottomleft", legend = unique_modules,
         fill = unique_modules, title = "Module",
         cex = 1.0, bty = "n")

  if (!is.null(output_file)) {
    grDevices::dev.off()
    message("Saved network plot to: ", output_file)
  }

  invisible(NULL)
}
