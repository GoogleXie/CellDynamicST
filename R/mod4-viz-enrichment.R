# =============================================================================
# Module 4d: Spatial GO Enrichment Visualization
# =============================================================================
# Refactored from: fixed_gsea_script_ucell.R + plot_spatial_go_enrichment_on_atlas.py
# Produces: Spatial scatter plots of GO enrichment scores overlaid on
#   brain coordinates, with bidirectional coloring (enriched/depleted)
# =============================================================================

#' Get Genes for a GO Term
#'
#' Retrieves gene symbols associated with a Gene Ontology term using
#' the \pkg{org.Mm.eg.db} (mouse) or \pkg{org.Hs.eg.db} (human) database.
#'
#' @param go_id Character string of the GO ID (e.g., "GO:0007268").
#' @param organism Character, either "mouse" or "human". Default "mouse".
#' @return Character vector of gene symbols.
#' @export
cdst_get_go_genes <- function(go_id, organism = "mouse") {
  if (organism == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      stop("Package 'org.Mm.eg.db' required. Install from Bioconductor.",
           call. = FALSE)
    }
    db <- org.Mm.eg.db::org.Mm.eg.db
  } else {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Package 'org.Hs.eg.db' required. Install from Bioconductor.",
           call. = FALSE)
    }
    db <- org.Hs.eg.db::org.Hs.eg.db
  }

  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' required.", call. = FALSE)
  }

  # Get Entrez IDs for the GO term
  entrez_ids <- tryCatch({
    AnnotationDbi::select(db, keys = go_id, keytype = "GOALL",
                          columns = "ENTREZID")$ENTREZID
  }, error = function(e) {
    warning("Could not retrieve genes for ", go_id, ": ", e$message)
    return(character(0))
  })

  entrez_ids <- unique(entrez_ids[!is.na(entrez_ids)])
  if (length(entrez_ids) == 0) return(character(0))

  # Convert to gene symbols
  symbols <- AnnotationDbi::select(db, keys = entrez_ids,
                                    keytype = "ENTREZID",
                                    columns = "SYMBOL")$SYMBOL
  unique(symbols[!is.na(symbols)])
}


#' Calculate UCell GO Enrichment Scores
#'
#' Computes per-cell GO pathway enrichment scores using UCell, then
#' z-scores relative to the control group distribution and identifies
#' significantly enriched/depleted cells.
#'
#' @param seurat_obj A Seurat object (subset to ROI and groups of interest).
#' @param go_genes Character vector of gene symbols for the GO term.
#' @param go_term Human-readable GO term name (for labeling).
#' @param control_group Name of the control group (used as null distribution).
#' @param treatment_group Name of the treatment group.
#' @param group_col Metadata column with group labels. Default "experiment_group".
#' @param upper_quantile Upper quantile threshold for significance. Default 0.95.
#' @param lower_quantile Lower quantile threshold for significance. Default 0.05.
#' @return A data frame with columns: cell_ID, AP_location, DV_location,
#'   experiment_group, enrichment_score, zscore, significant_enrichment.
#' @export
cdst_ucell_enrichment <- function(
    seurat_obj,
    go_genes,
    go_term = "GO_term",
    control_group,
    treatment_group,
    group_col = "experiment_group",
    upper_quantile = 0.95,
    lower_quantile = 0.05
) {
  if (!requireNamespace("UCell", quietly = TRUE)) {
    stop("Package 'UCell' required. Install from Bioconductor: ",
         "BiocManager::install('UCell')", call. = FALSE)
  }

  meta <- seurat_obj@meta.data
  available_genes <- intersect(go_genes, rownames(seurat_obj))
  if (length(available_genes) < 3) {
    warning("Fewer than 3 GO genes found in expression matrix for '",
            go_term, "'.")
    return(data.frame())
  }

  # Score with UCell
  gene_list <- list(go_genes)
  names(gene_list) <- gsub("[^A-Za-z0-9]", "_", go_term)
  seurat_obj <- UCell::AddModuleScore_UCell(seurat_obj, features = gene_list)

  score_col <- paste0(names(gene_list), "_UCell")
  if (!score_col %in% colnames(seurat_obj@meta.data)) {
    warning("UCell score column not found.")
    return(data.frame())
  }

  # Build result data frame
  df <- data.frame(
    cell_ID = colnames(seurat_obj),
    enrichment_score = seurat_obj@meta.data[[score_col]],
    experiment_group = meta[[group_col]],
    stringsAsFactors = FALSE
  )

  # Add spatial coordinates if available
  for (coord in c("AP_location", "DV_location")) {
    if (coord %in% colnames(meta)) {
      df[[coord]] <- meta[[coord]]
    }
  }

  # Z-score relative to control
  ctrl_scores <- df$enrichment_score[df$experiment_group == control_group]
  ctrl_scores <- ctrl_scores[!is.na(ctrl_scores)]

  if (length(ctrl_scores) == 0) {
    warning("No control cells found for group '", control_group, "'.")
    return(data.frame())
  }

  ctrl_mean <- mean(ctrl_scores, na.rm = TRUE)
  ctrl_sd <- stats::sd(ctrl_scores, na.rm = TRUE)

  if (!is.na(ctrl_sd) && ctrl_sd > 0) {
    df$zscore <- (df$enrichment_score - ctrl_mean) / ctrl_sd
  } else {
    df$zscore <- df$enrichment_score - ctrl_mean
  }

  # Significance thresholding
  upper_thresh <- stats::quantile(ctrl_scores, upper_quantile, na.rm = TRUE)
  lower_thresh <- stats::quantile(ctrl_scores, lower_quantile, na.rm = TRUE)

  df$significant_enrichment <- ifelse(
    df$enrichment_score > upper_thresh, df$zscore,
    ifelse(df$enrichment_score < lower_thresh, df$zscore, 0)
  )

  message("  Enriched: ", sum(df$significant_enrichment > 0, na.rm = TRUE),
          " | Depleted: ", sum(df$significant_enrichment < 0, na.rm = TRUE),
          " | Non-sig: ", sum(df$significant_enrichment == 0, na.rm = TRUE))

  df
}


#' Plot Spatial Enrichment Scores
#'
#' Creates a spatial scatter plot of GO enrichment scores on brain coordinates
#' with diverging colormap (enriched = warm, depleted = cold).
#'
#' @param enrichment_df Data frame from \code{\link{cdst_ucell_enrichment}}.
#' @param experiment_group Which group to plot (e.g., "AA MOR").
#' @param go_term GO term name for the title.
#' @param roi Brain region name for the title.
#' @param comparison Comparison label for the title.
#' @param score_col Column to use for coloring. Default "significant_enrichment".
#' @param x_col Column for x-axis. Default "DV_location".
#' @param y_col Column for y-axis. Default "AP_location".
#' @param low_color Color for depleted. Default "darkblue".
#' @param mid_color Color for non-significant. Default "darkgrey".
#' @param high_color Color for enriched. Default "darkorange".
#' @param point_size_range Numeric vector of length 2 for size range.
#'   Default c(0.6, 3).
#' @param output_file Optional file path to save the plot.
#' @param width Plot width in inches. Default 10.
#' @param height Plot height in inches. Default 8.
#' @return A ggplot object.
#' @export
cdst_plot_spatial_enrichment <- function(
    enrichment_df,
    experiment_group,
    go_term = "",
    roi = "",
    comparison = "",
    score_col = "significant_enrichment",
    x_col = "DV_location",
    y_col = "AP_location",
    low_color = "darkblue",
    mid_color = "darkgrey",
    high_color = "darkorange",
    point_size_range = c(0.6, 3),
    output_file = NULL,
    width = 10,
    height = 8
) {
  df <- enrichment_df[enrichment_df$experiment_group == experiment_group, ]
  if (nrow(df) == 0) {
    warning("No cells found for group '", experiment_group, "'.")
    return(ggplot2::ggplot())
  }

  # Calculate size and alpha from absolute score
  df$abs_score <- abs(df[[score_col]])
  df$plot_alpha <- ifelse(df$abs_score > 0,
                          pmin(sqrt(df$abs_score) / max(sqrt(df$abs_score + 1e-10)), 1),
                          0.3)

  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = x_col, y = y_col)) +
    ggplot2::geom_point(
      ggplot2::aes_string(
        color = score_col,
        size = "abs_score",
        alpha = "plot_alpha"
      )
    ) +
    ggplot2::scale_color_gradient2(
      low = low_color, mid = mid_color, high = high_color,
      midpoint = 0,
      name = "Enrichment\n(z-score)"
    ) +
    ggplot2::scale_size_continuous(
      range = point_size_range,
      trans = "sqrt",
      guide = "none"
    ) +
    ggplot2::scale_alpha_continuous(range = c(0.3, 1), guide = "none") +
    ggplot2::labs(
      title = paste0(go_term, " — ", experiment_group),
      subtitle = paste0("Region: ", roi, " | Comparison: ", comparison),
      x = "DV Location",
      y = "AP Location"
    ) +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 17),
      plot.subtitle = ggplot2::element_text(color = "grey40", size = 13),
      axis.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(face = "bold", size = 13),
      legend.text = ggplot2::element_text(size = 11),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::coord_fixed()

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height, dpi = 300)
    message("Saved spatial enrichment plot to: ", output_file)
    return(invisible(p))
  }

  p
}


#' Run Full Spatial Enrichment Pipeline
#'
#' End-to-end pipeline: retrieves GO genes, calculates UCell enrichment,
#' and generates spatial plots for specified ROIs and comparisons.
#'
#' @param seurat_obj A Seurat object with spatial coordinates in metadata.
#' @param go_terms Named character vector where names are GO IDs and values
#'   are human-readable term names. Example:
#'   \code{c("GO:0007268" = "chemical synaptic transmission")}.
#' @param regions Character vector of brain regions (level_3_parent_roi_name).
#' @param comparisons List of length-2 character vectors (control, treatment).
#' @param group_col Metadata column with group labels. Default "experiment_group".
#' @param region_col Metadata column with region labels.
#'   Default "level_3_parent_roi_name".
#' @param organism "mouse" or "human". Default "mouse".
#' @param output_dir Directory to save plots and CSV exports. Default "enrichment_output".
#' @return A list of enrichment data frames, keyed by "region_comparison_term".
#' @export
cdst_run_spatial_enrichment <- function(
    seurat_obj,
    go_terms,
    regions,
    comparisons,
    group_col = "experiment_group",
    region_col = "level_3_parent_roi_name",
    organism = "mouse",
    output_dir = "enrichment_output"
) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  csv_dir <- file.path(output_dir, "atlas_data")
  if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)

  meta <- seurat_obj@meta.data
  results <- list()

  for (go_id in names(go_terms)) {
    go_term <- go_terms[[go_id]]
    message("\n=== GO Term: ", go_term, " (", go_id, ") ===")

    go_genes <- cdst_get_go_genes(go_id, organism = organism)
    if (length(go_genes) == 0) {
      warning("No genes found for ", go_id, ". Skipping.")
      next
    }
    message("  Found ", length(go_genes), " genes")

    for (roi in regions) {
      for (comp in comparisons) {
        ctrl <- comp[1]
        treat <- comp[2]
        comp_label <- paste0(gsub(" ", "_", ctrl), "_vs_", gsub(" ", "_", treat))
        key <- paste(roi, comp_label, go_term, sep = "_")

        message("  Processing: ", roi, " | ", comp_label)

        # Subset to ROI and groups
        cells_use <- rownames(meta)[
          meta[[region_col]] == roi &
          meta[[group_col]] %in% c(ctrl, treat)
        ]

        if (length(cells_use) < 20) {
          message("    Skipping: only ", length(cells_use), " cells")
          next
        }

        sub_obj <- subset(seurat_obj, cells = cells_use)

        # Calculate enrichment
        enrich_df <- tryCatch(
          cdst_ucell_enrichment(
            sub_obj, go_genes, go_term,
            control_group = ctrl, treatment_group = treat,
            group_col = group_col
          ),
          error = function(e) {
            warning("Error: ", e$message)
            return(data.frame())
          }
        )

        if (nrow(enrich_df) == 0) next

        # Plot for treatment group
        safe_name <- gsub("[^A-Za-z0-9]", "_", key)
        out_file <- file.path(output_dir,
                              paste0(safe_name, "_", treat, ".png"))
        cdst_plot_spatial_enrichment(
          enrich_df,
          experiment_group = treat,
          go_term = go_term,
          roi = roi,
          comparison = comp_label,
          output_file = out_file
        )

        # Export CSV for Python atlas overlay
        export_df <- enrich_df[enrich_df$significant_enrichment != 0, ]
        if (nrow(export_df) > 0) {
          export_df$roi <- roi
          export_df$go_term <- go_term
          export_df$comparison <- comp_label
          csv_file <- file.path(csv_dir, paste0(safe_name, ".csv"))
          utils::write.csv(export_df, csv_file, row.names = FALSE)
          message("    Exported CSV: ", csv_file)
        }

        results[[key]] <- enrich_df
      }
    }
  }

  message("\n=== Spatial enrichment complete ===")
  message("Plots saved to: ", output_dir)
  message("Atlas CSV files saved to: ", csv_dir)

  results
}
