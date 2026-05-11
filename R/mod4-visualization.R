# =============================================================================
# Module 4: Visualization — Master Module
# =============================================================================
# This file provides R wrappers for all visualization functions, including
# those implemented in Python via reticulate. It also contains the core
# R-native plotting functions (volcano, spatial scatter, QC, proportions).
#
# Sub-modules (separate files):
#   mod4-viz-hierarchy.R     — Cell type hierarchy Sankey/dendrogram
#   mod4-viz-disproportion.R — Disproportion tSNE/UMAP
#   mod4-viz-enrichment.R    — Spatial GO enrichment
#   inst/python/celltype_distribution_heatmap.py — Multi-dimensional heatmap
#   inst/python/wgcna_network_plot.py — WGCNA network visualizations
# =============================================================================


# ---- Python Environment Setup ----

#' Set Up Python Environment for CellDynamicST Visualizations
#'
#' Ensures the required Python packages are available via reticulate.
#'
#' @param envname Name of the conda/virtualenv. Default "cdst".
#' @param method "conda" or "virtualenv". Default "conda".
#' @return Invisible NULL.
#' @export
cdst_setup_python <- function(envname = "cdst", method = "conda") {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' required. Install with: ",
         "install.packages('reticulate')", call. = FALSE)
  }
  required_pkgs <- c("numpy", "pandas", "matplotlib", "scipy", "networkx")
  if (method == "conda") {
    envs <- tryCatch(reticulate::conda_list(), error = function(e) data.frame())
    if (!envname %in% envs$name) {
      message("Creating conda environment '", envname, "'...")
      reticulate::conda_create(envname, packages = "python=3.11")
    }
    reticulate::use_condaenv(envname, required = FALSE)
    for (pkg in required_pkgs) {
      tryCatch(
        reticulate::conda_install(envname, packages = pkg, pip = TRUE),
        error = function(e) message("  Note: ", pkg, " install issue: ", e$message)
      )
    }
  } else {
    if (!reticulate::virtualenv_exists(envname)) {
      reticulate::virtualenv_create(envname)
    }
    reticulate::use_virtualenv(envname, required = FALSE)
    reticulate::virtualenv_install(envname, packages = required_pkgs)
  }
  message("Python environment '", envname, "' ready.")
  invisible(NULL)
}


# ---- R Wrapper for Python Heatmap ----

#' Plot Cell Type Distribution Heatmap (Multi-Dimensional)
#'
#' Creates the comprehensive multi-dimensional heatmap showing cell type
#' distribution across brain regions with circle encoding (size = count,
#' color = % change, border = significance, ring = variance, glyph = receptor).
#' Side annotations show stacked area charts of composition changes.
#'
#' Calls the Python implementation via reticulate.
#'
#' @param seurat_obj A Seurat object with annotated metadata, OR a path
#'   to a metadata CSV file.
#' @param genotype Genotype to analyze (e.g., "AA" or "GG").
#' @param control_suffix Control group suffix. Default "SAL".
#' @param treatment_suffix Treatment group suffix. Default "MOR".
#' @param output_dir Directory for output files. Default "heatmap_output".
#' @param nt_colors Named list of neurotransmitter colors (optional).
#' @param cell_group_colors Named list of cell group colors (optional).
#' @param meso_colors Named list of mesostructure colors (optional).
#' @return Invisible path to the generated heatmap PNG.
#' @export
cdst_plot_multidim_heatmap <- function(
    seurat_obj,
    genotype = "AA",
    control_suffix = "SAL",
    treatment_suffix = "MOR",
    output_dir = "heatmap_output",
    nt_colors = NULL,
    cell_group_colors = NULL,
    meso_colors = NULL
) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' required.", call. = FALSE)
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Export metadata to CSV if Seurat object provided
  if (inherits(seurat_obj, "Seurat")) {
    csv_path <- file.path(output_dir, "metadata_export.csv")
    utils::write.csv(seurat_obj@meta.data, csv_path)
    message("Exported metadata to: ", csv_path)
  } else if (is.character(seurat_obj) && file.exists(seurat_obj)) {
    csv_path <- seurat_obj
  } else {
    stop("Provide a Seurat object or path to metadata CSV.", call. = FALSE)
  }

  py_script <- system.file("python", "celltype_distribution_heatmap.py",
                            package = "CellDynamicST")
  if (py_script == "" || !file.exists(py_script)) {
    py_script <- file.path(find.package("CellDynamicST"),
                            "inst", "python", "celltype_distribution_heatmap.py")
  }
  if (!file.exists(py_script)) {
    stop("Python heatmap script not found.", call. = FALSE)
  }

  reticulate::source_python(py_script)
  meta_df <- reticulate::r_to_py(utils::read.csv(csv_path, row.names = 1))

  reticulate::py$process_and_visualize(
    meta_data = meta_df,
    output_dir = output_dir,
    genotype = genotype,
    control_suffix = control_suffix,
    treatment_suffix = treatment_suffix,
    nt_colors = nt_colors,
    cell_group_colors = cell_group_colors,
    meso_colors = meso_colors
  )

  out_file <- file.path(output_dir,
                        paste0(genotype, "_", treatment_suffix, "_vs_",
                               control_suffix, "_heatmap.png"))
  message("Heatmap saved to: ", out_file)
  invisible(out_file)
}


# ---- R Wrappers for WGCNA Network Plots ----

#' Plot WGCNA Network on Allen Brain Atlas
#'
#' Creates a network visualization of WGCNA modules overlaid on the Allen
#' Mouse Brain Atlas annotation volume. Users can specify the slice plane
#' (sagittal, coronal, or horizontal) and the exact slice coordinate.
#' Falls back to a simple silhouette if AllenSDK is not available.
#'
#' @param nodes_csv Path to CSV with columns: region, n_genes, module_color.
#' @param edges_csv Path to CSV with columns: source, target, correlation.
#' @param output_file Path for the output PNG.
#' @param title Plot title.
#' @param plane Slice plane: "sagittal", "coronal", or "horizontal".
#'   Default "sagittal".
#' @param slice Integer slice coordinate along the chosen axis.
#'   For sagittal: ML (z) coordinate. For coronal: AP (x) coordinate.
#'   For horizontal: DV (y) coordinate. Default NULL (midpoint).
#' @param resolution Allen annotation volume resolution in microns (10 or 25).
#'   Default 10.
#' @param use_simple Logical. Use simple silhouette fallback instead of
#'   AllenSDK. Default FALSE.
#' @param region_coords Named list of (x, y) coordinates per region (optional,
#'   only used with use_simple = TRUE).
#' @param min_edge_weight Minimum correlation to display. Default 0.3.
#' @return Invisible output file path.
#' @export
cdst_plot_wgcna_network <- function(
    nodes_csv, edges_csv,
    output_file = "wgcna_atlas_network.png",
    title = "WGCNA Module Network",
    plane = "sagittal",
    slice = NULL,
    resolution = 10,
    use_simple = FALSE,
    region_coords = NULL,
    min_edge_weight = 0.3
) {
  py_script <- system.file("python", "wgcna_network_plot.py",
                            package = "CellDynamicST")
  if (py_script == "" || !file.exists(py_script)) {
    py_script <- file.path(getwd(), "inst", "python",
                            "wgcna_network_plot.py")
  }
  if (!file.exists(py_script)) {
    stop("Cannot find wgcna_network_plot.py", call. = FALSE)
  }

  # Build CLI command
  cmd_parts <- c(
    "python3", shQuote(py_script), "atlas",
    "--nodes", shQuote(nodes_csv),
    "--edges", shQuote(edges_csv),
    "--output", shQuote(output_file),
    "--plane", plane,
    "--resolution", as.character(resolution)
  )
  if (!is.null(slice)) {
    cmd_parts <- c(cmd_parts, "--slice", as.character(as.integer(slice)))
  }
  if (!is.null(title)) {
    cmd_parts <- c(cmd_parts, "--title", shQuote(title))
  }
  if (use_simple) {
    cmd_parts <- c(cmd_parts, "--simple")
  }
  cmd <- paste(cmd_parts, collapse = " ")
  output <- system(paste(cmd, "2>&1"), intern = TRUE)
  for (line in output) message("  ", line)
  message("WGCNA atlas network saved to: ", output_file)
  invisible(output_file)
}


#' Plot WGCNA Module Enrichment Circos
#'
#' Creates a circos-style polar bar chart showing module enrichment across
#' multiple DEG comparisons.
#'
#' @param module_csv Path to CSV with module_color and NES_* columns.
#' @param output_file Path for the output PNG.
#' @param title Plot title.
#' @return Invisible output file path.
#' @export
cdst_plot_wgcna_circos <- function(
    module_csv,
    output_file = "wgcna_circos.png",
    title = "WGCNA Module Enrichment"
) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' required.", call. = FALSE)
  }
  py_script <- system.file("python", "wgcna_network_plot.py",
                            package = "CellDynamicST")
  if (py_script == "" || !file.exists(py_script)) {
    py_script <- file.path(find.package("CellDynamicST"),
                            "inst", "python", "wgcna_network_plot.py")
  }
  reticulate::source_python(py_script)
  module_df <- reticulate::r_to_py(utils::read.csv(module_csv))
  reticulate::py$plot_circos(
    module_df = module_df, output_path = output_file, title = title
  )
  message("WGCNA circos plot saved to: ", output_file)
  invisible(output_file)
}


#' Plot WGCNA Hub Gene Network
#'
#' Creates a force-directed network of hub genes colored by module.
#'
#' @param hub_csv Path to CSV with gene, module_color, kME columns.
#' @param edge_csv Path to CSV with gene1, gene2, weight columns.
#' @param output_file Path for the output PNG.
#' @param top_n Number of top hub genes to include. Default 50.
#' @param title Plot title.
#' @return Invisible output file path.
#' @export
cdst_plot_wgcna_hub_network <- function(
    hub_csv, edge_csv,
    output_file = "wgcna_hub_network.png",
    top_n = 50,
    title = "Hub Gene Network"
) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' required.", call. = FALSE)
  }
  py_script <- system.file("python", "wgcna_network_plot.py",
                            package = "CellDynamicST")
  if (py_script == "" || !file.exists(py_script)) {
    py_script <- file.path(find.package("CellDynamicST"),
                            "inst", "python", "wgcna_network_plot.py")
  }
  reticulate::source_python(py_script)
  hub_df <- reticulate::r_to_py(utils::read.csv(hub_csv))
  edge_df <- reticulate::r_to_py(utils::read.csv(edge_csv))
  reticulate::py$plot_hub_network(
    hub_df = hub_df, edge_df = edge_df,
    output_path = output_file, top_n = as.integer(top_n), title = title
  )
  message("WGCNA hub network saved to: ", output_file)
  invisible(output_file)
}


# ---- Core R-native Plotting Functions ----

#' Volcano Plot for Differential Expression
#'
#' Creates a publication-quality volcano plot from DEG results.
#'
#' @param deg_results Data frame with gene, avg_log2FC, p_val_adj columns.
#' @param contrast Character. Which contrast to plot. If NULL, uses all.
#' @param region Character. Which region to filter. If NULL, uses all.
#' @param logfc_threshold LogFC threshold. Default 0.5.
#' @param p_threshold Adjusted p-value threshold. Default 0.05.
#' @param n_label Number of top genes to label. Default 15.
#' @param title Plot title.
#' @param colors Named vector: "up", "down", "ns".
#' @param output_file Optional file path to save.
#' @param width Plot width. Default 10.
#' @param height Plot height. Default 8.
#' @return A ggplot object.
#' @export
cdst_plot_volcano <- function(
    deg_results,
    contrast = NULL,
    region = NULL,
    logfc_threshold = 0.5,
    p_threshold = 0.05,
    n_label = 15,
    title = NULL,
    colors = c(up = "#D62728", down = "#1F77B4", ns = "#CCCCCC"),
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

  df$neg_log10p <- -log10(pmax(df$p_val_adj, 1e-300))
  df$category <- "ns"
  df$category[df$avg_log2FC > logfc_threshold & df$p_val_adj < p_threshold] <- "up"
  df$category[df$avg_log2FC < -logfc_threshold & df$p_val_adj < p_threshold] <- "down"

  sig_df <- df[df$category != "ns", ]
  sig_df <- sig_df[order(-sig_df$neg_log10p), ]
  label_genes <- head(sig_df$gene, n_label)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = avg_log2FC, y = neg_log10p)) +
    ggplot2::geom_point(ggplot2::aes(color = category), size = 1.2, alpha = 0.7) +
    ggplot2::scale_color_manual(
      values = colors,
      labels = c(
        up = paste0("Up (", sum(df$category == "up"), ")"),
        down = paste0("Down (", sum(df$category == "down"), ")"),
        ns = "NS"
      )
    ) +
    ggplot2::geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
                        linetype = "dashed", color = "grey50", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(p_threshold),
                        linetype = "dashed", color = "grey50", linewidth = 0.4) +
    ggplot2::labs(
      title = title %||% "Volcano Plot",
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"Adjusted p-value"),
      color = "Category"
    ) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      axis.title = ggplot2::element_text(face = "bold", size = 15),
      axis.text = ggplot2::element_text(size = 13),
      legend.title = ggplot2::element_text(face = "bold", size = 14),
      legend.text = ggplot2::element_text(size = 13),
      legend.position = "right",
      panel.grid.minor = ggplot2::element_blank()
    )

  if (length(label_genes) > 0) {
    label_df <- df[df$gene %in% label_genes, ]
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = label_df, ggplot2::aes(label = gene),
        size = 3, max.overlaps = 20,
        segment.color = "grey60", segment.size = 0.3
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = label_df, ggplot2::aes(label = gene),
        size = 3, check_overlap = TRUE, nudge_y = 0.5
      )
    }
  }

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height, dpi = 300)
    message("Saved volcano plot to: ", output_file)
    return(invisible(p))
  }
  p
}


#' Spatial Feature Plot
#'
#' Plots gene expression or metadata values on spatial coordinates.
#'
#' @param seurat_obj A Seurat object with spatial coordinates in metadata.
#' @param feature Gene name or metadata column to plot.
#' @param x_col Column for x-axis. Default "DV_location".
#' @param y_col Column for y-axis. Default "AP_location".
#' @param point_size Point size. Default 0.8.
#' @param palette viridis palette name. Default "viridis".
#' @param title Plot title.
#' @param output_file Optional file path to save.
#' @param width Plot width. Default 10.
#' @param height Plot height. Default 8.
#' @return A ggplot object.
#' @export
cdst_plot_spatial_feature <- function(
    seurat_obj, feature,
    x_col = "DV_location", y_col = "AP_location",
    point_size = 0.8, palette = "viridis",
    title = NULL, output_file = NULL,
    width = 10, height = 8
) {
  meta <- seurat_obj@meta.data
  if (is.null(title)) title <- feature

  if (feature %in% colnames(meta)) {
    values <- meta[[feature]]
  } else if (feature %in% rownames(seurat_obj)) {
    values <- Seurat::GetAssayData(seurat_obj, layer = "data")[feature, ]
  } else {
    stop("Feature '", feature, "' not found.", call. = FALSE)
  }

  plot_df <- data.frame(x = meta[[x_col]], y = meta[[y_col]],
                        value = values, stringsAsFactors = FALSE)
  plot_df <- plot_df[order(plot_df$value), ]

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(color = value),
                        size = point_size, alpha = 0.8) +
    ggplot2::scale_color_viridis_c(option = palette, name = feature) +
    ggplot2::labs(title = title, x = "DV Location", y = "AP Location") +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 17),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::coord_fixed()

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height, dpi = 300)
    message("Saved spatial feature plot to: ", output_file)
    return(invisible(p))
  }
  p
}


#' Spatial Cell Type Map
#'
#' Generates a spatial scatter plot of cells colored by cell type.
#'
#' @param seurat_obj A Seurat object.
#' @param color_by Metadata column for coloring. Default "cell_type_annotation".
#' @param x_coord X coordinate column. Default "DV_location".
#' @param y_coord Y coordinate column. Default "AP_location".
#' @param facet_by Column for faceting. Default NULL.
#' @param point_size Point size. Default 0.3.
#' @param palette Named color vector. If NULL, auto-generated.
#' @param title Plot title.
#' @param show_legend Logical. Default TRUE.
#' @param output_file Optional file path to save.
#' @param width Plot width. Default 14.
#' @param height Plot height. Default 10.
#' @return A ggplot object.
#' @export
cdst_plot_spatial <- function(
    seurat_obj,
    color_by = "cell_type_annotation",
    x_coord = "DV_location",
    y_coord = "AP_location",
    facet_by = NULL,
    point_size = 0.3,
    palette = NULL,
    title = NULL,
    show_legend = TRUE,
    output_file = NULL,
    width = 14,
    height = 10
) {
  meta <- seurat_obj@meta.data
  p <- ggplot2::ggplot(meta, ggplot2::aes_string(
    x = x_coord, y = y_coord, color = color_by
  )) +
    ggplot2::geom_point(size = point_size, alpha = 0.7) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 15) +
    ggplot2::labs(
      title = title %||% paste("Spatial Map:", color_by),
      x = "DV", y = "AP", color = color_by
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 17),
      axis.title = ggplot2::element_text(face = "bold", size = 14),
      axis.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(face = "bold", size = 13),
      legend.text = ggplot2::element_text(size = 11),
      panel.grid = ggplot2::element_blank()
    )

  if (!is.null(palette)) {
    p <- p + ggplot2::scale_color_manual(values = palette)
  } else {
    n_colors <- length(unique(meta[[color_by]]))
    if (n_colors <= 30) {
      p <- p + ggplot2::scale_color_manual(
        values = grDevices::hcl.colors(n_colors, palette = "Set2")
      )
    }
  }

  if (!is.null(facet_by) && facet_by %in% colnames(meta)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)))
  }
  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height, dpi = 300)
    message("Saved spatial map to: ", output_file)
    return(invisible(p))
  }
  p
}


#' Cell Type Proportion Bar Plot
#'
#' Creates a grouped or stacked bar plot showing cell type proportions
#' across experimental conditions.
#'
#' @param seurat_obj A Seurat object.
#' @param celltype_col Metadata column with cell types.
#' @param group_col Metadata column with groups.
#' @param plot_type "stacked" or "grouped". Default "stacked".
#' @param colors Optional named color vector.
#' @param output_file Optional file path to save.
#' @param width Plot width. Default 12.
#' @param height Plot height. Default 8.
#' @return A ggplot object.
#' @export
cdst_plot_proportions <- function(
    seurat_obj,
    celltype_col = "cell_type_annotation",
    group_col = "experiment_group",
    plot_type = "stacked",
    colors = NULL,
    output_file = NULL,
    width = 12,
    height = 8
) {
  meta <- seurat_obj@meta.data
  counts <- as.data.frame(table(
    group = meta[[group_col]], celltype = meta[[celltype_col]]
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


#' WGCNA Module-Trait Correlation Heatmap
#'
#' Creates a heatmap of module-trait correlations with significance stars.
#'
#' @param trait_cor A list from cdst_wgcna_trait_cor with cor_* and p_* matrices.
#' @param trait_type "treatment", "genotype", or "region". Default "treatment".
#' @param sig_threshold P-value threshold for stars. Default 0.05.
#' @param title Plot title.
#' @param output_file Optional file path to save.
#' @param width Plot width. Default 10.
#' @param height Plot height. Default 8.
#' @return A pheatmap object.
#' @export
cdst_plot_trait_heatmap <- function(
    trait_cor,
    trait_type = "treatment",
    sig_threshold = 0.05,
    title = NULL,
    output_file = NULL,
    width = 10,
    height = 8
) {
  cor_name <- paste0("cor_", trait_type)
  p_name <- paste0("p_", trait_type)
  if (!cor_name %in% names(trait_cor)) {
    stop("Trait type '", trait_type, "' not found.", call. = FALSE)
  }
  cor_mat <- trait_cor[[cor_name]]
  p_mat <- trait_cor[[p_name]]

  sig_text <- matrix("", nrow = nrow(cor_mat), ncol = ncol(cor_mat))
  sig_text[p_mat < sig_threshold] <- "*"
  sig_text[p_mat < 0.01] <- "**"
  sig_text[p_mat < 0.001] <- "***"
  display_text <- matrix(
    paste0(round(cor_mat, 2), "\n", sig_text), nrow = nrow(cor_mat)
  )

  colors <- grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(100)
  if (is.null(title)) {
    title <- paste("Module-Trait Correlation:", tools::toTitleCase(trait_type))
  }

  hp <- pheatmap::pheatmap(
    cor_mat, color = colors, display_numbers = display_text,
    number_color = "black", fontsize_number = 9, main = title,
    border_color = "grey90", breaks = seq(-1, 1, length.out = 101),
    cluster_rows = nrow(cor_mat) >= 2,
    cluster_cols = ncol(cor_mat) >= 2
  )

  if (!is.null(output_file)) {
    grDevices::png(output_file, width = width, height = height,
                    units = "in", res = 300)
    print(hp)
    grDevices::dev.off()
    message("Saved trait heatmap to: ", output_file)
  }

  hp
}


#' Plot QC Summary
#'
#' Generates a multi-panel QC summary figure.
#'
#' @param seurat_obj A Seurat object (post-QC).
#' @param group_by Grouping column. Default "experiment_group".
#' @param output_file Optional file path to save.
#' @param width Plot width. Default 16.
#' @param height Plot height. Default 6.
#' @return A ggplot or patchwork object.
#' @export
cdst_plot_qc <- function(
    seurat_obj,
    group_by = "experiment_group",
    output_file = NULL,
    width = 16,
    height = 6
) {
  meta <- seurat_obj@meta.data

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
