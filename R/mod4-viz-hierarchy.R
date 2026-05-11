# =============================================================================
# CellDynamicST — Module 4a: Cell Type Hierarchy Visualization
# Faithfully generalized from: plot_cell_type_hierachy.R
# Produces: (1) interactive Sankey diagram  (2) ggalluvial static Sankey
#           (3) static dendrogram (legacy)
# =============================================================================

#' Build Cell Type Hierarchy Data
#'
#' Extracts hierarchical cell type relationships from metadata and prepares
#' node/link data structures for Sankey or dendrogram visualization.
#' Faithfully follows the logic in plot_cell_type_hierachy.R.
#'
#' @param data A Seurat object or a data.frame containing hierarchy columns.
#' @param levels Character vector of metadata column names ordered root -> leaf.
#' @param min_cells Minimum cell count for a link to be included (default 1).
#' @return A list with \code{nodes}, \code{links}, \code{hierarchy_df},
#'   and \code{root_mapping} (named vector: name -> root group).
#' @export
cdst_build_hierarchy <- function(
    data,
    levels = c("highest_level_cell_type",
               "high_level_cell_type",
               "combined_subtype",
               "cell_type_annotation"),
    min_cells = 1
) {
  if (inherits(data, "Seurat")) {
    meta <- data@meta.data
  } else if (is.data.frame(data)) {
    meta <- data
  } else {
    stop("`data` must be a Seurat object or a data frame.", call. = FALSE)
  }

  missing_cols <- setdiff(levels, colnames(meta))
  if (length(missing_cols) > 0)
    stop("Missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

  # ---- Unique hierarchy paths with cell counts ----
  hierarchy_df <- meta[, levels, drop = FALSE]
  hierarchy_df$cell_count <- 1L
  hierarchy_df <- stats::aggregate(cell_count ~ ., data = hierarchy_df, FUN = sum)
  hierarchy_df <- hierarchy_df[hierarchy_df$cell_count >= min_cells, , drop = FALSE]

  # ---- Build unique node list (all categories across all levels) ----
  unique_categories <- unique(unlist(lapply(levels, function(l) unique(as.character(hierarchy_df[[l]])))))
  nodes_df <- data.frame(name = unique_categories, stringsAsFactors = FALSE)
  category_to_id <- stats::setNames(seq_along(unique_categories) - 1L, unique_categories)

  # ---- Build links cell-by-cell then aggregate (matching original script) ----
  links_raw <- list()
  for (j in seq_len(length(levels) - 1)) {
    pairs <- hierarchy_df[, c(levels[j], levels[j + 1], "cell_count")]
    colnames(pairs) <- c("source_name", "target_name", "value")
    pairs$source <- category_to_id[as.character(pairs$source_name)]
    pairs$target <- category_to_id[as.character(pairs$target_name)]
    links_raw[[j]] <- pairs
  }
  links_df <- do.call(rbind, links_raw)
  links_df <- stats::aggregate(value ~ source + target + source_name + target_name,
                                data = links_df, FUN = sum)

  # ---- Root mapping: every name -> its root-level group ----
  root_col <- levels[1]
  root_mapping <- character()
  for (lv in levels) {
    pairs <- unique(hierarchy_df[, c(lv, root_col)])
    vec <- stats::setNames(as.character(pairs[[root_col]]), as.character(pairs[[lv]]))
    root_mapping <- c(root_mapping, vec)
  }
  root_mapping <- root_mapping[!duplicated(names(root_mapping))]

  nodes_df$group <- ifelse(nodes_df$name %in% names(root_mapping),
                           root_mapping[nodes_df$name], "Other")

  # Assign level index to each node
  node_level <- rep(NA_integer_, nrow(nodes_df))
  for (i in seq_along(levels)) {
    vals <- unique(as.character(hierarchy_df[[levels[i]]]))
    node_level[nodes_df$name %in% vals & is.na(node_level)] <- i
  }
  nodes_df$level <- node_level

  list(
    nodes = nodes_df,
    links = links_df,
    hierarchy_df = hierarchy_df,
    root_mapping = root_mapping,
    levels = levels
  )
}


#' Plot Cell Type Hierarchy as Interactive Sankey Diagram
#'
#' Creates a publication-quality interactive Sankey diagram using networkD3,
#' faithfully reproducing the visualization from plot_cell_type_hierachy.R.
#' Nodes are colored by their root-level group; link widths reflect cell counts.
#'
#' @param hierarchy_data Output from \code{cdst_build_hierarchy}.
#' @param color_map Optional named vector: group name -> hex color.
#' @param font_size Font size for node labels (default 12).
#' @param node_width Width of Sankey nodes in pixels (default 10).
#' @param node_padding Vertical padding between nodes (default 2).
#' @param width Widget width in pixels (default 1400).
#' @param height Widget height in pixels (default 900).
#' @param output_file Path to save HTML widget (NULL = return widget).
#' @return A networkD3 htmlwidget.
#' @export
cdst_plot_sankey <- function(
    hierarchy_data,
    color_map = NULL,
    font_size = 12,
    node_width = 10,
    node_padding = 2,
    width = 1400,
    height = 900,
    output_file = NULL
) {
  if (!requireNamespace("networkD3", quietly = TRUE))
    stop("Package 'networkD3' required. Install with install.packages('networkD3').")
  if (!requireNamespace("htmlwidgets", quietly = TRUE))
    stop("Package 'htmlwidgets' required.")

  nodes <- hierarchy_data$nodes
  links <- hierarchy_data$links

  # ---- Build color scale ----
  root_groups <- unique(nodes$group)
  if (is.null(color_map)) {
    pal <- if (length(root_groups) <= 3) {
      c("#1f77b4", "#ff7f0e", "#2ca02c")[seq_along(root_groups)]
    } else if (length(root_groups) <= 8) {
      RColorBrewer::brewer.pal(max(3, length(root_groups)), "Set2")[seq_along(root_groups)]
    } else {
      grDevices::hcl.colors(length(root_groups), palette = "Dynamic")
    }
    color_map <- stats::setNames(pal, root_groups)
  }

  # Build D3 colour scale JS
  domain_str <- paste0("'", paste(names(color_map), collapse = "','"), "'")
  range_str  <- paste0("'", paste(color_map, collapse = "','"), "'")
  colour_js  <- htmlwidgets::JS(
    sprintf("d3.scaleOrdinal().domain([%s]).range([%s])", domain_str, range_str)
  )

  sankey <- networkD3::sankeyNetwork(
    Links       = links,
    Nodes       = nodes,
    Source      = "source",
    Target      = "target",
    Value       = "value",
    NodeID      = "name",
    NodeGroup   = "group",
    colourScale = colour_js,
    sinksRight  = FALSE,
    nodeWidth   = node_width,
    nodePadding = node_padding,
    fontSize    = font_size,
    fontFamily  = "Arial",
    width       = width,
    height      = height
  )

  if (!is.null(output_file)) {
    htmlwidgets::saveWidget(sankey, output_file, selfcontained = FALSE)
    message("Sankey saved to: ", output_file)
    return(invisible(sankey))
  }
  sankey
}


#' Plot Cell Type Hierarchy as Static Alluvial/Sankey Diagram
#'
#' Creates a publication-quality static Sankey-style alluvial plot using
#' ggalluvial, closely matching the manuscript Fig 1A style. Each level
#' is a stratum axis; flows are proportional to cell counts; colors are
#' assigned by root-level group with lighter shades for sub-groups.
#'
#' @param hierarchy_data Output from \code{cdst_build_hierarchy}, OR a
#'   data.frame/Seurat object (in which case cdst_build_hierarchy is called).
#' @param levels Character vector of column names (only used if hierarchy_data
#'   is a data.frame/Seurat).
#' @param color_map Optional named vector: root group -> hex color.
#' @param level_labels Optional character vector of display labels for each
#'   level axis (same length as levels). Default uses the column names.
#' @param label_size Font size for stratum labels. Default 3.
#' @param title Plot title. Default "Cell Type Hierarchy".
#' @param width Plot width in inches (default 24).
#' @param height Plot height in inches (default 16).
#' @param output_file Path to save PNG (NULL = return ggplot).
#' @return A ggplot object.
#' @export
cdst_plot_dendrogram <- function(
    hierarchy_data,
    levels = NULL,
    color_map = NULL,
    level_labels = NULL,
    label_size = 3,
    title = "Cell Type Hierarchy",
    width = 24,
    height = 16,
    output_file = NULL
) {
  if (!requireNamespace("ggalluvial", quietly = TRUE))
    stop("Package 'ggalluvial' required. Install with install.packages('ggalluvial').")

  # Allow passing raw data directly
  if (inherits(hierarchy_data, "Seurat") || is.data.frame(hierarchy_data)) {
    if (is.null(levels))
      levels <- c("highest_level_cell_type", "high_level_cell_type",
                   "combined_subtype", "cell_type_annotation")
    hierarchy_data <- cdst_build_hierarchy(hierarchy_data, levels = levels)
  }

  hier_df <- hierarchy_data$hierarchy_df
  lvls <- hierarchy_data$levels
  root_col <- lvls[1]
  root_mapping <- hierarchy_data$root_mapping

  # ---- Prepare alluvial data (long format) ----
  # Each row = one unique hierarchy path with cell_count as frequency
  alluvial_df <- hier_df

  # ---- Color palette by root group ----
  root_groups <- unique(as.character(alluvial_df[[root_col]]))
  if (is.null(color_map)) {
    # Default palette matching the reference figure style
    default_colors <- c(
      "#4A90D9", "#E8734A", "#5CB85C", "#9B59B6",
      "#F5A623", "#1ABC9C", "#E74C3C", "#3498DB"
    )
    color_map <- stats::setNames(
      default_colors[seq_along(root_groups)], root_groups
    )
  }

  # Assign colors to each unique leaf-level category based on root group
  # with shade variations for visual distinction
  leaf_col <- lvls[length(lvls)]
  leaf_roots <- stats::setNames(
    as.character(alluvial_df[[root_col]]),
    as.character(alluvial_df[[leaf_col]])
  )
  leaf_roots <- leaf_roots[!duplicated(names(leaf_roots))]

  # Generate shade variations within each root group
  fill_colors <- character()
  for (grp in root_groups) {
    base_col <- color_map[grp]
    members <- names(leaf_roots)[leaf_roots == grp]
    n <- length(members)
    if (n == 0) next
    if (n == 1) {
      fill_colors[members] <- base_col
    } else {
      # Create lighter-to-darker shades
      base_rgb <- grDevices::col2rgb(base_col) / 255
      shades <- vapply(seq(0.3, 0.85, length.out = n), function(f) {
        r <- base_rgb[1, 1] * f + (1 - f)
        g <- base_rgb[2, 1] * f + (1 - f)
        b <- base_rgb[3, 1] * f + (1 - f)
        grDevices::rgb(min(r, 1), min(g, 1), min(b, 1))
      }, character(1))
      fill_colors[members] <- shades
    }
  }

  # ---- Build the alluvial plot ----
  # Rename columns to axis1, axis2, ... for ggalluvial
  plot_df <- alluvial_df
  axis_names <- paste0("axis", seq_along(lvls))
  for (i in seq_along(lvls)) {
    plot_df[[axis_names[i]]] <- factor(as.character(plot_df[[lvls[i]]]))
  }
  plot_df$Freq <- plot_df$cell_count

  # The fill variable is the leaf level (finest granularity)
  plot_df$fill_var <- as.character(plot_df[[leaf_col]])

  # Level labels for x-axis
  if (is.null(level_labels)) {
    level_labels <- gsub("_", " ", lvls)
    level_labels <- tools::toTitleCase(level_labels)
  }

  # Build the ggplot
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes_string(
      y = "Freq",
      axis1 = axis_names[1],
      axis2 = axis_names[2],
      axis3 = if (length(lvls) >= 3) axis_names[3] else NULL,
      axis4 = if (length(lvls) >= 4) axis_names[4] else NULL,
      axis5 = if (length(lvls) >= 5) axis_names[5] else NULL
    )
  ) +
    ggalluvial::geom_alluvium(
      ggplot2::aes_string(fill = "fill_var"),
      width = 1/5, alpha = 0.65, decreasing = FALSE
    ) +
    ggalluvial::geom_stratum(
      width = 1/5, fill = "grey95", color = "grey60",
      linewidth = 0.3, decreasing = FALSE
    ) +
    ggalluvial::geom_stratum(
      ggplot2::aes_string(fill = "fill_var"),
      width = 1/5, alpha = 0.85, color = "grey40",
      linewidth = 0.3, decreasing = FALSE
    ) +
    ggplot2::geom_text(
      stat = ggalluvial::StatStratum,
      ggplot2::aes(label = ggplot2::after_stat(stratum)),
      size = label_size, fontface = "bold",
      decreasing = FALSE
    ) +
    ggplot2::scale_x_discrete(
      limits = axis_names,
      labels = level_labels,
      expand = c(0.15, 0.05)
    ) +
    ggplot2::scale_fill_manual(
      values = fill_colors, guide = "none"
    ) +
    ggplot2::labs(
      title = title,
      y = "Number of Cells"
    ) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 20, hjust = 0.5
      ),
      axis.text.x = ggplot2::element_text(
        face = "bold", size = 16, color = "black"
      ),
      axis.text.y = ggplot2::element_text(size = 14),
      axis.title.y = ggplot2::element_text(face = "bold", size = 16),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(12, 24, 12, 24)
    )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, p, width = width, height = height, dpi = 300)
    message("Hierarchy Sankey saved to: ", output_file)
    return(invisible(p))
  }
  p
}


# Legacy aliases for backward compatibility
#' @rdname cdst_plot_sankey
#' @export
cdst_plot_hierarchy_sankey <- function(...) cdst_plot_sankey(...)

#' @rdname cdst_plot_dendrogram
#' @export
cdst_plot_hierarchy_tree <- function(...) cdst_plot_dendrogram(...)
