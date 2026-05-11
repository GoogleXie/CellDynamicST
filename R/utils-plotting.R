#' @title Plotting Utilities
#' @description Shared plotting theme, color palettes, and helper functions
#'   used across all CellDynamicST visualization modules.
#' @name plotting-utils
NULL

#' CellDynamicST Color Palette
#'
#' Returns a vector of colors for consistent visualization across all
#' CellDynamicST plots. Supports palettes for cell types, neurotransmitter
#' types, brain regions, and experimental groups.
#'
#' @param palette Character. Palette name. One of: "cell_type", "nt_type",
#'   "glia_type", "region", "group", "module", "diverging", "sequential".
#'   Default: "cell_type".
#' @param n Integer. Number of colors needed. If NULL, returns the full palette.
#'
#' @return A character vector of hex color codes.
#'
#' @examples
#' cdst_palette("nt_type")
#' cdst_palette("cell_type", n = 10)
#'
#' @export
cdst_palette <- function(palette = "cell_type", n = NULL) {
  palettes <- list(
    cell_type = c(
      "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
      "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
      "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
      "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA",
      "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
      "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1B9E77"
    ),
    nt_type = c(
      "Glutamatergic" = "#E41A1C",
      "GABAergic"     = "#377EB8",
      "Glycinergic"   = "#4DAF4A",
      "Cholinergic"   = "#984EA3",
      "Dopaminergic"  = "#FF7F00",
      "Serotonergic"  = "#FFFF33",
      "Noradrenergic" = "#A65628",
      "Histaminergic" = "#F781BF",
      "Mixed"         = "#999999",
      "Unknown"       = "#CCCCCC"
    ),
    glia_type = c(
      "Astrocyte"     = "#1B9E77",
      "Microglia"     = "#D95F02",
      "Oligodendrocyte" = "#7570B3",
      "OPC"           = "#E7298A",
      "Ependymal"     = "#66A61E",
      "Endothelial"   = "#E6AB02",
      "Pericyte"      = "#A6761D",
      "Unknown"       = "#CCCCCC"
    ),
    region = c(
      "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
      "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
      "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A"
    ),
    group = c(
      "#2166AC", "#B2182B", "#4393C3", "#D6604D",
      "#92C5DE", "#F4A582", "#053061", "#67001F"
    ),
    module = c(
      "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
      "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
      "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5"
    ),
    diverging = c(
      "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
      "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B",
      "#67001F"
    ),
    sequential = c(
      "#FFF7FB", "#ECE7F2", "#D0D1E6", "#A6BDDB", "#74A9CF",
      "#3690C0", "#0570B0", "#045A8D", "#023858"
    )
  )

  if (!palette %in% names(palettes)) {
    cli::cli_abort(c(
      "Unknown palette: {.val {palette}}",
      "i" = "Available palettes: {.val {names(palettes)}}"
    ))
  }

  pal <- palettes[[palette]]

  if (!is.null(n)) {
    if (n <= length(pal)) {
      pal <- pal[seq_len(n)]
    } else {
      pal <- grDevices::colorRampPalette(pal)(n)
    }
  }

  pal
}


#' CellDynamicST ggplot2 Theme
#'
#' A clean, publication-ready ggplot2 theme for all CellDynamicST
#' visualizations. Based on \code{theme_minimal} with customized typography
#' and spacing.
#'
#' @param base_size Numeric. Base font size. Default: 12.
#' @param base_family Character. Base font family. Default: "".
#'
#' @return A ggplot2 theme object.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) + geom_point() + cdst_plot_theme()
#'
#' @export
cdst_plot_theme <- function(base_size = 16, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size * 1.25,
                                          hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size * 0.9,
                                             color = "grey40", hjust = 0),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      axis.text = ggplot2::element_text(size = base_size * 0.85),
      legend.title = ggplot2::element_text(face = "bold", size = base_size * 0.9),
      legend.text = ggplot2::element_text(size = base_size * 0.85),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey92"),
      strip.text = ggplot2::element_text(face = "bold", size = base_size * 0.95),
      plot.margin = ggplot2::margin(12, 12, 12, 12)
    )
}


#' @rdname cdst_plot_theme
#' @export
cdst_theme_publication <- cdst_plot_theme
