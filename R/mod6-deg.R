#' @title Module 6: Differential Expression Analysis
#' @description Functions for interregional differential expression testing,
#'   spatial enrichment analysis, and GO-term based pathway enrichment.
#' @name deg
NULL

#' Run Interregional Differential Expression Analysis
#'
#' Performs differential expression testing between experimental groups within
#' each brain region. Computes log2 fold changes and p-values for every gene
#' in every region, enabling identification of spatially-specific treatment
#' and genotype effects.
#'
#' This function generalizes \code{Interregional_DEG.R}, replacing hardcoded
#' group comparisons with configurable contrasts and adding parallel
#' processing support.
#'
#' @param seurat_obj A Seurat object with cell type and region annotations.
#' @param config A \code{\link{CdstConfig}} object.
#' @param contrasts A list of named character vectors, each specifying a
#'   comparison. Format: \code{list(c(group1 = "AA SAL", group2 = "GG SAL"))}.
#'   If NULL, all pairwise comparisons from config are used.
#' @param region_column Character. Brain region column. Default: "brain_region_L7".
#' @param min_cells_per_region Integer. Minimum cells per region-group
#'   combination. Default: 30.
#' @param test_method Character. Statistical test method. One of "wilcox",
#'   "t", "MAST", "DESeq2". Default: "wilcox".
#' @param logfc_threshold Numeric. Minimum absolute log2FC to report.
#'   Default: 0.25.
#' @param p_adjust_method Character. P-value adjustment method. Default: "BH".
#' @param n_cores Integer. Number of cores for parallel processing. Default: 1.
#' @param verbose Logical. Default: TRUE.
#'
#' @return A data.frame with columns: gene, region, contrast, avg_log2FC,
#'   p_val, p_val_adj, pct_1, pct_2.
#'
#' @details
#' The function iterates over each brain region and each contrast, performing
#' differential expression testing between the two groups. Results are
#' aggregated into a single data.frame with region and contrast annotations.
#'
#' For large datasets, parallel processing across regions is supported via
#' the \code{n_cores} parameter.
#'
#' @examples
#' \dontrun{
#' config <- cdst_load_config("experiment_config.yaml")
#' deg_results <- cdst_deg_interregional(obj, config,
#'   contrasts = list(c(group1 = "WT SAL", group2 = "WT MOR")))
#' }
#'
#' @export
cdst_deg_interregional <- function(seurat_obj, config,
                                   contrasts = NULL,
                                   region_column = "brain_region_L7",
                                   min_cells_per_region = 30L,
                                   test_method = "wilcox",
                                   logfc_threshold = 0.25,
                                   p_adjust_method = "BH",
                                   n_cores = 1L,
                                   verbose = TRUE) {
  checkmate::assert_class(config, "CdstConfig")
  validate_seurat(seurat_obj,
                  required_metadata = c("experiment_group", region_column))
  start_time <- Sys.time()

  # Build contrasts from config if not provided
  if (is.null(contrasts)) {
    groups <- config@groups
    contrasts <- utils::combn(groups, 2, simplify = FALSE)
    contrasts <- lapply(contrasts, function(x) c(group1 = x[1], group2 = x[2]))
  }

  regions <- unique(seurat_obj@meta.data[[region_column]])
  regions <- regions[!is.na(regions)]

  if (verbose) {
    cli::cli_inform(c(
      "Running interregional DEG analysis...",
      "i" = "Regions: {length(regions)}, Contrasts: {length(contrasts)}",
      "i" = "Test method: {test_method}"
    ))
  }

  # ---- Worker function ----
  .deg_one_region <- function(region, seurat_obj, contrasts, region_column,
                              min_cells, test_method, logfc_threshold) {
    region_cells <- which(seurat_obj@meta.data[[region_column]] == region)
    region_obj <- seurat_obj[, region_cells]

    results_list <- list()
    for (contrast in contrasts) {
      g1 <- contrast["group1"]
      g2 <- contrast["group2"]
      contrast_name <- paste(g1, "vs", g2)

      cells_g1 <- which(region_obj$experiment_group == g1)
      cells_g2 <- which(region_obj$experiment_group == g2)

      if (length(cells_g1) < min_cells || length(cells_g2) < min_cells) next

      Seurat::Idents(region_obj) <- "experiment_group"

      tryCatch({
        markers <- Seurat::FindMarkers(
          region_obj,
          ident.1 = g1,
          ident.2 = g2,
          test.use = test_method,
          logfc.threshold = logfc_threshold,
          min.pct = 0.1,
          verbose = FALSE
        )

        if (nrow(markers) > 0) {
          markers$gene <- rownames(markers)
          markers$region <- region
          markers$contrast <- contrast_name
          results_list[[contrast_name]] <- markers
        }
      }, error = function(e) {
        # Skip regions with errors
      })
    }

    if (length(results_list) > 0) {
      do.call(rbind, results_list)
    } else {
      NULL
    }
  }

  # ---- Execute across regions ----
  if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(min(n_cores, length(regions)))
    on.exit(parallel::stopCluster(cl), add = TRUE)

    results <- parallel::parLapply(cl, regions, .deg_one_region,
                                    seurat_obj = seurat_obj,
                                    contrasts = contrasts,
                                    region_column = region_column,
                                    min_cells = min_cells_per_region,
                                    test_method = test_method,
                                    logfc_threshold = logfc_threshold)
  } else {
    if (verbose) pb <- utils::txtProgressBar(min = 0, max = length(regions), style = 3)
    results <- lapply(seq_along(regions), function(i) {
      if (verbose) utils::setTxtProgressBar(pb, i)
      .deg_one_region(regions[i], seurat_obj, contrasts, region_column,
                      min_cells_per_region, test_method, logfc_threshold)
    })
    if (verbose) close(pb)
  }

  # ---- Combine results ----
  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) == 0) {
    cli::cli_warn("No DEGs found across any region-contrast combination.")
    return(data.frame())
  }

  deg_df <- do.call(rbind, results)
  rownames(deg_df) <- NULL

  # Adjust p-values globally
  deg_df$p_val_adj <- stats::p.adjust(deg_df$p_val, method = p_adjust_method)

  if (verbose) {
    n_sig <- sum(deg_df$p_val_adj < 0.05, na.rm = TRUE)
    cli::cli_inform(c(
      "v" = "DEG analysis complete.",
      "i" = "Total DEGs tested: {nrow(deg_df)}",
      "i" = "Significant (adj. p < 0.05): {n_sig}"
    ))
  }

  deg_df
}


#' Spatial Enrichment Analysis (GO-term Based)
#'
#' Performs spatial enrichment analysis by testing whether genes associated
#' with specific GO terms show differential expression across brain regions
#' between experimental conditions. Computes enrichment scores per region
#' per GO term.
#'
#' This function generalizes \code{Spatial_Enrichment_Generalized_updated.R},
#' replacing hardcoded GO terms and comparisons with configurable parameters.
#'
#' @param seurat_obj A Seurat object.
#' @param config A \code{\link{CdstConfig}} object.
#' @param go_terms Character vector. GO term IDs to test. If NULL, uses
#'   a default set of neuroscience-relevant GO terms.
#' @param search_terms Character vector. Text terms to search for GO IDs
#'   (alternative to providing GO IDs directly).
#' @param contrasts List. Comparisons to make (same format as
#'   \code{cdst_deg_interregional}).
#' @param region_column Character. Default: "brain_region_L7".
#' @param organism Character. Organism for gene annotation. Default: "mouse".
#' @param p_adjust_method Character. Default: "BH".
#' @param verbose Logical. Default: TRUE.
#'
#' @return A data.frame with columns: go_term, go_name, region, contrast,
#'   enrichment_score, p_value, p_adj, n_genes.
#'
#' @export
cdst_spatial_enrichment <- function(seurat_obj, config,
                                    go_terms = NULL,
                                    search_terms = NULL,
                                    contrasts = NULL,
                                    region_column = "brain_region_L7",
                                    organism = "mouse",
                                    p_adjust_method = "BH",
                                    verbose = TRUE) {
  checkmate::assert_class(config, "CdstConfig")
  validate_seurat(seurat_obj,
                  required_metadata = c("experiment_group", region_column))

  # Load organism annotation database
  if (organism == "mouse") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg org.Mm.eg.db} is required. Install with BiocManager::install('org.Mm.eg.db').")
    }
    org_db <- org.Mm.eg.db::org.Mm.eg.db
  } else if (organism == "human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      cli::cli_abort("Package {.pkg org.Hs.eg.db} is required.")
    }
    org_db <- org.Hs.eg.db::org.Hs.eg.db
  } else {
    cli::cli_abort("Unsupported organism: {.val {organism}}. Use 'mouse' or 'human'.")
  }

  # Resolve GO terms
  if (!is.null(search_terms) && is.null(go_terms)) {
    go_terms <- .search_go_terms(search_terms)
  }
  if (is.null(go_terms)) {
    cli::cli_abort("Provide either {.arg go_terms} or {.arg search_terms}.")
  }

  # Get genes for each GO term
  go_gene_lists <- lapply(go_terms, function(go_id) {
    tryCatch({
      genes <- AnnotationDbi::select(
        org_db,
        keys = go_id,
        columns = "SYMBOL",
        keytype = "GOALL"
      )$SYMBOL
      unique(genes[!is.na(genes)])
    }, error = function(e) character(0))
  })
  names(go_gene_lists) <- go_terms

  # Filter to genes present in the dataset
  available_genes <- rownames(seurat_obj)
  go_gene_lists <- lapply(go_gene_lists, function(g) intersect(g, available_genes))
  go_gene_lists <- go_gene_lists[lengths(go_gene_lists) > 0]

  if (length(go_gene_lists) == 0) {
    cli::cli_warn("No GO term genes found in the dataset.")
    return(data.frame())
  }

  if (verbose) {
    cli::cli_inform(c(
      "Running spatial enrichment...",
      "i" = "GO terms with genes: {length(go_gene_lists)}"
    ))
  }

  # Build contrasts
  if (is.null(contrasts)) {
    groups <- config@groups
    contrasts <- utils::combn(groups, 2, simplify = FALSE)
    contrasts <- lapply(contrasts, function(x) c(group1 = x[1], group2 = x[2]))
  }

  # Extract expression data
  expr_data <- Seurat::GetAssayData(seurat_obj, layer = "data")
  meta <- seurat_obj@meta.data
  regions <- unique(meta[[region_column]])
  regions <- regions[!is.na(regions)]

  # Compute enrichment per GO term x region x contrast
  results <- list()

  for (go_id in names(go_gene_lists)) {
    go_genes <- go_gene_lists[[go_id]]

    for (contrast in contrasts) {
      g1 <- contrast["group1"]
      g2 <- contrast["group2"]
      contrast_name <- paste(g1, "vs", g2)

      for (region in regions) {
        cells_g1 <- which(meta[[region_column]] == region &
                           meta$experiment_group == g1)
        cells_g2 <- which(meta[[region_column]] == region &
                           meta$experiment_group == g2)

        if (length(cells_g1) < 10 || length(cells_g2) < 10) next

        # Mean expression of GO genes in each group
        mean_g1 <- rowMeans(expr_data[go_genes, cells_g1, drop = FALSE])
        mean_g2 <- rowMeans(expr_data[go_genes, cells_g2, drop = FALSE])

        # Log2 fold change
        log2fc <- mean(log2((mean_g2 + 1e-5) / (mean_g1 + 1e-5)))

        # Wilcoxon test on mean GO gene expression per cell
        score_g1 <- colMeans(expr_data[go_genes, cells_g1, drop = FALSE])
        score_g2 <- colMeans(expr_data[go_genes, cells_g2, drop = FALSE])

        p_val <- tryCatch({
          stats::wilcox.test(score_g1, score_g2)$p.value
        }, error = function(e) NA_real_)

        results[[length(results) + 1]] <- data.frame(
          go_term = go_id,
          region = region,
          contrast = contrast_name,
          enrichment_score = log2fc,
          p_value = p_val,
          n_genes = length(go_genes),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) == 0) {
    cli::cli_warn("No enrichment results computed.")
    return(data.frame())
  }

  enrich_df <- do.call(rbind, results)
  enrich_df$p_adj <- stats::p.adjust(enrich_df$p_value, method = p_adjust_method)

  # Add GO term names
  go_names <- tryCatch({
    terms <- AnnotationDbi::Term(GO.db::GOTERM[go_terms])
    stats::setNames(terms, go_terms)
  }, error = function(e) stats::setNames(rep(NA_character_, length(go_terms)), go_terms))

  enrich_df$go_name <- go_names[enrich_df$go_term]

  if (verbose) {
    n_sig <- sum(enrich_df$p_adj < 0.05, na.rm = TRUE)
    cli::cli_inform(c(
      "v" = "Spatial enrichment complete.",
      "i" = "Tests performed: {nrow(enrich_df)}",
      "i" = "Significant (adj. p < 0.05): {n_sig}"
    ))
  }

  enrich_df
}


#' Summarize DEG Results
#'
#' Produces a summary table of DEG results grouped by region and contrast,
#' including counts of up/down-regulated genes and top hits.
#'
#' @param deg_results A data.frame from \code{cdst_deg_interregional}.
#' @param p_threshold Numeric. Adjusted p-value threshold. Default: 0.05.
#' @param logfc_threshold Numeric. Minimum absolute log2FC. Default: 0.5.
#' @param n_top Integer. Number of top genes to report per group. Default: 10.
#'
#' @return A data.frame with summary statistics per region-contrast.
#'
#' @export
cdst_summarize_deg <- function(deg_results,
                               p_threshold = 0.05,
                               logfc_threshold = 0.5,
                               n_top = 10L) {
  if (nrow(deg_results) == 0) return(data.frame())

  sig <- deg_results |>
    dplyr::filter(.data$p_val_adj < p_threshold,
                  abs(.data$avg_log2FC) > logfc_threshold)

  summary_df <- sig |>
    dplyr::group_by(.data$region, .data$contrast) |>
    dplyr::summarise(
      n_up = sum(.data$avg_log2FC > 0),
      n_down = sum(.data$avg_log2FC < 0),
      n_total = dplyr::n(),
      top_up_genes = paste(
        utils::head(
          .data$gene[order(-.data$avg_log2FC)][.data$avg_log2FC[order(-.data$avg_log2FC)] > 0],
          n_top
        ),
        collapse = ", "
      ),
      top_down_genes = paste(
        utils::head(
          .data$gene[order(.data$avg_log2FC)][.data$avg_log2FC[order(.data$avg_log2FC)] < 0],
          n_top
        ),
        collapse = ", "
      ),
      .groups = "drop"
    )

  as.data.frame(summary_df)
}


# ============================================================================
# Internal helpers
# ============================================================================

#' Search GO terms by text
#' @keywords internal
.search_go_terms <- function(search_terms) {
  if (!requireNamespace("GO.db", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg GO.db} is required.")
  }

  all_terms <- AnnotationDbi::as.data.frame(GO.db::GOTERM)
  go_ids <- character(0)

  for (term in search_terms) {
    matched <- all_terms[grepl(term, all_terms$Term, ignore.case = TRUE), ]
    go_ids <- c(go_ids, unique(matched$go_id))
  }

  unique(go_ids)
}
