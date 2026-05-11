#' @title Module 7: WGCNA Network Analysis
#' @description Functions for interregional weighted gene co-expression
#'   network analysis, including module detection, trait correlation,
#'   hub gene identification, module preservation, and pathway enrichment.
#' @name wgcna
NULL

#' Build Region x Gene Expression Matrix for WGCNA
#'
#' Constructs the region-by-gene average expression matrix used as input for
#' WGCNA. Each row represents a unique region-group combination, and each
#' column represents a gene. Regions with fewer than \code{min_cells} are
#' excluded.
#'
#' This function generalizes the data preparation step from
#' \code{Interregional_WCGNA_module_detect.R}.
#'
#' @param seurat_obj A Seurat object.
#' @param config A \code{\link{CdstConfig}} object.
#' @param region_column Character. Default: "brain_region_L7".
#' @param min_cells Integer. Minimum cells per region-group. Default: 100.
#' @param verbose Logical. Default: TRUE.
#'
#' @return A list with components:
#'   \describe{
#'     \item{datExpr}{Numeric matrix (region-groups x genes)}
#'     \item{traits_df}{Data.frame with region, group, genotype, treatment}
#'   }
#'
#' @export
cdst_wgcna_prepare <- function(seurat_obj, config,
                                region_column = "brain_region_L7",
                                min_cells = 100L,
                                verbose = TRUE) {
  checkmate::assert_class(config, "CdstConfig")
  validate_seurat(seurat_obj,
                  required_metadata = c("experiment_group", region_column))

  if (verbose) cli::cli_inform("Building region x gene matrix for WGCNA...")

  expr_data <- Seurat::GetAssayData(seurat_obj, layer = "data")
  meta <- seurat_obj@meta.data
  genes <- rownames(expr_data)

  # Create region-group identifiers
  meta$region_group <- paste(meta[[region_column]],
                              meta$experiment_group, sep = "_")
  region_groups <- unique(meta$region_group)

  # Calculate average expression per region-group
  avg_expr <- matrix(NA_real_, nrow = length(genes), ncol = length(region_groups))
  rownames(avg_expr) <- genes
  colnames(avg_expr) <- region_groups

  if (verbose) pb <- utils::txtProgressBar(min = 0, max = length(region_groups), style = 3)

  for (i in seq_along(region_groups)) {
    rg <- region_groups[i]
    cells <- which(meta$region_group == rg)

    if (length(cells) >= min_cells) {
      valid_genes <- intersect(genes, rownames(expr_data))
      avg_expr[valid_genes, rg] <- Matrix::rowMeans(
        expr_data[valid_genes, cells, drop = FALSE]
      )
    }
    if (verbose) utils::setTxtProgressBar(pb, i)
  }
  if (verbose) close(pb)

  # Remove columns with NAs (regions below threshold)
  keep_cols <- colSums(is.na(avg_expr)) == 0
  avg_expr <- avg_expr[, keep_cols, drop = FALSE]

  # Transpose: rows = region-groups, cols = genes
  datExpr <- t(avg_expr)

  # Build traits data.frame
  traits_df <- data.frame(
    region_group = rownames(datExpr),
    stringsAsFactors = FALSE
  )

  # Parse region and group from region_group identifier
  split_parts <- strsplit(traits_df$region_group, "_")
  # Handle multi-word group names by finding the last match
  for (i in seq_len(nrow(traits_df))) {
    rg <- traits_df$region_group[i]
    # Find which group matches
    matched_group <- NA_character_
    for (g in config@groups) {
      if (grepl(paste0("_", g, "$"), rg)) {
        matched_group <- g
        break
      }
    }
    if (!is.na(matched_group)) {
      traits_df$region[i] <- sub(paste0("_", matched_group, "$"), "", rg)
      traits_df$experiment_group[i] <- matched_group
    } else {
      parts <- split_parts[[i]]
      traits_df$region[i] <- paste(parts[-length(parts)], collapse = "_")
      traits_df$experiment_group[i] <- parts[length(parts)]
    }
  }

  # Extract genotype and treatment from group names
  traits_df$genotype <- config@design$genotype_map[traits_df$experiment_group]
  traits_df$treatment <- config@design$treatment_map[traits_df$experiment_group]

  if (verbose) {
    cli::cli_inform(c(
      "v" = "WGCNA data prepared.",
      "i" = "Matrix: {nrow(datExpr)} region-groups x {ncol(datExpr)} genes"
    ))
  }

  list(datExpr = datExpr, traits_df = traits_df)
}


#' Run WGCNA Module Detection
#'
#' Performs weighted gene co-expression network analysis to detect gene
#' modules. Includes soft threshold selection, network construction, and
#' module identification via dynamic tree cutting.
#'
#' Generalizes \code{Interregional_WCGNA_module_detect.R}.
#'
#' @param wgcna_data A list from \code{cdst_wgcna_prepare} with datExpr and
#'   traits_df.
#' @param power Integer or NULL. Soft threshold power. If NULL, automatically
#'   selected. Default: NULL.
#' @param min_module_size Integer. Minimum genes per module. Default: 20.
#' @param deep_split Integer. Sensitivity of module detection (0-4). Default: 3.
#' @param merge_cut_height Numeric. Height for merging similar modules.
#'   Default: 0.10.
#' @param verbose Logical. Default: TRUE.
#'
#' @return A list with components:
#'   \describe{
#'     \item{net}{The blockwiseModules result object}
#'     \item{module_colors}{Named character vector of module assignments}
#'     \item{module_eigengenes}{Module eigengene matrix}
#'     \item{power}{Soft threshold power used}
#'     \item{sft}{Soft threshold analysis results}
#'     \item{datExpr}{The input expression matrix}
#'     \item{traits_df}{The traits data.frame}
#'   }
#'
#' @export
cdst_wgcna_detect <- function(wgcna_data,
                               power = NULL,
                               min_module_size = 20L,
                               deep_split = 3L,
                               merge_cut_height = 0.10,
                               verbose = TRUE) {
  datExpr <- wgcna_data$datExpr
  traits_df <- wgcna_data$traits_df

  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg WGCNA} is required. Install from CRAN.")
  }

  # Clean data
  gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    if (verbose) {
      cli::cli_inform("Removed {sum(!gsg$goodSamples)} bad samples and {sum(!gsg$goodGenes)} bad genes.")
    }
  }

  # Pick soft threshold
  if (verbose) cli::cli_inform("Selecting soft threshold power...")
  sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = 1:20, verbose = 0)

  if (is.null(power)) {
    power <- sft$powerEstimate
    if (is.na(power)) {
      power <- 6L
      if (verbose) cli::cli_warn("Could not auto-select power. Using default: {power}")
    }
  }

  if (verbose) cli::cli_inform("Using soft threshold power: {.val {power}}")

  # Module detection
  if (verbose) cli::cli_inform("Detecting modules...")
  net <- WGCNA::blockwiseModules(
    datExpr,
    power = power,
    minModuleSize = min_module_size,
    deepSplit = deep_split,
    pamRespectsDendro = FALSE,
    mergeCutHeight = merge_cut_height,
    numericLabels = FALSE,
    saveTOMs = FALSE,
    verbose = if (verbose) 3 else 0
  )

  module_colors <- net$colors
  n_modules <- length(unique(module_colors))

  # Module eigengenes
  MEs <- WGCNA::moduleEigengenes(datExpr, colors = module_colors)$eigengenes
  MEs <- WGCNA::orderMEs(MEs)

  if (verbose) {
    cli::cli_inform(c(
      "v" = "WGCNA module detection complete.",
      "i" = "Modules detected: {.val {n_modules}} (including grey)",
      "i" = "Module sizes: {paste(sort(table(module_colors), decreasing = TRUE)[1:min(5, n_modules)], collapse = ', ')}..."
    ))
  }

  list(
    net = net,
    module_colors = module_colors,
    module_eigengenes = MEs,
    power = power,
    sft = sft,
    datExpr = datExpr,
    traits_df = traits_df
  )
}


#' Compute Module-Trait Correlations
#'
#' Calculates Pearson correlations and p-values between module eigengenes
#' and experimental traits (genotype, treatment, brain region). Produces
#' the correlation matrices used for heatmap visualization.
#'
#' Generalizes \code{Interregional_WCGNA_all_trait_corr.R}.
#'
#' @param wgcna_result A list from \code{cdst_wgcna_detect}.
#' @param traits Character vector. Traits to correlate. Default: all available.
#' @param verbose Logical. Default: TRUE.
#'
#' @return A list with components:
#'   \describe{
#'     \item{cor_treatment}{Correlation matrix for treatment}
#'     \item{p_treatment}{P-value matrix for treatment}
#'     \item{cor_genotype}{Correlation matrix for genotype}
#'     \item{p_genotype}{P-value matrix for genotype}
#'     \item{cor_region}{Correlation matrix for regions}
#'     \item{p_region}{P-value matrix for regions}
#'   }
#'
#' @export
cdst_wgcna_trait_cor <- function(wgcna_result, traits = NULL, verbose = TRUE) {
  MEs <- wgcna_result$module_eigengenes
  traits_df <- wgcna_result$traits_df
  n_samples <- nrow(MEs)

  results <- list()

  # Treatment correlation
  if ("treatment" %in% colnames(traits_df)) {
    unique_treatments <- unique(traits_df$treatment)
    if (length(unique_treatments) == 2) {
      treatment_bin <- as.numeric(traits_df$treatment == unique_treatments[2])
      treatment_mat <- matrix(treatment_bin, ncol = 1)
      rownames(treatment_mat) <- rownames(MEs)
      colnames(treatment_mat) <- paste(unique_treatments[2], "vs", unique_treatments[1])

      results$cor_treatment <- stats::cor(MEs, treatment_mat, use = "pairwise.complete.obs")
      results$p_treatment <- WGCNA::corPvalueStudent(results$cor_treatment, n_samples)
    }
  }

  # Genotype correlation
  if ("genotype" %in% colnames(traits_df)) {
    unique_genotypes <- unique(traits_df$genotype)
    if (length(unique_genotypes) == 2) {
      genotype_bin <- as.numeric(traits_df$genotype == unique_genotypes[2])
      genotype_mat <- matrix(genotype_bin, ncol = 1)
      rownames(genotype_mat) <- rownames(MEs)
      colnames(genotype_mat) <- paste(unique_genotypes[2], "vs", unique_genotypes[1])

      results$cor_genotype <- stats::cor(MEs, genotype_mat, use = "pairwise.complete.obs")
      results$p_genotype <- WGCNA::corPvalueStudent(results$cor_genotype, n_samples)
    }
  }

  # Region correlations
  if ("region" %in% colnames(traits_df)) {
    unique_regions <- unique(traits_df$region)
    region_mat <- sapply(unique_regions, function(r) {
      as.numeric(traits_df$region == r)
    })
    rownames(region_mat) <- rownames(MEs)
    colnames(region_mat) <- unique_regions

    results$cor_region <- stats::cor(MEs, region_mat, use = "pairwise.complete.obs")
    results$p_region <- WGCNA::corPvalueStudent(results$cor_region, n_samples)
  }

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Module-trait correlations computed.",
      "i" = "Traits analyzed: {paste(names(results)[grepl('cor_', names(results))], collapse = ', ')}"
    ))
  }

  results
}


#' Identify Hub Genes per Module
#'
#' Identifies hub genes within each WGCNA module based on module membership
#' (kME) and intramodular connectivity. Hub genes are those with the highest
#' connectivity within their module.
#'
#' Generalizes \code{Interregional_WCGNA_Module_hub_gene.R}.
#'
#' @param wgcna_result A list from \code{cdst_wgcna_detect}.
#' @param n_hubs Integer. Number of hub genes per module. Default: 10.
#' @param kme_threshold Numeric. Minimum module membership (kME) to be
#'   considered a hub. Default: 0.7.
#' @param verbose Logical. Default: TRUE.
#'
#' @return A data.frame with columns: module, gene, kME, rank.
#'
#' @export
cdst_wgcna_hub_genes <- function(wgcna_result,
                                  n_hubs = 10L,
                                  kme_threshold = 0.7,
                                  verbose = TRUE) {
  datExpr <- wgcna_result$datExpr
  module_colors <- wgcna_result$module_colors
  MEs <- wgcna_result$module_eigengenes

  # Calculate module membership (kME)
  kME <- stats::cor(datExpr, MEs, use = "pairwise.complete.obs")

  modules <- unique(module_colors)
  modules <- modules[modules != "grey"]  # Exclude unassigned

  hub_list <- list()

  for (mod in modules) {
    mod_genes <- names(module_colors[module_colors == mod])
    me_col <- paste0("ME", mod)

    if (!me_col %in% colnames(kME)) next

    mod_kme <- kME[mod_genes, me_col]
    mod_kme <- sort(mod_kme, decreasing = TRUE)

    # Filter by kME threshold
    mod_kme <- mod_kme[mod_kme >= kme_threshold]

    # Take top n_hubs
    top_hubs <- utils::head(mod_kme, n_hubs)

    if (length(top_hubs) > 0) {
      hub_list[[mod]] <- data.frame(
        module = mod,
        gene = names(top_hubs),
        kME = unname(top_hubs),
        rank = seq_along(top_hubs),
        stringsAsFactors = FALSE
      )
    }
  }

  hub_df <- do.call(rbind, hub_list)
  rownames(hub_df) <- NULL

  if (verbose) {
    n_modules_with_hubs <- length(unique(hub_df$module))
    cli::cli_inform(c(
      "v" = "Hub gene identification complete.",
      "i" = "Modules with hubs: {n_modules_with_hubs}",
      "i" = "Total hub genes: {nrow(hub_df)}"
    ))
  }

  hub_df
}


#' Run Module Preservation Analysis
#'
#' Tests whether WGCNA modules detected in one experimental group are
#' preserved in other groups. Uses the WGCNA modulePreservation function.
#'
#' Generalizes \code{Interregional_WCGNA_preserve.R}.
#'
#' @param wgcna_result A list from \code{cdst_wgcna_detect}.
#' @param n_permutations Integer. Number of permutations. Default: 200.
#' @param verbose Logical. Default: TRUE.
#'
#' @return A list with preservation statistics per module per group comparison.
#'
#' @export
cdst_wgcna_preservation <- function(wgcna_result,
                                     n_permutations = 200L,
                                     verbose = TRUE) {
  datExpr <- wgcna_result$datExpr
  traits_df <- wgcna_result$traits_df
  module_colors <- wgcna_result$module_colors

  groups <- unique(traits_df$experiment_group)
  if (length(groups) < 2) {
    cli::cli_abort("At least 2 experimental groups are needed for preservation analysis.")
  }

  # Split data by group
  multiExpr <- lapply(groups, function(grp) {
    idx <- which(traits_df$experiment_group == grp)
    list(data = datExpr[idx, , drop = FALSE])
  })
  names(multiExpr) <- groups

  # Use first group as reference
  multiColor <- list(colors = module_colors)

  if (verbose) {
    cli::cli_inform(c(
      "Running module preservation analysis...",
      "i" = "Reference: {groups[1]}, Test: {paste(groups[-1], collapse = ', ')}",
      "i" = "Permutations: {n_permutations}"
    ))
  }

  preservation <- WGCNA::modulePreservation(
    multiData = multiExpr,
    multiColor = multiColor,
    referenceNetworks = 1,
    nPermutations = n_permutations,
    verbose = if (verbose) 3 else 0
  )

  if (verbose) cli::cli_inform(c("v" = "Module preservation analysis complete."))

  preservation
}


#' Run Pathway Enrichment on WGCNA Modules
#'
#' Performs GO and KEGG pathway enrichment analysis on genes within each
#' WGCNA module. Uses clusterProfiler for enrichment testing.
#'
#' Generalizes \code{Interregional_WCGNA_Pathway_Enrich.R}.
#'
#' @param wgcna_result A list from \code{cdst_wgcna_detect}.
#' @param organism Character. "mouse" or "human". Default: "mouse".
#' @param p_cutoff Numeric. P-value cutoff for enrichment. Default: 0.05.
#' @param q_cutoff Numeric. Q-value cutoff. Default: 0.1.
#' @param verbose Logical. Default: TRUE.
#'
#' @return A data.frame with enrichment results per module.
#'
#' @export
cdst_wgcna_pathway_enrich <- function(wgcna_result,
                                       organism = "mouse",
                                       p_cutoff = 0.05,
                                       q_cutoff = 0.1,
                                       verbose = TRUE) {
  module_colors <- wgcna_result$module_colors
  modules <- unique(module_colors)
  modules <- modules[modules != "grey"]

  if (organism == "mouse") {
    org_db <- "org.Mm.eg.db"
  } else {
    org_db <- "org.Hs.eg.db"
  }

  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg clusterProfiler} is required. Install from Bioconductor.")
  }

  results <- list()

  for (mod in modules) {
    mod_genes <- names(module_colors[module_colors == mod])

    tryCatch({
      # Convert gene symbols to Entrez IDs
      gene_ids <- AnnotationDbi::mapIds(
        get(org_db),
        keys = mod_genes,
        keytype = "SYMBOL",
        column = "ENTREZID"
      )
      gene_ids <- gene_ids[!is.na(gene_ids)]

      if (length(gene_ids) < 5) next

      # GO enrichment
      ego <- clusterProfiler::enrichGO(
        gene = gene_ids,
        OrgDb = get(org_db),
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = p_cutoff,
        qvalueCutoff = q_cutoff,
        readable = TRUE
      )

      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        ego_df <- as.data.frame(ego)
        ego_df$module <- mod
        ego_df$enrichment_type <- "GO_BP"
        results[[paste0(mod, "_GO")]] <- ego_df
      }
    }, error = function(e) {
      if (verbose) cli::cli_warn("Enrichment failed for module {mod}: {e$message}")
    })
  }

  if (length(results) == 0) {
    cli::cli_warn("No enrichment results found.")
    return(data.frame())
  }

  enrich_df <- do.call(rbind, results)
  rownames(enrich_df) <- NULL

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Pathway enrichment complete.",
      "i" = "Modules with enriched pathways: {length(unique(enrich_df$module))}"
    ))
  }

  enrich_df
}
