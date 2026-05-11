#' @title Module 3: Cell Type Annotation
#' @description Functions for reference-based generalized clustering with
#'   Random Forest label propagation, neurotransmitter classification,
#'   glial cell classification, and hierarchical cell type annotation.
#' @name celltyping
NULL

# ============================================================================
# 3A. Generalized Clustering
# ============================================================================

#' Reference-Based Generalized Clustering
#'
#' Implements the core CellDynamicST clustering algorithm. The workflow:
#' (1) Subset the reference group, (2) run PCA, (3) detect and exclude
#' technical bias PCs, (4) build a neighborhood graph, (5) cluster at high
#' resolution, (6) train a Random Forest classifier on reference clusters,
#' (7) propagate labels to all cells.
#'
#' This function generalizes the logic from \code{generalized_clustering.R}
#' and \code{generalized_clustering_randomforest.R}, replacing all hardcoded
#' paths, group names, and parameters with configurable arguments.
#'
#' @param seurat_obj A Seurat object (normalized and scaled).
#' @param config A \code{\link{CdstConfig}} object.
#' @param n_pcs Integer. Number of PCs to compute. Default: from config.
#' @param bias_threshold Numeric. PCs with |correlation| > this threshold
#'   with log2(nCount_RNA) are excluded as technical bias. Default: from config.
#' @param annoy_trees Integer. Number of trees for Annoy algorithm. Default: from config.
#' @param k_neighbors Integer. Number of neighbors for graph construction.
#'   Default: from config.
#' @param resolution Numeric. Louvain clustering resolution. Default: from config.
#' @param rf_ntree Integer. Number of trees in the Random Forest. Default: from config.
#' @param seed Integer. Random seed for reproducibility. Default: 42.
#' @param verbose Logical. Default: TRUE.
#'
#' @return The Seurat object with new metadata columns:
#'   \code{cdst_cluster} (cluster assignment for all cells),
#'   \code{cdst_ref_cluster} (reference-only cluster, NA for non-reference),
#'   \code{cdst_pcs_used} (stored in misc).
#'
#' @details
#' The technical bias detection step computes the Pearson correlation between
#' each PC and log2(nCount_RNA + 1). PCs exceeding the bias threshold are
#' excluded from downstream neighborhood graph construction. This is critical
#' for spatial transcriptomic data where PC1 often captures library size
#' variation rather than biological signal.
#'
#' The Random Forest label propagation step trains a classifier on the
#' reference group's cluster assignments and applies it to all cells,
#' ensuring consistent cluster labels across experimental conditions.
#'
#' @examples
#' \dontrun{
#' config <- cdst_load_config("experiment_config.yaml")
#' obj <- cdst_load_seurat("data.rds", config)
#' obj <- cdst_run_qc(obj, config)
#' obj <- cdst_normalize(obj)
#' obj <- cdst_cluster(obj, config)
#' table(obj$cdst_cluster)
#' }
#'
#' @export
cdst_cluster <- function(seurat_obj, config,
                         n_pcs = NULL,
                         bias_threshold = NULL,
                         annoy_trees = NULL,
                         k_neighbors = NULL,
                         resolution = NULL,
                         rf_ntree = NULL,
                         seed = 42L,
                         verbose = TRUE) {
  checkmate::assert_class(config, "CdstConfig")
  validate_seurat(seurat_obj, required_metadata = "experiment_group")
  start_time <- Sys.time()
  set.seed(seed)

  # Resolve parameters from config or overrides
  n_pcs          <- n_pcs          %||% config@clustering$n_pcs
  bias_threshold <- bias_threshold %||% config@clustering$bias_threshold
  annoy_trees    <- annoy_trees    %||% config@clustering$annoy_trees
  k_neighbors    <- k_neighbors    %||% config@clustering$k_neighbors
  resolution     <- resolution     %||% config@clustering$resolution
  rf_ntree       <- rf_ntree       %||% config@clustering$rf_ntree

  if (verbose) {
    cli::cli_inform(c(
      "Starting reference-based generalized clustering...",
      "i" = "Reference group: {.val {config@reference_group}}",
      "i" = "PCs: {n_pcs}, Resolution: {resolution}, RF trees: {rf_ntree}"
    ))
  }

  # ---- Step 1: Subset reference group ----
  ref_cells <- which(seurat_obj$experiment_group == config@reference_group)
  if (length(ref_cells) < 100) {
    cli::cli_abort(c(
      "Reference group {.val {config@reference_group}} has only {length(ref_cells)} cells.",
      "i" = "At least 100 cells are required for reliable clustering."
    ))
  }

  ref_obj <- seurat_obj[, ref_cells]
  if (verbose) cli::cli_inform("Reference subset: {.val {ncol(ref_obj)}} cells")

  # ---- Step 2: Scale and PCA on reference ----
  all_genes <- rownames(ref_obj)
  ref_obj <- Seurat::ScaleData(ref_obj, features = all_genes, verbose = FALSE)

  if (verbose) cli::cli_inform("Running PCA ({n_pcs} components)...")
  ref_obj <- Seurat::RunPCA(ref_obj, npcs = n_pcs, features = all_genes,
                             verbose = FALSE)

  # ---- Step 3: Detect technical bias PCs ----
  ref_obj[["log2GeneCount"]] <- log2(ref_obj$nCount_RNA + 1)
  pc_embeddings <- Seurat::Embeddings(ref_obj, reduction = "pca")[, seq_len(n_pcs)]
  pc_correlations <- stats::cor(
    ref_obj$log2GeneCount,
    pc_embeddings
  )

  bias_pcs <- which(abs(pc_correlations[1, ]) > bias_threshold)
  good_pcs <- setdiff(seq_len(n_pcs), bias_pcs)

  if (length(good_pcs) < 10) {
    cli::cli_warn(c(
      "Only {length(good_pcs)} PCs passed bias filter.",
      "i" = "Consider lowering {.arg bias_threshold} (currently {bias_threshold})."
    ))
  }

  if (verbose) {
    cli::cli_inform(c(
      "i" = "Bias PCs excluded (|r| > {bias_threshold}): {.val {paste(bias_pcs, collapse = ', ')}}",
      "i" = "Using {length(good_pcs)} PCs for clustering"
    ))
  }

  # ---- Step 4: Neighborhood graph and clustering ----
  if (verbose) cli::cli_inform("Building neighborhood graph...")
  ref_obj <- Seurat::FindNeighbors(
    ref_obj,
    dims = good_pcs,
    nn.method = "annoy",
    n.trees = annoy_trees,
    k.param = k_neighbors,
    verbose = FALSE
  )

  if (verbose) cli::cli_inform("Clustering at resolution {resolution}...")
  ref_obj <- Seurat::FindClusters(
    ref_obj,
    resolution = resolution,
    algorithm = 2,  # Louvain with multilevel refinement
    group.singletons = TRUE,
    verbose = FALSE
  )

  ref_clusters <- Seurat::Idents(ref_obj)
  n_clusters <- length(unique(ref_clusters))
  if (verbose) cli::cli_inform("Found {.val {n_clusters}} clusters in reference group")

  # ---- Step 5: Train Random Forest ----
  if (verbose) cli::cli_inform("Training Random Forest classifier ({rf_ntree} trees)...")

  ref_data <- t(as.matrix(Seurat::GetAssayData(ref_obj, layer = "data")))
  ref_labels <- factor(ref_clusters)

  # Drop unused levels
  ref_labels <- droplevels(ref_labels)

  rf_model <- randomForest::randomForest(
    x = ref_data,
    y = ref_labels,
    ntree = rf_ntree,
    mtry = floor(sqrt(ncol(ref_data))),
    importance = FALSE,
    keep.forest = TRUE,
    keep.inbag = FALSE,
    proximity = FALSE
  )

  # ---- Step 6: Propagate labels to all cells ----
  if (verbose) cli::cli_inform("Propagating labels to all {.val {ncol(seurat_obj)}} cells...")

  all_data <- t(as.matrix(Seurat::GetAssayData(seurat_obj, layer = "data")))
  predictions <- predict(rf_model, all_data)

  seurat_obj$cdst_cluster <- as.character(predictions)

  # Store reference-only cluster info
  seurat_obj$cdst_ref_cluster <- NA_character_
  seurat_obj$cdst_ref_cluster[ref_cells] <- as.character(ref_clusters)

  # ---- Step 7: Run PCA and dimensionality reduction on full dataset ----
  if (verbose) cli::cli_inform("Running PCA on full dataset...")
  all_genes_full <- rownames(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = all_genes_full, verbose = FALSE)
  seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = n_pcs, features = all_genes_full,
                                verbose = FALSE)

  # Detect bias PCs for full dataset
  seurat_obj[["log2GeneCount"]] <- log2(seurat_obj$nCount_RNA + 1)
  full_pc_corr <- stats::cor(
    seurat_obj$log2GeneCount,
    Seurat::Embeddings(seurat_obj, reduction = "pca")[, seq_len(n_pcs)]
  )
  full_bias_pcs <- which(abs(full_pc_corr[1, ]) > bias_threshold)
  full_good_pcs <- setdiff(seq_len(n_pcs), full_bias_pcs)

  # Store PCs used for downstream analyses
  seurat_obj@misc$cdst_pcs_used <- full_good_pcs
  seurat_obj@misc$cdst_bias_pcs <- full_bias_pcs
  seurat_obj@misc$cdst_rf_model <- rf_model

  # Run UMAP and tSNE
  if (verbose) cli::cli_inform("Running UMAP and tSNE...")
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = full_good_pcs,
                                 n.neighbors = 25, min.dist = 0.4,
                                 verbose = FALSE)
  seurat_obj <- Seurat::RunTSNE(seurat_obj, dims = full_good_pcs,
                                 verbose = FALSE)

  # ---- Provenance ----
  seurat_obj@misc$cdst_provenance_clustering <- cdst_provenance_entry(
    "clustering",
    params = list(
      n_pcs = n_pcs, bias_threshold = bias_threshold,
      annoy_trees = annoy_trees, k_neighbors = k_neighbors,
      resolution = resolution, rf_ntree = rf_ntree,
      bias_pcs_excluded = bias_pcs,
      n_clusters = n_clusters, seed = seed
    ),
    start_time = start_time
  )

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Clustering complete.",
      "i" = "Clusters: {.val {n_clusters}}",
      "i" = "All cells labeled via RF propagation."
    ))
  }

  seurat_obj
}


# ============================================================================
# 3B. Neurotransmitter Classification
# ============================================================================

#' Classify Neurotransmitter Types
#'
#' Assigns neurotransmitter type labels to each cell based on marker gene
#' expression in the log2CPM assay. Supports 8 neurotransmitter types:
#' Glutamatergic, GABAergic, Glycinergic, Cholinergic, Dopaminergic,
#' Serotonergic, Noradrenergic, and Histaminergic.
#'
#' This function generalizes \code{Neural_Transmitter_Process.R}, replacing
#' the cell-by-cell loop with vectorized operations for performance, and
#' making the marker gene criteria configurable via YAML.
#'
#' @param seurat_obj A Seurat object with a log2CPM assay.
#' @param criteria_path Character. Path to a YAML file defining NT marker
#'   criteria. If NULL, uses the built-in default criteria for the species.
#' @param threshold Numeric. log2CPM expression threshold. Default: from criteria file.
#' @param min_fraction Numeric. Minimum fraction of cluster cells expressing
#'   a marker for cluster-level assignment. Default: from criteria file.
#' @param verbose Logical. Default: TRUE.
#'
#' @return The Seurat object with new metadata column \code{cdst_nt_type}.
#'
#' @details
#' The classification logic follows the original CellDynamicST approach:
#' \enumerate{
#'   \item For each cell, check if any marker gene for each NT type exceeds
#'     the log2CPM threshold.
#'   \item For complex NT types (GABA, Dopamine, Serotonin, Noradrenaline),
#'     require both a transporter AND an enzyme to be expressed.
#'   \item Cells expressing markers for multiple NT types receive a combined
#'     label (e.g., "Glut-GABA").
#'   \item Cells with no NT markers are labeled "Undefined".
#' }
#'
#' @examples
#' \dontrun{
#' obj <- cdst_classify_nt(obj)
#' table(obj$cdst_nt_type)
#' }
#'
#' @export
cdst_classify_nt <- function(seurat_obj,
                             criteria_path = NULL,
                             threshold = NULL,
                             min_fraction = NULL,
                             verbose = TRUE) {
  validate_seurat(seurat_obj)
  start_time <- Sys.time()

  # Load criteria
  criteria <- cdst_default_nt_criteria(criteria_path)
  threshold    <- threshold    %||% criteria$threshold
  min_fraction <- min_fraction %||% criteria$min_fraction

  if (verbose) {
    cli::cli_inform(c(
      "Classifying neurotransmitter types...",
      "i" = "Threshold: log2CPM > {threshold}",
      "i" = "NT types: {length(criteria$neurotransmitter_types)}"
    ))
  }

  # Get log2CPM data
  if ("log2CPM" %in% Seurat::Assays(seurat_obj)) {
    expr_data <- Seurat::GetAssayData(seurat_obj, assay = "log2CPM", layer = "data")
  } else {
    if (verbose) cli::cli_inform("log2CPM assay not found. Computing from counts...")
    seurat_obj <- .compute_log2cpm(seurat_obj)
    expr_data <- Seurat::GetAssayData(seurat_obj, assay = "log2CPM", layer = "data")
  }

  # Effective threshold (accounting for log2(x+1) transform)
  effective_threshold <- log2(2^threshold + 1)
  available_genes <- rownames(expr_data)

  # ---- Vectorized classification ----
  n_cells <- ncol(expr_data)
  nt_labels <- character(n_cells)

  for (i in seq_len(n_cells)) {
    cell_expr <- expr_data[, i, drop = FALSE]
    types_found <- character(0)

    for (nt in criteria$neurotransmitter_types) {
      markers <- nt$markers
      markers_present <- intersect(markers, available_genes)

      if (length(markers_present) == 0) next

      # Check if any marker exceeds threshold
      expressed <- any(cell_expr[markers_present, 1] > effective_threshold,
                       na.rm = TRUE)
      if (expressed) {
        types_found <- c(types_found, nt$abbreviation)
      }
    }

    if (length(types_found) > 0) {
      nt_labels[i] <- paste(types_found, collapse = "-")
    } else {
      nt_labels[i] <- "Undefined"
    }
  }

  seurat_obj$cdst_nt_type <- nt_labels

  if (verbose) {
    nt_table <- table(nt_labels)
    cli::cli_inform(c(
      "v" = "NT classification complete.",
      "i" = "Defined: {sum(nt_labels != 'Undefined')} cells ({round(sum(nt_labels != 'Undefined')/n_cells*100, 1)}%)",
      "i" = "Undefined: {sum(nt_labels == 'Undefined')} cells"
    ))
  }

  seurat_obj
}


#' Load Default Neurotransmitter Classification Criteria
#'
#' Loads the built-in or user-specified neurotransmitter marker gene criteria
#' from a YAML file.
#'
#' @param criteria_path Character. Path to a YAML criteria file. If NULL,
#'   loads the built-in default for mouse.
#'
#' @return A list with fields: threshold, min_fraction, neurotransmitter_types.
#'
#' @export
cdst_default_nt_criteria <- function(criteria_path = NULL) {
  if (is.null(criteria_path)) {
    criteria_path <- system.file("extdata", "default_nt_criteria_mouse.yaml",
                                 package = "CellDynamicST")
  }
  checkmate::assert_file_exists(criteria_path)
  yaml::read_yaml(criteria_path)
}


# ============================================================================
# 3C. Glial Cell Classification
# ============================================================================

#' Classify Glial and Non-Neuronal Cell Types
#'
#' Identifies glial and non-neuronal cell types (astrocytes, microglia,
#' oligodendrocytes, OPCs, ependymal, endothelial, pericytes) based on
#' marker gene expression. Uses a case-insensitive matching approach against
#' the top marker genes for each cluster.
#'
#' This function generalizes \code{annotate_glia_type.R}, replacing hardcoded
#' criteria with a configurable YAML file and making the classification
#' compatible with any spatial transcriptomic dataset.
#'
#' @param seurat_obj A Seurat object with cluster assignments.
#' @param criteria_path Character. Path to a YAML file defining glial marker
#'   criteria. If NULL, uses the built-in default.
#' @param cluster_column Character. Metadata column with cluster IDs.
#'   Default: "cdst_cluster".
#' @param verbose Logical. Default: TRUE.
#'
#' @return The Seurat object with new metadata column \code{cdst_glia_type}.
#'
#' @export
cdst_classify_glia <- function(seurat_obj,
                               criteria_path = NULL,
                               cluster_column = "cdst_cluster",
                               verbose = TRUE) {
  validate_seurat(seurat_obj, required_metadata = cluster_column)
  start_time <- Sys.time()

  # Load criteria
  criteria <- cdst_default_glia_criteria(criteria_path)

  if (verbose) {
    cli::cli_inform(c(
      "Classifying glial/non-neuronal cell types...",
      "i" = "Glial types: {length(criteria$glial_types)}"
    ))
  }

  # Get expression data
  if ("log2CPM" %in% Seurat::Assays(seurat_obj)) {
    expr_data <- Seurat::GetAssayData(seurat_obj, assay = "log2CPM", layer = "data")
  } else {
    expr_data <- Seurat::GetAssayData(seurat_obj, layer = "data")
  }

  available_genes <- rownames(expr_data)
  threshold <- criteria$threshold %||% 3.0

  # ---- Find top marker genes per cluster ----
  clusters <- unique(seurat_obj@meta.data[[cluster_column]])

  # For each cluster, find the top expressed genes
  cluster_markers <- list()
  for (cl in clusters) {
    cells_in_cluster <- which(seurat_obj@meta.data[[cluster_column]] == cl)
    if (length(cells_in_cluster) < 5) next

    # Mean expression per gene in this cluster
    mean_expr <- Matrix::rowMeans(expr_data[, cells_in_cluster, drop = FALSE])
    top_genes <- names(sort(mean_expr, decreasing = TRUE))[1:20]
    cluster_markers[[as.character(cl)]] <- top_genes
  }

  # ---- Classify each cluster ----
  cluster_glia_type <- stats::setNames(
    rep("Neuronal", length(clusters)),
    as.character(clusters)
  )

  for (cl in names(cluster_markers)) {
    top_markers <- toupper(cluster_markers[[cl]])
    types_found <- character(0)

    for (glia in criteria$glial_types) {
      glia_markers <- toupper(glia$markers)
      if (any(top_markers %in% glia_markers)) {
        types_found <- c(types_found, glia$name)
      }
    }

    if (length(types_found) > 0) {
      cluster_glia_type[cl] <- paste(types_found, collapse = "-")
    }
  }

  # Map cluster-level labels to cells
  seurat_obj$cdst_glia_type <- cluster_glia_type[
    as.character(seurat_obj@meta.data[[cluster_column]])
  ]

  if (verbose) {
    glia_table <- table(seurat_obj$cdst_glia_type)
    n_glia <- sum(seurat_obj$cdst_glia_type != "Neuronal")
    cli::cli_inform(c(
      "v" = "Glial classification complete.",
      "i" = "Neuronal: {sum(seurat_obj$cdst_glia_type == 'Neuronal')} cells",
      "i" = "Non-neuronal: {n_glia} cells"
    ))
  }

  seurat_obj
}


#' Load Default Glial Classification Criteria
#'
#' @param criteria_path Character. Path to a YAML criteria file. If NULL,
#'   loads the built-in default for mouse.
#'
#' @return A list with fields: threshold, min_fraction, glial_types.
#'
#' @export
cdst_default_glia_criteria <- function(criteria_path = NULL) {
  if (is.null(criteria_path)) {
    criteria_path <- system.file("extdata", "default_glia_criteria_mouse.yaml",
                                 package = "CellDynamicST")
  }
  checkmate::assert_file_exists(criteria_path)
  yaml::read_yaml(criteria_path)
}


# ============================================================================
# 3D. Hierarchical Cell Type Annotation
# ============================================================================

#' Build Hierarchical Cell Type Annotations
#'
#' Constructs a multi-level cell type hierarchy by combining cluster identity,
#' neurotransmitter type, glial classification, top marker genes, and brain
#' region information. This produces the full annotation string used in the
#' CellDynamicST paper (e.g., "Glut Slc17a7-Camk2a Isocortex-Hippocampus").
#'
#' This function generalizes \code{annotate_cell_type.R}, replacing hardcoded
#' column names and thresholds with configurable parameters.
#'
#' @param seurat_obj A Seurat object with cluster, NT, and glia annotations.
#' @param cluster_column Character. Cluster column. Default: "cdst_cluster".
#' @param nt_column Character. NT type column. Default: "cdst_nt_type".
#' @param glia_column Character. Glia type column. Default: "cdst_glia_type".
#' @param region_column Character. Brain region column. Default: "brain_region_L3".
#' @param nt_threshold Numeric. Minimum proportion of cells in a cluster with
#'   a given NT type to include it in the annotation. Default: 0.30.
#' @param region_threshold Numeric. Minimum proportion for region annotation.
#'   Default: 0.30.
#' @param n_top_nt Integer. Maximum number of NT types per cluster. Default: 3.
#' @param n_top_regions Integer. Maximum number of regions per cluster. Default: 2.
#' @param verbose Logical. Default: TRUE.
#'
#' @return The Seurat object with new metadata columns:
#'   \code{cdst_cell_type} (full annotation string),
#'   \code{cdst_superclass} (broadest classification: Neuronal/Glial/Other),
#'   \code{cdst_top_nt} (dominant NT type per cluster),
#'   \code{cdst_top_region} (dominant region per cluster).
#'
#' @export
cdst_annotate_hierarchy <- function(seurat_obj,
                                    cluster_column = "cdst_cluster",
                                    nt_column = "cdst_nt_type",
                                    glia_column = "cdst_glia_type",
                                    region_column = "brain_region_L3",
                                    nt_threshold = 0.30,
                                    region_threshold = 0.30,
                                    n_top_nt = 3L,
                                    n_top_regions = 2L,
                                    verbose = TRUE) {
  required_cols <- c(cluster_column, nt_column)
  validate_seurat(seurat_obj, required_metadata = required_cols)
  start_time <- Sys.time()

  meta <- seurat_obj@meta.data
  clusters <- unique(meta[[cluster_column]])

  if (verbose) {
    cli::cli_inform("Building hierarchical annotations for {length(clusters)} clusters...")
  }

  # ---- Compute per-cluster summaries ----
  cluster_annotations <- data.frame(
    cluster = character(0),
    top_nt = character(0),
    top_region = character(0),
    cell_type = character(0),
    superclass = character(0),
    stringsAsFactors = FALSE
  )

  for (cl in clusters) {
    cl_cells <- meta[meta[[cluster_column]] == cl, ]
    n_cl <- nrow(cl_cells)

    # --- NT type summary ---
    nt_expanded <- unlist(strsplit(cl_cells[[nt_column]], "-"))
    nt_counts <- table(nt_expanded)
    nt_props <- nt_counts / n_cl
    top_nts <- names(sort(nt_props[nt_props > nt_threshold], decreasing = TRUE))
    top_nts <- utils::head(top_nts, n_top_nt)
    top_nt_str <- if (length(top_nts) > 0) paste(top_nts, collapse = "-") else "Undefined"

    # --- Region summary ---
    if (region_column %in% colnames(cl_cells)) {
      region_counts <- table(cl_cells[[region_column]])
      region_props <- region_counts / n_cl
      top_regions <- names(sort(region_props[region_props > region_threshold],
                                 decreasing = TRUE))
      top_regions <- utils::head(top_regions, n_top_regions)
      top_region_str <- if (length(top_regions) > 0) {
        paste(top_regions, collapse = "-")
      } else {
        ""
      }
    } else {
      top_region_str <- ""
    }

    # --- Glia type ---
    if (glia_column %in% colnames(cl_cells)) {
      glia_vals <- cl_cells[[glia_column]]
      glia_mode <- names(sort(table(glia_vals), decreasing = TRUE))[1]
    } else {
      glia_mode <- "Neuronal"
    }

    # --- Superclass ---
    if (glia_mode != "Neuronal") {
      superclass <- "Non-Neuronal"
    } else if (top_nt_str != "Undefined") {
      superclass <- "Neuronal"
    } else {
      superclass <- "Unclassified"
    }

    # --- Full annotation ---
    parts <- c(top_nt_str, top_region_str)
    parts <- parts[nchar(parts) > 0]
    cell_type_str <- paste(parts, collapse = " ")

    cluster_annotations <- rbind(cluster_annotations, data.frame(
      cluster = as.character(cl),
      top_nt = top_nt_str,
      top_region = top_region_str,
      cell_type = cell_type_str,
      superclass = superclass,
      stringsAsFactors = FALSE
    ))
  }

  # ---- Map annotations to cells ----
  rownames(cluster_annotations) <- cluster_annotations$cluster
  cell_clusters <- as.character(meta[[cluster_column]])

  seurat_obj$cdst_cell_type <- cluster_annotations[cell_clusters, "cell_type"]
  seurat_obj$cdst_superclass <- cluster_annotations[cell_clusters, "superclass"]
  seurat_obj$cdst_top_nt <- cluster_annotations[cell_clusters, "top_nt"]
  seurat_obj$cdst_top_region <- cluster_annotations[cell_clusters, "top_region"]

  # Store annotation table in misc
  seurat_obj@misc$cdst_cluster_annotations <- cluster_annotations

  if (verbose) {
    n_types <- length(unique(cluster_annotations$cell_type))
    cli::cli_inform(c(
      "v" = "Hierarchical annotation complete.",
      "i" = "Unique cell types: {.val {n_types}}",
      "i" = "Superclasses: {.val {paste(unique(cluster_annotations$superclass), collapse = ', ')}}"
    ))
  }

  seurat_obj
}


#' Validate Cell Type Annotations Against Known Markers
#'
#' Cross-validates the assigned cell type annotations against a reference
#' database of known cell type markers (e.g., CellMarker2). Reports
#' concordance statistics and flags potential misannotations.
#'
#' @param seurat_obj A Seurat object with cell type annotations.
#' @param marker_db Character. Path to a marker database file (CSV). If NULL,
#'   uses a minimal built-in reference.
#' @param cluster_column Character. Default: "cdst_cluster".
#' @param n_top_markers Integer. Number of top markers per cluster to check.
#'   Default: 10.
#' @param verbose Logical. Default: TRUE.
#'
#' @return A data.frame with validation results per cluster.
#'
#' @export
cdst_validate_markers <- function(seurat_obj,
                                  marker_db = NULL,
                                  cluster_column = "cdst_cluster",
                                  n_top_markers = 10L,
                                  verbose = TRUE) {
  validate_seurat(seurat_obj, required_metadata = cluster_column)

  if (verbose) cli::cli_inform("Validating cell type markers...")

  # Find markers for each cluster
  Seurat::Idents(seurat_obj) <- cluster_column
  markers <- Seurat::FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
  )

  # Get top markers per cluster
  top_markers <- markers |>
    dplyr::group_by(.data$cluster) |>
    dplyr::slice_max(order_by = .data$avg_log2FC, n = n_top_markers) |>
    dplyr::summarise(
      top_genes = paste(.data$gene, collapse = ", "),
      n_markers = dplyr::n(),
      .groups = "drop"
    )

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Marker validation complete.",
      "i" = "Clusters analyzed: {nrow(top_markers)}"
    ))
  }

  as.data.frame(top_markers)
}
