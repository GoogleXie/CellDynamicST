# ============================================================================
# Shared Test Fixtures for CellDynamicST
# ============================================================================
# These fixtures are automatically loaded by testthat before any test file.

#' Create a minimal Seurat object for testing
#' @return A Seurat object with 200 cells and 50 genes
create_test_seurat <- function(n_cells = 200, n_genes = 50, n_groups = 4) {
  set.seed(42)

  # Generate count matrix
  counts <- matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(
      paste0("Gene", seq_len(n_genes)),
      paste0("Cell", seq_len(n_cells))
    )
  )

  # Add some known marker genes
  marker_genes <- c("Slc17a7", "Gad1", "Aqp4", "Cx3cr1", "Mbp")
  available_slots <- min(length(marker_genes), n_genes)
  rownames(counts)[seq_len(available_slots)] <- marker_genes[seq_len(available_slots)]

  # Create Seurat object
  obj <- Seurat::CreateSeuratObject(counts = counts, project = "test")

  # Add metadata
  group_names <- c("WT SAL", "WT MOR", "KO SAL", "KO MOR")[seq_len(n_groups)]
  obj@meta.data$experiment_group <- sample(group_names, n_cells, replace = TRUE)
  obj@meta.data$genotype <- ifelse(
    grepl("WT", obj@meta.data$experiment_group), "WT", "KO"
  )
  obj@meta.data$treatment <- ifelse(
    grepl("SAL", obj@meta.data$experiment_group), "SAL", "MOR"
  )
  obj@meta.data$sample_id <- paste0("Sample_", sample(1:4, n_cells, replace = TRUE))
  obj@meta.data$x_coord <- runif(n_cells, 0, 10000)
  obj@meta.data$y_coord <- runif(n_cells, 0, 10000)
  obj@meta.data$fov_id <- paste0("FOV_", sample(1:10, n_cells, replace = TRUE))
  obj@meta.data$slide_id <- sample(c("Slide_1", "Slide_2"), n_cells, replace = TRUE)
  obj@meta.data$brain_region_L3 <- sample(
    c("Isocortex", "Hippocampus", "Thalamus", "Hypothalamus"),
    n_cells, replace = TRUE
  )
  obj@meta.data$brain_region_L7 <- sample(
    c("VISp", "CA1", "VPM", "LHA", "SSp", "DG", "ACA"),
    n_cells, replace = TRUE
  )

  # Normalize
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)

  obj
}


#' Create a test CdstConfig object
#' @return A CdstConfig object with 4-group design
create_test_config <- function() {
  new("CdstConfig",
    project_name    = "Test_Project",
    species         = "mouse",
    groups          = list(
      list(name = "WT SAL", genotype = "WT", treatment = "SAL", is_reference = TRUE),
      list(name = "WT MOR", genotype = "WT", treatment = "MOR", is_reference = FALSE),
      list(name = "KO SAL", genotype = "KO", treatment = "SAL", is_reference = FALSE),
      list(name = "KO MOR", genotype = "KO", treatment = "MOR", is_reference = FALSE)
    ),
    comparisons     = list(
      list(name = "treatment_WT", group1 = "WT SAL", group2 = "WT MOR"),
      list(name = "treatment_KO", group1 = "KO SAL", group2 = "KO MOR"),
      list(name = "genotype_baseline", group1 = "WT SAL", group2 = "KO SAL"),
      list(name = "genotype_treated", group1 = "WT MOR", group2 = "KO MOR")
    ),
    reference_group = "WT SAL",
    qc              = list(
      min_count_percentile = 0.05,
      min_feature_percentile = 0.05,
      min_cells_per_fov = 500,
      min_transcripts_per_cell = 100,
      max_neg_probe_per_cell = 1.0
    ),
    clustering      = list(
      n_pcs = 100L, bias_threshold = 0.7, annoy_trees = 100L,
      k_neighbors = 50L, resolution = 20, rf_ntree = 200L
    ),
    wgcna           = list(
      soft_power = NULL, min_module_size = 20L, deep_split = 3L,
      merge_cut_height = 0.10, min_cells_per_region = 100L
    ),
    atlas           = list(
      name = "Allen_CCFv3", plane = "sagittal", z_plane = 510,
      region_columns = list(L3 = "brain_region_L3", L5 = "brain_region_L5",
                            L7 = "brain_region_L7")
    ),
    output_dir      = tempdir()
  )
}
