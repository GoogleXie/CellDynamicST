# =============================================================================
# CellDynamicST: Sample Data Generator
# =============================================================================
# Generates realistic synthetic spatial transcriptomic data for testing
# all visualization and analysis functions. Mimics a 2x2 factorial design
# (genotype x treatment) with brain region annotations.
#
# Usage:
#   source("inst/scripts/generate_sample_data.R")
#   sample_data <- generate_cdst_sample_data()
# =============================================================================

#' Generate Sample CellDynamicST Data
#'
#' Creates a list of data structures that mimic real CellDynamicST pipeline
#' outputs, suitable for testing all visualization functions.
#'
#' @param n_cells Number of cells to simulate. Default 5000.
#' @param seed Random seed. Default 42.
#' @return A list with: metadata, deg_results, enrichment_results,
#'   wgcna_nodes, wgcna_edges, wgcna_modules, hub_genes, hub_edges.
#' @export
generate_cdst_sample_data <- function(n_cells = 5000, seed = 42) {
  set.seed(seed)

  # ---- Brain regions (Allen Brain Atlas abbreviations) ----
  regions <- c("CTXpl", "CTXsp", "STR", "PAL", "TH", "HY",
               "MB", "P", "MY", "CB", "HPF", "OLF")
  meso_structures <- c(
    "CTXpl" = "Cortex", "CTXsp" = "Cortex", "STR" = "Striatum",
    "PAL" = "Pallidum", "TH" = "Thalamus", "HY" = "Hypothalamus",
    "MB" = "Midbrain", "P" = "Pons", "MY" = "Medulla", "CB" = "Cerebellum",
    "HPF" = "Hippocampus", "OLF" = "Olfactory"
  )
  parent_structures <- c(
    "Cortex" = "Cerebrum", "Striatum" = "Cerebrum",
    "Pallidum" = "Cerebrum", "Hippocampus" = "Cerebrum",
    "Olfactory" = "Cerebrum", "Thalamus" = "Interbrain",
    "Hypothalamus" = "Interbrain", "Midbrain" = "Brainstem",
    "Pons" = "Brainstem", "Medulla" = "Brainstem",
    "Cerebellum" = "Cerebellum"
  )

  # ---- Cell types ----
  major_types <- c("Neuron", "Astrocyte", "Oligodendrocyte", "Microglia",
                   "OPC", "Endothelial", "Pericyte")
  nt_types <- c("Glutamatergic", "GABAergic", "Dopaminergic",
                "Serotonergic", "Cholinergic")
  glia_subtypes <- c("Protoplasmic", "Fibrous", "Reactive",
                     "Mature_OL", "Newly_formed_OL", "Homeostatic_MG",
                     "Activated_MG")

  # ---- Experimental groups (2x2 factorial) ----
  genotypes <- c("AA", "GG")
  treatments <- c("SAL", "MOR")
  groups <- paste(rep(genotypes, each = 2), rep(treatments, 2), sep = "_")

  # ---- Generate metadata ----
  cell_ids <- paste0("cell_", seq_len(n_cells))

  # Assign groups with slight imbalance
  group_probs <- c(0.27, 0.23, 0.26, 0.24)
  experiment_group <- sample(groups, n_cells, replace = TRUE, prob = group_probs)

  # Assign regions with realistic distribution
  region_probs <- c(0.15, 0.10, 0.12, 0.05, 0.10, 0.08,
                    0.08, 0.06, 0.05, 0.07, 0.09, 0.05)
  brain_region <- sample(regions, n_cells, replace = TRUE, prob = region_probs)

  # Assign cell types with region-dependent probabilities
  cell_type <- character(n_cells)
  nt_class <- character(n_cells)
  glia_subtype <- character(n_cells)

  for (i in seq_len(n_cells)) {
    reg <- brain_region[i]
    # Neuron-heavy in cortex, glia-heavy in white matter
    if (reg %in% c("CTXpl", "CTXsp", "HPF")) {
      type_probs <- c(0.55, 0.12, 0.10, 0.08, 0.05, 0.05, 0.05)
    } else if (reg %in% c("STR", "TH")) {
      type_probs <- c(0.45, 0.15, 0.12, 0.10, 0.06, 0.06, 0.06)
    } else {
      type_probs <- c(0.35, 0.18, 0.15, 0.12, 0.08, 0.06, 0.06)
    }
    cell_type[i] <- sample(major_types, 1, prob = type_probs)

    # NT classification for neurons
    if (cell_type[i] == "Neuron") {
      if (reg %in% c("CTXpl", "CTXsp", "HPF")) {
        nt_probs <- c(0.60, 0.30, 0.02, 0.02, 0.06)
      } else if (reg %in% c("STR", "PAL")) {
        nt_probs <- c(0.20, 0.60, 0.10, 0.02, 0.08)
      } else {
        nt_probs <- c(0.30, 0.40, 0.10, 0.10, 0.10)
      }
      nt_class[i] <- sample(nt_types, 1, prob = nt_probs)
    } else {
      nt_class[i] <- NA
    }

    # Glia subtype
    if (cell_type[i] == "Astrocyte") {
      glia_subtype[i] <- sample(c("Protoplasmic", "Fibrous", "Reactive"),
                                 1, prob = c(0.5, 0.3, 0.2))
    } else if (cell_type[i] == "Oligodendrocyte") {
      glia_subtype[i] <- sample(c("Mature_OL", "Newly_formed_OL"),
                                 1, prob = c(0.7, 0.3))
    } else if (cell_type[i] == "Microglia") {
      glia_subtype[i] <- sample(c("Homeostatic_MG", "Activated_MG"),
                                 1, prob = c(0.6, 0.4))
    } else {
      glia_subtype[i] <- NA
    }
  }

  # Build combined annotation
  cell_type_annotation <- ifelse(
    !is.na(nt_class), paste0(nt_class, "_Neuron"),
    ifelse(!is.na(glia_subtype), glia_subtype, cell_type)
  )

  # Spatial coordinates (AP and DV)
  ap_location <- runif(n_cells, 0, 13.2)  # Bregma range
  dv_location <- runif(n_cells, 0, 8.0)

  # QC metrics
  nCount_RNA <- rnbinom(n_cells, mu = 3000, size = 5)
  nFeature_RNA <- rnbinom(n_cells, mu = 800, size = 3)

  # Genotype and treatment
  genotype <- sub("_.*", "", experiment_group)
  treatment <- sub(".*_", "", experiment_group)

  # Meso structure
  meso_structure <- meso_structures[brain_region]

  metadata <- data.frame(
    cell_id = cell_ids,
    experiment_group = experiment_group,
    genotype = genotype,
    treatment = treatment,
    brain_region = brain_region,
    meso_structure = unname(meso_structure),
    parent_structure = unname(parent_structures[meso_structure]),
    cell_type = cell_type,
    nt_class = nt_class,
    glia_subtype = glia_subtype,
    cell_type_annotation = cell_type_annotation,
    AP_location = ap_location,
    DV_location = dv_location,
    nCount_RNA = nCount_RNA,
    nFeature_RNA = nFeature_RNA,
    stringsAsFactors = FALSE,
    row.names = cell_ids
  )

  # ---- Generate tSNE/UMAP coordinates ----
  # Cluster-aware embedding
  n_clusters <- length(unique(cell_type_annotation))
  cluster_centers <- matrix(rnorm(n_clusters * 2, sd = 15),
                            ncol = 2)
  rownames(cluster_centers) <- unique(cell_type_annotation)

  tsne_coords <- matrix(0, nrow = n_cells, ncol = 2)
  for (i in seq_len(n_cells)) {
    ct <- cell_type_annotation[i]
    center <- cluster_centers[ct, ]
    tsne_coords[i, ] <- center + rnorm(2, sd = 3)
  }
  metadata$tSNE_1 <- tsne_coords[, 1]
  metadata$tSNE_2 <- tsne_coords[, 2]
  metadata$UMAP_1 <- tsne_coords[, 1] * 0.8 + rnorm(n_cells, sd = 1)
  metadata$UMAP_2 <- tsne_coords[, 2] * 0.8 + rnorm(n_cells, sd = 1)

  # ---- Generate DEG results ----
  n_genes <- 500
  gene_names <- paste0("Gene", seq_len(n_genes))
  deg_list <- list()
  contrasts <- c("MOR_vs_SAL_AA", "MOR_vs_SAL_GG", "AA_vs_GG_SAL", "AA_vs_GG_MOR")

  for (contrast in contrasts) {
    for (reg in regions) {
      log2fc <- rnorm(n_genes, mean = 0, sd = 1.2)
      pvals <- 10^(-abs(log2fc) * runif(n_genes, 1, 5))
      pvals_adj <- p.adjust(pvals, method = "BH")

      deg_list[[length(deg_list) + 1]] <- data.frame(
        gene = gene_names,
        avg_log2FC = log2fc,
        p_val = pvals,
        p_val_adj = pvals_adj,
        pct.1 = runif(n_genes, 0.1, 0.9),
        pct.2 = runif(n_genes, 0.1, 0.9),
        contrast = contrast,
        region = reg,
        stringsAsFactors = FALSE
      )
    }
  }
  deg_results <- do.call(rbind, deg_list)

  # ---- Generate enrichment results ----
  go_terms <- c(
    "GO:0007268" = "chemical synaptic transmission",
    "GO:0007399" = "nervous system development",
    "GO:0006954" = "inflammatory response",
    "GO:0006915" = "apoptotic process",
    "GO:0007420" = "brain development",
    "GO:0050804" = "modulation of chemical synaptic transmission",
    "GO:0048699" = "generation of neurons",
    "GO:0030182" = "neuron differentiation",
    "GO:0042391" = "regulation of membrane potential",
    "GO:0007411" = "axon guidance"
  )

  enrich_list <- list()
  for (contrast in contrasts[1:2]) {
    for (reg in regions[1:6]) {
      for (go_id in names(go_terms)) {
        enrich_list[[length(enrich_list) + 1]] <- data.frame(
          go_id = go_id,
          go_name = go_terms[[go_id]],
          enrichment_score = rnorm(1, sd = 2),
          p_value = runif(1, 0.001, 0.1),
          n_genes = sample(10:100, 1),
          contrast = contrast,
          region = reg,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  enrichment_results <- do.call(rbind, enrich_list)

  # ---- Generate WGCNA data ----
  module_colors <- c("blue", "brown", "green", "red", "turquoise",
                     "yellow", "black", "pink", "magenta", "purple")

  # Nodes (one per region-module combination, simplified)
  wgcna_nodes <- data.frame(
    region = rep(regions, each = 3),
    module_color = sample(module_colors, length(regions) * 3, replace = TRUE),
    n_genes = sample(20:200, length(regions) * 3, replace = TRUE),
    stringsAsFactors = FALSE
  )

  # Edges (region-to-region correlations)
  region_pairs <- t(combn(regions, 2))
  wgcna_edges <- data.frame(
    source = region_pairs[, 1],
    target = region_pairs[, 2],
    correlation = rnorm(nrow(region_pairs), sd = 0.5),
    module_color = sample(module_colors, nrow(region_pairs), replace = TRUE),
    stringsAsFactors = FALSE
  )

  # Module enrichment for circos
  wgcna_modules <- data.frame(
    module_color = module_colors,
    n_genes = sample(50:500, length(module_colors)),
    NES_MOR_vs_SAL_AA = rnorm(length(module_colors), sd = 1.5),
    NES_MOR_vs_SAL_GG = rnorm(length(module_colors), sd = 1.5),
    NES_AA_vs_GG_SAL = rnorm(length(module_colors), sd = 1.5),
    NES_AA_vs_GG_MOR = rnorm(length(module_colors), sd = 1.5),
    stringsAsFactors = FALSE
  )

  # Hub genes
  hub_genes <- data.frame(
    gene = paste0("Hub_", seq_len(80)),
    module_color = sample(module_colors, 80, replace = TRUE),
    kME = runif(80, 0.5, 1.0),
    stringsAsFactors = FALSE
  )

  # Hub gene edges
  hub_pairs <- t(combn(head(hub_genes$gene, 30), 2))
  hub_edges <- data.frame(
    gene1 = hub_pairs[, 1],
    gene2 = hub_pairs[, 2],
    weight = runif(nrow(hub_pairs), 0.1, 1.0),
    stringsAsFactors = FALSE
  )

  # ---- Generate trait correlation matrices ----
  n_modules <- length(module_colors)
  cor_treatment <- matrix(rnorm(n_modules, sd = 0.4), ncol = 1)
  rownames(cor_treatment) <- paste0("ME", module_colors)
  colnames(cor_treatment) <- "treatment"
  p_treatment <- matrix(10^(-abs(cor_treatment) * runif(n_modules, 2, 8)),
                        ncol = 1)
  rownames(p_treatment) <- rownames(cor_treatment)
  colnames(p_treatment) <- "treatment"

  cor_genotype <- matrix(rnorm(n_modules, sd = 0.3), ncol = 1)
  rownames(cor_genotype) <- rownames(cor_treatment)
  colnames(cor_genotype) <- "genotype"
  p_genotype <- matrix(10^(-abs(cor_genotype) * runif(n_modules, 1, 6)),
                       ncol = 1)
  rownames(p_genotype) <- rownames(cor_treatment)
  colnames(p_genotype) <- "genotype"

  trait_cor <- list(
    cor_treatment = cor_treatment,
    p_treatment = p_treatment,
    cor_genotype = cor_genotype,
    p_genotype = p_genotype
  )

  # ---- Return everything ----
  list(
    metadata = metadata,
    deg_results = deg_results,
    enrichment_results = enrichment_results,
    wgcna_nodes = wgcna_nodes,
    wgcna_edges = wgcna_edges,
    wgcna_modules = wgcna_modules,
    hub_genes = hub_genes,
    hub_edges = hub_edges,
    trait_cor = trait_cor,
    regions = regions,
    meso_structures = meso_structures,
    parent_structures = parent_structures
  )
}
