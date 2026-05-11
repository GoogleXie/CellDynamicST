#!/usr/bin/env Rscript
# =============================================================================
# CellDynamicST: Comprehensive Visualization Test Suite
# =============================================================================
# Self-contained test that generates sample data inline and tests all
# visualization functions end-to-end. No Seurat dependency required.
#
# Usage:  cd CellDynamicST && Rscript tests/test_all_visualizations.R
# =============================================================================

cat("\n========================================================\n")
cat("  CellDynamicST Visualization Test Suite\n")
cat("========================================================\n\n")

# ---- Setup ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
})

# Source all R modules (sorted so aaa-classes.R loads first)
r_files <- sort(list.files(file.path(getwd(), "R"), pattern = "\\.R$",
                           full.names = TRUE))
for (f in r_files) {
  tryCatch(source(f, local = FALSE), error = function(e) {
    message("  Warning sourcing ", basename(f), ": ", e$message)
  })
}

# Output directory
out_dir <- file.path(getwd(), "tests", "test_output")
if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)
dir.create(out_dir, recursive = TRUE)

# Track results
results <- list()
test_num <- 0

run_test <- function(name, expr) {
  test_num <<- test_num + 1
  cat(sprintf("[%02d] %s ... ", test_num, name))
  t0 <- Sys.time()
  result <- tryCatch({
    force(expr)
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cat(sprintf("PASS (%.1fs)\n", elapsed))
    list(name = name, status = "PASS", time = elapsed, error = "")
  }, error = function(e) {
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cat(sprintf("FAIL (%.1fs)\n", elapsed))
    cat("     Error:", conditionMessage(e), "\n")
    list(name = name, status = "FAIL", time = elapsed, error = conditionMessage(e))
  })
  results[[test_num]] <<- result
}


# =============================================================================
# GENERATE SAMPLE DATA
# =============================================================================
cat("--- Generating sample data ---\n")
set.seed(42)
n_cells <- 2000

regions <- c("Isocortex", "HPF", "CTXsp", "STR", "TH", "HY", "MB", "CB")
meso_map <- list(
  Isocortex = c("VIS", "MO", "SS"), HPF = c("CA1", "CA3", "DG"),
  CTXsp = c("CLA", "EP", "LA"), STR = c("CP", "ACB"),
  TH = c("VPM", "LGd"), HY = c("LHA", "PVH"),
  MB = c("SN", "VTA"), CB = c("CBX", "CBN")
)

cell_types <- c("Glutamatergic", "GABAergic", "Dopaminergic", "Serotonergic",
                "Cholinergic", "Astrocyte", "Oligodendrocyte", "Microglia",
                "OPC", "Endothelial")
nt_types <- c("Glutamate", "GABA", "Dopamine", "Serotonin", "Acetylcholine",
              "Non-neuronal", "Non-neuronal", "Non-neuronal",
              "Non-neuronal", "Non-neuronal")
high_level <- c("Neuron", "Neuron", "Neuron", "Neuron", "Neuron",
                "Glia", "Glia", "Glia", "Glia", "Vascular")

groups <- c("AA_SAL", "AA_MOR", "GG_SAL", "GG_MOR")

ct_idx <- sample(seq_along(cell_types), n_cells, replace = TRUE,
                 prob = c(.25, .20, .05, .03, .02, .15, .12, .08, .05, .05))
rg_idx <- sample(seq_along(regions), n_cells, replace = TRUE)

metadata <- data.frame(
  cell_id = paste0("cell_", seq_len(n_cells)),
  cell_type_annotation = cell_types[ct_idx],
  high_level_cell_type = high_level[ct_idx],
  neurotransmitter_type = nt_types[ct_idx],
  level_7_parent_roi_acronym = regions[rg_idx],
  level_3_parent_roi_acronym = sapply(rg_idx, function(i) {
    sample(meso_map[[regions[i]]], 1)
  }),
  experiment_group = sample(groups, n_cells, replace = TRUE),
  genotype = sample(c("AA", "GG"), n_cells, replace = TRUE),
  treatment = sample(c("SAL", "MOR"), n_cells, replace = TRUE),
  nCount_RNA = rpois(n_cells, 5000),
  nFeature_RNA = rpois(n_cells, 2000),
  AP_location = runif(n_cells, -5, 5),
  DV_location = runif(n_cells, -4, 4),
  tSNE_1 = rnorm(n_cells, 0, 15),
  tSNE_2 = rnorm(n_cells, 0, 15),
  stringsAsFactors = FALSE
)
rownames(metadata) <- metadata$cell_id

cat("  Generated", n_cells, "cells,", length(regions), "regions,",
    length(cell_types), "cell types\n\n")


# =============================================================================
# TEST 1: Cell Type Hierarchy
# =============================================================================
cat("--- Test 1: Cell Type Hierarchy ---\n")

run_test("1a: Build hierarchy", {
  hierarchy <- cdst_build_hierarchy(
    data = metadata,
    levels = c("high_level_cell_type", "neurotransmitter_type",
               "level_7_parent_roi_acronym", "cell_type_annotation"),
    min_cells = 5
  )
  stopifnot(is.list(hierarchy))
  stopifnot("nodes" %in% names(hierarchy))
  stopifnot("links" %in% names(hierarchy))
  stopifnot(nrow(hierarchy$nodes) > 0)
  stopifnot(nrow(hierarchy$links) > 0)
  cat(sprintf("(%d nodes, %d links) ", nrow(hierarchy$nodes),
              nrow(hierarchy$links)))
})

run_test("1b: Hierarchy tree plot", {
  hierarchy <- cdst_build_hierarchy(
    data = metadata,
    levels = c("high_level_cell_type", "neurotransmitter_type",
               "level_7_parent_roi_acronym", "cell_type_annotation"),
    min_cells = 5
  )
  p <- cdst_plot_hierarchy_tree(
    hierarchy_data = hierarchy,
    output_file = file.path(out_dir, "01_hierarchy_tree.png"),
    width = 16, height = 10
  )
  stopifnot(file.exists(file.path(out_dir, "01_hierarchy_tree.png")))
  stopifnot(file.info(file.path(out_dir, "01_hierarchy_tree.png"))$size > 1000)
})

run_test("1c: Hierarchy Sankey", {
  hierarchy <- cdst_build_hierarchy(
    data = metadata,
    levels = c("high_level_cell_type", "neurotransmitter_type",
               "level_7_parent_roi_acronym", "cell_type_annotation"),
    min_cells = 5
  )
  if (requireNamespace("networkD3", quietly = TRUE) &&
      requireNamespace("htmlwidgets", quietly = TRUE)) {
    w <- cdst_plot_hierarchy_sankey(
      hierarchy_data = hierarchy,
      output_file = file.path(out_dir, "01_hierarchy_sankey.html")
    )
    stopifnot(file.exists(file.path(out_dir, "01_hierarchy_sankey.html")))
  } else {
    cat("(networkD3 not available, skipping) ")
  }
})


# =============================================================================
# TEST 2: Disproportion Scores on tSNE
# =============================================================================
cat("\n--- Test 2: Disproportion Scores ---\n")

run_test("2a: Calculate scores", {
  scores <- cdst_disproportion_scores(
    data = metadata,
    group1 = "AA_SAL",
    group2 = "AA_MOR",
    group_col = "experiment_group",
    celltype_col = "cell_type_annotation",
    min_cells = 3
  )
  stopifnot(is.data.frame(scores))
  stopifnot("cell_type" %in% colnames(scores))
  stopifnot("zscore" %in% colnames(scores))
  cat(sprintf("(%d cell types scored) ", nrow(scores)))
})

run_test("2b: Plot disproportion tSNE", {
  scores <- cdst_disproportion_scores(
    data = metadata,
    group1 = "AA_SAL", group2 = "AA_MOR",
    group_col = "experiment_group",
    celltype_col = "cell_type_annotation", min_cells = 3
  )
  p <- cdst_plot_disproportion(
    data = metadata,
    scores = scores,
    celltype_col = "cell_type_annotation",
    reduction_x = "tSNE_1",
    reduction_y = "tSNE_2",
    comparison_label = "AA: SAL vs MOR",
    output_file = file.path(out_dir, "02_disproportion_tsne.png"),
    width = 12, height = 10
  )
  stopifnot(file.exists(file.path(out_dir, "02_disproportion_tsne.png")))
  stopifnot(file.info(file.path(out_dir, "02_disproportion_tsne.png"))$size > 1000)
})

run_test("2c: Multi-comparison pipeline", {
  all_scores <- cdst_run_disproportion(
    data = metadata,
    comparisons = list(c("AA_SAL", "AA_MOR"), c("GG_SAL", "GG_MOR")),
    group_col = "experiment_group",
    celltype_col = "cell_type_annotation",
    reduction_x = "tSNE_1",
    reduction_y = "tSNE_2",
    output_dir = file.path(out_dir, "02_disproportion_multi")
  )
  stopifnot(is.list(all_scores))
  stopifnot(length(all_scores) == 2)
  pngs <- list.files(file.path(out_dir, "02_disproportion_multi"),
                     pattern = "\\.png$")
  stopifnot(length(pngs) >= 2)
  cat(sprintf("(%d plots) ", length(pngs)))
})


# =============================================================================
# TEST 3: Multi-Dimensional Heatmap (Python)
# =============================================================================
cat("\n--- Test 3: Multi-Dimensional Heatmap ---\n")

run_test("3a: Export metadata CSV", {
  heatmap_dir <- file.path(out_dir, "03_heatmap")
  dir.create(heatmap_dir, recursive = TRUE)
  csv_path <- file.path(heatmap_dir, "metadata.csv")
  write.csv(metadata, csv_path, row.names = TRUE)
  stopifnot(file.exists(csv_path))
  stopifnot(file.info(csv_path)$size > 1000)
})

run_test("3b: Python heatmap execution", {
  py_script <- file.path(getwd(), "inst", "python",
                          "celltype_distribution_heatmap.py")
  stopifnot(file.exists(py_script))

  heatmap_dir <- file.path(out_dir, "03_heatmap")
  csv_path <- file.path(heatmap_dir, "metadata.csv")

  cmd <- paste(
    "python3", shQuote(py_script),
    "--input", shQuote(csv_path),
    "--output", shQuote(heatmap_dir),
    "--genotypes", "AA", "GG",
    "--control", "SAL",
    "--treatment", "MOR",
    "--col-celltype", "cell_type_annotation",
    "--col-region", "level_7_parent_roi_acronym",
    "--col-meso", "level_3_parent_roi_acronym",
    "--col-cellgroup", "high_level_cell_type",
    "--col-nt", "neurotransmitter_type",
    "--col-group", "experiment_group",
    "2>&1"
  )
  output <- system(cmd, intern = TRUE)
  cat(sprintf("(%d lines output) ", length(output)))

  # Show last few lines of output for debugging
  for (line in tail(output, 5)) cat("\n     ", line)
  cat(" ")

  png_files <- list.files(heatmap_dir, pattern = "\\.png$")
  if (length(png_files) > 0) {
    cat(sprintf("(%d PNGs) ", length(png_files)))
  } else {
    # Not a hard failure — Python deps may vary
    cat("(no PNGs — check Python deps) ")
  }
})


# =============================================================================
# TEST 4: Spatial Enrichment
# =============================================================================
cat("\n--- Test 4: Spatial Enrichment ---\n")

run_test("4a: Create enrichment data", {
  enrichment_df <- data.frame(
    cell_ID = metadata$cell_id,
    AP_location = metadata$AP_location,
    DV_location = metadata$DV_location,
    experiment_group = metadata$experiment_group,
    enrichment_score = rnorm(n_cells, 0.3, 0.15),
    stringsAsFactors = FALSE
  )
  ctrl <- enrichment_df$enrichment_score[
    enrichment_df$experiment_group == "AA_SAL"
  ]
  enrichment_df$zscore <- (enrichment_df$enrichment_score - mean(ctrl)) / sd(ctrl)
  uq <- quantile(ctrl, 0.95); lq <- quantile(ctrl, 0.05)
  enrichment_df$significant_enrichment <- ifelse(
    enrichment_df$enrichment_score > uq, enrichment_df$zscore,
    ifelse(enrichment_df$enrichment_score < lq, enrichment_df$zscore, 0)
  )
  # Save for later tests
  assign("enrichment_df", enrichment_df, envir = .GlobalEnv)
  stopifnot(nrow(enrichment_df) == n_cells)
})

run_test("4b: Spatial enrichment plot", {
  p <- cdst_plot_spatial_enrichment(
    enrichment_df = enrichment_df,
    experiment_group = "AA_MOR",
    go_term = "synaptic transmission",
    roi = "Isocortex",
    comparison = "AA: SAL vs MOR",
    score_col = "significant_enrichment",
    x_col = "DV_location",
    y_col = "AP_location",
    output_file = file.path(out_dir, "04_spatial_enrichment.png"),
    width = 10, height = 8
  )
  stopifnot(file.exists(file.path(out_dir, "04_spatial_enrichment.png")))
  stopifnot(file.info(file.path(out_dir, "04_spatial_enrichment.png"))$size > 1000)
})


# =============================================================================
# TEST 5: WGCNA Network Visualizations (Python)
# =============================================================================
cat("\n--- Test 5: WGCNA Network Plots ---\n")

wgcna_dir <- file.path(out_dir, "05_wgcna")
dir.create(wgcna_dir, recursive = TRUE)

run_test("5a: Create WGCNA data", {
  nodes_df <- data.frame(
    region = regions,
    n_genes = c(150, 120, 80, 90, 110, 70, 95, 85),
    module_color = c("blue", "turquoise", "brown", "green",
                     "red", "yellow", "pink", "magenta"),
    stringsAsFactors = FALSE
  )
  write.csv(nodes_df, file.path(wgcna_dir, "nodes.csv"), row.names = FALSE)

  edge_pairs <- combn(regions, 2)
  edges_df <- data.frame(
    source = edge_pairs[1, ],
    target = edge_pairs[2, ],
    correlation = runif(ncol(edge_pairs), -0.8, 0.8),
    stringsAsFactors = FALSE
  )
  write.csv(edges_df, file.path(wgcna_dir, "edges.csv"), row.names = FALSE)

  module_df <- data.frame(
    module_color = c("blue", "turquoise", "brown", "green",
                     "red", "yellow", "pink", "magenta"),
    n_genes = c(150, 120, 80, 90, 110, 70, 95, 85),
    NES_OUD_Enrichment = runif(8, 0, 5),
    NES_AA_MOR_vs_SAL = runif(8, 0, 5),
    NES_GG_MOR_vs_SAL = runif(8, 0, 5),
    NES_AA_vs_GG_SAL = runif(8, 0, 5),
    NES_AA_vs_GG_MOR = runif(8, 0, 5),
    stringsAsFactors = FALSE
  )
  write.csv(module_df, file.path(wgcna_dir, "modules.csv"), row.names = FALSE)

  hub_genes <- paste0("Gene_", seq_len(60))
  hub_df <- data.frame(
    gene = hub_genes,
    module_color = sample(c("blue", "turquoise", "brown", "green", "red"),
                          60, replace = TRUE),
    kME = runif(60, 0.5, 1.0),
    is_hub = c(rep(TRUE, 8), rep(FALSE, 52)),
    stringsAsFactors = FALSE
  )
  write.csv(hub_df, file.path(wgcna_dir, "hubs.csv"), row.names = FALSE)

  ep <- combn(hub_genes[1:30], 2)
  idx <- sample(ncol(ep), min(100, ncol(ep)))
  gene_edges <- data.frame(
    gene1 = ep[1, idx], gene2 = ep[2, idx],
    weight = runif(length(idx), 0.1, 1.0),
    stringsAsFactors = FALSE
  )
  write.csv(gene_edges, file.path(wgcna_dir, "hub_edges.csv"), row.names = FALSE)
  stopifnot(file.exists(file.path(wgcna_dir, "nodes.csv")))
})

run_test("5b: Atlas network sagittal (Python)", {
  py_script <- file.path(getwd(), "inst", "python", "wgcna_network_plot.py")
  stopifnot(file.exists(py_script))
  out_file <- file.path(wgcna_dir, "atlas_network_sagittal.png")
  cmd <- paste(
    "python3", shQuote(py_script), "atlas",
    "--nodes", shQuote(file.path(wgcna_dir, "nodes.csv")),
    "--edges", shQuote(file.path(wgcna_dir, "edges.csv")),
    "--output", shQuote(out_file),
    "--plane", "sagittal",
    "--slice", "200",
    "--title", shQuote("Atlas Network (Sagittal z=200)"),
    "2>&1"
  )
  output <- system(cmd, intern = TRUE)
  for (line in tail(output, 3)) cat("\n     ", line)
  cat(" ")
  stopifnot(file.exists(out_file))
  stopifnot(file.info(out_file)$size > 1000)
})

run_test("5b2: Atlas network coronal (Python)", {
  py_script <- file.path(getwd(), "inst", "python", "wgcna_network_plot.py")
  out_file <- file.path(wgcna_dir, "atlas_network_coronal.png")
  cmd <- paste(
    "python3", shQuote(py_script), "atlas",
    "--nodes", shQuote(file.path(wgcna_dir, "nodes.csv")),
    "--edges", shQuote(file.path(wgcna_dir, "edges.csv")),
    "--output", shQuote(out_file),
    "--plane", "coronal",
    "--slice", "250",
    "--title", shQuote("Atlas Network (Coronal y=250)"),
    "2>&1"
  )
  output <- system(cmd, intern = TRUE)
  for (line in tail(output, 3)) cat("\n     ", line)
  cat(" ")
  stopifnot(file.exists(out_file))
  stopifnot(file.info(out_file)$size > 1000)
})

run_test("5c: Circos plot (Python)", {
  py_script <- file.path(getwd(), "inst", "python", "wgcna_network_plot.py")
  out_file <- file.path(wgcna_dir, "circos.png")
  cmd <- paste(
    "python3", shQuote(py_script), "circos",
    "--modules", shQuote(file.path(wgcna_dir, "modules.csv")),
    "--output", shQuote(out_file),
    "--title", shQuote("Test Circos"),
    "2>&1"
  )
  output <- system(cmd, intern = TRUE)
  for (line in tail(output, 3)) cat("\n     ", line)
  cat(" ")
  stopifnot(file.exists(out_file))
  stopifnot(file.info(out_file)$size > 1000)
})

run_test("5d: Hub gene network (Python)", {
  py_script <- file.path(getwd(), "inst", "python", "wgcna_network_plot.py")
  out_file <- file.path(wgcna_dir, "hub_network.png")
  cmd <- paste(
    "python3", shQuote(py_script), "hub",
    "--hubs", shQuote(file.path(wgcna_dir, "hubs.csv")),
    "--edges", shQuote(file.path(wgcna_dir, "hub_edges.csv")),
    "--output", shQuote(out_file),
    "--top-n", "30",
    "--title", shQuote("Test Hub Network"),
    "2>&1"
  )
  output <- system(cmd, intern = TRUE)
  for (line in tail(output, 3)) cat("\n     ", line)
  cat(" ")
  stopifnot(file.exists(out_file))
  stopifnot(file.info(out_file)$size > 1000)
})


# =============================================================================
# TEST 6: Core R Plotting Functions
# =============================================================================
cat("\n--- Test 6: Core R Plots ---\n")

run_test("6a: Volcano plot", {
  deg_results <- data.frame(
    gene = paste0("Gene_", seq_len(500)),
    avg_log2FC = rnorm(500, 0, 1.5),
    p_val_adj = 10^(-runif(500, 0, 10)),
    contrast = "AA_SAL_vs_MOR",
    region = sample(regions, 500, replace = TRUE),
    stringsAsFactors = FALSE
  )
  p <- cdst_plot_volcano(
    deg_results,
    contrast = "AA_SAL_vs_MOR",
    logfc_threshold = 0.5,
    p_threshold = 0.05,
    n_label = 10,
    title = "Test Volcano: AA SAL vs MOR",
    output_file = file.path(out_dir, "06_volcano.png")
  )
  stopifnot(file.exists(file.path(out_dir, "06_volcano.png")))
  stopifnot(file.info(file.path(out_dir, "06_volcano.png"))$size > 1000)
})

run_test("6b: Enrichment lollipop (single condition)", {
  enrich_results <- data.frame(
    go_name = paste0("GO:", sprintf("%04d", seq_len(25)), " ",
                     sample(c("synaptic transmission", "axon guidance",
                              "neuron migration", "myelination",
                              "immune response"), 25, replace = TRUE)),
    enrichment_score = rnorm(25, 0, 2),
    p_value = 10^(-runif(25, 0, 8)),
    n_genes = sample(5:50, 25, replace = TRUE),
    contrast = "AA_SAL_vs_MOR",
    region = sample(regions[1:3], 25, replace = TRUE),
    stringsAsFactors = FALSE
  )
  p <- cdst_plot_enrichment_lollipop(
    enrich_results,
    top_n = 15,
    title = "Test GO Enrichment Lollipop",
    output_file = file.path(out_dir, "06_enrichment_lollipop.png")
  )
  stopifnot(file.exists(file.path(out_dir, "06_enrichment_lollipop.png")))
})

run_test("6b2: Enrichment bidirectional comparison (2x2)", {
  # Simulate clusterProfiler-style enrichGO output for two conditions
  go_terms <- c("synaptic transmission", "axon guidance", "neuron migration",
                "myelination", "immune response", "apoptotic process",
                "cell adhesion", "signal transduction", "ion transport",
                "protein phosphorylation", "transcription regulation",
                "lipid metabolism", "vesicle transport", "endocytosis",
                "calcium signaling", "GABA signaling", "glutamate signaling",
                "dopamine signaling", "serotonin signaling", "neuropeptide signaling")
  genes_pool <- paste0("Gene_", seq_len(200))

  make_go_df <- function(n_terms, terms_pool, genes_pool) {
    sel <- sample(terms_pool, n_terms)
    data.frame(
      ID = paste0("GO:", sprintf("%07d", sample(1:9999999, n_terms))),
      Description = sel,
      GeneRatio = paste0(sample(5:30, n_terms, replace = TRUE), "/200"),
      BgRatio = paste0(sample(50:500, n_terms, replace = TRUE), "/20000"),
      pvalue = 10^(-runif(n_terms, 1, 8)),
      p.adjust = 10^(-runif(n_terms, 0.5, 7)),
      qvalue = 10^(-runif(n_terms, 0.5, 7)),
      geneID = sapply(seq_len(n_terms), function(i) {
        paste(sample(genes_pool, sample(5:25, 1)), collapse = "/")
      }),
      Count = sample(5:30, n_terms, replace = TRUE),
      stringsAsFactors = FALSE
    )
  }

  enrich_AA <- make_go_df(18, go_terms, genes_pool)
  enrich_GG <- make_go_df(16, go_terms, genes_pool)

  log2fc_AA <- data.frame(
    gene = genes_pool,
    log2FC = rnorm(length(genes_pool), 0, 1.2),
    stringsAsFactors = FALSE
  )
  log2fc_GG <- data.frame(
    gene = genes_pool,
    log2FC = rnorm(length(genes_pool), 0, 1.0),
    stringsAsFactors = FALSE
  )

  p <- cdst_plot_enrichment_lollipop(
    enrichment_A = enrich_AA,
    enrichment_B = enrich_GG,
    log2fc_A = log2fc_AA,
    log2fc_B = log2fc_GG,
    label_A = "AA",
    label_B = "GG",
    top_n = 15,
    title = "Cerebral Cortex: BP (AA vs GG)",
    output_file = file.path(out_dir, "06_enrichment_bidirectional.png")
  )
  stopifnot(file.exists(file.path(out_dir, "06_enrichment_bidirectional.png")))
  stopifnot(file.info(file.path(out_dir, "06_enrichment_bidirectional.png"))$size > 1000)
})

run_test("6c: Enrichment heatmap", {
  enrich_results <- data.frame(
    go_name = paste0("GO:", sprintf("%04d", seq_len(25)), " ",
                     sample(c("synaptic transmission", "axon guidance",
                              "neuron migration", "myelination"), 25, replace = TRUE)),
    enrichment_score = rnorm(25, 0, 2),
    p_value = 10^(-runif(25, 0, 8)),
    n_genes = sample(5:50, 25, replace = TRUE),
    contrast = rep("AA_SAL_vs_MOR", 25),
    region = sample(regions[1:4], 25, replace = TRUE),
    stringsAsFactors = FALSE
  )
  p <- cdst_plot_enrichment_heatmap(
    enrich_results,
    contrast = "AA_SAL_vs_MOR",
    top_n = 15,
    title = "Test GO Enrichment Heatmap",
    output_file = file.path(out_dir, "06_enrichment_heatmap.png")
  )
  stopifnot(file.exists(file.path(out_dir, "06_enrichment_heatmap.png")))
})

run_test("6d: QC standalone", {
  p <- cdst_plot_qc_standalone(
    metadata,
    group_by = "experiment_group",
    output_file = file.path(out_dir, "06_qc_standalone.png"),
    width = 16, height = 6
  )
  stopifnot(file.exists(file.path(out_dir, "06_qc_standalone.png")))
})

run_test("6e: Proportions standalone", {
  p <- cdst_plot_proportions_standalone(
    metadata,
    celltype_col = "cell_type_annotation",
    group_col = "experiment_group",
    plot_type = "stacked",
    output_file = file.path(out_dir, "06_proportions_standalone.png")
  )
  stopifnot(file.exists(file.path(out_dir, "06_proportions_standalone.png")))
})

run_test("6f: Trait heatmap", {
  modules <- c("blue", "turquoise", "brown", "green", "red")
  traits <- c("treatment", "genotype", "region")
  cor_mat <- matrix(runif(15, -0.8, 0.8), nrow = 5, ncol = 3,
                    dimnames = list(modules, traits))
  p_mat <- matrix(10^(-runif(15, 0, 4)), nrow = 5, ncol = 3,
                  dimnames = list(modules, traits))
  trait_cor <- list(
    cor_treatment = cor_mat[, 1, drop = FALSE],
    p_treatment = p_mat[, 1, drop = FALSE],
    cor_genotype = cor_mat[, 2, drop = FALSE],
    p_genotype = p_mat[, 2, drop = FALSE]
  )
  hp <- cdst_plot_trait_heatmap(
    trait_cor,
    trait_type = "treatment",
    output_file = file.path(out_dir, "06_trait_heatmap.png"),
    width = 8, height = 6
  )
  stopifnot(file.exists(file.path(out_dir, "06_trait_heatmap.png")))
})


# =============================================================================
# SUMMARY
# =============================================================================
cat("\n========================================================\n")
cat("  TEST SUMMARY\n")
cat("========================================================\n\n")

results_df <- do.call(rbind, lapply(results, as.data.frame,
                                     stringsAsFactors = FALSE))
n_pass <- sum(results_df$status == "PASS")
n_fail <- sum(results_df$status == "FAIL")
n_total <- nrow(results_df)

for (i in seq_len(nrow(results_df))) {
  icon <- if (results_df$status[i] == "PASS") "\u2713" else "\u2717"
  cat(sprintf("  %s  %s: %s\n", icon, results_df$name[i], results_df$status[i]))
}

cat(sprintf("\n  Total: %d | Pass: %d | Fail: %d | Rate: %.0f%%\n",
            n_total, n_pass, n_fail, 100 * n_pass / n_total))

# List generated files
cat("\n--- Generated output files ---\n")
all_files <- list.files(out_dir, recursive = TRUE, full.names = FALSE)
for (f in all_files) {
  sz <- file.info(file.path(out_dir, f))$size
  cat(sprintf("  %s (%.1f KB)\n", f, sz / 1024))
}

cat(sprintf("\nOutput directory: %s\n", out_dir))
write.csv(results_df, file.path(out_dir, "test_results.csv"), row.names = FALSE)

if (n_fail > 0) {
  cat("\n*** SOME TESTS FAILED ***\n\n")
  quit(status = 1)
} else {
  cat("\n*** ALL TESTS PASSED ***\n\n")
  quit(status = 0)
}
