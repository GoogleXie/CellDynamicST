#' @title Pipeline Orchestrator
#' @description High-level wrapper that runs the full CellDynamicST analysis
#'   pipeline from raw Seurat objects to final results, with provenance
#'   tracking and checkpoint support.
#' @name pipeline
NULL

#' Run the Full CellDynamicST Pipeline
#'
#' Executes the complete analysis pipeline in sequence: data registration,
#' preprocessing (QC + normalization), cell type classification (clustering,
#' NT classification, glial typing), and optionally downstream analyses
#' (dynamics, DEG, WGCNA). Each step is logged with provenance metadata and
#' intermediate results are saved to disk as checkpoints.
#'
#' @param config_path Character. Path to the experiment YAML configuration
#'   file. See \code{vignette("configuration")} for format details.
#' @param seurat_obj A Seurat object. If provided, the "register" step is
#'   skipped and this object is used directly. If NULL, samples are loaded
#'   from paths specified in the config.
#' @param steps Character vector. Pipeline steps to execute. Default: all.
#'   Options: "register", "preprocess", "cluster", "classify_nt",
#'   "classify_glia", "dynamics", "deg", "wgcna", "visualize".
#' @param output_dir Character. Directory for all outputs. Default: "cdst_output".
#' @param checkpoint Logical. Save intermediate checkpoints. Default: TRUE.
#' @param n_cores Integer. Cores for parallel steps. Default: 1.
#' @param verbose Logical. Default: TRUE.
#'
#' @return A \code{\link{CdstResult}} object containing all analysis outputs
#'   and provenance records.
#'
#' @examples
#' \dontrun{
#' # Minimal run — loads samples from config
#' result <- cdst_run("experiment_config.yaml")
#'
#' # Custom steps
#' result <- cdst_run("config.yaml",
#'   steps = c("register", "preprocess", "cluster"))
#'
#' # Resume from checkpoint
#' obj <- readRDS("cdst_output/checkpoint_03_clustered.rds")
#' result <- cdst_run("config.yaml", seurat_obj = obj,
#'   steps = c("classify_nt", "classify_glia", "deg", "wgcna"))
#' }
#'
#' @export
cdst_run <- function(config_path,
                     seurat_obj = NULL,
                     steps = c("register", "preprocess", "cluster",
                               "classify_nt", "classify_glia",
                               "dynamics", "deg", "wgcna", "visualize"),
                     output_dir = "cdst_output",
                     checkpoint = TRUE,
                     n_cores = 1L,
                     verbose = TRUE) {
  # ---- Setup ----
  start_time <- Sys.time()
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (verbose) {
    cli::cli_h1("CellDynamicST Pipeline")
    cli::cli_inform("Started: {format(start_time, '%Y-%m-%d %H:%M:%S')}")
  }

  # Load configuration
  config <- cdst_load_config(config_path)

  # Initialize result object
  result <- new("CdstResult",
    config = config,
    provenance = list(
      pipeline_start = as.character(start_time),
      package_version = as.character(utils::packageVersion("CellDynamicST")),
      r_version = R.version.string,
      steps_requested = steps
    )
  )

  # If a pre-built Seurat object is provided, use it directly
  if (!is.null(seurat_obj)) {
    result@seurat <- seurat_obj
  }

  # ---- Step 1: Register ----
  if ("register" %in% steps) {
    if (verbose) cli::cli_h2("Step 1/9: Data Registration")

    sample_paths <- vapply(config@samples, function(s) s$path, character(1))
    sample_ids   <- vapply(config@samples, function(s) s$id, character(1))
    sample_groups <- vapply(config@samples, function(s) s$group, character(1))

    obj_list <- lapply(seq_along(sample_paths), function(i) {
      if (verbose) cli::cli_inform("Loading sample {sample_ids[i]}...")
      obj <- cdst_load_seurat(sample_paths[i], config)
      obj$sample_id <- sample_ids[i]
      obj$experiment_group <- sample_groups[i]
      obj$genotype <- config@design$genotype_map[[sample_groups[i]]]
      obj$treatment <- config@design$treatment_map[[sample_groups[i]]]
      obj
    })
    names(obj_list) <- sample_ids

    # Merge
    if (length(obj_list) > 1) {
      merged <- obj_list[[1]]
      for (i in 2:length(obj_list)) {
        merged <- merge(merged, obj_list[[i]])
      }
    } else {
      merged <- obj_list[[1]]
    }

    result@seurat <- merged
    result@provenance$register <- list(
      timestamp = as.character(Sys.time()),
      n_samples = length(obj_list),
      n_cells = ncol(merged)
    )

    if (checkpoint) {
      saveRDS(merged, file.path(output_dir, "checkpoint_01_registered.rds"))
      if (verbose) cli::cli_inform("Checkpoint saved: checkpoint_01_registered.rds")
    }
  }

  # ---- Step 2: Preprocess ----
  if ("preprocess" %in% steps) {
    if (verbose) cli::cli_h2("Step 2/9: Preprocessing (QC + Normalization)")

    merged <- cdst_run_qc(result@seurat, config, verbose = verbose)
    merged <- cdst_normalize(merged, config, verbose = verbose)

    result@seurat <- merged
    result@provenance$preprocess <- list(
      timestamp = as.character(Sys.time()),
      n_cells_post_qc = ncol(merged)
    )

    if (checkpoint) {
      saveRDS(merged, file.path(output_dir, "checkpoint_02_preprocessed.rds"))
      if (verbose) cli::cli_inform("Checkpoint saved: checkpoint_02_preprocessed.rds")
    }
  }

  # ---- Step 3: Cluster ----
  if ("cluster" %in% steps) {
    if (verbose) cli::cli_h2("Step 3/9: Clustering")

    merged <- cdst_cluster(result@seurat, config, verbose = verbose)

    result@seurat <- merged
    result@provenance$cluster <- list(
      timestamp = as.character(Sys.time()),
      n_clusters = length(unique(merged$cdst_cluster))
    )

    if (checkpoint) {
      saveRDS(merged, file.path(output_dir, "checkpoint_03_clustered.rds"))
      if (verbose) cli::cli_inform("Checkpoint saved: checkpoint_03_clustered.rds")
    }
  }

  # ---- Step 4: Classify NT ----
  if ("classify_nt" %in% steps) {
    if (verbose) cli::cli_h2("Step 4/9: Neurotransmitter Classification")

    merged <- cdst_classify_nt(result@seurat, verbose = verbose)

    result@seurat <- merged
    result@provenance$classify_nt <- list(
      timestamp = as.character(Sys.time()),
      n_nt_types = length(unique(merged$cdst_nt_type))
    )

    if (checkpoint) {
      saveRDS(merged, file.path(output_dir, "checkpoint_04_nt_classified.rds"))
      if (verbose) cli::cli_inform("Checkpoint saved: checkpoint_04_nt_classified.rds")
    }
  }

  # ---- Step 5: Classify Glia ----
  if ("classify_glia" %in% steps) {
    if (verbose) cli::cli_h2("Step 5/9: Glial Subtype Classification")

    merged <- cdst_classify_glia(result@seurat, verbose = verbose)

    result@seurat <- merged
    result@provenance$classify_glia <- list(
      timestamp = as.character(Sys.time()),
      n_glia_types = length(unique(
        merged$cdst_glia_type[!is.na(merged$cdst_glia_type)]
      ))
    )

    if (checkpoint) {
      saveRDS(merged, file.path(output_dir, "checkpoint_05_glia_classified.rds"))
      if (verbose) cli::cli_inform("Checkpoint saved: checkpoint_05_glia_classified.rds")
    }
  }

  # ---- Step 6: Dynamics ----
  if ("dynamics" %in% steps) {
    if (verbose) cli::cli_h2("Step 6/9: Cell Dynamics Analysis")

    proportions <- cdst_cell_proportions(result@seurat, config, verbose = verbose)
    dist_matrix <- cdst_distribution_matrix(result@seurat, verbose = verbose)

    result@dynamics <- list(
      proportions = proportions,
      distribution_matrix = dist_matrix
    )
    result@provenance$dynamics <- list(
      timestamp = as.character(Sys.time())
    )

    if (checkpoint) {
      saveRDS(result@dynamics, file.path(output_dir, "checkpoint_06_dynamics.rds"))
      if (verbose) cli::cli_inform("Checkpoint saved: checkpoint_06_dynamics.rds")
    }
  }

  # ---- Step 7: DEG ----
  if ("deg" %in% steps) {
    if (verbose) cli::cli_h2("Step 7/9: Differential Expression Analysis")

    deg_results <- cdst_deg_interregional(
      result@seurat, config,
      n_cores = n_cores,
      verbose = verbose
    )

    result@deg <- deg_results
    result@provenance$deg <- list(
      timestamp = as.character(Sys.time()),
      n_total_degs = nrow(deg_results)
    )

    if (checkpoint) {
      saveRDS(deg_results, file.path(output_dir, "checkpoint_07_deg.rds"))
      if (verbose) cli::cli_inform("Checkpoint saved: checkpoint_07_deg.rds")
    }
  }

  # ---- Step 8: WGCNA ----
  if ("wgcna" %in% steps) {
    if (verbose) cli::cli_h2("Step 8/9: WGCNA Network Analysis")

    wgcna_data <- cdst_wgcna_prepare(result@seurat, config, verbose = verbose)
    wgcna_result <- cdst_wgcna_detect(wgcna_data, verbose = verbose)
    trait_cor <- cdst_wgcna_trait_cor(wgcna_result, verbose = verbose)
    hub_genes <- cdst_wgcna_hub_genes(wgcna_result, verbose = verbose)

    result@wgcna <- list(
      wgcna_result = wgcna_result,
      trait_correlations = trait_cor,
      hub_genes = hub_genes
    )
    result@provenance$wgcna <- list(
      timestamp = as.character(Sys.time()),
      n_modules = length(unique(wgcna_result$module_colors)),
      n_hub_genes = nrow(hub_genes)
    )

    if (checkpoint) {
      saveRDS(result@wgcna, file.path(output_dir, "checkpoint_08_wgcna.rds"))
      if (verbose) cli::cli_inform("Checkpoint saved: checkpoint_08_wgcna.rds")
    }
  }

  # ---- Step 9: Visualize ----
  if ("visualize" %in% steps) {
    if (verbose) cli::cli_h2("Step 9/9: Generating Visualizations")

    fig_dir <- file.path(output_dir, "figures")
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(fig_dir, "wgcna"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(fig_dir, "deg"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(fig_dir, "heatmaps"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(fig_dir, "hierarchy"), recursive = TRUE, showWarnings = FALSE)
    n_saved <- 0L

    # --- 1. QC Summary ---
    if (!is.null(result@seurat)) {
      tryCatch({
        p <- cdst_plot_qc(result@seurat)
        ggplot2::ggsave(file.path(fig_dir, "qc_summary.pdf"),
                         p, width = 15, height = 5)
        n_saved <- n_saved + 1L
        if (verbose) cli::cli_inform("Saved: qc_summary.pdf")
      }, error = function(e) {
        if (verbose) cli::cli_warn("QC plot failed: {e$message}")
      })
    }

    # --- 2. Spatial Cell Type Map ---
    if (!is.null(result@seurat) &&
        "cdst_cell_type" %in% colnames(result@seurat@meta.data)) {
      tryCatch({
        p <- cdst_plot_spatial(result@seurat)
        ggplot2::ggsave(file.path(fig_dir, "spatial_cell_types.pdf"),
                         p, width = 12, height = 10)
        n_saved <- n_saved + 1L
        if (verbose) cli::cli_inform("Saved: spatial_cell_types.pdf")
      }, error = function(e) {
        if (verbose) cli::cli_warn("Spatial plot failed: {e$message}")
      })
    }

    # --- 3. Distribution Heatmap ---
    if (!is.null(result@dynamics$distribution_matrix)) {
      tryCatch({
        cdst_plot_distribution_heatmap(
          result@dynamics$distribution_matrix,
          output_file = file.path(fig_dir, "heatmaps", "distribution_heatmap.png")
        )
        n_saved <- n_saved + 1L
        if (verbose) cli::cli_inform("Saved: heatmaps/distribution_heatmap.png")
      }, error = function(e) {
        if (verbose) cli::cli_warn("Distribution heatmap failed: {e$message}")
      })
    }

    # --- 4. Cell Proportions Bar Plot ---
    if (!is.null(result@dynamics$proportions)) {
      tryCatch({
        p <- cdst_plot_proportions(result@dynamics$proportions)
        ggplot2::ggsave(file.path(fig_dir, "cell_proportions.pdf"),
                         p, width = 12, height = 8)
        n_saved <- n_saved + 1L
        if (verbose) cli::cli_inform("Saved: cell_proportions.pdf")
      }, error = function(e) {
        if (verbose) cli::cli_warn("Proportions plot failed: {e$message}")
      })
    }

    # --- 5. Hierarchy Dendrogram ---
    if (!is.null(result@seurat) &&
        any(grepl("brain_region", colnames(result@seurat@meta.data)))) {
      tryCatch({
        hier <- cdst_build_hierarchy(result@seurat)
        cdst_plot_dendrogram(
          hier,
          output_file = file.path(fig_dir, "hierarchy", "hierarchy_tree.png")
        )
        n_saved <- n_saved + 1L
        if (verbose) cli::cli_inform("Saved: hierarchy/hierarchy_tree.png")
      }, error = function(e) {
        if (verbose) cli::cli_warn("Hierarchy dendrogram failed: {e$message}")
      })

      # Sankey diagram
      tryCatch({
        hier <- cdst_build_hierarchy(result@seurat)
        cdst_plot_sankey(
          hier,
          output_file = file.path(fig_dir, "hierarchy", "hierarchy_sankey.html")
        )
        n_saved <- n_saved + 1L
        if (verbose) cli::cli_inform("Saved: hierarchy/hierarchy_sankey.html")
      }, error = function(e) {
        if (verbose) cli::cli_warn("Hierarchy sankey failed: {e$message}")
      })
    }

    # --- 6. Volcano Plots (per contrast) ---
    if (is.data.frame(result@deg) && nrow(result@deg) > 0) {
      contrasts <- unique(result@deg$contrast)
      for (ct in contrasts) {
        tryCatch({
          p <- cdst_plot_volcano(result@deg, contrast = ct)
          safe_name <- gsub("[^A-Za-z0-9]", "_", ct)
          ggplot2::ggsave(
            file.path(fig_dir, "deg", paste0("volcano_", safe_name, ".pdf")),
            p, width = 8, height = 6
          )
          n_saved <- n_saved + 1L
        }, error = function(e) {
          if (verbose) cli::cli_warn("Volcano plot failed for {ct}: {e$message}")
        })
      }
      if (verbose) cli::cli_inform("Saved: {length(contrasts)} volcano plots")

      # DEG scatter (MA) plots
      for (ct in contrasts) {
        tryCatch({
          safe_name <- gsub("[^A-Za-z0-9]", "_", ct)
          cdst_plot_deg_scatter(
            result@deg, contrast = ct,
            output_file = file.path(fig_dir, "deg", paste0("ma_plot_", safe_name, ".png"))
          )
          n_saved <- n_saved + 1L
        }, error = function(e) {
          if (verbose) cli::cli_warn("MA plot failed for {ct}: {e$message}")
        })
      }
    }

    # --- 7. Multi-dimensional Heatmap (Python-based) ---
    if (is.data.frame(result@deg) && nrow(result@deg) > 0) {
      tryCatch({
        cdst_plot_multidim_heatmap(
          result@deg,
          output_dir = file.path(fig_dir, "heatmaps")
        )
        n_saved <- n_saved + 1L
        if (verbose) cli::cli_inform("Saved: heatmaps/ (multidim heatmaps)")
      }, error = function(e) {
        if (verbose) cli::cli_warn("Multidim heatmap failed: {e$message}")
      })
    }

    # --- 8. WGCNA Visualizations ---
    if (length(result@wgcna) > 0) {
      # Trait-module correlation heatmap
      if (!is.null(result@wgcna$trait_correlations)) {
        tryCatch({
          cdst_plot_trait_heatmap(
            result@wgcna$trait_correlations,
            output_file = file.path(fig_dir, "wgcna", "trait_heatmap.png")
          )
          n_saved <- n_saved + 1L
          if (verbose) cli::cli_inform("Saved: wgcna/trait_heatmap.png")
        }, error = function(e) {
          if (verbose) cli::cli_warn("WGCNA trait heatmap failed: {e$message}")
        })
      }

      # Hub gene network (R-native via igraph)
      if (!is.null(result@wgcna$hub_genes)) {
        tryCatch({
          cdst_plot_network(
            result@wgcna$hub_genes,
            wgcna_result = result@wgcna$wgcna_result,
            output_file = file.path(fig_dir, "wgcna", "hub_network.png")
          )
          n_saved <- n_saved + 1L
          if (verbose) cli::cli_inform("Saved: wgcna/hub_network.png")
        }, error = function(e) {
          if (verbose) cli::cli_warn("Hub network failed: {e$message}")
        })
      }

      # Python-based circos and hub network (requires reticulate + Python)
      # These use CSV files as input; write temp CSVs from result objects
      if (!is.null(result@wgcna$hub_genes)) {
        tryCatch({
          tmp_hub <- tempfile(fileext = ".csv")
          utils::write.csv(result@wgcna$hub_genes, tmp_hub, row.names = FALSE)
          cdst_plot_wgcna_circos(
            tmp_hub,
            output_file = file.path(fig_dir, "wgcna", "circos.png")
          )
          n_saved <- n_saved + 1L
          if (verbose) cli::cli_inform("Saved: wgcna/circos.png")
        }, error = function(e) {
          if (verbose) cli::cli_warn("Circos plot failed: {e$message}")
        })
      }
    }

    if (verbose) cli::cli_inform(c(
      "v" = "Visualization complete: {n_saved} figures saved to {fig_dir}"
    ))

    result@provenance$visualize <- list(
      timestamp = as.character(Sys.time()),
      output_dir = fig_dir,
      n_figures = n_saved
    )
  }

  # ---- Finalize ----
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")

  result@provenance$pipeline_end <- as.character(end_time)
  result@provenance$elapsed_minutes <- as.numeric(elapsed)

  # Save final result
  saveRDS(result, file.path(output_dir, "cdst_result.rds"))

  # Save provenance log as YAML
  provenance_yaml <- yaml::as.yaml(result@provenance)
  writeLines(provenance_yaml, file.path(output_dir, "provenance.yaml"))

  if (verbose) {
    cli::cli_h1("Pipeline Complete")
    cli::cli_inform(c(
      "v" = "Elapsed: {round(as.numeric(elapsed), 1)} minutes",
      "i" = "Results saved to: {output_dir}",
      "i" = "Provenance log: {file.path(output_dir, 'provenance.yaml')}"
    ))
  }

  result
}


#' Print Pipeline Summary
#'
#' Prints a human-readable summary of a CdstResult object, showing which
#' pipeline steps were completed and key statistics from each step.
#'
#' @param result A \code{\link{CdstResult}} object.
#'
#' @export
cdst_summary <- function(result) {
  checkmate::assert_class(result, "CdstResult")

  cli::cli_h1("CellDynamicST Analysis Summary")

  prov <- result@provenance

  cli::cli_inform(c(
    "i" = "Package version: {prov$package_version %||% 'unknown'}",
    "i" = "R version: {prov$r_version %||% 'unknown'}",
    "i" = "Started: {prov$pipeline_start %||% 'N/A'}",
    "i" = "Elapsed: {round(prov$elapsed_minutes %||% 0, 1)} minutes"
  ))

  cli::cli_h2("Steps Completed")

  step_names <- c("register", "preprocess", "cluster",
                   "classify_nt", "classify_glia", "dynamics", "deg",
                   "wgcna", "visualize")

  for (step in step_names) {
    if (!is.null(prov[[step]])) {
      info <- prov[[step]]
      details <- setdiff(names(info), "timestamp")
      if (length(details) > 0) {
        detail_str <- paste(
          vapply(details, function(d) paste0(d, "=", info[[d]]), character(1)),
          collapse = ", "
        )
        cli::cli_inform(c("v" = "{step}: {detail_str}"))
      } else {
        cli::cli_inform(c("v" = "{step}: completed"))
      }
    }
  }

  invisible(result)
}
