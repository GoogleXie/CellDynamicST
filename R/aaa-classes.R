#' @title CellDynamicST S4 Class Definitions
#' @description Core S4 classes used throughout the CellDynamicST package.
#' @name CellDynamicST-classes
NULL

# ============================================================================
# CdstConfig: Configuration Object
# ============================================================================

#' CellDynamicST Configuration Object
#'
#' An S4 class that encapsulates all configuration parameters for a
#' CellDynamicST analysis pipeline run. Created by \code{\link{cdst_load_config}}
#' from a YAML configuration file, or programmatically via
#' \code{\link{cdst_default_config}}.
#'
#' @slot project_name Character. Project identifier used in output filenames.
#' @slot species Character. Species identifier ("mouse" or "human").
#' @slot groups List. Each element is a named list with fields: name, genotype,
#'   treatment, is_reference (logical).
#' @slot comparisons List. Each element is a named list with fields: name,
#'   group1, group2 defining pairwise comparisons for DEG analysis.
#' @slot reference_group Character. Name of the reference/baseline group used
#'   for clustering. Must match one of the group names.
#' @slot samples List. Each element is a named list with fields: id, path,
#'   group. Used by \code{cdst_run()} to load data automatically.
#' @slot design List. Experimental design metadata including genotype_map and
#'   treatment_map (named character vectors mapping group names to genotypes
#'   and treatments respectively).
#' @slot qc List. QC threshold parameters: min_count_percentile,
#'   min_feature_percentile, min_cells_per_fov, etc.
#' @slot clustering List. Clustering parameters: n_pcs, bias_threshold,
#'   annoy_trees, k_neighbors, resolution, rf_ntree.
#' @slot wgcna List. WGCNA parameters: soft_power, min_module_size,
#'   deep_split, merge_cut_height, min_cells_per_region.
#' @slot atlas List. Brain atlas configuration: name, plane, z_plane,
#'   region_columns.
#' @slot output_dir Character. Base output directory for results.
#'
#' @export
#' @examples
#' # Load from YAML
#' config <- cdst_load_config(system.file("templates", "experiment_config.yaml",
#'                                         package = "CellDynamicST"))
#' config@project_name
setClass("CdstConfig",
  representation(
    project_name    = "character",
    species         = "character",
    groups          = "list",
    comparisons     = "list",
    reference_group = "character",
    samples         = "list",
    design          = "list",
    qc              = "list",
    clustering      = "list",
    wgcna           = "list",
    atlas           = "list",
    output_dir      = "character"
  ),
  prototype(
    project_name    = "CellDynamicST_Project",
    species         = "mouse",
    groups          = list(),
    comparisons     = list(),
    reference_group = character(0),
    samples         = list(),
    design          = list(genotype_map = character(0),
                           treatment_map = character(0)),
    qc              = list(),
    clustering      = list(),
    wgcna           = list(),
    atlas           = list(),
    output_dir      = "results"
  ),
  validity = function(object) {
    errors <- character()

    if (length(object@project_name) != 1 || nchar(object@project_name) == 0) {
      errors <- c(errors, "project_name must be a non-empty single string")
    }
    if (!object@species %in% c("mouse", "human")) {
      errors <- c(errors, "species must be 'mouse' or 'human'")
    }
    if (length(object@groups) == 0) {
      errors <- c(errors, "At least one experimental group must be defined")
    }
    if (length(object@reference_group) == 1) {
      group_names <- vapply(object@groups, function(g) g$name, character(1))
      if (!object@reference_group %in% group_names) {
        errors <- c(errors, paste0(
          "reference_group '", object@reference_group,
          "' not found in groups. Available: ",
          paste(group_names, collapse = ", ")
        ))
      }
    }

    if (length(errors) == 0) TRUE else errors
  }
)

#' Show method for CdstConfig
#' @param object A CdstConfig object
#' @keywords internal
setMethod("show", "CdstConfig", function(object) {
  cat("CellDynamicST Configuration\n")
  cat("===========================\n")
  cat("Project:    ", object@project_name, "\n")
  cat("Species:    ", object@species, "\n")
  cat("Groups:     ", length(object@groups), "\n")
  cat("Comparisons:", length(object@comparisons), "\n")
  cat("Reference:  ", object@reference_group, "\n")
  cat("Samples:    ", length(object@samples), "\n")
  cat("Output dir: ", object@output_dir, "\n")
})


# ============================================================================
# CdstResult: Pipeline Result Object
# ============================================================================

#' CellDynamicST Pipeline Result Object
#'
#' An S4 class that encapsulates all results from a CellDynamicST pipeline run,
#' including the annotated Seurat object, configuration, and all downstream
#' analysis results. The provenance slot records computational provenance for
#' FAIR4RS compliance.
#'
#' @slot seurat A Seurat object with all annotations added by the pipeline.
#' @slot config The CdstConfig object used for the analysis.
#' @slot modules_run Character vector of module names that were executed.
#' @slot dynamics List. Cell dynamics results (proportions, distribution_matrix).
#' @slot deg Data.frame. Differential expression results.
#' @slot wgcna List. WGCNA results (wgcna_result, trait_correlations, hub_genes).
#' @slot enrichment List. Enrichment analysis results.
#' @slot provenance List recording package version, R session info, timestamps,
#'   and parameter values for each module execution.
#'
#' @export
setClass("CdstResult",
  representation(
    seurat     = "ANY",
    config     = "CdstConfig",
    modules_run = "character",
    dynamics   = "list",
    deg        = "ANY",
    wgcna      = "list",
    enrichment = "list",
    provenance = "list"
  ),
  prototype(
    seurat     = NULL,
    modules_run = character(0),
    dynamics   = list(),
    deg        = data.frame(),
    wgcna      = list(),
    enrichment = list(),
    provenance = list()
  )
)

#' Show method for CdstResult
#' @param object A CdstResult object
#' @keywords internal
setMethod("show", "CdstResult", function(object) {
  cat("CellDynamicST Result\n")
  cat("====================\n")
  cat("Project:      ", object@config@project_name, "\n")
  cat("Modules run:  ", paste(object@modules_run, collapse = ", "), "\n")
  if (!is.null(object@seurat)) {
    cat("Cells:        ", ncol(object@seurat), "\n")
    cat("Genes:        ", nrow(object@seurat), "\n")
  }
  if (is.data.frame(object@deg) && nrow(object@deg) > 0) {
    cat("DEG results:  ", nrow(object@deg), " genes tested\n")
  }
  if (length(object@wgcna) > 0) {
    cat("WGCNA:         completed\n")
  }
  if (length(object@dynamics) > 0) {
    cat("Dynamics:      completed\n")
  }
})
