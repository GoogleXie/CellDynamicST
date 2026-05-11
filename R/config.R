#' @title Module 0: Configuration
#' @description Functions for loading, validating, and creating CellDynamicST
#'   configuration objects from YAML files.
#' @name config
NULL

# ============================================================================
# Default parameter sets
# ============================================================================

.default_qc <- list(

  min_count_percentile   = 0.05,
  min_feature_percentile = 0.05,
  min_cells_per_fov      = 500,
  min_transcripts_per_cell = 100,
  max_neg_probe_per_cell = 1.0
)

.default_clustering <- list(
  n_pcs                  = 100L,
  bias_threshold         = 0.7,
  annoy_trees            = 100L,
  k_neighbors            = 50L,
  resolution             = 20,
  rf_ntree               = 200L
)

.default_wgcna <- list(
  soft_power             = NULL,
  min_module_size        = 20L,
  deep_split             = 3L,
  merge_cut_height       = 0.10,
  min_cells_per_region   = 100L
)

.default_atlas <- list(
  name            = "Allen_CCFv3",
  plane           = "sagittal",
  z_plane         = 510,
  region_columns  = list(
    L3 = "brain_region_L3",
    L5 = "brain_region_L5",
    L7 = "brain_region_L7"
  )
)


# ============================================================================
# Public API
# ============================================================================

#' Load CellDynamicST Configuration from YAML
#'
#' Reads a YAML configuration file and returns a validated \code{CdstConfig}
#' object. The YAML file specifies experimental groups, pairwise comparisons,
#' QC thresholds, clustering parameters, WGCNA parameters, and atlas settings.
#' Any parameters not specified in the YAML file are filled with sensible
#' defaults derived from the original CellDynamicST paper.
#'
#' @param config_path Character. Path to the YAML configuration file.
#' @param verbose Logical. Print informational messages. Default: TRUE.
#'
#' @return A validated \code{\link{CdstConfig}} object.
#'
#' @details
#' The YAML file must contain at minimum:
#' \itemize{
#'   \item \code{project_name}: A string identifier for the project.
#'   \item \code{groups}: A list of experimental groups, each with \code{name},
#'     \code{genotype}, and \code{treatment} fields. Exactly one group should
#'     have \code{is_reference: true}.
#' }
#'
#' All other sections (\code{samples}, \code{comparisons}, \code{qc},
#' \code{clustering}, \code{wgcna}, \code{atlas}) are optional and will be
#' filled with defaults.
#'
#' @examples
#' config_file <- system.file("templates", "experiment_config.yaml",
#'                            package = "CellDynamicST")
#' config <- cdst_load_config(config_file)
#' config
#'
#' @export
cdst_load_config <- function(config_path, verbose = TRUE) {
  checkmate::assert_file_exists(config_path, extension = c("yaml", "yml"))

  if (verbose) cli::cli_inform("Loading configuration from {.file {config_path}}")

  raw <- yaml::read_yaml(config_path)

  # Validate required top-level fields
  if (is.null(raw$project_name)) {
    cli::cli_abort("Configuration file must contain a {.field project_name} field.")
  }
  if (is.null(raw$groups) || length(raw$groups) == 0) {
    cli::cli_abort("Configuration file must contain at least one group in {.field groups}.")
  }

  # Parse groups
  groups <- lapply(raw$groups, function(g) {
    list(
      name         = as.character(g$name),
      genotype     = as.character(g$genotype %||% ""),
      treatment    = as.character(g$treatment %||% ""),
      is_reference = isTRUE(g$is_reference)
    )
  })

  # Determine reference group
  ref_groups <- Filter(function(g) g$is_reference, groups)
  if (length(ref_groups) == 0) {
    cli::cli_warn("No group marked as {.field is_reference: true}. Using first group as reference.")
    groups[[1]]$is_reference <- TRUE
    ref_group <- groups[[1]]$name
  } else if (length(ref_groups) > 1) {
    cli::cli_warn("Multiple groups marked as reference. Using the first one: {.val {ref_groups[[1]]$name}}")
    ref_group <- ref_groups[[1]]$name
  } else {
    ref_group <- ref_groups[[1]]$name
  }

  # Parse comparisons (or auto-generate)
  if (!is.null(raw$comparisons) && length(raw$comparisons) > 0) {
    comparisons <- lapply(raw$comparisons, function(c) {
      list(
        name   = as.character(c$name),
        group1 = as.character(c$group1),
        group2 = as.character(c$group2)
      )
    })
  } else {
    comparisons <- .auto_generate_comparisons(groups)
    if (verbose) {
      cli::cli_inform("Auto-generated {length(comparisons)} pairwise comparisons.")
    }
  }

  # Parse samples (optional — user may provide seurat_obj directly)
  samples <- list()
  if (!is.null(raw$samples) && length(raw$samples) > 0) {
    samples <- lapply(raw$samples, function(s) {
      list(
        id    = as.character(s$id %||% s$name %||% ""),
        path  = as.character(s$path %||% ""),
        group = as.character(s$group %||% "")
      )
    })
    if (verbose) {
      cli::cli_inform("Samples defined: {length(samples)}")
    }
  }

  # Build design maps from groups
  group_names <- vapply(groups, function(g) g$name, character(1))
  genotype_map <- stats::setNames(
    vapply(groups, function(g) g$genotype, character(1)),
    group_names
  )
  treatment_map <- stats::setNames(
    vapply(groups, function(g) g$treatment, character(1)),
    group_names
  )
  design <- list(
    genotype_map  = genotype_map,
    treatment_map = treatment_map
  )

  # Merge user parameters with defaults
  qc_params         <- .merge_params(raw$qc, .default_qc)
  clustering_params <- .merge_params(raw$clustering, .default_clustering)
  wgcna_params      <- .merge_params(raw$wgcna, .default_wgcna)
  atlas_params      <- .merge_params(raw$atlas, .default_atlas)

  # Build config object
  config <- new("CdstConfig",
    project_name    = as.character(raw$project_name),
    species         = as.character(raw$species %||% "mouse"),
    groups          = groups,
    comparisons     = comparisons,
    reference_group = ref_group,
    samples         = samples,
    design          = design,
    qc              = qc_params,
    clustering      = clustering_params,
    wgcna           = wgcna_params,
    atlas           = atlas_params,
    output_dir      = as.character(raw$output_dir %||% "results")
  )

  # Validate
  validObject(config)

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Configuration loaded successfully.",
      "i" = "Project: {.val {config@project_name}}",
      "i" = "Species: {.val {config@species}}",
      "i" = "Groups: {length(config@groups)} ({.val {ref_group}} is reference)",
      "i" = "Comparisons: {length(config@comparisons)}",
      "i" = "Samples: {length(config@samples)}"
    ))
  }

  config
}


#' Validate a CdstConfig Object
#'
#' Performs comprehensive validation of a \code{CdstConfig} object, checking
#' that all required fields are present and internally consistent.
#'
#' @param config A \code{\link{CdstConfig}} object.
#'
#' @return Invisibly returns TRUE if validation passes. Throws an error with
#'   a detailed message if validation fails.
#'
#' @export
cdst_validate_config <- function(config) {
  checkmate::assert_class(config, "CdstConfig")

  errors <- character()

  # Validate groups
  group_names <- vapply(config@groups, function(g) g$name, character(1))
  if (any(duplicated(group_names))) {
    errors <- c(errors, paste0(
      "Duplicate group names found: ",
      paste(group_names[duplicated(group_names)], collapse = ", ")
    ))
  }

  # Validate comparisons reference existing groups
  for (comp in config@comparisons) {
    if (!comp$group1 %in% group_names) {
      errors <- c(errors, paste0(
        "Comparison '", comp$name, "' references unknown group1: '", comp$group1, "'"
      ))
    }
    if (!comp$group2 %in% group_names) {
      errors <- c(errors, paste0(
        "Comparison '", comp$name, "' references unknown group2: '", comp$group2, "'"
      ))
    }
  }

  # Validate samples reference existing groups
  for (samp in config@samples) {
    if (nchar(samp$group) > 0 && !samp$group %in% group_names) {
      errors <- c(errors, paste0(
        "Sample '", samp$id, "' references unknown group: '", samp$group, "'"
      ))
    }
  }

  # Validate QC parameters
  if (config@qc$min_count_percentile < 0 || config@qc$min_count_percentile > 1) {
    errors <- c(errors, "qc$min_count_percentile must be between 0 and 1")
  }

  # Validate clustering parameters
  if (config@clustering$n_pcs < 10 || config@clustering$n_pcs > 500) {
    errors <- c(errors, "clustering$n_pcs must be between 10 and 500")
  }

  if (length(errors) > 0) {
    cli::cli_abort(c(
      "Configuration validation failed:",
      set_names(errors, rep("x", length(errors)))
    ))
  }

  invisible(TRUE)
}


#' Create a Default Configuration Object
#'
#' Creates a \code{CdstConfig} object with default parameters and a minimal
#' two-group experimental design. Useful for testing and as a starting point
#' for customization.
#'
#' @param project_name Character. Project name. Default: "CellDynamicST_Project".
#' @param species Character. "mouse" or "human". Default: "mouse".
#'
#' @return A \code{\link{CdstConfig}} object with default parameters.
#'
#' @examples
#' config <- cdst_default_config()
#' config
#'
#' @export
cdst_default_config <- function(project_name = "CellDynamicST_Project",
                                species = "mouse") {
  groups <- list(
    list(name = "Control", genotype = "WT", treatment = "VEH", is_reference = TRUE),
    list(name = "Treatment", genotype = "WT", treatment = "TRT", is_reference = FALSE)
  )

  group_names <- vapply(groups, function(g) g$name, character(1))
  genotype_map <- stats::setNames(
    vapply(groups, function(g) g$genotype, character(1)),
    group_names
  )
  treatment_map <- stats::setNames(
    vapply(groups, function(g) g$treatment, character(1)),
    group_names
  )

  new("CdstConfig",
    project_name    = project_name,
    species         = species,
    groups          = groups,
    comparisons     = list(
      list(name = "treatment_effect", group1 = "Control", group2 = "Treatment")
    ),
    reference_group = "Control",
    samples         = list(),
    design          = list(genotype_map = genotype_map,
                           treatment_map = treatment_map),
    qc              = .default_qc,
    clustering      = .default_clustering,
    wgcna           = .default_wgcna,
    atlas           = .default_atlas,
    output_dir      = "results"
  )
}


#' Create a Configuration Template File
#'
#' Copies the package's built-in YAML configuration template to a specified
#' directory. The template contains all available parameters with documentation
#' comments, making it easy for users to customize for their experiment.
#'
#' @param output_path Character. Path where the template should be written.
#'   Default: "experiment_config.yaml" in the current working directory.
#' @param overwrite Logical. Overwrite existing file. Default: FALSE.
#'
#' @return Invisibly returns the output path.
#'
#' @examples
#' \dontrun{
#' cdst_create_config_template("my_experiment.yaml")
#' }
#'
#' @export
cdst_create_config_template <- function(output_path = "experiment_config.yaml",
                                        overwrite = FALSE) {
  template <- system.file("templates", "experiment_config.yaml",
                          package = "CellDynamicST")

  if (!file.exists(template)) {
    cli::cli_abort("Package template file not found. Reinstall CellDynamicST.")
  }

  if (file.exists(output_path) && !overwrite) {
    cli::cli_abort(c(
      "File already exists: {.file {output_path}}",
      "i" = "Use {.code overwrite = TRUE} to replace it."
    ))
  }

  file.copy(template, output_path, overwrite = overwrite)
  cli::cli_inform(c(
    "v" = "Configuration template created: {.file {output_path}}",
    "i" = "Edit this file to match your experimental design, then load with {.fn cdst_load_config}."
  ))

  invisible(output_path)
}


# ============================================================================
# Internal helpers
# ============================================================================

#' Merge user parameters with defaults
#' @keywords internal
.merge_params <- function(user, defaults) {
  if (is.null(user)) return(defaults)
  result <- defaults
  for (name in names(user)) {
    result[[name]] <- user[[name]]
  }
  result
}

#' Auto-generate pairwise comparisons from groups
#' @keywords internal
.auto_generate_comparisons <- function(groups) {
  group_names <- vapply(groups, function(g) g$name, character(1))
  n <- length(group_names)
  if (n < 2) return(list())

  comparisons <- list()
  idx <- 1
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      comparisons[[idx]] <- list(
        name   = paste0(group_names[i], "_vs_", group_names[j]),
        group1 = group_names[i],
        group2 = group_names[j]
      )
      idx <- idx + 1
    }
  }
  comparisons
}

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x
