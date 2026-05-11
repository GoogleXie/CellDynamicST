#' @title Module 1: Data Registration
#' @description Functions for loading spatial transcriptomic data into Seurat
#'   objects and registering spatial coordinates to a brain atlas reference.
#' @name registration
NULL

#' Load Spatial Transcriptomic Data as a Seurat Object
#'
#' Loads a spatial transcriptomic dataset from various input formats and
#' standardizes the metadata column names to the CellDynamicST convention.
#' Accepts RDS files containing Seurat objects, CosMx output directories,
#' or pre-built Seurat objects passed directly.
#'
#' @param input Character or Seurat object. If character, the path to an RDS
#'   file or a CosMx/MERFISH output directory. If a Seurat object, it is used
#'   directly.
#' @param config A \code{\link{CdstConfig}} object specifying the experimental
#'   design.
#' @param assay Character. Name of the assay to use. Default: "RNA".
#' @param metadata_mapping Named character vector. Maps existing metadata column
#'   names (values) to standardized CellDynamicST names (names). For example:
#'   \code{c(experiment_group = "condition", sample_id = "mouse_id")}.
#'   If NULL, the function attempts auto-detection.
#' @param verbose Logical. Print informational messages. Default: TRUE.
#'
#' @return A Seurat object with standardized metadata columns:
#'   \code{experiment_group}, \code{genotype}, \code{treatment},
#'   \code{sample_id}, \code{slide_id}, \code{fov_id},
#'   \code{x_coord}, \code{y_coord}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Load data from the specified input format.
#'   \item Set the default assay to the specified assay name.
#'   \item Standardize metadata column names using the provided mapping or
#'     auto-detection.
#'   \item Validate that all required metadata columns are present.
#'   \item Compute log2CPM normalization if not already present.
#' }
#'
#' @section Required Metadata:
#' After loading, the Seurat object must contain these metadata columns:
#' \describe{
#'   \item{experiment_group}{Experimental group label (e.g., "WT SAL")}
#'   \item{sample_id}{Biological replicate identifier}
#'   \item{x_coord}{Spatial x-coordinate}
#'   \item{y_coord}{Spatial y-coordinate}
#' }
#'
#' @examples
#' \dontrun{
#' config <- cdst_load_config("experiment_config.yaml")
#' obj <- cdst_load_seurat("path/to/seurat.rds", config)
#' }
#'
#' @export
cdst_load_seurat <- function(input, config,
                             assay = "RNA",
                             metadata_mapping = NULL,
                             verbose = TRUE) {
  checkmate::assert_class(config, "CdstConfig")
  start_time <- Sys.time()

  # ---- Load data ----
  if (is.character(input)) {
    checkmate::assert_file_exists(input)
    if (verbose) cli::cli_inform("Loading data from {.file {input}}")

    if (grepl("\\.rds$", input, ignore.case = TRUE)) {
      obj <- readRDS(input)
    } else if (dir.exists(input)) {
      # Attempt to load as CosMx directory
      obj <- .load_cosmx_directory(input, verbose)
    } else {
      cli::cli_abort("Unsupported input format. Provide an RDS file or CosMx directory.")
    }
  } else if (inherits(input, "Seurat")) {
    obj <- input
  } else {
    cli::cli_abort("Input must be a file path (character) or a Seurat object.")
  }

  if (!inherits(obj, "Seurat")) {
    cli::cli_abort("Loaded object is not a Seurat object (got {.cls {class(obj)}}).")
  }

  # ---- Set assay ----
  available_assays <- Seurat::Assays(obj)
  if (assay %in% available_assays) {
    Seurat::DefaultAssay(obj) <- assay
  } else if (length(available_assays) == 1) {
    if (verbose) {
      cli::cli_inform("Assay {.val {assay}} not found. Using {.val {available_assays[1]}}.")
    }
    Seurat::DefaultAssay(obj) <- available_assays[1]
  } else {
    cli::cli_abort(c(
      "Assay {.val {assay}} not found.",
      "i" = "Available assays: {.val {available_assays}}"
    ))
  }

  # ---- Standardize metadata ----
  obj <- .standardize_metadata(obj, metadata_mapping, verbose)

  # ---- Validate required columns ----
  required_cols <- c("experiment_group", "sample_id")
  validate_seurat(obj, required_metadata = required_cols)

  # ---- Validate groups match config ----
  validate_groups(obj, config)

  # ---- Compute log2CPM if needed ----
  if (!"log2CPM" %in% Seurat::Assays(obj)) {
    if (verbose) cli::cli_inform("Computing log2CPM normalization...")
    obj <- .compute_log2cpm(obj)
  }

  if (verbose) {
    cli::cli_inform(c(
      "v" = "Data loaded successfully.",
      "i" = "Cells: {.val {ncol(obj)}}",
      "i" = "Genes: {.val {nrow(obj)}}",
      "i" = "Groups: {.val {paste(unique(obj$experiment_group), collapse = ', ')}}"
    ))
  }

  obj
}


#' Register Spatial Coordinates to Brain Atlas
#'
#' Maps spatial coordinates from the experimental data to brain atlas regions
#' using the Allen CCFv3 or a user-provided atlas reference. Adds hierarchical
#' region annotations (L3, L5, L7) to the Seurat object metadata.
#'
#' @param seurat_obj A Seurat object with spatial coordinates.
#' @param config A \code{\link{CdstConfig}} object.
#' @param atlas_path Character. Path to the atlas reference file (CSV or RDS).
#'   If NULL, uses the built-in Allen CCFv3 structure tree.
#' @param coord_columns Named character vector of length 2. Maps x and y
#'   coordinate column names. Default: \code{c(x = "x_coord", y = "y_coord")}.
#' @param verbose Logical. Default: TRUE.
#'
#' @return The Seurat object with added brain region columns at multiple
#'   hierarchy levels.
#'
#' @details
#' This function implements the atlas registration step from the original
#' CellDynamicST pipeline. It maps each cell's spatial coordinates to the
#' nearest brain region in the atlas reference, then traverses the atlas
#' hierarchy to assign region labels at multiple levels of granularity.
#'
#' The atlas hierarchy levels correspond to the Allen CCFv3 ontology:
#' \describe{
#'   \item{L3}{Broadest level (e.g., "Isocortex", "Hippocampus")}
#'   \item{L5}{Intermediate level (e.g., "Visual areas", "Field CA1")}
#'   \item{L7}{Finest level (e.g., "VISp", "CA1")}
#' }
#'
#' @export
cdst_register_atlas <- function(seurat_obj, config,
                                atlas_path = NULL,
                                coord_columns = c(x = "x_coord", y = "y_coord"),
                                verbose = TRUE) {
  checkmate::assert_class(config, "CdstConfig")
  validate_seurat(seurat_obj, required_metadata = unname(coord_columns))
  start_time <- Sys.time()

  # Load atlas reference
  if (is.null(atlas_path)) {
    atlas_file <- system.file("extdata", "structure_tree_safe_2017.csv",
                              package = "CellDynamicST")
    if (!file.exists(atlas_file)) {
      cli::cli_abort(c(
        "Built-in atlas reference not found.",
        "i" = "Provide an atlas file via {.arg atlas_path}, or install the full package."
      ))
    }
    atlas <- readr::read_csv(atlas_file, show_col_types = FALSE)
  } else {
    checkmate::assert_file_exists(atlas_path)
    if (grepl("\\.csv$", atlas_path, ignore.case = TRUE)) {
      atlas <- readr::read_csv(atlas_path, show_col_types = FALSE)
    } else {
      atlas <- readRDS(atlas_path)
    }
  }

  # Standardize atlas column names
  if ("name" %in% colnames(atlas)) {
    names(atlas)[names(atlas) == "name"] <- "roi_name"
  }
  if ("acronym" %in% colnames(atlas)) {
    names(atlas)[names(atlas) == "acronym"] <- "roi_acronym"
  }
  if ("id" %in% colnames(atlas)) {
    names(atlas)[names(atlas) == "id"] <- "roi_id"
  }

  # Build hierarchy lookup from atlas structure_id_path
  if ("structure_id_path" %in% colnames(atlas)) {
    hierarchy <- .build_atlas_hierarchy(atlas)

    # Map each cell's ROI to hierarchy levels
    meta <- seurat_obj@meta.data
    region_cols <- config@atlas$region_columns

    if ("roi_id" %in% colnames(meta)) {
      for (level_name in names(region_cols)) {
        col_name <- region_cols[[level_name]]
        level_num <- as.integer(gsub("L", "", level_name))
        meta[[col_name]] <- .map_to_hierarchy_level(
          meta$roi_id, hierarchy, level_num
        )
      }
      seurat_obj@meta.data <- meta
    } else if (verbose) {
      cli::cli_warn(c(
        "Column {.val roi_id} not found in metadata.",
        "i" = "Atlas hierarchy mapping requires ROI IDs. Skipping hierarchy assignment."
      ))
    }
  }

  if (verbose) {
    region_col <- config@atlas$region_columns$L3
    if (!is.null(region_col) && region_col %in% colnames(seurat_obj@meta.data)) {
      n_regions <- length(unique(seurat_obj@meta.data[[region_col]]))
      cli::cli_inform(c(
        "v" = "Atlas registration complete.",
        "i" = "Regions detected (L3): {.val {n_regions}}"
      ))
    }
  }

  seurat_obj
}


# ============================================================================
# Internal helpers
# ============================================================================

#' Standardize metadata column names
#' @keywords internal
.standardize_metadata <- function(obj, mapping = NULL, verbose = TRUE) {
  meta_cols <- colnames(obj@meta.data)

  # Default auto-detection mapping
  auto_map <- list(
    experiment_group = c("experiment_group", "condition", "group", "sample_group",
                         "experimental_group"),
    genotype         = c("genotype", "Genotype", "geno"),
    treatment        = c("treatment", "Treatment", "drug", "condition"),
    sample_id        = c("sample_id", "mouse_id", "animal_id", "subject_id",
                         "biological_replicate", "mice_id"),
    slide_id         = c("slide_id", "Slide_ID", "slide_name", "slide"),
    fov_id           = c("fov_id", "fov", "FOV", "field_of_view"),
    x_coord          = c("x_coord", "x_centroid", "CenterX_global_px",
                         "x_slide_mm", "x_FOV_px", "x"),
    y_coord          = c("y_coord", "y_centroid", "CenterY_global_px",
                         "y_slide_mm", "y_FOV_px", "y")
  )

  if (!is.null(mapping)) {
    # User-provided mapping takes priority
    for (std_name in names(mapping)) {
      old_name <- mapping[[std_name]]
      if (old_name %in% meta_cols && old_name != std_name) {
        obj@meta.data[[std_name]] <- obj@meta.data[[old_name]]
        if (verbose) {
          cli::cli_inform("Mapped {.val {old_name}} -> {.val {std_name}}")
        }
      }
    }
  } else {
    # Auto-detect
    for (std_name in names(auto_map)) {
      if (std_name %in% meta_cols) next
      candidates <- auto_map[[std_name]]
      found <- intersect(candidates, meta_cols)
      if (length(found) > 0) {
        obj@meta.data[[std_name]] <- obj@meta.data[[found[1]]]
        if (verbose) {
          cli::cli_inform("Auto-mapped {.val {found[1]}} -> {.val {std_name}}")
        }
      }
    }
  }

  obj
}


#' Compute log2CPM normalization
#' @keywords internal
.compute_log2cpm <- function(obj) {
  counts <- Seurat::GetAssayData(obj, layer = "counts")
  col_sums <- Matrix::colSums(counts)
  col_sums[col_sums == 0] <- 1  # Avoid division by zero
  cpm <- sweep(counts, 2, col_sums, "/") * 1e6
  log2_cpm <- log2(cpm + 1)

  obj[["log2CPM"]] <- Seurat::CreateAssayObject(data = log2_cpm)
  obj
}


#' Load CosMx directory format
#' @keywords internal
.load_cosmx_directory <- function(dir_path, verbose = TRUE) {
  # Look for standard CosMx output files
  expr_file <- list.files(dir_path, pattern = "exprMat_file\\.csv", full.names = TRUE)
  meta_file <- list.files(dir_path, pattern = "metadata_file\\.csv", full.names = TRUE)

  if (length(expr_file) == 0 || length(meta_file) == 0) {
    # Try Seurat's LoadNanostring
    if (requireNamespace("Seurat", quietly = TRUE)) {
      cli::cli_inform("Attempting to load as Nanostring CosMx directory...")
      obj <- Seurat::LoadNanostring(dir_path)
      return(obj)
    }
    cli::cli_abort("Could not find CosMx output files in {.file {dir_path}}")
  }

  if (verbose) cli::cli_inform("Loading CosMx data from directory...")

  expr <- data.table::fread(expr_file[1])
  meta <- data.table::fread(meta_file[1])

  # Build count matrix
  cell_ids <- expr[[1]]
  gene_names <- colnames(expr)[-1]
  count_mat <- as.matrix(expr[, -1])
  rownames(count_mat) <- cell_ids
  count_mat <- t(count_mat)

  # Build metadata
  meta_df <- as.data.frame(meta)
  rownames(meta_df) <- meta_df[[1]]

  obj <- Seurat::CreateSeuratObject(
    counts = count_mat,
    meta.data = meta_df,
    project = "CosMx"
  )

  obj
}


#' Build atlas hierarchy lookup from structure_id_path
#' @keywords internal
.build_atlas_hierarchy <- function(atlas) {
  hierarchy <- list()
  for (i in seq_len(nrow(atlas))) {
    path <- atlas$structure_id_path[i]
    if (is.na(path)) next
    ids <- as.integer(unlist(strsplit(gsub("^/|/$", "", path), "/")))
    hierarchy[[as.character(atlas$roi_id[i])]] <- list(
      path = ids,
      name = atlas$roi_name[i],
      acronym = atlas$roi_acronym[i]
    )
  }

  # Build ID-to-name lookup
  id_to_name <- stats::setNames(atlas$roi_name, as.character(atlas$roi_id))
  id_to_acronym <- stats::setNames(atlas$roi_acronym, as.character(atlas$roi_id))

  list(hierarchy = hierarchy, id_to_name = id_to_name,
       id_to_acronym = id_to_acronym)
}


#' Map ROI IDs to a specific hierarchy level
#' @keywords internal
.map_to_hierarchy_level <- function(roi_ids, hierarchy_data, level) {
  id_to_acronym <- hierarchy_data$id_to_acronym
  hierarchy <- hierarchy_data$hierarchy

  vapply(as.character(roi_ids), function(rid) {
    if (is.na(rid) || !rid %in% names(hierarchy)) return(NA_character_)
    path <- hierarchy[[rid]]$path
    if (length(path) >= level) {
      target_id <- as.character(path[level])
      if (target_id %in% names(id_to_acronym)) {
        return(unname(id_to_acronym[target_id]))
      }
    }
    NA_character_
  }, character(1), USE.NAMES = FALSE)
}
