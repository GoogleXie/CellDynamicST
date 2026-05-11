# ============================================================================
# Tests for Input Validation Utilities
# ============================================================================

test_that("validate_seurat rejects non-Seurat objects", {
  expect_error(
    validate_seurat(data.frame(x = 1)),
    "Expected a Seurat object"
  )
})

test_that("validate_seurat accepts valid Seurat objects", {
  obj <- create_test_seurat()
  expect_silent(validate_seurat(obj))
})

test_that("validate_seurat catches missing metadata columns", {
  obj <- create_test_seurat()
  expect_error(
    validate_seurat(obj, required_metadata = c("nonexistent_column")),
    "missing"
  )
})

test_that("validate_seurat checks for required assay", {
  obj <- create_test_seurat()
  expect_error(
    validate_seurat(obj, required_assay = "SCT"),
    "not found"
  )
})

test_that("validate_groups catches missing groups", {
  obj <- create_test_seurat()
  config <- create_test_config()

  # Remove one group from the data
  obj@meta.data$experiment_group[obj@meta.data$experiment_group == "KO MOR"] <- "WT SAL"

  expect_error(
    validate_groups(obj, config),
    "not found"
  )
})

test_that("require_package gives informative error for missing packages", {
  expect_error(
    require_package("nonexistent_package_xyz123"),
    "not installed"
  )
})
