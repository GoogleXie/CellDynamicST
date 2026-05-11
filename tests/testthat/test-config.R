# ============================================================================
# Tests for Configuration Loading and Validation
# ============================================================================

test_that("cdst_load_config loads from YAML template", {
  config_file <- system.file("templates", "experiment_config.yaml",
                             package = "CellDynamicST")
  skip_if(!file.exists(config_file), "Template config file not found")

  config <- cdst_load_config(config_file, verbose = FALSE)

  expect_s4_class(config, "CdstConfig")
  expect_equal(config@species, "mouse")
  expect_true(length(config@groups) >= 1)
  expect_true(nchar(config@reference_group) > 0)
})

test_that("cdst_default_config creates valid config", {
  config <- cdst_default_config()

  expect_s4_class(config, "CdstConfig")
  expect_equal(config@project_name, "CellDynamicST_Project")
  expect_equal(config@species, "mouse")
  expect_equal(length(config@groups), 2)
  expect_equal(config@reference_group, "Control")
  expect_equal(length(config@comparisons), 1)
})

test_that("cdst_validate_config catches invalid configs", {
  config <- create_test_config()

  # Valid config should pass

  expect_silent(cdst_validate_config(config))

  # Invalid: non-existent comparison group
  bad_config <- config
  bad_config@comparisons <- list(
    list(name = "bad", group1 = "NONEXISTENT", group2 = "WT SAL")
  )
  expect_error(cdst_validate_config(bad_config), "unknown group1")
})

test_that("CdstConfig show method works", {
  config <- create_test_config()
  expect_output(show(config), "CellDynamicST Configuration")
})

test_that("CdstConfig validation catches bad species", {
  expect_error(
    new("CdstConfig",
      project_name = "test",
      species = "fish",
      groups = list(list(name = "A", genotype = "WT", treatment = "SAL",
                         is_reference = TRUE)),
      reference_group = "A",
      comparisons = list(),
      qc = list(), clustering = list(), wgcna = list(),
      atlas = list(), output_dir = "results"
    ),
    "mouse.*human"
  )
})

test_that("auto-generated comparisons work for multiple groups", {
  yaml_content <- "
project_name: test
groups:
  - name: A
    genotype: WT
    treatment: SAL
    is_reference: true
  - name: B
    genotype: WT
    treatment: MOR
  - name: C
    genotype: KO
    treatment: SAL
"
  tmp <- tempfile(fileext = ".yaml")
  writeLines(yaml_content, tmp)

  config <- cdst_load_config(tmp, verbose = FALSE)

  # 3 groups should produce 3 pairwise comparisons
  expect_equal(length(config@comparisons), 3)

  unlink(tmp)
})
