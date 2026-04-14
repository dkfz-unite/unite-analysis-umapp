library(testthat)
library(readr)
library(jsonlite)

# Helper to run the full pipeline script
run_umap_pipeline <- function(data, metadata, options) {
  # Create temporary directory
  tmpdir <- tempdir()
  
  # Write test data files
  data_file <- file.path(tmpdir, "data.tsv")
  metadata_file <- file.path(tmpdir, "metadata.tsv")
  options_file <- file.path(tmpdir, "options.json")
  results_file <- file.path(tmpdir, "results.tsv")
  
  write_tsv(data, data_file)
  write_tsv(metadata, metadata_file)
  write(toJSON(options, null = "null", auto_unbox = TRUE), options_file)
  
  # Execute run.R as subprocess
  cmd <- sprintf("cd ../run && Rscript run.R %s", tmpdir)
  output <- system(cmd, intern = TRUE)
  exit_code <- attr(output, "status") %||% 0  # get exit code from output
  if (exit_code != 0) {
    stop(sprintf("run.R exited with code %d", exit_code))
  }
  
  # Check that output file was created
  if (!file.exists(results_file)) {
    stop(sprintf("results.tsv not created at %s", results_file))
  }
  
  # Read and return results
  read_tsv(results_file)
}

# Default options fixture
default_options <- function() {
  list(
    normalization_method = "median",
    normalization_log_offset = 1,
    imputation_method = "mindet",
    stratify_imputation_by_batch = FALSE,
    batch_correction_method = NULL,
    min_non_missing_fraction = 0.5,
    require_min_fraction_one_class = FALSE,
    feature_selection_method = NULL,
    feature_selection_n_features = 1000,
    umap_n_neighbors = 5,
    umap_metric = "euclidean",
    umap_random_state = 42,
    umap_min_dist = 0.1,
    umap_n_principal_components = 10
  )
}

# Test data fixtures
create_test_data <- function(n_samples = 10, n_features = 50) {
  data <- data.frame(
    feature = paste0("protein", 1:n_features),
    matrix(runif(n_samples * n_features, min = 1, max = 100), 
           nrow = n_features, ncol = n_samples)
  )
  colnames(data)[2:(n_samples + 1)] <- paste0("sample", 1:n_samples)
  data
}

create_test_metadata <- function(n_samples = 10) {
  data.frame(
    sample = paste0("sample", 1:n_samples),
    batch = rep(c(1, 2), length.out = n_samples),
    condition = rep(c("A", "B"), length.out = n_samples),
    # arbitrary metadata columns can be added here if needed but are ignored
    metadata1 = runif(n_samples, min = 0, max = 1),
    metadata2 = sample(c("A", "B"), n_samples, replace = TRUE)
  )
}

# Tests
test_that("n_neighbors greater than n_samples returns 2D array with no NAs", {
  options <- default_options()
  options$umap_n_neighbors <- 15  # > 5 samples
  
  data <- create_test_data(n_samples = 5, n_features = 20)
  metadata <- create_test_metadata(n_samples = 5)
  
  umap_coords <- run_umap_pipeline(data, metadata, options)
  
  expect_equal(nrow(umap_coords), 5)
  expect_equal(ncol(umap_coords), 3)
  expect_false(anyNA(umap_coords))
})

test_that("n_principal_components < n_variables returns 2D array with no NAs", {
  options <- default_options()
  options$umap_n_principal_components <- 10  # < 50 features
  
  data <- create_test_data(n_samples = 8, n_features = 50)
  metadata <- create_test_metadata(n_samples = 8)
  
  umap_coords <- run_umap_pipeline(data, metadata, options)
  
  expect_equal(nrow(umap_coords), 8)
  expect_equal(ncol(umap_coords), 3)
  expect_false(anyNA(umap_coords))
})

test_that("n_principal_components > n_variables returns 2D array with no NAs", {
  options <- default_options()
  options$umap_n_principal_components <- 50  # > 15 features
  
  data <- create_test_data(n_samples = 8, n_features = 15)
  metadata <- create_test_metadata(n_samples = 8)
  
  umap_coords <- run_umap_pipeline(data, metadata, options)
  
  expect_equal(nrow(umap_coords), 8)
  expect_equal(ncol(umap_coords), 3)
  expect_false(anyNA(umap_coords))
})

test_that("test edge case where feature_selection_n_features == 1", {
  options <- default_options()
  options$feature_selection_method <- "variance"
  options$feature_selection_n_features <- 1
  
  data <- create_test_data(n_samples = 8, n_features = 10)
  metadata <- create_test_metadata(n_samples = 8)
  
  umap_coords <- run_umap_pipeline(data, metadata, options)
  
  expect_equal(nrow(umap_coords), 8)
  expect_equal(ncol(umap_coords), 3)
  expect_false(anyNA(umap_coords))
})

# test an error is raised when any elements of condition column are missing when require_min_fraction_one_class is TRUE
test_that("error is raised when condition contains missing values and require_min_fraction_one_class is TRUE", {
  options <- default_options()
  options$require_min_fraction_one_class <- TRUE
  
  data <- create_test_data(n_samples = 8, n_features = 10)
  metadata <- create_test_metadata(n_samples = 8)
  metadata$condition[1] <- NA  # introduce missing value
  
  expect_error(run_umap_pipeline(data, metadata, options), "run.R exited with code 1")
})

# test that no error is raised when require_min_fraction_one_class is FALSE
test_that("no error is raised when condition contains missing values and require_min_fraction_one_class is FALSE", {
  options <- default_options()
  options$require_min_fraction_one_class <- FALSE
  
  data <- create_test_data(n_samples = 8, n_features = 10)
  metadata <- create_test_metadata(n_samples = 8)
  metadata$condition[1] <- NA  # introduce missing value
  # check it runs without error
  umap_coords <- run_umap_pipeline(data, metadata, options)
  expect_equal(nrow(umap_coords), 8)
  expect_equal(ncol(umap_coords), 3)
  expect_false(anyNA(umap_coords))
})

# test no error is raised when all values in condition column are missing and require_min_fraction_one_class is FALSE
test_that("no error is raised when all values in condition column are missing and require_min_fraction_one_class is FALSE", {
  options <- default_options()
  options$require_min_fraction_one_class <- FALSE

  data <- create_test_data(n_samples = 8, n_features = 10)
  metadata <- create_test_metadata(n_samples = 8)
  metadata$condition <- NA  # all values missing
  
  umap_coords <- run_umap_pipeline(data, metadata, options)
  expect_equal(nrow(umap_coords), 8)
  expect_equal(ncol(umap_coords), 3)
  expect_false(anyNA(umap_coords))
})

