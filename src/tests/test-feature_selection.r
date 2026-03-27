library(testthat)
source(file.path(getwd(), "helpers", "feature_selection.r"))

test_that("returns all data when n_features >= ncol(data)", {
  data <- data.frame(a = 1:5, b = 6:10, c = 11:15)
  result <- feature_selection(data, "variance", 5)
  expect_equal(result, data)
})

test_that("returns all data for method 'none'", {
  data <- data.frame(a = 1:5, b = 6:10, c = 11:15)
  result <- feature_selection(data, "none", 2)
  expect_equal(result, data)
})

test_that("selects top n_features by variance", {
  data <- data.frame(
    high_var = c(1, 100, 1, 100, 1),
    low_var = c(5, 5, 5, 5, 5),
    med_var = c(1, 5, 1, 5, 1)
  )
  result <- feature_selection(data, "variance", 2)
  expect_equal(ncol(result), 2)
  expect_true("high_var" %in% colnames(result))
  expect_true("med_var" %in% colnames(result))
  expect_false("low_var" %in% colnames(result))
})

test_that("raises error for unsupported method", {
  data <- data.frame(a = 1:5, b = 6:10)
  expect_error(feature_selection(data, "unsupported", 1))
})

test_that("preserves rownames and colnames", {
  data <- data.frame(a = 1:3, b = 4:6, c = 7:9)
  rownames(data) <- c("s1", "s2", "s3")
  result <- feature_selection(data, "variance", 2)
  expect_equal(rownames(result), c("s1", "s2", "s3"))
  expect_equal(length(colnames(result)), 2)
})