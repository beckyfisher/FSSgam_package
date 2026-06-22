test_that("check_non_linear_correlations builds a matrix with a unit diagonal", {
  data(case_study1)
  cm <- check_non_linear_correlations(case_study1[, c("depth", "complexity", "ZONE")])

  expect_true(is.matrix(cm))
  expect_equal(dim(cm), c(3, 3))
  expect_equal(rownames(cm), c("depth", "complexity", "ZONE"))
  expect_equal(colnames(cm), c("depth", "complexity", "ZONE"))
  expect_equal(unname(diag(cm)), rep(1, 3))
})

test_that("check_non_linear_correlations is not necessarily symmetric", {
  # row = response, column = predictor, so off-diagonal entries can differ
  data(case_study1)
  cm <- check_non_linear_correlations(case_study1[, c("depth", "complexity")])
  expect_true(all(cm >= 0 & cm <= 1))
})

test_that("check_non_linear_correlations errors on an unsupported column class", {
  data(case_study1)
  bad.dat <- case_study1[, c("depth", "complexity")]
  bad.dat$listcol <- as.list(seq_len(nrow(bad.dat)))
  expect_error(check_non_linear_correlations(bad.dat), "not supported")
})
