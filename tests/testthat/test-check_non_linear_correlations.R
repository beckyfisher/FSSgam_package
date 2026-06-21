test_that("check.non.linear.correlations builds a matrix with a unit diagonal", {
  data(case_study1)
  cm <- check.non.linear.correlations(case_study1[, c("depth", "complexity", "ZONE")])

  expect_true(is.matrix(cm))
  expect_equal(dim(cm), c(3, 3))
  expect_equal(rownames(cm), c("depth", "complexity", "ZONE"))
  expect_equal(colnames(cm), c("depth", "complexity", "ZONE"))
  expect_equal(unname(diag(cm)), rep(1, 3))
})

test_that("check.non.linear.correlations is not necessarily symmetric", {
  # row = response, column = predictor, so off-diagonal entries can differ
  data(case_study1)
  cm <- check.non.linear.correlations(case_study1[, c("depth", "complexity")])
  expect_true(all(cm >= 0 & cm <= 1))
})

test_that("check.non.linear.correlations errors on an unsupported column class", {
  data(case_study1)
  bad.dat <- case_study1[, c("depth", "complexity")]
  bad.dat$listcol <- as.list(seq_len(nrow(bad.dat)))
  expect_error(check.non.linear.correlations(bad.dat), "not supported")
})
