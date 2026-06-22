test_that("check_correlations builds a symmetric matrix for continuous and factor predictors", {
  data(case_study1)
  cm <- check_correlations(case_study1[, c("depth", "complexity", "ZONE")])

  expect_true(is.matrix(cm))
  expect_equal(dim(cm), c(3, 3))
  expect_equal(rownames(cm), c("depth", "complexity", "ZONE"))
  expect_equal(colnames(cm), c("depth", "complexity", "ZONE"))
  expect_equal(cm, t(cm), tolerance = 1e-6)
  # continuous-continuous correlation should match a plain cor() call
  expect_equal(
    unname(cm["depth", "complexity"]),
    unname(stats::cor(case_study1$depth, case_study1$complexity)),
    tolerance = 1e-6
  )
})

test_that("check_correlations handles a single continuous predictor", {
  data(case_study1)
  cm <- check_correlations(case_study1[, "depth", drop = FALSE])
  expect_equal(dim(cm), c(1, 1))
  expect_equal(unname(cm[1, 1]), 1)
})

test_that("check_correlations errors on an unsupported column class", {
  data(case_study1)
  bad.dat <- case_study1[, c("depth", "complexity")]
  bad.dat$listcol <- as.list(seq_len(nrow(bad.dat)))
  expect_error(check_correlations(bad.dat), "not supported")
})
