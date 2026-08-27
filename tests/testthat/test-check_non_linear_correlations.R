test_that("check_non_linear_correlations builds a matrix with a unit diagonal", {
  cm <- check_non_linear_correlations(FSSgam::case_study1[, c("depth", "complexity", "ZONE")])

  expect_true(is.matrix(cm))
  expect_equal(dim(cm), c(3, 3))
  expect_equal(rownames(cm), c("depth", "complexity", "ZONE"))
  expect_equal(colnames(cm), c("depth", "complexity", "ZONE"))
  expect_equal(unname(diag(cm)), rep(1, 3))
})

test_that("check_non_linear_correlations is not necessarily symmetric", {
  # row = response, column = predictor, so off-diagonal entries can differ
  cm <- check_non_linear_correlations(FSSgam::case_study1[, c("depth", "complexity")])
  expect_true(all(cm >= 0 & cm <= 1))
})

test_that("check_non_linear_correlations errors on an unsupported column class", {
  bad.dat <- FSSgam::case_study1[, c("depth", "complexity")]
  bad.dat$listcol <- as.list(seq_len(nrow(bad.dat)))
  expect_error(check_non_linear_correlations(bad.dat), "not supported")
})

test_that("check_non_linear_correlations returns the unit matrix for a single column", {
  # regression test: a one-column input produced a zero-row pair grid and failed
  # in .rowNamesDF<- rather than returning the 1x1 matrix check_correlations()
  # returns for the same input
  cm <- check_non_linear_correlations(FSSgam::case_study1[, "depth", drop = FALSE])

  expect_equal(dim(cm), c(1, 1))
  expect_equal(rownames(cm), "depth")
  expect_equal(unname(cm[1, 1]), 1)
})

test_that("continuous-continuous entries are the sqrt of a gam R squared", {
  dat <- FSSgam::case_study1[, c("depth", "complexity")]
  cm <- check_non_linear_correlations(dat)

  # row = response, column = predictor
  fit <- mgcv::gam(depth ~ s(complexity, k = 4), data = dat)
  expect_equal(unname(cm["depth", "complexity"]), sqrt(summary(fit)$r.sq),
                tolerance = 1e-6)

  # and the reverse direction is a different fit, hence a different value
  fit.rev <- mgcv::gam(complexity ~ s(depth, k = 4), data = dat)
  expect_equal(unname(cm["complexity", "depth"]), sqrt(summary(fit.rev)$r.sq),
                tolerance = 1e-6)
})

test_that("a continuous response on a factor predictor uses lm, not gam", {
  dat <- fixture_cs1_data()[, c("depth", "ZONE")]
  cm <- check_non_linear_correlations(dat)

  expected <- sqrt(summary(stats::lm(depth ~ ZONE, data = dat))$r.sq)
  expect_equal(unname(cm["depth", "ZONE"]), expected, tolerance = 1e-6)
})

test_that("a factor response uses a multinomial model", {
  dat <- fixture_cs1_data()[, c("ZONE", "depth")]
  cm <- check_non_linear_correlations(dat)

  fit <- summary(nnet::multinom(ZONE ~ depth, data = dat, trace = FALSE))$deviance
  null.fit <- summary(nnet::multinom(ZONE ~ 1, data = dat, trace = FALSE))$deviance
  expect_equal(unname(cm["ZONE", "depth"]), sqrt(1 - (fit / null.fit)),
                tolerance = 1e-6)
})

test_that("the gam basis dimension is reduced for predictors with few unique values", {
  # k.use starts at 4 and drops to the number of unique predictor values when
  # that is smaller, so a coarsely discretised predictor still yields an estimate
  # rather than failing on an over-specified basis
  dat <- fixture_cs1_data()[, c("depth", "SCORE2")]
  dat$coarse <- as.numeric(cut(dat$depth, breaks = 3, labels = FALSE))
  expect_length(unique(dat$coarse), 3)

  cm <- check_non_linear_correlations(dat[, c("SCORE2", "coarse")])

  expect_equal(dim(cm), c(2, 2))
  expect_false(anyNA(cm))
  expect_true(all(cm >= 0 & cm <= 1))

  expected <- sqrt(summary(mgcv::gam(SCORE2 ~ s(coarse, k = 3), data = dat))$r.sq)
  expect_equal(unname(cm["SCORE2", "coarse"]), expected, tolerance = 1e-6)
})

test_that("a factor response with an uninformative predictor is reported as uncorrelated", {
  # The multinomial branch estimates sqrt(1 - fit/null). When the predictor
  # carries no information the two deviances are identical, and the guarded
  # branch returns exactly 0 rather than a square root of a quantity that is
  # zero only to the optimiser's tolerance. predf splits each level of respf
  # exactly in half, so the two models are the same model.
  dat <- data.frame(
    respf = factor(rep(c("a", "b"), each = 34)),
    predf = factor(rep(c("x", "y"), times = 34))
  )
  cm <- check_non_linear_correlations(dat)

  expect_equal(unname(cm["respf", "predf"]), 0)
  expect_equal(unname(cm["predf", "respf"]), 0)
})
