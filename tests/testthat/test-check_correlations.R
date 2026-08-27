test_that("check_correlations builds a symmetric matrix for continuous and factor predictors", {
  cm <- check_correlations(FSSgam::case_study1[, c("depth", "complexity", "ZONE")])

  expect_true(is.matrix(cm))
  expect_equal(dim(cm), c(3, 3))
  expect_equal(rownames(cm), c("depth", "complexity", "ZONE"))
  expect_equal(colnames(cm), c("depth", "complexity", "ZONE"))
  expect_equal(cm, t(cm), tolerance = 1e-6)
  # continuous-continuous correlation should match a plain cor() call
  expect_equal(
    unname(cm["depth", "complexity"]),
    unname(stats::cor(FSSgam::case_study1$depth, FSSgam::case_study1$complexity)),
    tolerance = 1e-6
  )
})

test_that("check_correlations handles a single continuous predictor", {
  cm <- check_correlations(FSSgam::case_study1[, "depth", drop = FALSE])
  expect_equal(dim(cm), c(1, 1))
  expect_equal(unname(cm[1, 1]), 1)
})

test_that("check_correlations errors on an unsupported column class", {
  bad.dat <- FSSgam::case_study1[, c("depth", "complexity")]
  bad.dat$listcol <- as.list(seq_len(nrow(bad.dat)))
  expect_error(check_correlations(bad.dat), "not supported")
})

test_that("check_correlations estimates factor-factor correlations via multinom", {
  dat <- fixture_cs1_data()
  dat$ZONE2 <- factor(
    ifelse(dat$SCORE2 > stats::median(dat$SCORE2), "high", "low")
  )
  cm <- check_correlations(dat[, c("depth", "ZONE", "ZONE2")])

  # the factor-factor block is the sqrt of the pseudo-R2 of a multinomial fit
  fit <- summary(nnet::multinom(dat$ZONE ~ dat$ZONE2, trace = FALSE))$deviance
  null.fit <- summary(nnet::multinom(dat$ZONE ~ 1, trace = FALSE))$deviance
  expected <- sqrt(1 - (fit / null.fit))

  expect_equal(unname(cm["ZONE", "ZONE2"]), expected, tolerance = 1e-6)
  # multinom() is fitted separately in each direction, so the two off-diagonal
  # entries agree only to the optimiser's convergence tolerance
  expect_equal(unname(cm["ZONE2", "ZONE"]), unname(cm["ZONE", "ZONE2"]), tolerance = 1e-3)
  # The factor-factor diagonal is a fitted multinomial pseudo-R2 rather than a
  # set constant, so it is 1 only to the optimiser's convergence tolerance.
  # Asserting equality to 1 would pin nnet's optimiser, not this package.
  expect_true(all(diag(cm)[c("ZONE", "ZONE2")] > 0.99))
})

test_that("check_correlations handles factor predictors with no continuous predictors", {
  dat <- fixture_cs1_data()
  dat$ZONE2 <- factor(
    ifelse(dat$SCORE2 > stats::median(dat$SCORE2), "high", "low")
  )
  cm <- check_correlations(dat[, c("ZONE", "ZONE2")])

  expect_equal(dim(cm), c(2, 2))
  expect_equal(rownames(cm), c("ZONE", "ZONE2"))
  expect_true(all(diag(cm) > 0.99))
})

test_that("check_correlations estimates factor-continuous correlations via lm", {
  dat <- fixture_cs1_data()
  cm <- check_correlations(dat[, c("depth", "ZONE")])

  expected <- sqrt(summary(stats::lm(dat$depth ~ factor(dat$ZONE)))$r.sq)
  expect_equal(unname(cm["depth", "ZONE"]), expected, tolerance = 1e-6)
  expect_equal(unname(cm["ZONE", "depth"]), expected, tolerance = 1e-6)
})

test_that("check_correlations handles a single factor predictor", {
  dat <- fixture_cs1_data()
  cm <- check_correlations(dat[, "ZONE", drop = FALSE])

  expect_equal(dim(cm), c(1, 1))
  expect_true(cm[1, 1] > 0.99)
})

test_that("check_correlations parallel = TRUE matches the sequential result", {
  skip_on_cran()
  skip_unless_parallel_opt_in()

  dat <- fixture_cs1_data()
  dat$ZONE2 <- factor(
    ifelse(dat$SCORE2 > stats::median(dat$SCORE2), "high", "low")
  )
  use.dat <- dat[, c("depth", "ZONE", "ZONE2")]

  sequential <- check_correlations(use.dat)
  parallel <- check_correlations(use.dat, parallel = TRUE, n.cores = 2)

  expect_equal(parallel, sequential, tolerance = 1e-8)
})

test_that("two exactly balanced factors are reported as uncorrelated", {
  # The factor-factor estimate is sqrt(1 - fit/null). When the predictor carries
  # no information at all the two deviances are identical, and the guarded
  # branch sets the estimate to exactly 0 rather than evaluating a square root
  # of a quantity that is zero only to the optimiser's tolerance.
  #
  # predf splits each level of respf exactly in half, so the fitted multinomial
  # model and the intercept-only model are the same model.
  dat <- data.frame(
    respf = factor(rep(c("a", "b"), each = 34)),
    predf = factor(rep(c("x", "y"), times = 34))
  )
  cm <- check_correlations(dat)

  expect_equal(unname(cm["respf", "predf"]), 0)
  expect_equal(unname(cm["predf", "respf"]), 0)
})
