# Structural tests for the bundled datasets. Each dimension asserted here is the
# one recorded in that dataset's @format entry in R/data.R, so these tests fail
# if data-raw/case_study_datasets.R is re-run against changed source CSVs
# without the documentation being updated to match.

test_that("case_study1 matches its documented format", {
  dat <- FSSgam::case_study1

  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(68L, 27L))
  expect_true(all(c(
    "site", "ZONE", "depth", "complexity", "SCORE2",
    "Herbivore.abundance", "log.Herbivore.biomass"
  ) %in% names(dat)))
  expect_type(dat$depth, "double")
  expect_false(anyNA(dat$depth))
})

test_that("case_study2 matches its documented format", {
  dat <- FSSgam::case_study2

  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(285L, 21L))
  expect_true(all(c(
    "Location", "Status", "Site", "Distance", "depth",
    "snapper", "lobster", "Taxa", "Abundance"
  ) %in% names(dat)))
  # Distance is the vignette's linear.vars predictor, so it must stay numeric
  expect_type(dat$Distance, "integer")
})

test_that("case_study3 matches its documented format", {
  dat <- FSSgam::case_study3

  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(1110L, 12L))
  expect_true(all(c(
    "GSI", "Species", "Sex", "year", "month", "lunar.date"
  ) %in% names(dat)))
  # month and lunar.date are the cyclic predictors and must be numeric, not
  # factors, for bs = 'cc' to be applicable
  expect_type(dat$month, "integer")
  expect_type(dat$lunar.date, "integer")
  # GSI is modelled with Gamma(), which requires strictly positive values
  expect_true(all(dat$GSI > 0))
})

test_that("coral_data matches its documented format", {
  dat <- FSSgam::coral_data

  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(350L, 73L))
  expect_true(all(c(
    "Site", "Survey", "av.wave", "Depth", "allcoral", "totalpoints",
    "bleach.pres", "dredge.pres"
  ) %in% names(dat)))
  # the binomial fixture builds cbind(successes, failures) from these two
  expect_true(all(dat$allcoral <= dat$totalpoints))
})

test_that("every bundled dataset has columns of classes the correlation code supports", {
  # classify_correlation_predictors() accepts factor, character, integer and
  # numeric only, and errors on anything else. A dataset carrying, say, a list
  # or Date column would make check_correlations() unusable on it.
  supported <- c("factor", "character", "integer", "numeric")

  for (nm in c("case_study1", "case_study2", "case_study3", "coral_data")) {
    dat <- get(nm, envir = asNamespace("FSSgam"))
    classes <- vapply(dat, function(x) class(x)[1], character(1))
    expect_true(all(classes %in% supported), info = nm)
  }
})
