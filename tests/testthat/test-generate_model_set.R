test_that("generate_model_set builds a model set for a Gaussian response", {
  model.set <- fixture_cs1_model_set()

  expect_named(
    model.set,
    c("n.mods", "predictor.correlations", "mod.formula", "used.data", "test.fit", "included.vars")
  )
  expect_equal(model.set$n.mods, length(model.set$mod.formula))
  expect_true(model.set$n.mods > 0)
  expect_s3_class(model.set$mod.formula[[1]], "formula")
  expect_true("null" %in% names(model.set$mod.formula))
  expect_equal(model.set$included.vars, c("complexity", "depth", "ZONE"))
})

test_that("generate_model_set errors when use.dat is not a data.frame", {
  expect_error(
    generate_model_set(use.dat = list(a = 1), test.fit = NULL, pred.vars.cont = "a"),
    "data\\.frame"
  )
})

test_that("generate_model_set errors clearly when null.terms is not a single character string", {
  fit <- fixture_cs1_gaussian()

  bad.null.terms <- list(
    NA,
    NULL,
    123,
    TRUE,
    c("s(site,bs='re')", "s(other,bs='re')"),
    factor("s(site,bs='re')")
  )

  for (null.terms in bad.null.terms) {
    expect_error(
      fixture_cs1_model_set(fit = fit, null.terms = null.terms),
      "null.terms must be a single character string"
    )
  }
})

test_that("generate_model_set errors when predictors contain NA", {
  fit <- fixture_cs1_gaussian()
  fit$use.dat$depth[1] <- NA

  expect_error(
    fixture_cs1_model_set(fit = fit),
    "Predictor variables contain NA"
  )
})

test_that("generate_model_set reports a failed null model with the underlying error", {
  fit <- fixture_cs1_gaussian()

  expect_error(
    fixture_cs1_model_set(fit = fit, null.terms = "s(not_a_column,bs='re')"),
    "Null model not successfully fitted"
  )
})

test_that("generate_model_set errors when max.predictors exceeds the number of predictors", {
  fit <- fixture_cs1_gaussian()

  expect_error(
    fixture_cs1_model_set(
      fit = fit,
      pred.vars.cont = c("complexity", "depth"),
      pred.vars.fact = NA,
      max.predictors = 3
    ),
    "max.predictors is greater than the number of predictors"
  )
})

# ---- single-predictor model sets (regression) -------------------------------

test_that("generate_model_set builds a model set from a single predictor", {
  # Regression test: use.dat[, all.predictors] dropped to a bare vector when
  # there was exactly one predictor, so check_correlations() received something
  # that was not a data.frame and failed with "argument is of length zero"
  # rather than building a 1x1 correlation matrix.
  fit <- fixture_cs1_gaussian()

  cont.only <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = NA, max.predictors = 1
  )
  expect_equal(names(cont.only$mod.formula), c("null", "depth"))
  expect_equal(dim(cont.only$predictor.correlations), c(1, 1))
  expect_equal(unname(cont.only$predictor.correlations[1, 1]), 1)

  fact.only <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = NA, pred.vars.fact = "ZONE", max.predictors = 1
  )
  expect_equal(names(fact.only$mod.formula), c("null", "ZONE"))

  linear.only <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = NA, pred.vars.fact = NA,
    linear.vars = "depth", max.predictors = 1
  )
  expect_equal(names(linear.only$mod.formula), c("null", "depth"))
})

test_that("generate_model_set builds a single-predictor set with non.linear.correlations", {
  # The same regression, through check_non_linear_correlations(), which
  # additionally had no zero-pair case: a one-column input produced a zero-row
  # pair grid and failed in .rowNamesDF<-.
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = "depth", pred.vars.fact = NA, max.predictors = 1,
    non.linear.correlations = TRUE
  )

  expect_equal(names(model.set$mod.formula), c("null", "depth"))
  expect_equal(dim(model.set$predictor.correlations), c(1, 1))
})

test_that("generate_model_set builds no phantom by-terms when pred.vars.cont is NA", {
  # Regression test: pred.vars.cont = NA is the documented way to run without
  # smooth predictors, but expand.grid(NA, factors) paired the NA itself with
  # each factor, producing a phantom "NA.by.<factor>" candidate term. It was
  # discarded again before the model set was returned, but only after
  # enumerate_candidate_models() had taken max() of an empty correlation
  # sub-matrix and warned twice.
  fit <- fixture_cs1_gaussian()

  expect_no_warning(
    model.set <- fixture_cs1_model_set(
      fit = fit, pred.vars.cont = NA, pred.vars.fact = "ZONE", max.predictors = 1
    )
  )
  expect_false(any(grepl("NA.by.", names(model.set$mod.formula), fixed = TRUE)))
})
