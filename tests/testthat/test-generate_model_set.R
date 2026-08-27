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

test_that("generate_model_set returns exactly the elements it documents", {
  # Regression test: the @return block listed used.data, predictor.correlations
  # and a "generated.models" element that has never existed, and omitted three
  # that do. The same class of defect as the Phase 7 full_subsets_gam() one.
  model.set <- fixture_cs1_model_set()

  expect_named(
    model.set,
    c("n.mods", "predictor.correlations", "mod.formula", "used.data",
      "test.fit", "included.vars")
  )
  expect_type(model.set$n.mods, "integer")
  expect_true(is.matrix(model.set$predictor.correlations))
  expect_type(model.set$mod.formula, "list")
  expect_s3_class(model.set$used.data, "data.frame")
  expect_s3_class(model.set$test.fit, "gam")
  expect_type(model.set$included.vars, "character")
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

  # expect_no_warning() is the load-bearing assertion here: the phantom term was
  # always discarded before the model set was returned, so the expect_false()
  # below passed with the defect present too. Do not remove the warning check.
  expect_no_warning(
    model.set <- fixture_cs1_model_set(
      fit = fit, pred.vars.cont = NA, pred.vars.fact = "ZONE", max.predictors = 1
    )
  )
  expect_false(any(grepl("NA.by.", names(model.set$mod.formula), fixed = TRUE)))
})

# ---- cor.matrix -------------------------------------------------------------

test_that("generate_model_set uses a user-supplied cor.matrix verbatim", {
  fit <- fixture_cs1_gaussian()
  default.set <- fixture_cs1_model_set(fit = fit)

  # complexity and depth are correlated above the default cov.cutoff in
  # case_study1, so they never share a model. Supply a zero matrix instead and
  # the combined model must appear.
  zero.cor <- default.set$predictor.correlations
  zero.cor[] <- 0
  diag(zero.cor) <- 1

  supplied <- fixture_cs1_model_set(fit = fit, cor.matrix = zero.cor)

  expect_equal(supplied$predictor.correlations, zero.cor)
  expect_false("complexity+depth" %in% names(default.set$mod.formula))
  expect_true("complexity+depth" %in% names(supplied$mod.formula))
})

test_that("generate_model_set errors when a supplied cor.matrix omits a predictor", {
  fit <- fixture_cs1_gaussian()
  full.cor <- fixture_cs1_model_set(fit = fit)$predictor.correlations
  partial.cor <- full.cor[c("complexity", "depth"), c("complexity", "depth")]

  expect_error(
    fixture_cs1_model_set(fit = fit, cor.matrix = partial.cor),
    "Supplied cor.matrix is missing required predictors: ZONE"
  )
})

# ---- factor.smooth.interactions ---------------------------------------------

test_that("factor.smooth.interactions = NA suppresses all by-terms", {
  model.set <- fixture_cs1_model_set(factor.smooth.interactions = NA)

  expect_false(any(grepl(".by.", names(model.set$mod.formula), fixed = TRUE)))
  # Candidate names are built by sorting the term names, so they follow the
  # active collation. testthat forces LC_COLLATE=C for reproducible output, which
  # is why the factor sorts ahead of the lowercase continuous predictors here.
  expect_setequal(
    names(model.set$mod.formula),
    c("null", "complexity", "depth", "ZONE", "ZONE+complexity", "ZONE+depth")
  )
})

test_that("factor.smooth.interactions as a character vector selects which factors interact", {
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  model.set <- fixture_cs1_model_set(
    fit = fit,
    pred.vars.cont = "depth",
    pred.vars.fact = c("ZONE", "ZONE2"),
    factor.smooth.interactions = "ZONE"
  )

  by.terms <- grep(".by.", names(model.set$mod.formula), fixed = TRUE, value = TRUE)
  expect_true(length(by.terms) > 0)
  expect_true(all(grepl(".by.ZONE$", by.terms)))
  expect_false(any(grepl(".by.ZONE2", by.terms, fixed = TRUE)))
})

test_that("factor.smooth.interactions as a list gives per-predictor-type control", {
  # The list form names which factors, continuous predictors and linear
  # predictors take part, independently of pred.vars.fact/cont and linear.vars.
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = "ZONE",
    linear.vars = "SCORE2",
    factor.smooth.interactions = list(
      fact.vars = "ZONE", cont.vars = "complexity", linear.vars = "SCORE2"
    )
  )
  mod.names <- names(model.set$mod.formula)

  # complexity was named as a cont.var, so its by-term is present; depth was not
  expect_true("ZONE+complexity.by.ZONE" %in% mod.names)
  expect_false(any(grepl("depth.by.", mod.names, fixed = TRUE)))
  # SCORE2 was named as a linear.var, so it gets a .t. (product) interaction
  expect_true("SCORE2.t.ZONE+ZONE" %in% mod.names)
  expect_match(
    deparse1(model.set$mod.formula[["SCORE2.t.ZONE+ZONE"]]), "SCORE2 * ZONE",
    fixed = TRUE
  )
})

test_that("factor.smooth.interactions list errors on a variable that is not a predictor", {
  expect_error(
    fixture_cs1_model_set(
      factor.smooth.interactions = list(fact.vars = "ZONE", cont.vars = "not_a_column")
    ),
    "not_a_column"
  )
})

# ---- factor.factor.interactions ---------------------------------------------

test_that("factor.factor.interactions errors when fewer than two factors are named", {
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  expect_error(
    fixture_cs1_model_set(
      fit = fit, pred.vars.fact = c("ZONE", "ZONE2"),
      factor.factor.interactions = "ZONE"
    ),
    "less than 2 factors"
  )
})

test_that("factor.factor.interactions errors when a named factor is not in use.dat", {
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  expect_error(
    fixture_cs1_model_set(
      fit = fit, pred.vars.fact = c("ZONE", "ZONE2"),
      factor.factor.interactions = c("ZONE", "not_a_column")
    ),
    "Not all specified factor.factor.interactions are supplied in use.dat"
  )
})

# ---- smooth.smooth.interactions ---------------------------------------------

test_that("smooth.smooth.interactions errors when fewer than two predictors are named", {
  expect_error(
    fixture_cs1_model_set(smooth.smooth.interactions = "depth"),
    "less than 2 variables as smooth.smooth.interactions"
  )
})

test_that("smooth.smooth.interactions errors when a named predictor is not in use.dat", {
  expect_error(
    fixture_cs1_model_set(
      smooth.smooth.interactions = c("depth", "not_a_column")
    ),
    "Not all specified smooth.smooth.interactions are supplied in use.dat"
  )
})

test_that("smooth.smooth.interactions errors when TRUE with fewer than two continuous predictors", {
  expect_error(
    fixture_cs1_model_set(
      pred.vars.cont = "depth", pred.vars.fact = "ZONE",
      smooth.smooth.interactions = TRUE, max.predictors = 2
    ),
    "less than 2 continuous predictors"
  )
})
