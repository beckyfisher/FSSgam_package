test_that("generate_model_set builds factor-factor interaction terms when requested", {
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  use.dat$ZONE2 <- factor(ifelse(use.dat$depth > median(use.dat$depth), "deep", "shallow"))
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  model.set <- generate_model_set(
    use.dat = use.dat,
    test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = c("ZONE", "ZONE2"),
    factor.factor.interactions = TRUE,
    null.terms = "s(site,bs='re')",
    max.predictors = 2,
    k = 3
  )

  expect_true("ZONE.I.ZONE2" %in% colnames(model.set$used.data))
  expect_true("ZONE.I.ZONE2" %in% names(model.set$mod.formula))
  expect_true(model.set$n.mods > 0)
})

test_that("generate_model_set builds a te() smooth-smooth interaction for uncorrelated predictors", {
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  # depth and SCORE2 have near-zero correlation in case_study1, so the te()
  # interaction survives the cov.cutoff exclusion (default 0.28)
  model.set <- generate_model_set(
    use.dat = use.dat,
    test.fit = test.fit,
    pred.vars.cont = c("depth", "SCORE2"),
    smooth.smooth.interactions = TRUE,
    null.terms = "s(site,bs='re')",
    max.predictors = 2,
    k = 3
  )

  expect_true("depth.te.SCORE2" %in% names(model.set$mod.formula))
  te.formula <- model.set$mod.formula[["depth.te.SCORE2"]]
  expect_s3_class(te.formula, "formula")
  expect_match(as.character(te.formula)[2], "^te\\(")
})

test_that("smooth-smooth interactions are excluded when predictors exceed cov.cutoff", {
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  # complexity and depth are correlated above the default cov.cutoff (0.28),
  # so the te() interaction should be dropped from the candidate set
  model.set <- generate_model_set(
    use.dat = use.dat,
    test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = "ZONE",
    smooth.smooth.interactions = TRUE,
    null.terms = "s(site,bs='re')",
    max.predictors = 2,
    k = 3
  )

  expect_false("complexity.te.depth" %in% names(model.set$mod.formula))
})

test_that("generate_model_set uses check_non_linear_correlations when requested", {
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  model.set <- generate_model_set(
    use.dat = use.dat,
    test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = "ZONE",
    non.linear.correlations = TRUE,
    null.terms = "s(site,bs='re')",
    max.predictors = 2,
    k = 3
  )

  expected <- check_non_linear_correlations(use.dat[, c("complexity", "depth", "ZONE")])
  expect_equal(model.set$predictor.correlations, expected)
  # the non-linear correlation matrix is asymmetric, unlike the default
  # check_correlations() matrix
  expect_false(isTRUE(all.equal(
    model.set$predictor.correlations,
    t(model.set$predictor.correlations)
  )))
})

test_that("factor-factor interactions are excluded when the factors exceed cov.cutoff", {
  fit <- fixture_cs1_gaussian()
  # ZONE.copy is a relabelling of ZONE, so the two are perfectly correlated and
  # their interaction is dropped, leaving nothing to interact
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))

  expect_warning(
    model.set <- fixture_cs1_model_set(
      fit = fit,
      pred.vars.cont = "complexity",
      pred.vars.fact = c("ZONE", "ZONE.copy"),
      factor.factor.interactions = TRUE
    ),
    "no\\s+factors to interaction"
  )

  expect_false(any(grepl(".I.", names(model.set$mod.formula), fixed = TRUE)))
  # and the two collinear factors never appear in the same model either
  expect_false(any(grepl("ZONE\\+ZONE.copy", names(model.set$mod.formula))))
})

test_that("factor.factor.interactions as a character vector selects which factors interact", {
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )
  fit$use.dat$ZONE3 <- factor(
    ifelse(fit$use.dat$complexity > stats::median(fit$use.dat$complexity), "hi", "lo")
  )

  model.set <- fixture_cs1_model_set(
    fit = fit,
    pred.vars.cont = "depth",
    pred.vars.fact = c("ZONE", "ZONE2", "ZONE3"),
    factor.factor.interactions = c("ZONE", "ZONE2"),
    max.predictors = 2
  )

  expect_true("ZONE.I.ZONE2" %in% colnames(model.set$used.data))
  expect_true("ZONE.I.ZONE2" %in% names(model.set$mod.formula))
  # ZONE3 was not named, so it forms no interaction
  expect_false(any(grepl("ZONE3.I.", names(model.set$mod.formula), fixed = TRUE)))
  expect_false(any(grepl(".I.ZONE3", names(model.set$mod.formula), fixed = TRUE)))
})

test_that("factor.factor.interactions respects max.predictors when building combinations", {
  # with three factors and max.predictors = 3 the two- and three-way pasted
  # interactions are both generated
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )
  fit$use.dat$ZONE3 <- factor(
    ifelse(fit$use.dat$complexity > stats::median(fit$use.dat$complexity), "hi", "lo")
  )

  # suppressWarnings: a pasted interaction column is perfectly predicted by its
  # own components, and the multinom() summary that check_correlations() takes
  # the deviance from warns "NaNs produced" while computing standard errors it
  # then discards (issue #10). The correlations themselves are correct.
  model.set <- suppressWarnings(fixture_cs1_model_set(
    fit = fit,
    pred.vars.cont = "depth",
    pred.vars.fact = c("ZONE", "ZONE2", "ZONE3"),
    factor.factor.interactions = TRUE,
    max.predictors = 3
  ))

  interaction.cols <- grep(".I.", colnames(model.set$used.data), fixed = TRUE, value = TRUE)
  expect_true("ZONE.I.ZONE2" %in% interaction.cols)
  expect_true("ZONE.I.ZONE2.I.ZONE3" %in% interaction.cols)
})

test_that("smooth.smooth.interactions as a character vector selects which predictors interact", {
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("depth", "SCORE2", "complexity"),
    pred.vars.fact = NA,
    smooth.smooth.interactions = c("depth", "SCORE2"),
    max.predictors = 2
  )

  expect_true("depth.te.SCORE2" %in% names(model.set$mod.formula))
  # complexity was not named, so it forms no te() term
  expect_false(any(grepl("complexity.te.", names(model.set$mod.formula), fixed = TRUE)))
  expect_false(any(grepl(".te.complexity", names(model.set$mod.formula), fixed = TRUE)))
})

test_that("smooth.smooth.interactions uses non-linear correlations when asked", {
  # the te() candidates are screened against check_non_linear_correlations()
  # rather than check_correlations() when non.linear.correlations = TRUE
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("depth", "SCORE2"),
    pred.vars.fact = NA,
    smooth.smooth.interactions = TRUE,
    non.linear.correlations = TRUE,
    max.predictors = 2
  )

  expect_equal(
    model.set$predictor.correlations,
    check_non_linear_correlations(fixture_cs1_data()[, c("depth", "SCORE2")])
  )
  expect_true("depth.te.SCORE2" %in% names(model.set$mod.formula))
})

test_that("named smooth.smooth.interactions are screened with non-linear correlations too", {
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("depth", "SCORE2", "complexity"),
    pred.vars.fact = NA,
    smooth.smooth.interactions = c("depth", "SCORE2"),
    non.linear.correlations = TRUE,
    max.predictors = 2
  )

  expect_true("depth.te.SCORE2" %in% names(model.set$mod.formula))
})
