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
