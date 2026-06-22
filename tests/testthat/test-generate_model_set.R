test_that("generate_model_set builds a model set for a Gaussian response", {
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
    null.terms = "s(site,bs='re')",
    max.predictors = 2,
    k = 3
  )

  expect_named(
    model.set,
    c("n.mods", "predictor.correlations", "mod.formula", "used.data", "test.fit", "included.vars")
  )
  expect_equal(model.set$n.mods, length(model.set$mod.formula))
  expect_true(model.set$n.mods > 0)
  expect_s3_class(model.set$mod.formula[[1]], "formula")
  expect_true("null" %in% names(model.set$mod.formula))
})

test_that("generate_model_set errors when use.dat is not a data.frame", {
  expect_error(
    generate_model_set(use.dat = list(a = 1), test.fit = NULL, pred.vars.cont = "a"),
    "data\\.frame"
  )
})

test_that("generate_model_set errors when predictors contain NA", {
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  use.dat$depth[1] <- NA
  expect_error(
    generate_model_set(
      use.dat = use.dat,
      test.fit = test.fit,
      pred.vars.cont = "depth",
      null.terms = "s(site,bs='re')",
      max.predictors = 1,
      k = 3
    ),
    "NA"
  )
})
