test_that("full.subsets.gam returns a fully populated output list", {
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  out <- full.subsets.gam(
    use.dat = use.dat,
    test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = "ZONE",
    null.terms = "s(site,bs='re')",
    max.predictors = 2,
    k = 3
  )

  expect_named(
    out,
    c("mod.data.out", "used.data", "predictor.correlations", "failed.models",
      "success.models", "variable.importance")
  )
  # regression test: used.data/predictor.correlations were previously NULL
  # due to a field-name mismatch (model.set$use.dat/$cor.matrix instead of
  # $used.data/$predictor.correlations)
  expect_false(is.null(out$used.data))
  expect_equal(dim(out$used.data), dim(use.dat))
  expect_false(is.null(out$predictor.correlations))
  expect_equal(
    out$predictor.correlations,
    check.correlations(use.dat[, c("complexity", "depth", "ZONE")])
  )
})

test_that("full.subsets.gam matches a separate generate_model_set + fit_model_set call", {
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  args <- list(
    use.dat = use.dat, test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  combined <- full.subsets.gam(
    use.dat = args$use.dat, test.fit = args$test.fit,
    pred.vars.cont = args$pred.vars.cont, pred.vars.fact = args$pred.vars.fact,
    null.terms = args$null.terms, max.predictors = args$max.predictors, k = args$k
  )
  model.set <- generate_model_set(
    use.dat = args$use.dat, test.fit = args$test.fit,
    pred.vars.cont = args$pred.vars.cont, pred.vars.fact = args$pred.vars.fact,
    null.terms = args$null.terms, max.predictors = args$max.predictors, k = args$k
  )
  separate <- fit_model_set(model.set, parallel = FALSE)

  expect_equal(combined$mod.data.out$AICc, separate$mod.data.out$AICc)
})

test_that("full.subsets.gam's legacy factor.interactions argument forwards with a warning", {
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  expect_warning(
    legacy <- full.subsets.gam(
      use.dat = use.dat, test.fit = test.fit,
      pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
      factor.interactions = TRUE,
      null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
    ),
    "factor.factor.interactions"
  )
  current <- full.subsets.gam(
    use.dat = use.dat, test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    factor.factor.interactions = TRUE,
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  expect_equal(legacy$mod.data.out$AICc, current$mod.data.out$AICc)
})

test_that("full.subsets.gam's legacy smooth.interactions argument forwards with a warning", {
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  expect_warning(
    legacy <- full.subsets.gam(
      use.dat = use.dat, test.fit = test.fit,
      pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
      smooth.interactions = NA,
      null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
    ),
    "factor.smooth.interactions"
  )
  current <- full.subsets.gam(
    use.dat = use.dat, test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    factor.smooth.interactions = NA,
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  expect_equal(legacy$mod.data.out$AICc, current$mod.data.out$AICc)
  # confirm the .by. interaction models are actually suppressed
  expect_false(any(grepl("\\.by\\.", legacy$mod.data.out$modname)))
})

test_that("full.subsets.gam's legacy size argument constrains max.predictors", {
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  # regression test: size= used to be a no-op (it warned but never actually
  # fed its value into max.predictors)
  expect_warning(
    legacy <- full.subsets.gam(
      use.dat = use.dat, test.fit = test.fit,
      pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
      size = 1,
      null.terms = "s(site,bs='re')", k = 3
    ),
    "max.predictors"
  )
  current <- full.subsets.gam(
    use.dat = use.dat, test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    max.predictors = 1,
    null.terms = "s(site,bs='re')", k = 3
  )

  expect_equal(sort(legacy$mod.data.out$modname), sort(current$mod.data.out$modname))
  # only single-predictor models (plus null) should be present
  expect_false(any(grepl("\\+", legacy$mod.data.out$modname)))
})
