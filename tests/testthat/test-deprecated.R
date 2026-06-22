test_that("generate.model.set warns and forwards to generate_model_set", {
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  args <- list(
    use.dat = use.dat,
    test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = "ZONE",
    null.terms = "s(site,bs='re')",
    max.predictors = 2,
    k = 3
  )

  expect_warning(old <- do.call(generate.model.set, args), class = "deprecatedWarning")
  new <- suppressWarnings(do.call(generate_model_set, args))

  expect_identical(names(old$mod.formula), names(new$mod.formula))
  expect_equal(old$n.mods, new$n.mods)
})

test_that("fit.model.set warns and forwards to fit_model_set", {
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

  expect_warning(old <- fit.model.set(model.set, parallel = FALSE), class = "deprecatedWarning")
  new <- fit_model_set(model.set, parallel = FALSE)

  expect_equal(old$mod.data.out$AICc, new$mod.data.out$AICc)
})

test_that("full.subsets.gam does not itself emit deprecation warnings", {
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  expect_no_warning(
    full.subsets.gam(
      use.dat = use.dat,
      test.fit = test.fit,
      pred.vars.cont = c("complexity", "depth"),
      pred.vars.fact = "ZONE",
      null.terms = "s(site,bs='re')",
      max.predictors = 2,
      k = 3
    )
  )
})
