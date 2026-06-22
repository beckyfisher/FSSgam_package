test_that("fit_model_set fits and ranks a small Gaussian candidate set", {
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

  out <- fit_model_set(model.set, parallel = FALSE)

  expect_named(out, c("mod.data.out", "failed.models", "success.models", "variable.importance"))
  expect_equal(nrow(out$mod.data.out), model.set$n.mods)
  expect_true(all(c("AICc", "BIC", "wi.AICc", "wi.BIC") %in% colnames(out$mod.data.out)))
  expect_length(out$failed.models, 0)
  expect_true(all(is.finite(out$mod.data.out$AICc)))
  # AICc weights for a model set should sum to 1
  expect_equal(sum(out$mod.data.out$wi.AICc), 1, tolerance = 1e-2)
})

test_that("fit_model_set works with a Tweedie family (non-Gaussian)", {
  library(mgcv)
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    family = tw(),
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

  out <- fit_model_set(model.set, parallel = FALSE)

  expect_equal(nrow(out$mod.data.out), model.set$n.mods)
  expect_length(out$failed.models, 0)
  expect_true(all(is.finite(out$mod.data.out$AICc)))
})

test_that("fit_model_set works with a negative binomial family (extended family with estimated theta)", {
  library(mgcv)
  data(case_study1)
  use.dat <- case_study1
  use.dat$site <- as.factor(use.dat$site)
  use.dat$Herbivore.abundance <- round(use.dat$Herbivore.abundance)
  test.fit <- gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    family = nb(),
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

  out <- fit_model_set(model.set, parallel = FALSE)

  expect_equal(nrow(out$mod.data.out), model.set$n.mods)
  expect_length(out$failed.models, 0)
  expect_true(all(is.finite(out$mod.data.out$AICc)))
})
