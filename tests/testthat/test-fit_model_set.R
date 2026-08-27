test_that("fit_model_set fits and ranks a small Gaussian candidate set", {
  model.set <- fixture_cs1_model_set()

  out <- fit_quietly(model.set, parallel = FALSE)

  expect_named(out, c("mod.data.out", "failed.models", "success.models", "variable.importance"))
  expect_equal(nrow(out$mod.data.out), model.set$n.mods)
  expect_true(all(c("AICc", "BIC", "wi.AICc", "wi.BIC") %in% colnames(out$mod.data.out)))
  expect_length(out$failed.models, 0)
  expect_true(all(is.finite(out$mod.data.out$AICc)))
  # AICc weights for a model set should sum to 1
  expect_equal(sum(out$mod.data.out$wi.AICc), 1, tolerance = 1e-2)
})

test_that("fit_model_set works with a Tweedie family (non-Gaussian)", {
  fit <- fixture_cs1_tweedie()
  model.set <- fixture_cs1_model_set(fit = fit)

  out <- fit_quietly(model.set, parallel = FALSE)

  expect_equal(nrow(out$mod.data.out), model.set$n.mods)
  expect_length(out$failed.models, 0)
  expect_true(all(is.finite(out$mod.data.out$AICc)))
})

test_that("fit_model_set works when family is supplied as a list element rather than a literal call", {
  # Regression test for GitHub issue #10. Passing family = family.opts[[2]]
  # (rather than family = tw() written out directly) used to make every
  # candidate model fail under parallel = TRUE, because fit_mod_l() relied
  # on update() re-evaluating that *expression* wherever the refit itself
  # ran -- fine on the calling process, but a doSNOW worker never has
  # family.opts in scope. The parallel = TRUE crash itself isn't exercised
  # here (committed tests don't spin up a cluster, see CLAUDE.md Section 6
  # Phase 6b); this checks the same underlying fix sequentially: that
  # resolve_candidate_family() (R/utils.R) can resolve and independently
  # clone a family supplied this way at all, for every candidate, without
  # error.
  #
  # family.opts is assigned to the global environment (via <<-) rather
  # than left as a local variable of this test_that() block, to mirror how
  # a real top-level user script would have it in scope. mgcv::gam()
  # normalises away the environment a formula/call was originally written
  # in, so re-evaluating test.fit's family expression always resolves
  # variables via the global environment, never an intermediate calling
  # function's local frame -- true of both the pre-fix and post-fix
  # mechanisms, and not something this fix changes or needs to change.
  library(mgcv)
  use.dat <- FSSgam::case_study1
  use.dat$site <- as.factor(use.dat$site)
  family.opts <<- list(gaussian(), tw())
  on.exit(rm("family.opts", envir = globalenv()), add = TRUE)
  test.fit <- gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    family = family.opts[[2]],
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
  use.dat <- FSSgam::case_study1
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

test_that("fit_model_set fits a gamm4/uGamm model set with a random effect", {
  fit <- fixture_coral_ugamm()

  model.set <- generate_model_set(
    use.dat = fit$use.dat,
    test.fit = fit$test.fit,
    pred.vars.cont = c("av.wave", "Depth"),
    pred.vars.fact = "bleach.pres",
    max.predictors = 2,
    k = 4
  )
  out <- fit_quietly(model.set, parallel = FALSE)

  expect_equal(nrow(out$mod.data.out), model.set$n.mods)
  expect_length(out$failed.models, 0)
  expect_true(all(is.finite(out$mod.data.out$AICc)))
  expect_s3_class(out$success.models[[1]], "gamm4")
  # the random effect is supplied outside the formula in a uGamm call, so it is
  # carried into every candidate without appearing in null.terms. The grouping
  # factors are read from the merMod's flist slot rather than through
  # lme4::ranef(), to avoid a `::` call into a package this one does not declare
  # (lme4 reaches the tests only as a dependency of gamm4).
  expect_true(all(vapply(
    out$success.models, function(x) "Site" %in% names(x$mer@flist), logical(1)
  )))
})

test_that("fit_model_set fits a factor-heavy FSSgam::case_study2 model set with nested random effects", {
  fit <- fixture_cs2_tweedie()

  model.set <- generate_model_set(
    use.dat = fit$use.dat,
    test.fit = fit$test.fit,
    pred.vars.cont = c("lobster", "snapper"),
    pred.vars.fact = "Status",
    linear.vars = "Distance",
    null.terms = "s(Location,Site,bs='re')",
    max.predictors = 2,
    k = 3
  )
  out <- fit_quietly(model.set, parallel = FALSE)

  expect_length(out$failed.models, 0)
  expect_equal(nrow(out$mod.data.out), model.set$n.mods)
  expect_setequal(model.set$included.vars, c("lobster", "snapper", "Status", "Distance"))
  expect_true(all(is.finite(out$mod.data.out$AICc)))
})

test_that("fit_model_set orders mod.data.out to match the candidate set, not by AICc", {
  # mod.data.out is returned in candidate order; callers sort it themselves
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE)

  expect_equal(out$mod.data.out$modname, names(model.set$mod.formula))
  expect_equal(out$mod.data.out$modname[1], "null")
})

test_that("fit_model_set's inclusion columns mark exactly the predictors in each model", {
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE)
  mdo <- out$mod.data.out

  expect_equal(unname(unlist(mdo[mdo$modname == "null", model.set$included.vars])),
                c(0, 0, 0))
  expect_equal(unname(unlist(mdo[mdo$modname == "ZONE+depth", model.set$included.vars])),
                c(0, 1, 1))
  # a by-term counts as including both of its predictors
  expect_equal(
    unname(unlist(mdo[mdo$modname == "ZONE+depth.by.ZONE", model.set$included.vars])),
    c(0, 1, 1)
  )
})
