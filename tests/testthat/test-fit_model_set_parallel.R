# parallel = TRUE for both fitting paths.
#
# The suite has deliberately avoided spinning up a real doSNOW cluster because
# of CRAN check time and flakiness, so everything here is skipped on CRAN. It is
# also skipped whenever FSSgam has been loaded by pkgload rather than installed:
# a doSNOW worker is a fresh R process that loads FSSgam with library() from the
# installed library path, so testing this path against a pkgload copy tests
# nothing (see CLAUDE.md Section 6, Phase 12). R CMD check installs the package
# first, which is where these tests are meant to run.
#
# skip_on_ci() as well: doSNOW clusters have been observed to hang indefinitely
# on a loaded machine (a trivial foreach() with .packages= hangs the same way, so
# it is not specific to this package), and a hang in CI costs the whole job where
# a skip costs nothing. Run these deliberately, against an installed package:
#
#   R CMD INSTALL . && NOT_CRAN=true Rscript -e 'testthat::test_local()'
#
# n.cores = 2 throughout, per CRAN's limit on parallel workers in checks.

test_that("parallel = TRUE reproduces the sequential fit with saved model fits", {
  skip_on_cran()
  skip_on_ci()
  skip_if_dev_loaded()

  model.set <- fixture_cs1_model_set()

  sequential <- fit_quietly(model.set, parallel = FALSE)
  parallel <- fit_quietly(model.set, parallel = TRUE, n.cores = 2)

  expect_equal(parallel$mod.data.out$modname, sequential$mod.data.out$modname)
  expect_equal(parallel$mod.data.out$AICc, sequential$mod.data.out$AICc, tolerance = 1e-8)
  expect_equal(parallel$mod.data.out$r2.vals, sequential$mod.data.out$r2.vals,
                tolerance = 1e-8)
  expect_equal(parallel$variable.importance, sequential$variable.importance,
                tolerance = 1e-8)
  expect_s3_class(parallel$success.models[[1]], "gam")
})

test_that("parallel = TRUE reproduces the sequential fit with save.model.fits = FALSE", {
  skip_on_cran()
  skip_on_ci()
  skip_if_dev_loaded()

  model.set <- fixture_cs1_model_set()

  sequential <- fit_quietly(model.set, parallel = FALSE, save.model.fits = FALSE)
  parallel <- fit_quietly(model.set, parallel = TRUE, n.cores = 2,
                           save.model.fits = FALSE)

  expect_equal(parallel$mod.data.out$AICc, sequential$mod.data.out$AICc, tolerance = 1e-8)
  expect_equal(parallel$mod.data.out$edf, sequential$mod.data.out$edf, tolerance = 1e-8)
  expect_s3_class(parallel$success.models[[1]], "formula")
})

test_that("parallel = TRUE keeps each candidate's extended-family state independent", {
  # The parallel half of the issue #10/#12 conflict. Every candidate's family is
  # resolved on the calling process, before makeCluster(), precisely so that a
  # worker never has to re-evaluate test.fit's family expression itself. A
  # Tweedie test.fit exercises that: its estimated power parameter lives in the
  # family object's own mutable environment.
  skip_on_cran()
  skip_on_ci()
  skip_if_dev_loaded()

  fit <- fixture_cs1_tweedie()
  model.set <- fixture_cs1_model_set(fit = fit)

  sequential <- fit_quietly(model.set, parallel = FALSE)
  parallel <- fit_quietly(model.set, parallel = TRUE, n.cores = 2)

  expect_length(parallel$failed.models, 0)
  expect_equal(parallel$mod.data.out$AICc, sequential$mod.data.out$AICc, tolerance = 1e-6)
})

test_that("parallel = TRUE works when family was supplied as a list element", {
  # Issue #10 proper: under parallel = TRUE every candidate used to fail with
  # "object 'family.opts' not found", because update() re-evaluated test.fit's
  # family expression on a worker that had never seen family.opts.
  #
  # family.opts is assigned to the global environment (via <<-) rather than left
  # local to this block, mirroring a real top-level user script -- mgcv::gam()
  # normalises away the environment a formula was written in, so the expression
  # always resolves through the global environment. See the matching comment in
  # test-fit_model_set.R.
  skip_on_cran()
  skip_on_ci()
  skip_if_dev_loaded()

  use.dat <- fixture_cs1_data()
  family.opts <<- list(gaussian(), tw())
  on.exit(rm("family.opts", envir = globalenv()), add = TRUE)
  test.fit <- gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    family = family.opts[[2]],
    data = use.dat
  )

  model.set <- fixture_cs1_model_set(
    fit = list(use.dat = use.dat, test.fit = test.fit)
  )
  out <- fit_quietly(model.set, parallel = TRUE, n.cores = 2)

  expect_length(out$failed.models, 0)
  expect_equal(nrow(out$mod.data.out), model.set$n.mods)
  expect_true(all(is.finite(out$mod.data.out$AICc)))
})

test_that("full_subsets_gam forwards parallel and n.cores", {
  skip_on_cran()
  skip_on_ci()
  skip_if_dev_loaded()

  fit <- fixture_cs1_gaussian()
  args <- list(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  sequential <- do.call(full_subsets_quietly, args)
  parallel <- do.call(full_subsets_quietly, c(args, list(parallel = TRUE, n.cores = 2)))

  expect_equal(parallel$mod.data.out$AICc, sequential$mod.data.out$AICc, tolerance = 1e-8)
})
