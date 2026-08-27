# Numerical regression tests.
#
# Everything else in the suite asserts on structure, relationships between
# outputs, or values recomputed from the same fitted object, so a refactor that
# quietly changed a fitted value would pass. These tests pin the actual numbers
# for three representative end-to-end scenarios so that any such change shows up
# as a reviewable diff in tests/testthat/_snaps/.
#
# Tolerance: mgcv's smoothing-parameter optimisation is not bit-identical across
# BLAS implementations, so the comparison is made at 1e-3 relative rather than
# exactly. That is also loose enough to absorb a last-digit flip in the
# delta/weight columns, which the package itself rounds to three decimal places.
# expect_snapshot_value() defaults to cran = FALSE, so these do not run on CRAN.
#
# Nothing snapshotted here embeds a path, a timestamp or an environment address.

# Reduces a fit_model_set() result to a stable, comparable summary.
snapshot_summary <- function(model.set, out) {
  mdo <- out$mod.data.out[order(out$mod.data.out$modname), ]
  list(
    formulas = vapply(
      model.set$mod.formula[sort(names(model.set$mod.formula))], deparse1, character(1)
    ),
    modname = mdo$modname,
    AICc = mdo$AICc,
    BIC = mdo$BIC,
    delta.AICc = mdo$delta.AICc,
    delta.BIC = mdo$delta.BIC,
    wi.AICc = mdo$wi.AICc,
    wi.BIC = mdo$wi.BIC,
    r2.vals = mdo$r2.vals,
    edf = mdo$edf,
    vi.aic = out$variable.importance$aic$variable.weights.raw,
    vi.bic = out$variable.importance$bic$variable.weights.raw
  )
}

test_that("the Gaussian case_study1 model set is numerically unchanged", {
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE, report.unique.r2 = TRUE)

  expect_snapshot_value(
    snapshot_summary(model.set, out), style = "json2", tolerance = 1e-3
  )
})

test_that("variable importance under VI.mods = 'all' is numerically unchanged", {
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE, VI.mods = "all")

  expect_snapshot_value(
    list(
      aic = out$variable.importance$aic$variable.weights.raw,
      bic = out$variable.importance$bic$variable.weights.raw
    ),
    style = "json2", tolerance = 1e-3
  )
})

test_that("the negative binomial case_study1 model set is numerically unchanged", {
  # an extended family, so this also pins the behaviour of the per-candidate
  # family resolution that issues #10 and #12 turned on
  fit <- fixture_cs1_nb()
  model.set <- fixture_cs1_model_set(fit = fit)
  out <- fit_quietly(model.set, parallel = FALSE)

  expect_snapshot_value(
    snapshot_summary(model.set, out), style = "json2", tolerance = 1e-3
  )
})

test_that("the binomial gamm4/uGamm model set is numerically unchanged", {
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

  expect_snapshot_value(
    snapshot_summary(model.set, out), style = "json2", tolerance = 1e-3
  )
})

test_that("the cyclic case_study3 model set is numerically unchanged", {
  fit <- fixture_cs3_cyclic()
  model.set <- generate_model_set(
    use.dat = fit$use.dat,
    test.fit = fit$test.fit,
    pred.vars.cont = c("lunar.date", "month"),
    pred.vars.fact = "Sex",
    cyclic.vars = c("lunar.date", "month"),
    max.predictors = 2,
    k = 5
  )
  out <- fit_quietly(model.set, parallel = FALSE)

  expect_snapshot_value(
    snapshot_summary(model.set, out), style = "json2", tolerance = 1e-3
  )
})
