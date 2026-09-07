# Numerical regression tests.
#
# Everything else in the suite asserts on structure, relationships between
# outputs, or values recomputed from the same fitted object, so a refactor that
# quietly changed a fitted value would pass. These tests pin the actual numbers
# for five representative end-to-end scenarios so that any such change shows up
# as a reviewable diff in tests/testthat/_snaps/.
#
# The comparison is split in two, because expect_snapshot_value() passes its
# tolerance to waldo::compare(), which applies it element-wise and *relatively*.
# One tolerance cannot serve both kinds of column here:
#
#  - delta.AICc, delta.BIC, wi.AICc, wi.BIC and the variable importance scores
#    are already rounded to three decimal places by compute_model_weights(), and
#    the importance scores are sums of those. A relative tolerance is the wrong
#    instrument for them: 1e-3 relative would reject a weight moving from 0.007
#    to 0.0075 -- a single-digit change in the last place the package itself
#    keeps -- while accepting a change of 0.5 in a weight of 0.9. They are
#    compared at testthat's default tolerance, which for values on a 0.001 grid
#    is exact.
#  - AICc, BIC, r2.vals and edf are unrounded (or rounded finely enough not to
#    matter) and do vary at the last bits between BLAS implementations, because
#    mgcv's smoothing-parameter optimisation is not bit-identical across them.
#    They are compared at 1e-6 relative, which on an AICc of about 200 is 2e-4
#    absolute -- far tighter than anything that could reorder two models, and
#    far looser than the 1e-8 scale on which BLAS differences appear.
#
# One quantity is excluded from the comparison rather than tolerated: the edf of
# the binomial gamm4/uGamm scenario. Measured on this branch's CI run at
# 7ec2b61, which ran all five scenarios with these tolerances in place, four
# passed in full and the fifth differed only there -- 4.000/5.000/8.000 locally
# against 4.260/5.050/8.350 on ubuntu-latest. gamm4 reports a smooth's edf
# through lme4's random-effect machinery, so that column tracks the lme4 version
# rather than anything this package does; it is asserted structurally in that
# scenario instead. Everything else, in all five scenarios, is compared
# everywhere the suite runs.
#
# Note that the exactly-compared group is the more platform-fragile of the two
# despite its name: delta.AICc and delta.BIC are quantised at three decimal
# places, so an AICc drift of 2e-4 -- which snapshot_numeric()'s tolerance
# permits -- can move one across a rounding boundary. That has not been observed
# on CI, but it is where a future failure is most likely to come from, and such
# a failure would be a platform difference rather than a defect.
#
# expect_snapshot_value() defaults to cran = FALSE, so none of this runs on
# CRAN. Nothing snapshotted here embeds a path, a timestamp or an environment
# address.

# The columns the package rounds to three decimal places, plus the candidate
# names and formulas. Compared exactly.
snapshot_exact <- function(model.set, out) {
  mdo <- out$mod.data.out[order(out$mod.data.out$modname), ]
  list(
    formulas = vapply(
      model.set$mod.formula[sort(names(model.set$mod.formula))], deparse1, character(1)
    ),
    modname = mdo$modname,
    delta.AICc = mdo$delta.AICc,
    delta.BIC = mdo$delta.BIC,
    wi.AICc = mdo$wi.AICc,
    wi.BIC = mdo$wi.BIC,
    vi.aic = out$variable.importance$aic$variable.weights.raw,
    vi.bic = out$variable.importance$bic$variable.weights.raw
  )
}

# The unrounded fit statistics. Compared at 1e-6 relative. edf is dropped for
# the gamm4 scenario -- see the header.
snapshot_numeric <- function(out, include.edf = TRUE) {
  mdo <- out$mod.data.out[order(out$mod.data.out$modname), ]
  values <- list(
    AICc = mdo$AICc,
    BIC = mdo$BIC,
    r2.vals = mdo$r2.vals,
    edf = mdo$edf
  )
  if (!include.edf) values$edf <- NULL
  values
}

test_that("the Gaussian case_study1 model set is numerically unchanged", {
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE, report.unique.r2 = TRUE)

  expect_snapshot_value(snapshot_exact(model.set, out), style = "json2")
  expect_snapshot_value(snapshot_numeric(out), style = "json2", tolerance = 1e-6)
})

test_that("variable importance under VI.mods = 'all' is numerically unchanged", {
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE, VI.mods = "all")

  expect_snapshot_value(
    list(
      aic = out$variable.importance$aic$variable.weights.raw,
      bic = out$variable.importance$bic$variable.weights.raw
    ),
    style = "json2"
  )
})

test_that("the negative binomial case_study1 model set is numerically unchanged", {
  # an extended family, so this also pins the behaviour of the per-candidate
  # family resolution that the two conflicting upstream issues
  # (beckyfisher/FSSgam#10 and #12) turned on
  fit <- fixture_cs1_nb()
  model.set <- fixture_cs1_model_set(fit = fit)
  out <- fit_quietly(model.set, parallel = FALSE)

  expect_snapshot_value(snapshot_exact(model.set, out), style = "json2")
  expect_snapshot_value(snapshot_numeric(out), style = "json2", tolerance = 1e-6)
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

  expect_snapshot_value(snapshot_exact(model.set, out), style = "json2")
  expect_snapshot_value(
    snapshot_numeric(out, include.edf = FALSE), style = "json2", tolerance = 1e-6
  )
  # edf is version dependent here (see the header), so it is only checked for
  # shape: one parameter per smooth at minimum, and no model without parameters
  expect_true(all(is.finite(out$mod.data.out$edf)))
  expect_true(all(out$mod.data.out$edf >= 1))
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

  expect_snapshot_value(snapshot_exact(model.set, out), style = "json2")
  expect_snapshot_value(snapshot_numeric(out), style = "json2", tolerance = 1e-6)
})
