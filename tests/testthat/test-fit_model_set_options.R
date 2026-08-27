# Non-default arguments of fit_model_set(). The parallel = TRUE variants of the
# same paths live in test-fit_model_set_parallel.R.

# ---- save.model.fits --------------------------------------------------------

test_that("save.model.fits = FALSE reproduces the saved path's model table", {
  # fit_and_summarise_unsaved_models() is a separate implementation, kept
  # separate so fitted objects are never all held in memory at once. It has to
  # agree with the saved path on the same model set.
  model.set <- fixture_cs1_model_set()

  saved <- fit_quietly(model.set, parallel = FALSE, save.model.fits = TRUE)
  unsaved <- fit_quietly(model.set, parallel = FALSE, save.model.fits = FALSE)

  expect_equal(colnames(unsaved$mod.data.out), colnames(saved$mod.data.out))
  expect_equal(unsaved$mod.data.out$modname, saved$mod.data.out$modname)
  expect_equal(unsaved$mod.data.out$AICc, saved$mod.data.out$AICc, tolerance = 1e-8)
  expect_equal(unsaved$mod.data.out$BIC, saved$mod.data.out$BIC, tolerance = 1e-8)
  expect_equal(unsaved$mod.data.out$r2.vals, saved$mod.data.out$r2.vals, tolerance = 1e-8)
  expect_equal(unsaved$mod.data.out$edf, saved$mod.data.out$edf, tolerance = 1e-8)
  expect_equal(unsaved$mod.data.out$wi.AICc, saved$mod.data.out$wi.AICc, tolerance = 1e-8)
  expect_equal(unsaved$variable.importance, saved$variable.importance, tolerance = 1e-8)
})

test_that("save.model.fits = FALSE returns formulas rather than fitted models", {
  model.set <- fixture_cs1_model_set()

  saved <- fit_quietly(model.set, parallel = FALSE, save.model.fits = TRUE)
  unsaved <- fit_quietly(model.set, parallel = FALSE, save.model.fits = FALSE)

  expect_s3_class(saved$success.models[[1]], "gam")
  expect_s3_class(unsaved$success.models[[1]], "formula")
  expect_named(unsaved$success.models, names(saved$success.models))
})

test_that("save.model.fits = FALSE reports failed models by formula", {
  model.set <- break_one_candidate(fixture_cs1_model_set())

  unsaved <- fit_quietly(model.set, parallel = FALSE, save.model.fits = FALSE)

  # in the unsaved path a failed model keeps its row, with NA statistics
  failed <- is.na(unsaved$mod.data.out$AICc)
  expect_equal(sum(failed), 1)
  expect_length(unsaved$failed.models, 1)
  expect_length(unsaved$success.models, model.set$n.mods - 1)
  expect_equal(nrow(unsaved$mod.data.out), model.set$n.mods)
})

# ---- max.models -------------------------------------------------------------

test_that("exceeding max.models warns and forces save.model.fits = FALSE", {
  model.set <- fixture_cs1_model_set()
  expect_true(model.set$n.mods > 2)

  expect_warning(
    capped <- fit_quietly(model.set, max.models = 2, parallel = FALSE),
    "Individual models fits will not be saved"
  )

  # the downgrade is silent in the return value except through success.models,
  # which now holds formulas rather than fitted model objects
  expect_s3_class(capped$success.models[[1]], "formula")
  uncapped <- fit_quietly(model.set, parallel = FALSE, save.model.fits = FALSE)
  expect_equal(capped$mod.data.out$AICc, uncapped$mod.data.out$AICc, tolerance = 1e-8)
})

test_that("a model set at or below max.models keeps its fitted models", {
  model.set <- fixture_cs1_model_set()

  expect_no_warning(
    out <- fit_quietly(model.set, max.models = model.set$n.mods, parallel = FALSE)
  )
  expect_s3_class(out$success.models[[1]], "gam")
})

# ---- r2.type ----------------------------------------------------------------

test_that("r2.type = 'r2.lm.est' is the default and reaches the model table", {
  # What r2.lm.est actually is -- the R squared of the response against the
  # link-scale prediction, which is neither of the quantities mgcv reports -- is
  # pinned in test-functions_supporting.R against an independently derived
  # value. This checks only that fit_model_set() defaults to it and puts it in
  # the table, which is the part that belongs here.
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE)
  named <- fit_quietly(model.set, parallel = FALSE, r2.type = "r2.lm.est")

  expect_equal(out$mod.data.out$r2.vals, named$mod.data.out$r2.vals)
  expect_false(anyNA(out$mod.data.out$r2.vals))
  # and it is not either of the values the other two r2.type settings give
  expect_false(isTRUE(all.equal(
    out$mod.data.out$r2.vals,
    fit_quietly(model.set, parallel = FALSE, r2.type = "r2")$mod.data.out$r2.vals
  )))
  expect_false(isTRUE(all.equal(
    out$mod.data.out$r2.vals,
    fit_quietly(model.set, parallel = FALSE, r2.type = "dev")$mod.data.out$r2.vals
  )))
})

test_that("r2.type = 'r2' reports mgcv's adjusted R squared", {
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE, r2.type = "r2")

  expected <- vapply(
    out$success.models, function(x) round(summary(x)$r.sq, 5), numeric(1)
  )
  expect_equal(out$mod.data.out$r2.vals, unname(expected))
})

test_that("r2.type = 'dev' reports mgcv's deviance explained", {
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE, r2.type = "dev")

  expected <- vapply(
    out$success.models, function(x) round(summary(x)$dev.expl, 5), numeric(1)
  )
  expect_equal(out$mod.data.out$r2.vals, unname(expected))
})

# ---- report.unique.r2 -------------------------------------------------------

test_that("report.unique.r2 = TRUE adds r2.vals.unique as r2 less the null model r2", {
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE, report.unique.r2 = TRUE)

  null.r2 <- out$mod.data.out$r2.vals[out$mod.data.out$modname == "null"]
  expect_true(null.r2 > 0)
  expect_equal(out$mod.data.out$r2.vals.unique, out$mod.data.out$r2.vals - null.r2)
  expect_equal(
    out$mod.data.out$r2.vals.unique[out$mod.data.out$modname == "null"], 0
  )
})

test_that("report.unique.r2 = FALSE leaves r2.vals.unique as the all-NA placeholder", {
  # extract_mod_dat() always returns an r2.vals.unique slot; it stays NA unless
  # compute_model_weights() is asked to fill it in
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE, report.unique.r2 = FALSE)

  expect_true("r2.vals.unique" %in% colnames(out$mod.data.out))
  expect_true(all(is.na(out$mod.data.out$r2.vals.unique)))
})

# ---- VI.mods ----------------------------------------------------------------

test_that("VI.mods = 'all' sums weights over every model containing a predictor", {
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE, VI.mods = "all")

  mdo <- out$mod.data.out
  for (v in model.set$included.vars) {
    expect_equal(
      unname(out$variable.importance$aic$variable.weights.raw[v]),
      sum(mdo$wi.AICc[mdo[, v] == 1]),
      tolerance = 1e-8,
      info = v
    )
    expect_equal(
      unname(out$variable.importance$bic$variable.weights.raw[v]),
      sum(mdo$wi.BIC[mdo[, v] == 1]),
      tolerance = 1e-8,
      info = v
    )
  }
})

test_that("VI.mods = 'min.n' (the default) sums only the best n models per predictor", {
  model.set <- fixture_cs1_model_set()
  out <- fit_quietly(model.set, parallel = FALSE)

  mdo <- out$mod.data.out
  min.mods <- min(colSums(mdo[, model.set$included.vars]))
  for (v in model.set$included.vars) {
    expect_equal(
      unname(out$variable.importance$aic$variable.weights.raw[v]),
      sum(sort(mdo$wi.AICc[mdo[, v] == 1], decreasing = TRUE)[seq_len(min.mods)]),
      tolerance = 1e-8,
      info = v
    )
  }
})

# ---- failure handling -------------------------------------------------------

test_that("fit_model_set errors when no model in the set can be fitted", {
  model.set <- fixture_cs1_model_set()
  # an all-NA response makes every candidate, including the null model, fail
  model.set$used.data$log.Herbivore.biomass <- NA_real_

  expect_error(
    fit_quietly(model.set, parallel = FALSE),
    "None of your models fitted successfully"
  )
})

test_that("fit_model_set reports partially failing model sets in failed.models", {
  model.set <- break_one_candidate(fixture_cs1_model_set())

  out <- fit_quietly(model.set, parallel = FALSE)

  expect_length(out$failed.models, 1)
  expect_named(out$failed.models, "ZONE+depth")
  expect_equal(
    length(out$failed.models) + length(out$success.models), model.set$n.mods
  )
  expect_true(all(vapply(out$failed.models, inherits, logical(1), "formula")))
  # in the saved path the failed model is absent from the table, not NA in it
  expect_equal(nrow(out$mod.data.out), length(out$success.models))
  expect_false("ZONE+depth" %in% out$mod.data.out$modname)
})
