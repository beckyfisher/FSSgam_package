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

  # the whole table, not a selection of columns: the delta, weight, r2.vals.unique,
  # edf.less.1 and variable-inclusion columns are the ones a refactor of the
  # separate unsaved implementation is most likely to get wrong
  expect_equal(unsaved$mod.data.out, saved$mod.data.out, tolerance = 1e-8)
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

test_that("an unrecognised r2.type is an error, not a silent column of NA", {
  # extract_mod_dat() matches r2.type. against three literals and leaves the
  # value at NA when none matches, so r2.vals came back all-NA with no message.
  model.set <- fixture_cs1_model_set()

  expect_error(
    fit_quietly(model.set, parallel = FALSE, r2.type = "rsq"),
    "'arg' should be one of"
  )

  fit <- fixture_cs1_gaussian()
  expect_error(
    full_subsets_quietly(
      use.dat = fit$use.dat, test.fit = fit$test.fit,
      pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
      null.terms = "s(site,bs='re')", max.predictors = 2, k = 3,
      r2.type = "rsq"
    ),
    "'arg' should be one of"
  )
})

test_that("each documented r2.type is accepted and reaches the model table", {
  # "r2" is an exact match against the second element, so match.arg()'s partial
  # matching does not make it ambiguous with "r2.lm.est".
  model.set <- fixture_cs1_model_set()
  for (rt in c("r2.lm.est", "r2", "dev")) {
    out <- fit_quietly(model.set, parallel = FALSE, r2.type = rt)
    expect_true(all(is.finite(out$mod.data.out$r2.vals)), info = rt)
  }
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

test_that("an unrecognised VI.mods is an error, not a missing object", {
  # compute_variable_importance() had two if blocks and no else, so a value
  # matching neither assigned nothing and the failure surfaced as
  # "object 'aic.var.weights' not found".
  model.set <- fixture_cs1_model_set()

  expect_error(
    fit_quietly(model.set, parallel = FALSE, VI.mods = "alll"),
    "'arg' should be one of"
  )
})

# ---- failure handling -------------------------------------------------------

test_that("normalise_mod_dat_rows replaces anything that is not a summary row", {
  # The unsaved parallel path runs foreach() with .errorhandling = 'pass', so a
  # failing element comes back as the condition object rather than as the named
  # vector unlist(extract_mod_dat(...)) returns. rbind() over the mixed list
  # would produce a malformed table.
  good <- unlist(list(AICc = 1, BIC = 2, r2.vals = 0.5, r2.vals.unique = NA,
                      edf = 3, edf.less.1 = 0))
  na.row <- unlist(list(AICc = NA, BIC = NA, r2.vals = NA, r2.vals.unique = NA,
                        edf = NA, edf.less.1 = NA))

  out <- normalise_mod_dat_rows(list(
    a = good,
    b = simpleError("worker failed"),
    c = na.row,
    d = NULL,
    e = 1:3
  ))

  expect_named(out, c("a", "b", "c", "d", "e"))
  expect_identical(out$a, good)
  # already an all-NA row: passed through, not replaced
  expect_identical(out$c, na.row)
  for (nm in c("b", "d", "e")) {
    expect_identical(out[[nm]], na.row, info = nm)
  }

  # the point of the replacement: rbind() now yields one row per element
  tab <- do.call("rbind", out)
  expect_equal(dim(tab), c(5L, 6L))
  expect_true(is.numeric(tab))
  expect_equal(unname(tab[1, "AICc"]), 1)
  expect_true(all(is.na(tab[c(2, 3, 4, 5), "AICc"])))
})

test_that("normalise_mod_dat_rows leaves a fully successful list untouched", {
  rows <- replicate(3, unlist(list(AICc = 1, BIC = 2, r2.vals = 0.5,
                                   r2.vals.unique = NA, edf = 3,
                                   edf.less.1 = 0)), simplify = FALSE)
  expect_identical(normalise_mod_dat_rows(rows), rows)
  expect_identical(normalise_mod_dat_rows(list()), list())
})

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

# ---- single-predictor model sets --------------------------------------------

test_that("a single-predictor model set reaches variable importance", {
  # compute_variable_importance() indexed mod.data.out with a single column
  # name and no drop = FALSE, so colSums() received a vector and errored.
  # generate_model_set() was fixed for the single-predictor case in Phase 13;
  # fit_model_set() was not.
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = "depth", pred.vars.fact = NA, max.predictors = 1
  )
  expect_equal(model.set$included.vars, "depth")

  for (vi in c("min.n", "all")) {
    out <- fit_quietly(model.set, parallel = FALSE, VI.mods = vi)
    for (ic in c("aic", "bic")) {
      w <- out$variable.importance[[ic]]$variable.weights.raw
      expect_length(w, 1)
      expect_named(w, "depth")
      expect_true(is.finite(w), info = paste(vi, ic))
    }
  }
})

test_that("a predictor present in no surviving model contributes no weight", {
  # min.mods is 0 in that case, and 1:min.mods is c(1, 0), which selects the
  # single largest weight instead of none. seq_len() is empty at 0.
  model.set <- fixture_cs1_model_set(max.predictors = 1)
  model.set <- break_one_candidate(model.set, modname = "ZONE")

  out <- fit_quietly(model.set, parallel = FALSE)

  expect_false("ZONE" %in% out$mod.data.out$modname)
  expect_equal(min(colSums(out$mod.data.out[, model.set$included.vars])), 0)
  expect_equal(
    unname(out$variable.importance$aic$variable.weights.raw),
    rep(0, length(model.set$included.vars))
  )
  expect_equal(
    unname(out$variable.importance$bic$variable.weights.raw),
    rep(0, length(model.set$included.vars))
  )
})

# ---- progress ---------------------------------------------------------------

test_that("progress = FALSE writes nothing to stdout", {
  # FSSgam_package#9: the txtProgressBar was unconditional, so every call wrote
  # to stdout and the only way to keep it out of a report or a test reporter
  # was to wrap the call in capture.output().
  model.set <- fixture_cs1_model_set()

  for (save in c(TRUE, FALSE)) {
    out <- utils::capture.output(
      res <- fit_model_set(model.set, parallel = FALSE, progress = FALSE,
                           save.model.fits = save)
    )
    expect_equal(out, character(0), info = paste("save.model.fits =", save))
    # the fit itself is unaffected
    expect_equal(nrow(res$mod.data.out), model.set$n.mods, info = paste(save))
  }
})

test_that("progress = TRUE writes a progress bar to stdout", {
  model.set <- fixture_cs1_model_set()

  for (save in c(TRUE, FALSE)) {
    out <- utils::capture.output(
      fit_model_set(model.set, parallel = FALSE, progress = TRUE,
                    save.model.fits = save)
    )
    expect_true(any(grepl("=", out, fixed = TRUE)),
                info = paste("save.model.fits =", save))
  }
})

test_that("progress defaults to interactive()", {
  # so the bar appears at the console but not in scripts, reports or checks
  expect_identical(formals(fit_model_set)$progress, quote(interactive()))
  expect_identical(formals(full_subsets_gam)$progress, quote(interactive()))
})

test_that("full_subsets_gam forwards progress", {
  fit <- fixture_cs1_gaussian()
  args <- list(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  quiet <- utils::capture.output(
    res <- do.call(full_subsets_gam, c(args, list(progress = FALSE)))
  )
  noisy <- utils::capture.output(
    do.call(full_subsets_gam, c(args, list(progress = TRUE)))
  )

  expect_equal(quiet, character(0))
  expect_true(any(grepl("=", noisy, fixed = TRUE)))
  expect_equal(nrow(res$mod.data.out), do.call(generate_model_set, args)$n.mods)
})
