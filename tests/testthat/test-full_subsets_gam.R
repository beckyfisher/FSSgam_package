test_that("full_subsets_gam returns a fully populated output list", {
  use.dat <- FSSgam::case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  out <- full_subsets_gam(
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
    check_correlations(use.dat[, c("complexity", "depth", "ZONE")])
  )
})

test_that("full_subsets_gam matches a separate generate_model_set + fit_model_set call", {
  use.dat <- FSSgam::case_study1
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

  combined <- full_subsets_quietly(
    use.dat = args$use.dat, test.fit = args$test.fit,
    pred.vars.cont = args$pred.vars.cont, pred.vars.fact = args$pred.vars.fact,
    null.terms = args$null.terms, max.predictors = args$max.predictors, k = args$k
  )
  model.set <- generate_model_set(
    use.dat = args$use.dat, test.fit = args$test.fit,
    pred.vars.cont = args$pred.vars.cont, pred.vars.fact = args$pred.vars.fact,
    null.terms = args$null.terms, max.predictors = args$max.predictors, k = args$k
  )
  separate <- fit_quietly(model.set, parallel = FALSE)

  expect_equal(combined$mod.data.out$AICc, separate$mod.data.out$AICc)
})

test_that("full_subsets_gam's legacy factor.interactions argument forwards with a warning", {
  use.dat <- FSSgam::case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  expect_warning(
    legacy <- full_subsets_quietly(
      use.dat = use.dat, test.fit = test.fit,
      pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
      factor.interactions = TRUE,
      null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
    ),
    "factor.factor.interactions"
  )
  current <- full_subsets_quietly(
    use.dat = use.dat, test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    factor.factor.interactions = TRUE,
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  expect_equal(legacy$mod.data.out$AICc, current$mod.data.out$AICc)
})

test_that("full_subsets_gam forwards factor.factor.interactions", {
  # The legacy-argument test above cannot show this: it compares two calls both
  # made through full_subsets_gam(), so a change inside the wrapper moves both
  # sides together, and it supplies a single factor, which the factor-factor
  # block skips entirely. This needs two factors and a direct assertion on the
  # interaction column the argument is supposed to produce.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  out <- suppress_nnet_nans(full_subsets_quietly(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = "depth", pred.vars.fact = c("ZONE", "ZONE2"),
    factor.factor.interactions = TRUE,
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  ))

  expect_true("ZONE.I.ZONE2" %in% colnames(out$used.data))
  expect_true("ZONE.I.ZONE2" %in% out$mod.data.out$modname)
})

test_that("full_subsets_gam's legacy smooth.interactions argument forwards with a warning", {
  use.dat <- FSSgam::case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  expect_warning(
    legacy <- full_subsets_quietly(
      use.dat = use.dat, test.fit = test.fit,
      pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
      smooth.interactions = NA,
      null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
    ),
    "factor.smooth.interactions"
  )
  current <- full_subsets_quietly(
    use.dat = use.dat, test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    factor.smooth.interactions = NA,
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  expect_equal(legacy$mod.data.out$AICc, current$mod.data.out$AICc)
  # confirm the .by. interaction models are actually suppressed
  expect_false(any(grepl("\\.by\\.", legacy$mod.data.out$modname)))
})

test_that("full_subsets_gam's legacy size argument constrains max.predictors", {
  use.dat <- FSSgam::case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  # regression test: size= used to be a no-op (it warned but never actually
  # fed its value into max.predictors)
  expect_warning(
    legacy <- full_subsets_gam(
      use.dat = use.dat, test.fit = test.fit,
      pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
      size = 1,
      null.terms = "s(site,bs='re')", k = 3
    ),
    "max.predictors"
  )
  current <- full_subsets_gam(
    use.dat = use.dat, test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    max.predictors = 1,
    null.terms = "s(site,bs='re')", k = 3
  )

  expect_equal(sort(legacy$mod.data.out$modname), sort(current$mod.data.out$modname))
  # only single-predictor models (plus null) should be present
  expect_false(any(grepl("\\+", legacy$mod.data.out$modname)))
})

test_that("full_subsets_gam forwards its model-set arguments unchanged", {
  fit <- fixture_cs1_gaussian()
  args <- list(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    linear.vars = "SCORE2", null.terms = "s(site,bs='re')",
    max.predictors = 2, k = 3, bs.arg = "'ts'", cov.cutoff = 0.5
  )

  combined <- do.call(full_subsets_quietly, args)
  model.set <- do.call(generate_model_set, args)

  expect_equal(combined$mod.data.out$modname, names(model.set$mod.formula))
  expect_equal(combined$predictor.correlations, model.set$predictor.correlations)
  expect_equal(dim(combined$used.data), dim(model.set$used.data))
  # bs.arg has to be asserted on the fitted formulas: comparing the two calls'
  # model sets to each other only shows they agree, not that either used the
  # value supplied. The candidate names carry no basis information at all.
  smooths <- Filter(
    function(x) grepl("s(complexity", x, fixed = TRUE) || grepl("s(depth", x, fixed = TRUE),
    combined$mod.data.out$formula
  )
  expect_true(length(smooths) > 0)
  expect_true(all(grepl("bs = \"ts\"", smooths, fixed = TRUE)))
  expect_false(any(grepl("bs = \"cr\"", smooths, fixed = TRUE)))
})

test_that("full_subsets_gam forwards its fitting arguments unchanged", {
  fit <- fixture_cs1_gaussian()
  set.args <- list(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  combined <- do.call(
    full_subsets_quietly,
    c(set.args, list(r2.type = "dev", report.unique.r2 = TRUE, save.model.fits = FALSE))
  )
  separate <- fit_quietly(
    do.call(generate_model_set, set.args),
    r2.type = "dev", report.unique.r2 = TRUE, save.model.fits = FALSE
  )

  expect_equal(combined$mod.data.out$r2.vals, separate$mod.data.out$r2.vals)
  expect_equal(combined$mod.data.out$r2.vals.unique, separate$mod.data.out$r2.vals.unique)
  expect_s3_class(combined$success.models[[1]], "formula")
  expect_equal(combined$variable.importance, separate$variable.importance)
})

test_that("full_subsets_gam warns and downgrades when max.models is exceeded", {
  fit <- fixture_cs1_gaussian()

  expect_warning(
    out <- full_subsets_quietly(
      use.dat = fit$use.dat, test.fit = fit$test.fit,
      pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
      null.terms = "s(site,bs='re')", max.predictors = 2, k = 3,
      max.models = 2
    ),
    "Individual models fits will not be saved"
  )
  expect_s3_class(out$success.models[[1]], "formula")
})

test_that("full_subsets_gam forwards cyclic.vars and smooth.smooth.interactions", {
  fit <- fixture_cs3_cyclic()

  out <- full_subsets_quietly(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("lunar.date", "month"),
    cyclic.vars = c("lunar.date", "month"),
    smooth.smooth.interactions = TRUE,
    max.predictors = 2, k = 5
  )

  expect_true("lunar.date.te.month" %in% out$mod.data.out$modname)
  expect_length(out$failed.models, 0)
  expect_true(all(is.finite(out$mod.data.out$AICc)))
  # as with bs.arg above, cyclic.vars leaves no trace in the candidate names --
  # only in the basis each smooth was given, so that is what has to be asserted
  main.effects <- out$mod.data.out$formula[
    out$mod.data.out$modname %in% c("lunar.date", "month")
  ]
  expect_length(main.effects, 2)
  expect_true(all(grepl("bs = \"cc\"", main.effects, fixed = TRUE)))
  expect_equal(
    out$mod.data.out$formula[out$mod.data.out$modname == "lunar.date.te.month"],
    "te(lunar.date, month, k = 5, bs = c(\"cc\", \"cc\"))"
  )
})

test_that("full_subsets_gam forwards a user-supplied cor.matrix and non.linear.correlations", {
  fit <- fixture_cs1_gaussian()
  base.args <- list(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  nonlinear <- do.call(full_subsets_quietly, c(base.args, list(non.linear.correlations = TRUE)))
  expect_equal(
    nonlinear$predictor.correlations,
    check_non_linear_correlations(fit$use.dat[, c("complexity", "depth", "ZONE")])
  )

  zero.cor <- nonlinear$predictor.correlations
  zero.cor[] <- 0
  diag(zero.cor) <- 1
  supplied <- do.call(full_subsets_quietly, c(base.args, list(cor.matrix = zero.cor)))
  expect_equal(supplied$predictor.correlations, zero.cor)
  expect_true("complexity+depth" %in% supplied$mod.data.out$modname)
})

test_that("full_subsets_gam has no VI.mods argument and always uses the min.n default", {
  # Recorded as current behaviour: fit_model_set() takes VI.mods but
  # full_subsets_gam() does not forward it, so the all-models variant is
  # unreachable through the wrapper. Raised as FSSgam_package#7.
  expect_false("VI.mods" %in% names(formals(full_subsets_gam)))

  fit <- fixture_cs1_gaussian()
  set.args <- list(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  combined <- do.call(full_subsets_quietly, set.args)
  separate <- fit_quietly(do.call(generate_model_set, set.args), VI.mods = "min.n")
  expect_equal(combined$variable.importance, separate$variable.importance)
})

test_that("full_subsets_gam's legacy smooth.interactions argument accepts a character value", {
  fit <- fixture_cs1_gaussian()

  expect_warning(
    legacy <- full_subsets_quietly(
      use.dat = fit$use.dat, test.fit = fit$test.fit,
      pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
      smooth.interactions = "ZONE",
      null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
    ),
    "factor.smooth.interactions"
  )
  current <- full_subsets_quietly(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    factor.smooth.interactions = "ZONE",
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  expect_equal(legacy$mod.data.out$AICc, current$mod.data.out$AICc)
  expect_true(any(grepl(".by.ZONE", legacy$mod.data.out$modname, fixed = TRUE)))
})

# ---- the deprecated arguments on valid but awkward input --------------------
#
# These three were declared with a "previous.arg" sentinel default and tested
# with != or is.na(). Both comparisons fail on input the arguments are
# documented to accept, so passing a perfectly valid value raised an error
# instead of forwarding it. missing() replaced the sentinel.

test_that("a legacy smooth.interactions of length greater than one forwards", {
  # is.na(smooth.interactions) on a character vector of length 2 gave
  # "the condition has length > 1".
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  args <- list(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = c("ZONE", "ZONE2"),
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  expect_warning(
    legacy <- do.call(
      full_subsets_quietly, c(args, list(smooth.interactions = c("ZONE", "ZONE2")))
    ),
    "factor.smooth.interactions"
  )
  current <- do.call(
    full_subsets_quietly,
    c(args, list(factor.smooth.interactions = c("ZONE", "ZONE2")))
  )

  expect_equal(legacy$mod.data.out$AICc, current$mod.data.out$AICc)
  expect_true(any(grepl(".by.ZONE2", legacy$mod.data.out$modname, fixed = TRUE)))
})

test_that("a legacy factor.interactions of NA forwards", {
  # factor.interactions != "previous.arg" on NA gave "missing value where
  # TRUE/FALSE needed".
  fit <- fixture_cs1_gaussian()
  args <- list(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  expect_warning(
    legacy <- do.call(full_subsets_quietly, c(args, list(factor.interactions = NA))),
    "factor.factor.interactions"
  )
  current <- do.call(
    full_subsets_quietly, c(args, list(factor.factor.interactions = NA))
  )
  expect_equal(legacy$mod.data.out$AICc, current$mod.data.out$AICc)
})

test_that("a legacy smooth.interactions of NA still warns", {
  # NA used to take a branch of its own, which assigned the value but reached
  # the warning only through that branch. Under missing() it is assigned and
  # warned about like any other value.
  fit <- fixture_cs1_gaussian()

  expect_warning(
    legacy <- full_subsets_quietly(
      use.dat = fit$use.dat, test.fit = fit$test.fit,
      pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
      smooth.interactions = NA,
      null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
    ),
    "factor.smooth.interactions"
  )
  expect_false(any(grepl("\\.by\\.", legacy$mod.data.out$modname)))
})

test_that("the deprecated arguments stay absent when not supplied", {
  # missing() has to keep reporting correctly through full.subsets.gam(), which
  # forwards everything through ... rather than naming the arguments.
  fit <- fixture_cs1_gaussian()
  args <- list(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    null.terms = "s(site,bs='re')", max.predictors = 2, k = 3
  )

  expect_no_warning(do.call(full_subsets_quietly, args))

  # through the deprecated alias: the .Deprecated() warning fires, but none of
  # the three legacy-argument warnings do
  legacy.warnings <- character(0)
  withCallingHandlers(
    utils::capture.output(out <- do.call(full.subsets.gam, args)),
    warning = function(w) {
      legacy.warnings <<- c(legacy.warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("full_subsets_gam", legacy.warnings)))
  expect_false(any(grepl("has been replaced with", legacy.warnings)))

  # and supplying one through the alias is still detected
  expect_warning(
    withCallingHandlers(
      utils::capture.output(
        do.call(full.subsets.gam, c(args[names(args) != "max.predictors"],
                                    list(size = 1)))
      ),
      warning = function(w) {
        if (grepl("full_subsets_gam", conditionMessage(w))) invokeRestart("muffleWarning")
      }
    ),
    "max.predictors"
  )
})

test_that("max.models has the same default in the wrapper and in fit_model_set", {
  # The two defaults were 500 and 200, undocumented in either, so a candidate
  # set of 300 saved its fits through full_subsets_gam() and not through
  # generate_model_set() + fit_model_set().
  expect_equal(
    eval(formals(full_subsets_gam)$max.models),
    eval(formals(fit_model_set)$max.models)
  )
  expect_equal(eval(formals(full_subsets_gam)$max.models), 200)
})
