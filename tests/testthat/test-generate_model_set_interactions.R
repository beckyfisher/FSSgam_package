test_that("generate_model_set builds factor-factor interaction terms when requested", {
  use.dat <- FSSgam::case_study1
  use.dat$site <- as.factor(use.dat$site)
  use.dat$ZONE2 <- factor(ifelse(use.dat$depth > median(use.dat$depth), "deep", "shallow"))
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  # suppressWarnings: check_correlations() surfaces an "NaNs produced" warning
  # out of nnet whenever a pasted interaction column is perfectly predicted by
  # its own components, which is by construction (FSSgam_package#10). Whether
  # the optimiser reaches that state depends on the nnet version, so it warns on
  # some platforms and not others.
  model.set <- suppressWarnings(generate_model_set(
    use.dat = use.dat,
    test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = c("ZONE", "ZONE2"),
    factor.factor.interactions = TRUE,
    null.terms = "s(site,bs='re')",
    max.predictors = 2,
    k = 3
  ))

  expect_true("ZONE.I.ZONE2" %in% colnames(model.set$used.data))
  expect_true("ZONE.I.ZONE2" %in% names(model.set$mod.formula))
  expect_true(model.set$n.mods > 0)
})

test_that("generate_model_set builds a te() smooth-smooth interaction for uncorrelated predictors", {
  use.dat <- FSSgam::case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  # depth and SCORE2 have near-zero correlation in FSSgam::case_study1, so the te()
  # interaction survives the cov.cutoff exclusion (default 0.28)
  model.set <- generate_model_set(
    use.dat = use.dat,
    test.fit = test.fit,
    pred.vars.cont = c("depth", "SCORE2"),
    smooth.smooth.interactions = TRUE,
    null.terms = "s(site,bs='re')",
    max.predictors = 2,
    k = 3
  )

  expect_true("depth.te.SCORE2" %in% names(model.set$mod.formula))
  te.formula <- model.set$mod.formula[["depth.te.SCORE2"]]
  expect_s3_class(te.formula, "formula")
  expect_match(as.character(te.formula)[2], "^te\\(")
})

test_that("smooth-smooth interactions are excluded when predictors exceed cov.cutoff", {
  use.dat <- FSSgam::case_study1
  use.dat$site <- as.factor(use.dat$site)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  # complexity and depth are correlated above the default cov.cutoff (0.28),
  # so the te() interaction should be dropped from the candidate set
  model.set <- generate_model_set(
    use.dat = use.dat,
    test.fit = test.fit,
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = "ZONE",
    smooth.smooth.interactions = TRUE,
    null.terms = "s(site,bs='re')",
    max.predictors = 2,
    k = 3
  )

  expect_false("complexity.te.depth" %in% names(model.set$mod.formula))
})

test_that("generate_model_set uses check_non_linear_correlations when requested", {
  use.dat <- FSSgam::case_study1
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
    non.linear.correlations = TRUE,
    null.terms = "s(site,bs='re')",
    max.predictors = 2,
    k = 3
  )

  expected <- check_non_linear_correlations(use.dat[, c("complexity", "depth", "ZONE")])
  expect_equal(model.set$predictor.correlations, expected)
  # the non-linear correlation matrix is asymmetric, unlike the default
  # check_correlations() matrix
  expect_false(isTRUE(all.equal(
    model.set$predictor.correlations,
    t(model.set$predictor.correlations)
  )))
})

test_that("factor-factor interactions are excluded when the factors exceed cov.cutoff", {
  fit <- fixture_cs1_gaussian()
  # ZONE.copy is a relabelling of ZONE, so the two are perfectly correlated and
  # their interaction is dropped, leaving nothing to interact
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))

  expect_warning(
    model.set <- fixture_cs1_model_set(
      fit = fit,
      pred.vars.cont = "complexity",
      pred.vars.fact = c("ZONE", "ZONE.copy"),
      factor.factor.interactions = TRUE
    ),
    "no\\s+factors to interaction"
  )

  expect_false(any(grepl(".I.", names(model.set$mod.formula), fixed = TRUE)))
  # and the two collinear factors never appear in the same model either
  expect_false(any(grepl("ZONE\\+ZONE.copy", names(model.set$mod.formula))))
})

test_that("factor.factor.interactions as a character vector selects which factors interact", {
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )
  fit$use.dat$ZONE3 <- factor(
    ifelse(fit$use.dat$complexity > stats::median(fit$use.dat$complexity), "hi", "lo")
  )

  # suppressWarnings for the same reason as the three-factor test below: a
  # pasted interaction column is perfectly predicted by its own components, and
  # the multinom() summary check_correlations() takes its deviance from warns
  # "NaNs produced" while computing standard errors it discards
  # (FSSgam_package#10). Whether the optimiser reaches that state depends on the
  # nnet version, so this warns on some platforms and not others.
  model.set <- suppressWarnings(fixture_cs1_model_set(
    fit = fit,
    pred.vars.cont = "depth",
    pred.vars.fact = c("ZONE", "ZONE2", "ZONE3"),
    factor.factor.interactions = c("ZONE", "ZONE2"),
    max.predictors = 2
  ))

  expect_true("ZONE.I.ZONE2" %in% colnames(model.set$used.data))
  expect_true("ZONE.I.ZONE2" %in% names(model.set$mod.formula))
  # ZONE3 was not named, so it forms no interaction
  expect_false(any(grepl("ZONE3.I.", names(model.set$mod.formula), fixed = TRUE)))
  expect_false(any(grepl(".I.ZONE3", names(model.set$mod.formula), fixed = TRUE)))
})

test_that("factor.factor.interactions respects max.predictors when building combinations", {
  # with three factors and max.predictors = 3 the two- and three-way pasted
  # interactions are both generated
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )
  fit$use.dat$ZONE3 <- factor(
    ifelse(fit$use.dat$complexity > stats::median(fit$use.dat$complexity), "hi", "lo")
  )

  # suppressWarnings: a pasted interaction column is perfectly predicted by its
  # own components, and the multinom() summary that check_correlations() takes
  # the deviance from warns "NaNs produced" while computing standard errors it
  # then discards (FSSgam_package#10). The correlations themselves are correct.
  model.set <- suppressWarnings(fixture_cs1_model_set(
    fit = fit,
    pred.vars.cont = "depth",
    pred.vars.fact = c("ZONE", "ZONE2", "ZONE3"),
    factor.factor.interactions = TRUE,
    max.predictors = 3
  ))

  interaction.cols <- grep(".I.", colnames(model.set$used.data), fixed = TRUE, value = TRUE)
  expect_true("ZONE.I.ZONE2" %in% interaction.cols)
  expect_true("ZONE.I.ZONE2.I.ZONE3" %in% interaction.cols)
})

test_that("smooth.smooth.interactions as a character vector selects which predictors interact", {
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("depth", "SCORE2", "complexity"),
    pred.vars.fact = NA,
    smooth.smooth.interactions = c("depth", "SCORE2"),
    max.predictors = 2
  )

  expect_true("depth.te.SCORE2" %in% names(model.set$mod.formula))
  # complexity was not named, so it forms no te() term
  expect_false(any(grepl("complexity.te.", names(model.set$mod.formula), fixed = TRUE)))
  expect_false(any(grepl(".te.complexity", names(model.set$mod.formula), fixed = TRUE)))
})

test_that("smooth.smooth.interactions uses non-linear correlations when asked", {
  # the te() candidates are screened against check_non_linear_correlations()
  # rather than check_correlations() when non.linear.correlations = TRUE
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("depth", "SCORE2"),
    pred.vars.fact = NA,
    smooth.smooth.interactions = TRUE,
    non.linear.correlations = TRUE,
    max.predictors = 2
  )

  expect_equal(
    model.set$predictor.correlations,
    check_non_linear_correlations(fixture_cs1_data()[, c("depth", "SCORE2")])
  )
  expect_true("depth.te.SCORE2" %in% names(model.set$mod.formula))
})

test_that("named smooth.smooth.interactions are screened with non-linear correlations too", {
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("depth", "SCORE2", "complexity"),
    pred.vars.fact = NA,
    smooth.smooth.interactions = c("depth", "SCORE2"),
    non.linear.correlations = TRUE,
    max.predictors = 2
  )

  expect_true("depth.te.SCORE2" %in% names(model.set$mod.formula))
})

test_that("named factor.factor.interactions are still screened against cov.cutoff", {
  # The character-vector form of factor.factor.interactions has its own
  # correlation-exclusion pass, separate from the logical TRUE form tested above.
  # ZONE.copy is a relabelling of ZONE, so that pair must be dropped while the
  # independent ZONE/ZONE2 pair survives.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))

  model.set <- suppressWarnings(fixture_cs1_model_set(
    fit = fit,
    pred.vars.cont = "depth",
    pred.vars.fact = c("ZONE", "ZONE2", "ZONE.copy"),
    factor.factor.interactions = c("ZONE", "ZONE2", "ZONE.copy"),
    max.predictors = 2
  ))
  interaction.cols <- grep(".I.", colnames(model.set$used.data), fixed = TRUE, value = TRUE)

  expect_true("ZONE.I.ZONE2" %in% interaction.cols)
  expect_false("ZONE.I.ZONE.copy" %in% interaction.cols)
})

test_that("named smooth.smooth.interactions are still screened against cov.cutoff", {
  # complexity and depth correlate above the default cov.cutoff in FSSgam::case_study1,
  # so naming them explicitly must not bypass the exclusion the logical TRUE
  # form applies.
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = NA,
    smooth.smooth.interactions = c("complexity", "depth"),
    max.predictors = 2
  )

  expect_false(any(grepl(".te.", names(model.set$mod.formula), fixed = TRUE)))
})

test_that("cov.cutoff decides which predictors may share a model", {
  # complexity and depth correlate at about 0.61 in FSSgam::case_study1: excluded at the
  # default 0.28, admitted once cov.cutoff is raised past their correlation.
  fit <- fixture_cs1_gaussian()
  observed.cor <- abs(stats::cor(fit$use.dat$complexity, fit$use.dat$depth))
  expect_true(observed.cor > 0.28 && observed.cor < 0.9)

  strict <- fixture_cs1_model_set(fit = fit, pred.vars.fact = NA, cov.cutoff = 0.28)
  relaxed <- fixture_cs1_model_set(fit = fit, pred.vars.fact = NA, cov.cutoff = 0.9)

  expect_false("complexity+depth" %in% names(strict$mod.formula))
  expect_true("complexity+depth" %in% names(relaxed$mod.formula))
  # and the single-predictor models are unaffected either way
  expect_true(all(c("complexity", "depth") %in% names(strict$mod.formula)))
})

test_that("interaction enumeration is capped by the number of predictors, not max.predictors", {
  # Both interaction builders clamp their combn() depth to the number of
  # predictors available to them when max.predictors is larger, which is what
  # keeps combn() from being asked for combinations that cannot exist.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  factors <- suppressWarnings(fixture_cs1_model_set(
    fit = fit,
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = c("ZONE", "ZONE2"),
    factor.factor.interactions = TRUE,
    max.predictors = 3
  ))
  # only the two-way interaction exists; a three-way one would need a third factor
  interaction.cols <- grep(".I.", colnames(factors$used.data), fixed = TRUE, value = TRUE)
  expect_equal(interaction.cols, "ZONE.I.ZONE2")

  # the character-vector form clamps against the number of factors it was given,
  # not against the whole of pred.vars.fact
  named <- suppressWarnings(fixture_cs1_model_set(
    fit = fit,
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = c("ZONE", "ZONE2"),
    factor.factor.interactions = c("ZONE", "ZONE2"),
    max.predictors = 3
  ))
  expect_equal(
    grep(".I.", colnames(named$used.data), fixed = TRUE, value = TRUE), "ZONE.I.ZONE2"
  )

  smooths <- fixture_cs1_model_set(
    fit = fit,
    pred.vars.cont = c("depth", "SCORE2", "complexity"),
    pred.vars.fact = NA,
    smooth.smooth.interactions = c("depth", "SCORE2"),
    max.predictors = 3
  )
  te.terms <- grep(".te.", names(smooths$mod.formula), fixed = TRUE, value = TRUE)
  expect_true(all(grepl("^depth\\.te\\.SCORE2", te.terms)))
})
