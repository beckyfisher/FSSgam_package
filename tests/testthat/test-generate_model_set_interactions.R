test_that("generate_model_set builds factor-factor interaction terms when requested", {
  use.dat <- FSSgam::case_study1
  use.dat$site <- as.factor(use.dat$site)
  use.dat$ZONE2 <- factor(ifelse(use.dat$depth > median(use.dat$depth), "deep", "shallow"))
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  # suppress_nnet_nans: check_correlations() surfaces an "NaNs produced" warning
  # out of nnet whenever a pasted interaction column is perfectly predicted by
  # its own components, which is by construction (FSSgam_package#10). Whether
  # the optimiser reaches that state depends on the nnet version, so it warns on
  # some platforms and not others.
  model.set <- suppress_nnet_nans(generate_model_set(
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

  # The generated interaction column has to reach all.predictors before the
  # correlation matrix is built, or it is never screened for collinearity. It is
  # perfectly predicted by its own components by construction, so no candidate
  # may pair it with either of them.
  expect_true("ZONE.I.ZONE2" %in% rownames(model.set$predictor.correlations))
  pairs.own.components <- vapply(
    strsplit(names(model.set$mod.formula), "+", fixed = TRUE),
    function(terms.) {
      any(terms. == "ZONE.I.ZONE2") && any(terms. %in% c("ZONE", "ZONE2"))
    },
    logical(1)
  )
  expect_false(any(pairs.own.components))
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
    "exceeded cov.cutoff"
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

  # suppress_nnet_nans for the same reason as the three-factor test below: a
  # pasted interaction column is perfectly predicted by its own components, and
  # the multinom() summary check_correlations() takes its deviance from warns
  # "NaNs produced" while computing standard errors it discards
  # (FSSgam_package#10). Whether the optimiser reaches that state depends on the
  # nnet version, so this warns on some platforms and not others.
  model.set <- suppress_nnet_nans(fixture_cs1_model_set(
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

  # suppress_nnet_nans: a pasted interaction column is perfectly predicted by its
  # own components, and the multinom() summary that check_correlations() takes
  # the deviance from warns "NaNs produced" while computing standard errors it
  # then discards (FSSgam_package#10). The correlations themselves are correct.
  model.set <- suppress_nnet_nans(fixture_cs1_model_set(
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

  model.set <- suppress_nnet_nans(fixture_cs1_model_set(
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

  factors <- suppress_nnet_nans(fixture_cs1_model_set(
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
  named <- suppress_nnet_nans(fixture_cs1_model_set(
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

test_that("a pasted interaction counts as its component variables against max.predictors", {
  # The .I. counterpart of the te() case: a factor-factor interaction column is
  # one term but two variables, and max.predictors counts variables. Without the
  # .I. split in that count, a two-predictor model set admits candidates carrying
  # three.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )
  fit$use.dat$ZONE3 <- factor(
    ifelse(fit$use.dat$complexity > stats::median(fit$use.dat$complexity), "hi", "lo")
  )

  model.set <- suppress_nnet_nans(fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth",
    pred.vars.fact = c("ZONE", "ZONE2", "ZONE3"),
    factor.factor.interactions = TRUE, cov.cutoff = 0.95, max.predictors = 2
  ))
  mod.names <- names(model.set$mod.formula)

  # the interaction alone already uses both predictors the set allows
  expect_true("ZONE.I.ZONE2" %in% mod.names)
  expect_false(any(grepl("^depth\\+ZONE\\.I\\.|\\.I\\..*\\+depth$", mod.names)))
  # stated generally: no candidate may carry more than max.predictors variables
  n.vars <- vapply(mod.names, function(nm) {
    if (nm == "null") return(0L)
    terms. <- unlist(strsplit(nm, "+", fixed = TRUE))
    for (sep in c(".by.", ".I.", ".t.", ".te.")) {
      terms. <- unlist(strsplit(terms., sep, fixed = TRUE))
    }
    length(unique(terms.))
  }, integer(1))
  expect_true(all(n.vars <= 2))
})

test_that("collinearity is screened in both triangles of the correlation matrix", {
  # Factor-factor correlations are asymmetric: check_correlations() estimates
  # each ordered pair with its own multinom() fit. A nested factor predicts its
  # parent exactly while the parent predicts it only partly, so the exceedance
  # sits in one triangle and not the other, and which one depends purely on the
  # order the predictors were named. Testing only the upper triangle would make
  # the collinearity screen sensitive to argument order.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$fine <- factor(paste(fit$use.dat$ZONE, fit$use.dat$site))

  model.set <- suppress_nnet_nans(fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth",
    pred.vars.fact = c("fine", "ZONE"), cov.cutoff = 0.8, max.predictors = 2
  ))
  cm <- model.set$predictor.correlations

  # the fixture must be asymmetric across the cutoff, or this pins nothing
  expect_true(cm["ZONE", "fine"] > 0.8)
  expect_true(cm["fine", "ZONE"] < 0.8)
  # so the pair must still be excluded, from whichever triangle carries it
  expect_false("fine+ZONE" %in% names(model.set$mod.formula))
  expect_false("ZONE+fine" %in% names(model.set$mod.formula))
  # and each factor on its own is unaffected
  expect_true(all(c("fine", "ZONE") %in% names(model.set$mod.formula)))
})

test_that("factor-factor interactions are screened in both triangles too", {
  # resolve_factor_interactions() has its own collinearity screen, separate from
  # the enumeration-stage one tested above, and it is the screen that decides
  # which pasted interaction columns get built at all. It needs the same
  # both-triangles treatment for the same reason: with fine nested in ZONE the
  # exceedance sits in the lower triangle only.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$fine <- factor(paste(fit$use.dat$ZONE, fit$use.dat$site))

  expect_warning(
    model.set <- suppress_nnet_nans(fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "depth",
      pred.vars.fact = c("fine", "ZONE"), factor.factor.interactions = TRUE,
      cov.cutoff = 0.8, max.predictors = 2
    )),
    "exceeded cov.cutoff"
  )

  # the only available pair is collinear, so no interaction column is built
  expect_false(any(grepl(".I.", colnames(model.set$used.data), fixed = TRUE)))
  expect_false(any(grepl(".I.", names(model.set$mod.formula), fixed = TRUE)))
})

test_that("max.predictors counts variables, not terms", {
  # A te() term carries two variables. Counting terms rather than variables
  # would let a two-predictor model set admit a candidate with three.
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("depth", "SCORE2"), pred.vars.fact = "ZONE",
    smooth.smooth.interactions = TRUE, max.predictors = 2
  )

  expect_true("depth.te.SCORE2" %in% names(model.set$mod.formula))
  # the te() term already uses both of the two predictors allowed, so it cannot
  # be combined with a third
  expect_false("ZONE+depth.te.SCORE2" %in% names(model.set$mod.formula))
  expect_false(any(grepl(".te.", grep("+", names(model.set$mod.formula), fixed = TRUE,
                                       value = TRUE), fixed = TRUE)))
})

test_that("smooth-smooth interactions are pairwise, whatever max.predictors allows", {
  # cont.cmbns.max.predictors is fixed at 2 rather than following max.predictors:
  # te() terms are bivariate by decision, not by accident, and the source carries
  # a commented-out max.predictors on that line to say so.
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("depth", "SCORE2", "complexity"), pred.vars.fact = NA,
    smooth.smooth.interactions = TRUE, cov.cutoff = 0.9, max.predictors = 3
  )
  te.terms <- grep(".te.", names(model.set$mod.formula), fixed = TRUE, value = TRUE)

  expect_true(length(te.terms) > 0)
  # every te() term joins exactly two predictors, even though three would fit
  # within max.predictors
  n.joined <- vapply(strsplit(te.terms, "+", fixed = TRUE), function(terms.) {
    max(lengths(strsplit(terms., ".te.", fixed = TRUE)))
  }, integer(1))
  expect_true(all(n.joined == 2))
})

test_that("max.predictors = 1 builds no interaction columns and warns why", {
  # A3. The enumeration was written as for(i in 2:n), and at n = 1 that counts
  # backwards to c(2, 1), so a size-1 "combination" was enumerated too. Its
  # pasted name is just the original variable name, so cbind() appended a
  # duplicate column, and taking max() of the empty off-diagonal of a 1x1
  # sub-matrix warned "no non-missing arguments to max" once per variable.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  expect_warning(
    model.set <- fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "depth",
      pred.vars.fact = c("ZONE", "ZONE2"),
      factor.factor.interactions = TRUE, max.predictors = 1
    ),
    "max.predictors is 1"
  )

  expect_equal(anyDuplicated(colnames(model.set$used.data)), 0)
  expect_false(any(grepl(".I.", colnames(model.set$used.data), fixed = TRUE)))
  expect_false(any(grepl(".I.", names(model.set$mod.formula), fixed = TRUE)))
  # the correlation matrix covers the three real predictors, not a duplicated set
  expect_equal(dim(model.set$predictor.correlations), c(3, 3))
})

test_that("max.predictors = 1 emits no 'no non-missing arguments to max' warning", {
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  warnings.seen <- character(0)
  withCallingHandlers(
    model.set <- fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "depth",
      pred.vars.fact = c("ZONE", "ZONE2"),
      factor.factor.interactions = TRUE, max.predictors = 1
    ),
    warning = function(w) {
      warnings.seen <<- c(warnings.seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_false(any(grepl("no non-missing arguments to max", warnings.seen)))
})

test_that("a character factor.factor.interactions warns when the screen leaves nothing", {
  # The character branch used to reach cbind() with a zero-column data.frame and
  # fail with "arguments imply differing number of rows", naming neither the
  # argument nor the cutoff, where the logical branch warned. Both now warn.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  expect_warning(
    model.set <- fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "depth",
      pred.vars.fact = c("ZONE", "ZONE2"),
      factor.factor.interactions = c("ZONE", "ZONE2"), max.predictors = 1
    ),
    "max.predictors is 1"
  )

  expect_false(any(grepl(".I.", colnames(model.set$used.data), fixed = TRUE)))
})

test_that("both factor.factor.interactions forms screen both triangles alike", {
  # The logical branch screened upper.tri and lower.tri; the character branch
  # screened upper.tri alone, so it admitted a pair whose correlation exceeded
  # the cutoff in the other direction. check_correlations() fits multinom()
  # separately in each direction, so its factor block is not symmetric.
  #
  # ZONE3 and ZONE5 below measure about 0.5425 one way and 0.5434 the other, so
  # cov.cutoff = 0.543 falls between them. The third factor is needed because
  # the character branch is otherwise skipped by its own guard.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE3 <- factor(
    ifelse(fit$use.dat$complexity > stats::median(fit$use.dat$complexity), "hi", "lo")
  )
  fit$use.dat$ZONE5 <- factor(
    ifelse(fit$use.dat$rugosity > stats::median(fit$use.dat$rugosity), "r1", "r2")
  )
  facts <- c("ZONE", "ZONE3", "ZONE5")

  # the asymmetry the test depends on is real, and the cutoff sits inside it
  cm <- suppress_nnet_nans(check_correlations(fit$use.dat[, facts]))
  expect_lt(cm["ZONE3", "ZONE5"], 0.543)
  expect_gt(cm["ZONE5", "ZONE3"], 0.543)

  args <- list(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = facts,
    cov.cutoff = 0.543, max.predictors = 2
  )
  logical.set <- suppress_nnet_nans(do.call(
    fixture_cs1_model_set, c(args, list(factor.factor.interactions = TRUE))
  ))
  character.set <- suppress_nnet_nans(do.call(
    fixture_cs1_model_set, c(args, list(factor.factor.interactions = facts))
  ))

  # the straddling pair is rejected by both
  expect_false("ZONE3.I.ZONE5" %in% colnames(logical.set$used.data))
  expect_false("ZONE3.I.ZONE5" %in% colnames(character.set$used.data))
  # the uncorrelated pairs are still built
  expect_true("ZONE.I.ZONE3" %in% colnames(logical.set$used.data))
  expect_true("ZONE.I.ZONE3" %in% colnames(character.set$used.data))
  # and the two forms now agree entirely
  expect_equal(names(character.set$mod.formula), names(logical.set$mod.formula))
})

test_that("a character factor.factor.interactions reports a fully screened pair", {
  # An outer guard, length(which(factor.correlations < cov.cutoff)) > 1, used to
  # skip the whole block in silence. It counted matrix cells rather than pairs,
  # so two factors whose two off-diagonal estimates straddle cov.cutoff gave a
  # count of 1: no interaction was built and nothing was said. The logical
  # branch had no such guard.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE3 <- factor(
    ifelse(fit$use.dat$complexity > stats::median(fit$use.dat$complexity), "hi", "lo")
  )
  fit$use.dat$ZONE5 <- factor(
    ifelse(fit$use.dat$rugosity > stats::median(fit$use.dat$rugosity), "r1", "r2")
  )

  expect_warning(
    model.set <- suppress_nnet_nans(fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "depth",
      pred.vars.fact = c("ZONE3", "ZONE5"),
      factor.factor.interactions = c("ZONE3", "ZONE5"),
      cov.cutoff = 0.543, max.predictors = 2
    )),
    "exceeded cov.cutoff"
  )

  expect_false(any(grepl(".I.", colnames(model.set$used.data), fixed = TRUE)))
})

# ---- a supplied cor.matrix governs every stage ------------------------------

test_that("a supplied cor.matrix decides which te() terms are built", {
  # FSSgam_package#13. resolve_smooth_smooth_interactions() computed its own
  # correlations from use.dat and ignored cor.matrix entirely, so a user who
  # supplied a matrix saying two predictors were uncorrelated still found their
  # te() term absent.
  fit <- fixture_cs1_gaussian()
  cont <- c("depth", "complexity", "SCORE2")

  # depth and complexity correlate at about 0.33 in the data, above the default
  # cov.cutoff, so their te() term is excluded when the data decide
  from.data <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = cont, pred.vars.fact = "ZONE",
    smooth.smooth.interactions = TRUE, max.predictors = 2
  )
  expect_false("depth.te.complexity" %in% names(from.data$mod.formula))

  vars <- c(cont, "ZONE")
  cm <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  diag(cm) <- 1
  from.matrix <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = cont, pred.vars.fact = "ZONE",
    smooth.smooth.interactions = TRUE, max.predictors = 2, cor.matrix = cm
  )
  expect_true("depth.te.complexity" %in% names(from.matrix$mod.formula))
})

test_that("a supplied cor.matrix decides which factor interaction columns are built", {
  # The third stage that screens on correlation. ZONE.copy is a relabelling of
  # ZONE, so the data say they are perfectly correlated and no interaction is
  # built; a supplied matrix saying otherwise now governs, and must then also
  # carry the interaction column it causes to exist.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))
  facts <- c("ZONE", "ZONE.copy")

  expect_warning(
    from.data <- fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "depth", pred.vars.fact = facts,
      factor.factor.interactions = TRUE, max.predictors = 2
    ),
    "exceeded cov.cutoff"
  )
  expect_false("ZONE.I.ZONE.copy" %in% colnames(from.data$used.data))

  vars <- c("depth", facts, "ZONE.I.ZONE.copy")
  cm <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  diag(cm) <- 1
  from.matrix <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = facts,
    factor.factor.interactions = TRUE, max.predictors = 2, cor.matrix = cm
  )
  expect_true("ZONE.I.ZONE.copy" %in% colnames(from.matrix$used.data))
  expect_true("ZONE.I.ZONE.copy" %in% names(from.matrix$mod.formula))
})

test_that("a supplied cor.matrix missing a factor is reported by name", {
  fit <- fixture_cs1_gaussian()
  cm <- matrix(c(1, 0, 0, 1), 2, 2,
               dimnames = list(c("depth", "complexity"), c("depth", "complexity")))

  expect_error(
    fixture_cs1_model_set(
      fit = fit, pred.vars.cont = c("depth", "complexity"),
      pred.vars.fact = c("ZONE", "ZONE"), factor.factor.interactions = TRUE,
      max.predictors = 2, cor.matrix = cm
    ),
    "missing required predictors"
  )
})

test_that("supplying no cor.matrix leaves te() selection unchanged", {
  # The reorder subsets the resolved matrix instead of computing a fresh one
  # over the continuous predictors alone. The two agree, which is the thing
  # that makes the reorder safe for everyone not supplying a matrix.
  fit <- fixture_cs1_gaussian()
  cont <- c("depth", "complexity", "SCORE2")

  for (nlc in c(FALSE, TRUE)) {
    model.set <- suppress_nnet_nans(fixture_cs1_model_set(
      fit = fit, pred.vars.cont = cont, pred.vars.fact = "ZONE",
      smooth.smooth.interactions = TRUE, non.linear.correlations = nlc,
      max.predictors = 2
    ))
    te.terms <- grep(".te.", names(model.set$mod.formula), fixed = TRUE, value = TRUE)
    # computed directly from the continuous predictors, the way the function
    # used to do it
    cc <- if (nlc) {
      check_non_linear_correlations(fit$use.dat[, cont])
    } else {
      check_correlations(fit$use.dat[, cont])
    }
    expected <- combn(cont, 2, simplify = FALSE)
    expected <- expected[vapply(expected, function(x) {
      m <- cc[x, x]
      max(abs(m[upper.tri(m)])) <= 0.28 && max(abs(m[lower.tri(m)])) <= 0.28
    }, logical(1))]
    expected <- vapply(expected, paste, character(1), collapse = ".te.")
    expect_setequal(te.terms, expected)
  }
})

test_that("a supplied cor.matrix replaces the automatic estimate rather than adding to it", {
  # FSSgam_package#13. build_predictor_correlation_matrix() used to compute
  # check_correlations() either way and discard it when a matrix was supplied,
  # so supplying one did not avoid the multinom()/gam() fits and did not let a
  # user past check_correlations()'s restriction on column classes.
  fit <- fixture_cs1_gaussian()
  preds <- c("complexity", "depth", "ZONE")
  cm <- check_correlations(fit$use.dat[, preds])

  calls <- 0
  ns <- asNamespace("FSSgam")
  real <- get("check_correlations", envir = ns)
  # the binding is locked in an installed package and unlocked under pkgload,
  # so its original state is recorded rather than assumed
  was.locked <- bindingIsLocked("check_correlations", ns)
  if (was.locked) unlockBinding("check_correlations", ns)
  assign("check_correlations",
         function(...) { calls <<- calls + 1; real(...) }, envir = ns)
  on.exit({
    assign("check_correlations", real, envir = ns)
    if (was.locked) lockBinding("check_correlations", ns)
  }, add = TRUE)

  fixture_cs1_model_set(fit = fit, cor.matrix = cm)
  expect_identical(calls, 0)

  # and the same call with no matrix supplied does reach it, so the counter is
  # measuring something
  calls <- 0
  fixture_cs1_model_set(fit = fit)
  expect_gt(calls, 0)
})

test_that("a supplied cor.matrix covers a predictor check_correlations cannot classify", {
  # a Date column is not one of the four classes classify_correlation_predictors()
  # accepts, so the automatic estimate cannot be formed for it -- which is one
  # reason to supply a matrix in the first place
  fit <- fixture_cs1_gaussian()
  use.dat <- fit$use.dat
  use.dat$when <- as.Date("2020-01-01") + seq_len(nrow(use.dat))
  preds <- c("depth", "when")
  cm <- matrix(0, 2, 2, dimnames = list(preds, preds))
  diag(cm) <- 1

  model.set <- generate_model_set(
    use.dat = use.dat, test.fit = fit$test.fit,
    pred.vars.cont = "depth", linear.vars = "when",
    null.terms = "s(site,bs='re')", cor.matrix = cm, max.predictors = 2, k = 3
  )

  expect_true("depth+when" %in% names(model.set$mod.formula))
  expect_identical(model.set$predictor.correlations, cm)
})

# ---- a supplied cor.matrix and the interaction columns it causes -------------

test_that("a supplied cor.matrix need not name the interaction columns it causes", {
  # FSSgam_package#15. With factor.factor.interactions set, the hard coded
  # interaction columns are named by resolve_factor_interactions() from the
  # combinations that survive the screen, and that screen reads the supplied
  # matrix. A user therefore cannot know which of them to supply, and the call
  # used to stop partway through construction demanding them by name.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))
  facts <- c("ZONE", "ZONE.copy")

  vars <- c("depth", facts)
  cm <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  diag(cm) <- 1

  model.set <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = facts,
    factor.factor.interactions = TRUE, max.predictors = 2, cor.matrix = cm
  )
  pc <- model.set$predictor.correlations

  expect_false("ZONE.I.ZONE.copy" %in% colnames(fit$use.dat))
  expect_true("ZONE.I.ZONE.copy" %in% colnames(model.set$used.data))
  expect_true("ZONE.I.ZONE.copy" %in% rownames(pc))

  # every cell the user supplied is kept as supplied
  expect_equal(pc[vars, vars], cm)

  # the derived row and column are what check_correlations gives for that column
  computed <- check_correlations(
    model.set$used.data[, c("depth", facts, "ZONE.I.ZONE.copy")]
  )
  expect_equal(
    unname(pc["ZONE.I.ZONE.copy", vars]),
    unname(computed["ZONE.I.ZONE.copy", vars]),
    tolerance = 1e-6
  )
  expect_equal(
    unname(pc[vars, "ZONE.I.ZONE.copy"]),
    unname(computed[vars, "ZONE.I.ZONE.copy"]),
    tolerance = 1e-6
  )
})

test_that("a supplied value for an interaction column is used, not recomputed", {
  # Supplying the derived names is still honoured: the pasted column is
  # perfectly predicted by each of its components, so computing would give 1,
  # and a supplied 0 has to survive. This passes against the pre-fix code as
  # well; it pins behaviour the fix must not take away, rather than the fix.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))
  facts <- c("ZONE", "ZONE.copy")

  vars <- c("depth", facts, "ZONE.I.ZONE.copy")
  cm <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  diag(cm) <- 1

  model.set <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = facts,
    factor.factor.interactions = TRUE, max.predictors = 2, cor.matrix = cm
  )
  expect_equal(model.set$predictor.correlations, cm)
  # and with the components declared uncorrelated with the interaction, the
  # candidate set contains a model holding both, which the computed matrix
  # would have excluded
  expect_true("ZONE+ZONE.I.ZONE.copy" %in% names(model.set$mod.formula))
})

test_that("only the predictors the user named are reported as missing", {
  # A misspelled or forgotten predictor of the user's own is still an error.
  # The interaction column absent from the same matrix is not named in it,
  # because supplying it was never possible. This also passes against the
  # pre-fix code, where resolve_factor_interactions() rejects the named factor
  # before any interaction column exists; it is here so that the fallback
  # cannot start swallowing a user's own missing predictor.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))
  cm <- matrix(c(1, 0, 0, 1), 2, 2,
               dimnames = list(c("depth", "ZONE"), c("depth", "ZONE")))

  err <- expect_error(
    fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "depth",
      pred.vars.fact = c("ZONE", "ZONE.copy"),
      factor.factor.interactions = TRUE, max.predictors = 2, cor.matrix = cm
    )
  )
  expect_match(conditionMessage(err), "ZONE\\.copy")
  expect_false(grepl("ZONE.I.ZONE.copy", conditionMessage(err), fixed = TRUE))
})

test_that("the interaction rows are computed under non.linear.correlations too", {
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))
  facts <- c("ZONE", "ZONE.copy")

  vars <- c("depth", facts)
  cm <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  diag(cm) <- 1

  model.set <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = facts,
    factor.factor.interactions = TRUE, max.predictors = 2,
    non.linear.correlations = TRUE, cor.matrix = cm
  )
  pc <- model.set$predictor.correlations
  expect_true("ZONE.I.ZONE.copy" %in% rownames(pc))
  expect_equal(pc[vars, vars], cm)

  computed <- check_non_linear_correlations(
    model.set$used.data[, c("depth", facts, "ZONE.I.ZONE.copy")]
  )
  expect_equal(
    unname(pc["ZONE.I.ZONE.copy", vars]),
    unname(computed["ZONE.I.ZONE.copy", vars]),
    tolerance = 1e-6
  )
})

test_that("an unclassifiable predictor and a missing interaction column report both", {
  # A supplied cor.matrix lets a predictor of a class check_correlations rejects
  # be used (FSSgam_package#13), but computing rows for an interaction column
  # needs every predictor's data, so the two cannot be combined. The error says
  # so rather than reporting the class on its own.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))
  fit$use.dat$when <- as.Date("2020-01-01") + seq_len(nrow(fit$use.dat))
  facts <- c("ZONE", "ZONE.copy")

  vars <- c("depth", "when", facts)
  cm <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  diag(cm) <- 1

  expect_error(
    fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "depth", pred.vars.fact = facts,
      linear.vars = "when", factor.factor.interactions = TRUE,
      max.predictors = 2, cor.matrix = cm
    ),
    "could not be computed from use.dat"
  )
})

test_that("a cor.matrix supplied as a data.frame is augmented as well", {
  # Assigning a data.frame into a sub-block of a matrix replaces the matrix with
  # a plain list, so the augmentation failed with "incorrect number of
  # subscripts on matrix", naming neither the argument nor the cause. A
  # data.frame is accepted everywhere else, including by this branch when the
  # matrix is complete.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))
  facts <- c("ZONE", "ZONE.copy")

  vars <- c("depth", facts)
  cm <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  diag(cm) <- 1

  from.matrix <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = facts,
    factor.factor.interactions = TRUE, max.predictors = 2, cor.matrix = cm
  )
  from.frame <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = facts,
    factor.factor.interactions = TRUE, max.predictors = 2,
    cor.matrix = as.data.frame(cm)
  )
  expect_equal(from.frame$predictor.correlations,
               from.matrix$predictor.correlations)
  expect_identical(names(from.frame$mod.formula), names(from.matrix$mod.formula))
})

test_that("a name in the supplied matrix that is not a predictor gets no NA cell", {
  # The branch that computes the whole matrix replaces every NA with zero. The
  # augmented one has to as well: the character form of
  # smooth.smooth.interactions validates against the rownames of the resolved
  # matrix rather than against the predictors, so a name the user supplied that
  # is not a predictor reaches combine_uncorrelated(), where max() of an all-NA
  # sub-matrix returns NA and `if (NA)` raises an error naming nothing.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))
  facts <- c("ZONE", "ZONE.copy")

  vars <- c("depth", "complexity", facts, "junk")
  cm <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  diag(cm) <- 1

  model.set <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = c("depth", "complexity"),
    pred.vars.fact = facts, factor.factor.interactions = TRUE,
    smooth.smooth.interactions = c("junk", "ZONE.I.ZONE.copy"),
    max.predictors = 2, cor.matrix = cm
  )
  pc <- model.set$predictor.correlations
  expect_false(anyNA(pc["ZONE.I.ZONE.copy", ]))
  expect_false(anyNA(pc[, "ZONE.I.ZONE.copy"]))
  expect_identical(unname(pc["junk", "ZONE.I.ZONE.copy"]), 0)
})

test_that("a supplied matrix keeps its own row and column names", {
  # The supplied names are kept in their own order and the derived ones
  # appended to each dimension separately, so a matrix whose two dimnames
  # differ is not reshaped into a square one with an invented all-NA row.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))
  facts <- c("ZONE", "ZONE.copy")

  vars <- c("depth", facts)
  cm <- matrix(0, length(vars), length(vars) + 1,
               dimnames = list(vars, c(vars, "extra.col")))
  cm[cbind(vars, vars)] <- 1

  model.set <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = facts,
    factor.factor.interactions = TRUE, max.predictors = 2, cor.matrix = cm
  )
  pc <- model.set$predictor.correlations
  expect_identical(rownames(pc), c(vars, "ZONE.I.ZONE.copy"))
  expect_identical(colnames(pc), c(vars, "extra.col", "ZONE.I.ZONE.copy"))
  expect_false("extra.col" %in% rownames(pc))
})

test_that("a continuous predictor omitted from cor.matrix is still reported", {
  # The named/derived split. The blocks above omit a factor, which
  # resolve_factor_interactions() rejects several stages earlier, so they pass
  # whatever the split does. A continuous predictor is not screened there, so
  # this is the case that reaches build_predictor_correlation_matrix() with
  # both a named and a derived predictor missing.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))
  facts <- c("ZONE", "ZONE.copy")

  vars <- c("depth", facts)
  cm <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  diag(cm) <- 1

  err <- expect_error(
    fixture_cs1_model_set(
      fit = fit, pred.vars.cont = c("depth", "complexity"),
      pred.vars.fact = facts, factor.factor.interactions = TRUE,
      max.predictors = 2, cor.matrix = cm
    )
  )
  expect_match(conditionMessage(err), "complexity")
  expect_false(grepl("ZONE.I.ZONE.copy", conditionMessage(err), fixed = TRUE))
})

test_that("a derived name supplied in one dimension keeps the values given there", {
  # The missing names are tracked per dimension. Pooling them meant that a
  # derived column named as a column and not as a row was treated as missing
  # from both, and the computed values overwrote the column the user did supply.
  #
  # What the user sees changed is predictor.correlations. The candidate set is
  # not, in this configuration: enumerate_candidate_models() screens both
  # triangles, and the row was not supplied, so the computed 1 excludes
  # ZONE+ZONE.I.ZONE.copy whichever way the column is resolved.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))
  facts <- c("ZONE", "ZONE.copy")

  rows <- c("depth", facts)
  cols <- c(rows, "ZONE.I.ZONE.copy")
  cm <- matrix(0, length(rows), length(cols), dimnames = list(rows, cols))
  cm[cbind(rows, rows)] <- 1
  cm[, "ZONE.I.ZONE.copy"] <- c(0.11, 0.22, 0.33)

  model.set <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = facts,
    factor.factor.interactions = TRUE, max.predictors = 2, cor.matrix = cm
  )
  pc <- model.set$predictor.correlations

  expect_equal(unname(pc[rows, "ZONE.I.ZONE.copy"]), c(0.11, 0.22, 0.33))
  # the row was not supplied, so it is computed
  expect_equal(unname(pc["ZONE.I.ZONE.copy", "ZONE"]), 1, tolerance = 1e-6)
  expect_identical(rownames(pc), c(rows, "ZONE.I.ZONE.copy"))
  expect_identical(colnames(pc), cols)
})

test_that("a supplied cor.matrix with a duplicated name is reported", {
  # unique() over the dimnames would otherwise keep one of the two silently, and
  # which of them governs a screen depends on how the sub-matrix is indexed.
  fit <- fixture_cs1_gaussian()
  vars <- c("depth", "complexity", "complexity")
  cm <- matrix(0, 3, 3, dimnames = list(vars, vars))
  diag(cm) <- 1

  expect_error(
    fixture_cs1_model_set(
      fit = fit, pred.vars.cont = c("depth", "complexity"),
      pred.vars.fact = NA, cor.matrix = cm
    ),
    "duplicated names"
  )
})

test_that("where a supplied derived name sits in its dimension does not change the model set", {
  # The derived names are appended to the end of each dimension, so a user who
  # supplied one as a column anywhere other than last gets a resolved matrix
  # whose two dimnames hold the same set in different orders. The screens used
  # to index the sub-matrix by position in each dimension separately, so the
  # cells they compared were not the pairs they intended and the candidate set
  # moved. Both are now indexed by name in one order.
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE.copy <- factor(paste0(fit$use.dat$ZONE, ".copy"))
  facts <- c("ZONE", "ZONE.copy")
  rows <- c("depth", "complexity", facts)

  make.cm <- function(cols) {
    cm <- matrix(0, length(rows), length(cols), dimnames = list(rows, cols))
    cm[cbind(rows, rows)] <- 1
    cm[, "ZONE.I.ZONE.copy"] <- c(0.11, 0.12, 0.22, 0.33)
    cm
  }
  run <- function(cm) {
    fixture_cs1_model_set(
      fit = fit, pred.vars.cont = c("depth", "complexity"),
      pred.vars.fact = facts, factor.factor.interactions = TRUE,
      factor.smooth.interactions = NA, max.predictors = 3, cor.matrix = cm
    )
  }
  last <- run(make.cm(c(rows, "ZONE.I.ZONE.copy")))
  first <- run(make.cm(c("ZONE.I.ZONE.copy", rows)))

  expect_identical(names(last$mod.formula), names(first$mod.formula))
})

test_that("a complete supplied matrix screens the same however its columns are ordered", {
  # The same defect without any interaction column: a square matrix whose
  # colnames are a permutation of its rownames used to drop candidates whose
  # correlation was well below cov.cutoff.
  fit <- fixture_cs1_gaussian()
  vars <- c("depth", "complexity", "ZONE")
  cm <- matrix(0, 3, 3, dimnames = list(vars, vars))
  diag(cm) <- 1
  permuted <- cm[, c("ZONE", "depth", "complexity")]

  plain <- fixture_cs1_model_set(fit = fit, cor.matrix = cm)
  shuffled <- fixture_cs1_model_set(fit = fit, cor.matrix = permuted)

  expect_identical(names(shuffled$mod.formula), names(plain$mod.formula))
})

test_that("a duplicated name that is not a predictor is accepted", {
  # Only a duplicated name the model set indexes is ambiguous. One nothing reads
  # is left alone, as it is when no matrix is supplied at all.
  fit <- fixture_cs1_gaussian()
  vars <- c("depth", "complexity", "ZONE", "spare", "spare")
  cm <- matrix(0, 5, 5, dimnames = list(vars, vars))
  diag(cm) <- 1

  model.set <- fixture_cs1_model_set(fit = fit, cor.matrix = cm)
  expect_gt(model.set$n.mods, 1)
})

test_that("a factor named twice yields no interaction of itself", {
  # A side effect of indexing the sub-matrix by name. Where the same factor was
  # named twice, selecting rows by position returned a 1x1 sub-matrix, both
  # triangles of which are empty, so max() warned "no non-missing arguments to
  # max" and returned -Inf, the combination survived, and a ZONE.I.ZONE column
  # was built. Selecting by name returns the 2x2 sub-matrix whose off-diagonal
  # is the correlation of the factor with itself, which is 1, so the
  # combination is dropped as perfectly collinear.
  #
  # enumerate_candidate_models() deduplicates its terms before indexing, so it
  # still reaches a 1x1 sub-matrix and still admits a nonsense ZONE+ZONE
  # candidate with two "no non-missing arguments to max" warnings. That is
  # pre-existing and unchanged, and is FSSgam_package#28.
  fit <- fixture_cs1_gaussian()

  model.set <- suppressWarnings(fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "depth",
      pred.vars.fact = c("ZONE", "ZONE"),
      factor.factor.interactions = TRUE, max.predictors = 2
  ))
  # both assertions below are negative, so a degenerate return would satisfy
  # them; this one fails if the model set collapses for an unrelated reason
  expect_true("ZONE+depth" %in% names(model.set$mod.formula))
  expect_false("ZONE.I.ZONE" %in% colnames(model.set$used.data))
  expect_false(any(grepl("ZONE.I.ZONE", names(model.set$mod.formula), fixed = TRUE)))
})
