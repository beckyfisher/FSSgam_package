test_that("generate_model_set builds a model set for a Gaussian response", {
  model.set <- fixture_cs1_model_set()

  expect_named(
    model.set,
    c("n.mods", "predictor.correlations", "mod.formula", "used.data", "test.fit", "included.vars")
  )
  expect_equal(model.set$n.mods, length(model.set$mod.formula))
  expect_true(model.set$n.mods > 0)
  expect_s3_class(model.set$mod.formula[[1]], "formula")
  expect_true("null" %in% names(model.set$mod.formula))
  expect_equal(model.set$included.vars, c("complexity", "depth", "ZONE"))
})

test_that("generate_model_set returns exactly the elements it documents", {
  # Regression test: the @return block listed used.data, predictor.correlations
  # and a "generated.models" element that has never existed, and omitted three
  # that do. The same class of defect as the Phase 7 full_subsets_gam() one.
  # The element names are asserted in the block above; this one pins their types.
  model.set <- fixture_cs1_model_set()

  expect_type(model.set$n.mods, "integer")
  expect_true(is.matrix(model.set$predictor.correlations))
  expect_type(model.set$mod.formula, "list")
  expect_s3_class(model.set$used.data, "data.frame")
  expect_s3_class(model.set$test.fit, "gam")
  expect_type(model.set$included.vars, "character")
})

test_that("generate_model_set errors when use.dat is not a data.frame", {
  expect_error(
    generate_model_set(use.dat = list(a = 1), test.fit = NULL, pred.vars.cont = "a"),
    "data\\.frame"
  )
})

test_that("generate_model_set errors clearly when null.terms is not a single character string", {
  fit <- fixture_cs1_gaussian()

  bad.null.terms <- list(
    NA,
    NULL,
    123,
    TRUE,
    c("s(site,bs='re')", "s(other,bs='re')"),
    factor("s(site,bs='re')")
  )

  for (null.terms in bad.null.terms) {
    expect_error(
      fixture_cs1_model_set(fit = fit, null.terms = null.terms),
      "null.terms must be a single character string"
    )
  }
})

test_that("generate_model_set errors when predictors contain NA", {
  fit <- fixture_cs1_gaussian()
  fit$use.dat$depth[1] <- NA

  expect_error(
    fixture_cs1_model_set(fit = fit),
    "Predictor variables contain NA"
  )
})

test_that("generate_model_set reports a failed null model with the underlying error", {
  fit <- fixture_cs1_gaussian()

  expect_error(
    fixture_cs1_model_set(fit = fit, null.terms = "s(not_a_column,bs='re')"),
    "Null model not successfully fitted"
  )
})

test_that("generate_model_set errors when max.predictors exceeds the number of predictors", {
  fit <- fixture_cs1_gaussian()

  expect_error(
    fixture_cs1_model_set(
      fit = fit,
      pred.vars.cont = c("complexity", "depth"),
      pred.vars.fact = NA,
      max.predictors = 3
    ),
    "max.predictors is greater than the number of predictors"
  )
})

# ---- single-predictor model sets (regression) -------------------------------

test_that("generate_model_set builds a model set from a single predictor", {
  # Regression test: use.dat[, all.predictors] dropped to a bare vector when
  # there was exactly one predictor, so check_correlations() received something
  # that was not a data.frame and failed with "argument is of length zero"
  # rather than building a 1x1 correlation matrix.
  fit <- fixture_cs1_gaussian()

  cont.only <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = NA, max.predictors = 1
  )
  expect_equal(names(cont.only$mod.formula), c("null", "depth"))
  expect_equal(dim(cont.only$predictor.correlations), c(1, 1))
  expect_equal(unname(cont.only$predictor.correlations[1, 1]), 1)

  fact.only <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = NA, pred.vars.fact = "ZONE", max.predictors = 1
  )
  expect_equal(names(fact.only$mod.formula), c("null", "ZONE"))

  linear.only <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = NA, pred.vars.fact = NA,
    linear.vars = "depth", max.predictors = 1
  )
  expect_equal(names(linear.only$mod.formula), c("null", "depth"))
})

test_that("generate_model_set builds a single-predictor set with non.linear.correlations", {
  # The same regression, through check_non_linear_correlations(), which
  # additionally had no zero-pair case: a one-column input produced a zero-row
  # pair grid and failed in .rowNamesDF<-.
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = "depth", pred.vars.fact = NA, max.predictors = 1,
    non.linear.correlations = TRUE
  )

  expect_equal(names(model.set$mod.formula), c("null", "depth"))
  expect_equal(dim(model.set$predictor.correlations), c(1, 1))
})

test_that("generate_model_set builds no phantom by-terms when pred.vars.cont is NA", {
  # Regression test: pred.vars.cont = NA is the documented way to run without
  # smooth predictors, but expand.grid(NA, factors) paired the NA itself with
  # each factor, producing a phantom "NA.by.<factor>" candidate term.
  #
  # At max.predictors = 1 the term was discarded again before the model set was
  # returned, so the only visible symptom was enumerate_candidate_models()
  # taking max() of an empty correlation sub-matrix and warning twice -- which
  # is why expect_no_warning() below is the load-bearing assertion in that case.
  # From max.predictors = 2 the term survived into the returned set, giving a
  # candidate whose formula smoothed the literal NA. The second block covers
  # that.
  fit <- fixture_cs1_gaussian()

  # expect_no_warning() is the load-bearing assertion here: the phantom term was
  # always discarded before the model set was returned, so the expect_false()
  # below passed with the defect present too. Do not remove the warning check.
  expect_no_warning(
    model.set <- fixture_cs1_model_set(
      fit = fit, pred.vars.cont = NA, pred.vars.fact = "ZONE", max.predictors = 1
    )
  )
  expect_false(any(grepl("NA.by.", names(model.set$mod.formula), fixed = TRUE)))
})

test_that("no phantom by-term survives into the model set at higher max.predictors", {
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  expect_no_warning(
    model.set <- fixture_cs1_model_set(
      fit = fit, pred.vars.cont = NA, pred.vars.fact = c("ZONE", "ZONE2"),
      max.predictors = 2
    )
  )

  expect_false(any(grepl("NA.by.", names(model.set$mod.formula), fixed = TRUE)))
  formulas <- vapply(model.set$mod.formula, deparse1, character(1))
  expect_false(any(grepl("s(NA", formulas, fixed = TRUE)))
  expect_setequal(
    names(model.set$mod.formula),
    c("null", "ZONE", "ZONE2", "ZONE+ZONE2")
  )
})

# ---- cor.matrix -------------------------------------------------------------

test_that("generate_model_set uses a user-supplied cor.matrix verbatim", {
  fit <- fixture_cs1_gaussian()
  default.set <- fixture_cs1_model_set(fit = fit)

  # complexity and depth are correlated above the default cov.cutoff in
  # case_study1, so they never share a model. Supply a zero matrix instead and
  # the combined model must appear.
  zero.cor <- default.set$predictor.correlations
  zero.cor[] <- 0
  diag(zero.cor) <- 1

  supplied <- fixture_cs1_model_set(fit = fit, cor.matrix = zero.cor)

  expect_equal(supplied$predictor.correlations, zero.cor)
  expect_false("complexity+depth" %in% names(default.set$mod.formula))
  expect_true("complexity+depth" %in% names(supplied$mod.formula))
})

test_that("generate_model_set errors when a supplied cor.matrix omits a predictor", {
  fit <- fixture_cs1_gaussian()
  full.cor <- fixture_cs1_model_set(fit = fit)$predictor.correlations
  partial.cor <- full.cor[c("complexity", "depth"), c("complexity", "depth")]

  expect_error(
    fixture_cs1_model_set(fit = fit, cor.matrix = partial.cor),
    "Supplied cor.matrix is missing required predictors: ZONE"
  )
})

test_that("a 1x1 supplied cor.matrix is used rather than recomputed", {
  # The only case in which a supplied matrix is 1x1 is a single-predictor model
  # set, and a 1x1 matrix has length 1, which was the sentinel for "nothing
  # supplied" -- so the matrix was silently discarded and computed from use.dat
  # (FSSgam_package#26). The supplied value is deliberately not what
  # check_correlations() would give, which for a single predictor is 1.
  fit <- fixture_cs1_gaussian()
  supplied.cor <- matrix(0.5, 1, 1, dimnames = list("ZONE", "ZONE"))

  model.set <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = NA, pred.vars.fact = "ZONE",
    max.predictors = 1, cor.matrix = supplied.cor
  )

  expect_equal(model.set$predictor.correlations, supplied.cor)
})

test_that("cor.matrix rejects a value that is neither a matrix nor the NA default", {
  # Under the length-based sentinel what happened depended on the length, and
  # neither outcome was right. A length-one value of any class satisfied the
  # sentinel and was treated as nothing supplied, so the model set was built
  # from a computed matrix as though the argument had not been given. A value
  # of any other length fell through to the supplied branch and failed against
  # the missing-predictor check, reporting predictors rather than the argument
  # (FSSgam_package#26).
  fit <- fixture_cs1_gaussian()

  # length one, silently ignored before
  expect_error(
    fixture_cs1_model_set(fit = fit, cor.matrix = "oops"),
    "cor.matrix must be a two-dimensional matrix or data.frame"
  )
  expect_error(
    fixture_cs1_model_set(fit = fit, cor.matrix = list(matrix(1, 1, 1))),
    "cor.matrix must be a two-dimensional matrix or data.frame"
  )
  # other lengths, reported as a missing predictor before
  expect_error(
    fixture_cs1_model_set(fit = fit, cor.matrix = NULL),
    "cor.matrix must be a two-dimensional matrix or data.frame"
  )
  expect_error(
    fixture_cs1_model_set(fit = fit, cor.matrix = c(NA, NA)),
    "cor.matrix must be a two-dimensional matrix or data.frame"
  )
})

test_that("a supplied cor.matrix is validated before the factor interaction screen", {
  # resolve_factor_interactions() screens factor pairs against the raw supplied
  # matrix, and runs before build_predictor_correlation_matrix(). Checking only
  # in the latter left this screen reached with an NA still in the matrix, so
  # factor.factor.interactions still failed with "missing value where
  # TRUE/FALSE needed" (FSSgam_package#27).
  set.seed(1)
  use.dat <- fixture_cs1_data()
  use.dat$fa <- factor(sample(c("p", "q"), nrow(use.dat), TRUE))
  use.dat$fb <- factor(sample(c("x", "y"), nrow(use.dat), TRUE))
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  args <- list(
    use.dat = use.dat, test.fit = test.fit, k = 3, pred.vars.cont = "depth",
    pred.vars.fact = c("fa", "fb"), factor.factor.interactions = TRUE,
    max.predictors = 2
  )
  na.cor <- do.call(generate_model_set, args)$predictor.correlations
  na.cor["fa", "fb"] <- NA

  expect_error(
    do.call(generate_model_set, c(args, list(cor.matrix = na.cor))),
    "The correlation matrix has NA between terms"
  )
})

test_that("a name screened but not a predictor is now rejected outright", {
  # This case used to reach the screen. The character form of
  # smooth.smooth.interactions was validated against rownames(cor.matrix), so a
  # supplied matrix carrying an extra name admitted a non-predictor, and the NA
  # check had to be widened to cover it (FSSgam_package#27). That widening is
  # now belt and braces: the argument is validated against the predictor lists
  # before either screen runs (FSSgam_package#37), so the name never reaches
  # them.
  use.dat <- fixture_cs1_data()
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  cor.mat <- check_correlations(use.dat[, c("depth", "complexity", "macro")])

  expect_error(
    generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("depth", "complexity"),
      smooth.smooth.interactions = c("macro", "depth"), max.predictors = 2,
      cor.matrix = cor.mat
    ),
    "macro named in smooth.smooth.interactions are not predictors"
  )
})

test_that("an NA is reported on a pair crossing a derived column and a te() term", {
  # enumerate_candidate_models() indexes all.predictors, the hard coded .I.
  # columns, and the names given in the character form of
  # smooth.smooth.interactions. A pair crossing the last two is screened and was
  # covered by neither of two earlier placements of this check
  # (FSSgam_package#27).
  #
  # zz is a real predictor here. It was a non-predictor when this was written,
  # which is now rejected before either screen (FSSgam_package#37).
  set.seed(1)
  use.dat <- fixture_cs1_data()
  use.dat$fa <- factor(sample(c("p", "q"), nrow(use.dat), TRUE))
  use.dat$fb <- factor(sample(c("x", "y"), nrow(use.dat), TRUE))
  use.dat$zz <- rnorm(nrow(use.dat))
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  args <- list(
    use.dat = use.dat, test.fit = test.fit, k = 3,
    pred.vars.cont = c("depth", "zz"), pred.vars.fact = c("fa", "fb"),
    factor.factor.interactions = TRUE,
    smooth.smooth.interactions = c("zz", "depth"), max.predictors = 3
  )
  full.cor <- do.call(generate_model_set, args)$predictor.correlations
  expect_true("fa.I.fb" %in% rownames(full.cor))

  na.cor <- full.cor
  na.cor["fa.I.fb", "zz"] <- NA
  na.cor["zz", "fa.I.fb"] <- NA

  expect_error(
    do.call(generate_model_set, c(args, list(cor.matrix = na.cor))),
    "The correlation matrix has NA between terms"
  )
})

test_that("a single-level factor predictor is rejected, and the computed matrix is zero-filled", {
  # Two changes meet here. check_correlations() returns NA for a pair involving
  # a single-level factor, and factor_correlations() computed its own matrix
  # without the zero-fill build_predictor_correlation_matrix() applies, so a
  # call supplying no cor.matrix reached the screen with NA in it and stopped
  # (FSSgam_package#27). The zero-fill made it build instead -- a model set in
  # which fa.I.const is paste(fa, "a"), the same model as fa under two names.
  #
  # generate_model_set() now rejects such a factor by name before either
  # (FSSgam_package#33), so the user-visible behaviour is the rejection. The
  # zero-fill is asserted directly on the helper, since nothing public reaches
  # it with an NA any more.
  set.seed(1)
  use.dat <- fixture_cs1_data()
  use.dat$fa <- factor(sample(c("p", "q"), nrow(use.dat), TRUE))
  use.dat$const <- factor("a")
  test.fit <- mgcv::gam(log.Herbivore.biomass ~ s(site, bs = "re"), data = use.dat)

  expect_error(
    generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3, pred.vars.cont = NA,
      pred.vars.fact = c("fa", "const"), factor.factor.interactions = TRUE,
      max.predictors = 2
    ),
    "fewer than two observed levels: const"
  )

  # check_correlations() itself degrades rather than aborting, and the value
  # the factor-factor screen would see is zero rather than NA.
  expect_true(is.na(check_correlations(use.dat[, c("fa", "const")])["fa", "const"]))
})

test_that("an NA in a derived column supplied in one dimension only is reported", {
  # A hard coded interaction name supplied as a row and not as a column is
  # absent from one set of dimnames, so a check taken over the intersection of
  # the two missed it and its NA reached enumerate_candidate_models(). The
  # check runs after augment_supplied_correlation_matrix(), which zero-fills
  # every derived row and column it adds, so the only NA left is a user's
  # (FSSgam_package#27).
  set.seed(1)
  use.dat <- fixture_cs1_data()
  use.dat$fa <- factor(sample(c("p", "q"), nrow(use.dat), TRUE))
  use.dat$fb <- factor(sample(c("x", "y"), nrow(use.dat), TRUE))
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  args <- list(
    use.dat = use.dat, test.fit = test.fit, k = 3, pred.vars.cont = "depth",
    pred.vars.fact = c("fa", "fb"), factor.factor.interactions = TRUE,
    max.predictors = 2
  )
  full.cor <- do.call(generate_model_set, args)$predictor.correlations
  expect_true("fa.I.fb" %in% rownames(full.cor))

  one.dim <- full.cor[, setdiff(colnames(full.cor), "fa.I.fb"), drop = FALSE]
  one.dim["fa.I.fb", "depth"] <- NA

  expect_error(
    do.call(generate_model_set, c(args, list(cor.matrix = one.dim))),
    "The correlation matrix has NA between terms"
  )
})

test_that("an NA between two predictors of a supplied cor.matrix is reported by name", {
  # Reaching combine_uncorrelated() with the NA gives "missing value where
  # TRUE/FALSE needed", which names neither the matrix, the argument, nor the
  # pair (FSSgam_package#27).
  fit <- fixture_cs1_gaussian()
  na.cor <- fixture_cs1_model_set(fit = fit)$predictor.correlations
  na.cor["complexity", "depth"] <- NA

  expect_error(
    fixture_cs1_model_set(fit = fit, cor.matrix = na.cor),
    "The correlation matrix has NA between terms.*complexity/depth"
  )
})

test_that("an NA is reported only on pairs some screen actually compares", {
  # A deliberate departure from FSSgam_package#27, which asked for the report
  # not to depend on max.predictors. Reporting a cell that nothing reads is a
  # false failure, so the report covers exactly the pairs compared, and which
  # pairs those are depends on max.predictors and on the interaction arguments.
  #
  # Not "at max.predictors = 1 nothing is compared". A .by. term is one term of
  # a candidate but splits into two for the screen, so a continuous/factor pair
  # is still compared here. Only complexity/depth is not.
  fit <- fixture_cs1_gaussian()
  base.cor <- fixture_cs1_model_set(fit = fit, max.predictors = 1)$predictor.correlations

  reported <- function(pair) {
    na.cor <- base.cor
    na.cor[pair[1], pair[2]] <- NA
    res <- try(
      fixture_cs1_model_set(fit = fit, cor.matrix = na.cor, max.predictors = 1),
      silent = TRUE
    )
    inherits(res, "try-error")
  }

  expect_false(reported(c("complexity", "depth")))
  expect_true(reported(c("depth", "ZONE")))
  expect_true(reported(c("complexity", "ZONE")))
})

test_that("an NA on a pair no screen compares is accepted", {
  # Companion to the block above, reached differently. At max.predictors = 1 no
  # candidate holds two continuous predictors, so complexity and depth are never
  # compared and an NA between them is inert. Rejecting it would be a false
  # failure (FSSgam_package#27).
  #
  # This was written against a factor named in factor.factor.interactions but
  # not in pred.vars.fact, which is now rejected outright
  # (FSSgam_package#37), so the property is asserted this way instead.
  fit <- fixture_cs1_gaussian()
  na.cor <- fixture_cs1_model_set(fit = fit)$predictor.correlations
  na.cor["complexity", "depth"] <- NA

  expect_no_error(
    fixture_cs1_model_set(fit = fit, cor.matrix = na.cor, max.predictors = 1)
  )
})

test_that("an S4 Matrix is accepted as a cor.matrix", {
  # cor_matrix_supplied() tests dimensionality rather than a list of accepted
  # classes. A whitelist of matrix and data.frame rejected the S4 matrix
  # Matrix::Matrix() returns, which worked before this
  # argument was validated at all.
  skip_if_not_installed("Matrix")
  fit <- fixture_cs1_gaussian()
  base.cor <- fixture_cs1_model_set(fit = fit)$predictor.correlations
  s4.cor <- Matrix::Matrix(base.cor, dimnames = dimnames(base.cor))

  expect_true(isS4(s4.cor))
  expect_no_error(fixture_cs1_model_set(fit = fit, cor.matrix = s4.cor))
})

test_that("an NA on the diagonal of a supplied cor.matrix is accepted", {
  # Both screens take max(abs(.)) over upper.tri() and lower.tri(), which
  # exclude the diagonal, so a diagonal NA is never read and must not be
  # reported as though it were.
  fit <- fixture_cs1_gaussian()
  diag.na.cor <- fixture_cs1_model_set(fit = fit)$predictor.correlations
  diag(diag.na.cor) <- NA

  expect_no_error(fixture_cs1_model_set(fit = fit, cor.matrix = diag.na.cor))
})

# ---- factor.smooth.interactions ---------------------------------------------

test_that("factor.smooth.interactions = NA suppresses all by-terms", {
  model.set <- fixture_cs1_model_set(factor.smooth.interactions = NA)

  expect_false(any(grepl(".by.", names(model.set$mod.formula), fixed = TRUE)))
  # Candidate names are built by sorting the term names, so they follow the
  # active collation. testthat forces LC_COLLATE=C for reproducible output, which
  # is why the factor sorts ahead of the lowercase continuous predictors here.
  expect_setequal(
    names(model.set$mod.formula),
    c("null", "complexity", "depth", "ZONE", "ZONE+complexity", "ZONE+depth")
  )
})

test_that("factor.smooth.interactions as a character vector selects which factors interact", {
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  model.set <- fixture_cs1_model_set(
    fit = fit,
    pred.vars.cont = "depth",
    pred.vars.fact = c("ZONE", "ZONE2"),
    factor.smooth.interactions = "ZONE"
  )

  by.terms <- grep(".by.", names(model.set$mod.formula), fixed = TRUE, value = TRUE)
  expect_true(length(by.terms) > 0)
  expect_true(all(grepl(".by.ZONE$", by.terms)))
  expect_false(any(grepl(".by.ZONE2", by.terms, fixed = TRUE)))
})

test_that("factor.smooth.interactions as a list gives per-predictor-type control", {
  # The list form names which factors, continuous predictors and linear
  # predictors take part, independently of pred.vars.fact/cont and linear.vars.
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = "ZONE",
    linear.vars = "SCORE2",
    factor.smooth.interactions = list(
      fact.vars = "ZONE", cont.vars = "complexity", linear.vars = "SCORE2"
    )
  )
  mod.names <- names(model.set$mod.formula)

  # complexity was named as a cont.var, so its by-term is present; depth was not
  expect_true("ZONE+complexity.by.ZONE" %in% mod.names)
  expect_false(any(grepl("depth.by.", mod.names, fixed = TRUE)))
  # SCORE2 was named as a linear.var, so it gets a .t. (product) interaction
  expect_true("SCORE2.t.ZONE+ZONE" %in% mod.names)
  expect_match(
    deparse1(model.set$mod.formula[["SCORE2.t.ZONE+ZONE"]]), "SCORE2 * ZONE",
    fixed = TRUE
  )
})

test_that("factor.smooth.interactions list errors on a variable that is not a predictor", {
  expect_error(
    fixture_cs1_model_set(
      factor.smooth.interactions = list(fact.vars = "ZONE", cont.vars = "not_a_column")
    ),
    "not_a_column"
  )
})

# ---- factor.factor.interactions ---------------------------------------------

test_that("factor.factor.interactions errors when fewer than two factors are named", {
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  expect_error(
    fixture_cs1_model_set(
      fit = fit, pred.vars.fact = c("ZONE", "ZONE2"),
      factor.factor.interactions = "ZONE"
    ),
    "less than 2 factors"
  )
})

test_that("factor.factor.interactions errors when a named factor is not a predictor", {
  # The check was against colnames(use.dat), so it caught only a name absent
  # from the data and accepted any column, predictor or not. It is now against
  # pred.vars.fact, which is what the argument documents, so the message
  # changed with it (FSSgam_package#37).
  fit <- fixture_cs1_gaussian()
  fit$use.dat$ZONE2 <- factor(
    ifelse(fit$use.dat$SCORE2 > stats::median(fit$use.dat$SCORE2), "high", "low")
  )

  expect_error(
    fixture_cs1_model_set(
      fit = fit, pred.vars.fact = c("ZONE", "ZONE2"),
      factor.factor.interactions = c("ZONE", "not_a_column")
    ),
    "not_a_column named in factor.factor.interactions are not predictors"
  )
})

# ---- smooth.smooth.interactions ---------------------------------------------

test_that("smooth.smooth.interactions errors when fewer than two predictors are named", {
  expect_error(
    fixture_cs1_model_set(smooth.smooth.interactions = "depth"),
    "less than 2 variables as smooth.smooth.interactions"
  )
})

test_that("smooth.smooth.interactions errors when a named variable is not a predictor", {
  expect_error(
    fixture_cs1_model_set(
      smooth.smooth.interactions = c("depth", "not_a_column")
    ),
    "not_a_column"
  )

  # A column that is in use.dat but not in pred.vars.cont has no row in the
  # correlation matrix, so it cannot be screened. It used to be admitted --
  # combine_uncorrelated() took max() of an empty sub-matrix, warned "no
  # non-missing arguments to max" and returned -Inf, which is below any
  # cov.cutoff -- and the te() term was built regardless of its correlations.
  expect_error(
    fixture_cs1_model_set(
      smooth.smooth.interactions = c("depth", "SCORE2")
    ),
    "SCORE2.*are not predictors"
  )
})

test_that("smooth.smooth.interactions errors when TRUE with fewer than two continuous predictors", {
  expect_error(
    fixture_cs1_model_set(
      pred.vars.cont = "depth", pred.vars.fact = "ZONE",
      smooth.smooth.interactions = TRUE, max.predictors = 2
    ),
    "less than 2 continuous predictors"
  )
})

# ---- candidate names are locale-independent ---------------------------------

test_that("candidate model names do not depend on the collation locale", {
  # FSSgam_package#8. Names are built by sorting term names, and sort() follows
  # LC_COLLATE, so the same analysis produced "complexity+ZONE" in an
  # en_US.UTF-8 session and "ZONE+complexity" in a C-locale one, with
  # mod.data.out in a correspondingly different row order.
  #
  # testthat forces LC_COLLATE=C, so this has to set the locale itself to test
  # anything at all: under C collation the old and new code agree.
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)
  set.ok <- suppressWarnings(Sys.setlocale("LC_COLLATE", "en_US.UTF-8"))
  skip_if(
    !identical(Sys.getlocale("LC_COLLATE"), "en_US.UTF-8"),
    "en_US.UTF-8 collation not available"
  )
  # guard against the locale silently not taking effect
  expect_false(identical(sort(c("ZONE", "complexity")),
                         sort(c("ZONE", "complexity"), method = "radix")))

  model.set <- fixture_cs1_model_set()

  expect_equal(
    names(model.set$mod.formula),
    c("null", "complexity", "depth", "ZONE", "ZONE+complexity", "ZONE+depth",
      "ZONE+complexity.by.ZONE", "ZONE+depth.by.ZONE")
  )
})

# ---- predictor validation ----------------------------------------------------

test_that("check_correlations degrades on a single-level factor rather than aborting", {
  # The lm() in build_factor_continuous_skeleton() was the only fit in either
  # correlation function not wrapped in try(), so a single-level factor stopped
  # the whole call with "contrasts can be applied only to factors with 2 or more
  # levels" -- naming neither the function, the argument, nor the predictor
  # (FSSgam_package#33).
  dat <- data.frame(x = rnorm(20), f = factor(rep("a", 20)))

  expect_no_error(check_correlations(dat))
  expect_no_error(check_non_linear_correlations(dat))
  expect_true(is.na(check_correlations(dat)["x", "f"]))
  expect_true(is.na(check_non_linear_correlations(dat)["x", "f"]))
})

test_that("a factor whose extra levels are unobserved is rejected", {
  # droplevels() first: a factor declaring levels no row takes has one level in
  # the data and is the same problem (FSSgam_package#33).
  fit <- fixture_cs1_gaussian()
  fit$use.dat$decl <- factor(rep("a", nrow(fit$use.dat)), levels = c("a", "b"))

  expect_error(
    fixture_cs1_model_set(fit = fit, pred.vars.fact = c("ZONE", "decl")),
    "fewer than two observed levels: decl"
  )
})

test_that("a factor interaction column colliding with an existing column is declined", {
  # cbind() accepts the duplicate, and every later stage selects predictors by
  # name: use.dat[, "fa.I.fb"] returns the user's own column while the formulas
  # and the correlation matrix were built for the generated one. Declining
  # rather than overwriting, because overwriting silently changes a predictor
  # the user named (FSSgam_package#22).
  set.seed(1)
  use.dat <- fixture_cs1_data()
  use.dat$fa <- factor(sample(c("p", "q"), nrow(use.dat), TRUE))
  use.dat$fb <- factor(sample(c("x", "y"), nrow(use.dat), TRUE))
  use.dat$fa.I.fb <- rnorm(nrow(use.dat))
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  expect_warning(
    model.set <- generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3, pred.vars.cont = "depth",
      pred.vars.fact = c("fa", "fb"), factor.factor.interactions = TRUE,
      max.predictors = 2
    ),
    "already has a column of the same name: fa.I.fb"
  )

  expect_equal(sum(colnames(model.set$used.data) == "fa.I.fb"), 1L)
  expect_true(is.numeric(model.set$used.data$fa.I.fb))
  expect_false(any(grepl("fa.I.fb", names(model.set$mod.formula), fixed = TRUE)))
})

test_that("passing used.data back in declines every generated column", {
  # The shape arises without contrivance: running generate_model_set() twice and
  # passing the used.data of the first call as the use.dat of the second
  # presents every generated name as an existing column (FSSgam_package#22).
  set.seed(1)
  use.dat <- fixture_cs1_data()
  use.dat$fa <- factor(sample(c("p", "q"), nrow(use.dat), TRUE))
  use.dat$fb <- factor(sample(c("x", "y"), nrow(use.dat), TRUE))
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  args <- list(
    test.fit = test.fit, k = 3, pred.vars.cont = "depth",
    pred.vars.fact = c("fa", "fb"), factor.factor.interactions = TRUE,
    max.predictors = 2
  )

  first <- do.call(generate_model_set, c(list(use.dat = use.dat), args))
  expect_true("fa.I.fb" %in% colnames(first$used.data))

  expect_warning(
    second <- do.call(generate_model_set, c(list(use.dat = first$used.data), args)),
    "already has a column of the same name"
  )
  expect_equal(sum(colnames(second$used.data) == "fa.I.fb"), 1L)
})
