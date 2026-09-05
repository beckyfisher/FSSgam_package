test_that("generate_model_set builds a model set for a Gaussian response", {
  model.set <- fixture_cs1_model_set()

  expect_named(
    model.set,
    c("n.mods", "predictor.correlations", "null.term.correlations", "mod.formula",
      "used.data", "test.fit", "included.vars")
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

test_that("a single-level factor predictor is rejected", {
  # check_correlations() returns NA for a pair involving a single-level factor,
  # and factor_correlations() computed its own matrix without the zero-fill
  # build_predictor_correlation_matrix() applies, so a call supplying no
  # cor.matrix reached the screen with NA in it and stopped
  # (FSSgam_package#27). generate_model_set() now rejects such a factor by name
  # before either (FSSgam_package#33).
  #
  # The zero-fill added for FSSgam_package#27 is retained and is still
  # reachable, which an earlier version of this comment denied. A single-level
  # factor is only one way check_correlations() returns NA; two high-cardinality
  # factors exceed multinom()'s MaxNWts and do the same. It is asserted in its
  # own block below.
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

  # the NA that made the zero-fill necessary, at its source
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

test_that("the computed factor-factor matrix is zero-filled", {
  # factor_correlations() computes its own matrix inside
  # resolve_factor_interactions() and must zero-fill it, as
  # build_predictor_correlation_matrix() does, or the factor-factor screen is
  # reached with NA in it and the call stops (FSSgam_package#27).
  #
  # Two high-cardinality factors, not a single-level one: multinom() exceeds its
  # MaxNWts and check_correlations() returns NA, with nothing a validation
  # rejects. A single-level factor also produces the NA but is now refused
  # earlier (FSSgam_package#33), which is why that construction cannot be used
  # here and why an earlier version of this suite wrongly concluded the zero-fill
  # was unreachable.
  set.seed(1)
  n <- 200
  use.dat <- data.frame(
    y = rnorm(n),
    f1 = factor(sample(paste0("a", 1:40), n, TRUE)),
    f2 = factor(sample(paste0("b", 1:40), n, TRUE)),
    x = rnorm(n)
  )
  expect_true(is.na(suppressWarnings(check_correlations(use.dat[, c("f1", "f2")]))["f1", "f2"]))

  test.fit <- mgcv::gam(y ~ s(x, k = 3, bs = "cr"), data = use.dat)

  expect_no_error(suppressWarnings(generate_model_set(
    use.dat = use.dat, test.fit = test.fit, k = 3, pred.vars.cont = "x",
    pred.vars.fact = c("f1", "f2"), factor.factor.interactions = TRUE,
    max.predictors = 2
  )))
})

test_that("a supplied matrix augmented with derived columns is zero-filled", {
  # augment_supplied_correlation_matrix() splices computed rows and columns for
  # the hard coded interaction names a user cannot know to supply
  # (FSSgam_package#15), and zero-fills what the computation leaves NA. Without
  # that, those cells reach the candidate screen as NA.
  #
  # The supplied matrix must carry a name that is not a predictor, which is what
  # leaves an NA cell against the derived column. FSSgam_package#37 rejects a
  # non-predictor named in an interaction argument; it does not stop a supplied
  # matrix from carrying extra names, and this is the case the zero-fill exists
  # for. A version of this block without such a name could not fail: with every
  # name a predictor, the computed cells are never NA and deleting the zero-fill
  # returns a byte-identical matrix.
  set.seed(1)
  use.dat <- fixture_cs1_data()
  use.dat$fa <- factor(sample(c("p", "q"), nrow(use.dat), TRUE))
  use.dat$fb <- factor(sample(c("x", "y"), nrow(use.dat), TRUE))
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  nms <- c("depth", "fa", "fb", "junk")
  cor.mat <- matrix(0, 4, 4, dimnames = list(nms, nms))
  diag(cor.mat) <- 1

  model.set <- generate_model_set(
    use.dat = use.dat, test.fit = test.fit, k = 3, pred.vars.cont = "depth",
    pred.vars.fact = c("fa", "fb"), factor.factor.interactions = TRUE,
    max.predictors = 2, cor.matrix = cor.mat
  )

  expect_true("fa.I.fb" %in% rownames(model.set$predictor.correlations))
  expect_true("junk" %in% rownames(model.set$predictor.correlations))
  expect_false(anyNA(model.set$predictor.correlations))
  expect_identical(unname(model.set$predictor.correlations["junk", "fa.I.fb"]), 0)
})

test_that("a predictor named null is rejected", {
  # "null" is the name given to the null model's candidate, so a predictor of
  # that name produced two candidates called null: mod.formula[["null"]]
  # returned the first, and build_inclusion_mat() matched the wrong row, zeroing
  # the variable importance table (FSSgam_package#39).
  fit <- fixture_cs1_gaussian()
  fit$use.dat$null <- rnorm(nrow(fit$use.dat))

  expect_error(
    fixture_cs1_model_set(fit = fit, pred.vars.cont = c("depth", "null")),
    "may not be named"
  )
})

test_that("a predictor name containing + or * is rejected", {
  # "+" joins the terms of a candidate name and "*" writes a linear interaction,
  # so either is parsed as more than one term (FSSgam_package#39).
  fit <- fixture_cs1_gaussian()
  fit$use.dat[["a+b"]] <- rnorm(nrow(fit$use.dat))

  expect_error(
    fixture_cs1_model_set(fit = fit, pred.vars.cont = c("depth", "a+b")),
    "reserved separator"
  )
})

test_that("factor.smooth.interactions in its character form must name factors", {
  # The third interaction argument. Its list form was validated where it is
  # consumed; the character form was not, and names factors exactly as
  # factor.factor.interactions does (FSSgam_package#37).
  fit <- fixture_cs1_gaussian()

  expect_error(
    fixture_cs1_model_set(fit = fit, factor.smooth.interactions = c("ZONE", "site")),
    "site named in factor.smooth.interactions are not predictors"
  )
  expect_error(
    fixture_cs1_model_set(fit = fit, factor.smooth.interactions = c("ZONE", "ZONE")),
    "Each variable may be named once in factor.smooth.interactions"
  )
})

# ---- null.terms screening -----------------------------------------------------

test_that("a predictor correlated with a null.terms variable is dropped and named", {
  # null.terms forces a term into every candidate, and those terms are outside
  # cov.cutoff's matrix, which covers pred.vars.cont, pred.vars.fact and
  # linear.vars only. So a candidate could be arbitrarily strongly correlated
  # with a forced term and still appear in every model in the set, inflating the
  # variance of that term's estimate with nothing in the output to indicate it
  # (FSSgam_package#23).
  set.seed(1)
  use.dat <- fixture_cs1_data()
  use.dat$depth.copy <- use.dat$depth + rnorm(nrow(use.dat), 0, 0.01)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  args <- list(
    use.dat = use.dat, test.fit = test.fit, k = 3,
    pred.vars.cont = c("complexity", "depth.copy"), pred.vars.fact = "ZONE",
    null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')", max.predictors = 2
  )

  expect_warning(
    model.set <- do.call(generate_model_set, args),
    "depth.copy \\(depth, max 1\\)"
  )
  expect_false("depth.copy" %in% model.set$included.vars)
  expect_false(any(grepl("depth.copy", names(model.set$mod.formula), fixed = TRUE)))

  # and raising the cutoff keeps it
  kept <- do.call(generate_model_set, c(args, list(null.cov.cutoff = 1)))
  expect_true("depth.copy" %in% kept$included.vars)
})

test_that("null.term.correlations is returned whether or not anything is dropped", {
  # The point of returning them: the screen is silent when nothing exceeds the
  # cutoff, and a user comparing a forced term against candidates needs the
  # numbers either way (FSSgam_package#23).
  #
  # A fixed forced term, not the fixture's s(site,bs='re'): a random-effect
  # grouping factor is excluded from the screen entirely, so it contributes no
  # row.
  model.set <- fixture_cs1_model_set(
    pred.vars.cont = "complexity", pred.vars.fact = "ZONE",
    null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')"
  )

  expect_true(is.matrix(model.set$null.term.correlations))
  expect_identical(rownames(model.set$null.term.correlations), "depth")
  expect_true(all(c("complexity", "ZONE") %in%
                    colnames(model.set$null.term.correlations)))
  expect_false(anyNA(model.set$null.term.correlations))
})

test_that("a random-effect grouping factor is excluded from the screen", {
  # A grouping factor is correlated with the predictors measured within it by
  # construction. In the companion repository's case study 2, null.terms is
  # s(Location,Site,bs='re') and Status is nested in Location, so their
  # correlation is 1 -- the design of the study, not collinearity to screen out.
  # Screening it dropped Status, which is wrong (FSSgam_package#23).
  fit <- fixture_cs1_gaussian()
  fit$use.dat$fine <- factor(paste(fit$use.dat$ZONE, fit$use.dat$site))

  model.set <- suppress_nnet_nans(fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = "fine",
    null.terms = "s(site,bs='re')", max.predictors = 2
  ))

  # fine is built from site and so is perfectly correlated with it
  expect_true("fine" %in% model.set$included.vars)
  expect_null(model.set$null.term.correlations)
})

test_that("correlations among the null.terms variables are neither computed nor screened", {
  # Those terms are forced in by the user's decision, and dropping one is what
  # must not happen. Two perfectly correlated forced terms are both retained.
  set.seed(1)
  use.dat <- fixture_cs1_data()
  use.dat$depth.copy <- use.dat$depth + rnorm(nrow(use.dat), 0, 0.01)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )

  expect_no_warning(model.set <- generate_model_set(
    use.dat = use.dat, test.fit = test.fit, k = 3,
    pred.vars.cont = "complexity", pred.vars.fact = "ZONE",
    null.terms = "s(depth,k=3,bs='cr')+s(depth.copy,k=3,bs='cr')+s(site,bs='re')",
    max.predictors = 2
  ))

  # both forced terms appear in every candidate, correlated 1.0 with each other
  expect_true(all(grepl("depth.copy", vapply(model.set$mod.formula, deparse1, ""),
                        fixed = TRUE)))
  # and neither appears as a row against the other
  expect_false("depth.copy" %in% colnames(model.set$null.term.correlations))
})

test_that("a supplied cor.matrix is used for the forced terms it names", {
  # A supplied matrix is indexed by predictor and a forced term is not one, so
  # requiring it to be named there would disable this screen for every caller
  # who supplies a matrix -- reinstating the silence FSSgam_package#23 exists to
  # end. A forced term the matrix does name is read from it verbatim; the rest of
  # the block is computed.
  fit <- fixture_cs1_gaussian()
  base.cor <- fixture_cs1_model_set(fit = fit)$predictor.correlations
  args <- list(
    fit = fit, pred.vars.cont = "complexity", pred.vars.fact = "ZONE",
    null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')", max.predictors = 1
  )

  supplied <- base.cor
  supplied["depth", "complexity"] <- 0.95
  supplied["complexity", "depth"] <- 0.95

  expect_warning(
    screened <- do.call(fixture_cs1_model_set, c(args, list(cor.matrix = supplied))),
    "complexity \\(depth, max 0.95\\)"
  )
  expect_false("complexity" %in% screened$included.vars)
})

test_that("a supplied cor.matrix that names no forced term does not disable the screen", {
  # The ordinary case: the matrix covers the predictors and says nothing about
  # the forced term, so that term's row is computed from use.dat. An earlier
  # version returned NULL here, leaving every cor.matrix caller unscreened
  # (FSSgam_package#23).
  set.seed(1)
  use.dat <- fixture_cs1_data()
  use.dat$depth.copy <- use.dat$depth + rnorm(nrow(use.dat), 0, 0.01)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  args <- list(
    use.dat = use.dat, test.fit = test.fit, k = 3,
    pred.vars.cont = c("complexity", "depth.copy"), pred.vars.fact = "ZONE",
    null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')", max.predictors = 1
  )
  predictors.only <- suppressWarnings(
    do.call(generate_model_set, c(args, list(null.cov.cutoff = 1)))
  )$predictor.correlations
  expect_false("depth" %in% rownames(predictors.only))

  expect_warning(
    screened <- do.call(generate_model_set,
                        c(args, list(cor.matrix = predictors.only))),
    "depth.copy \\(depth"
  )
  expect_false("depth.copy" %in% screened$included.vars)
})

test_that("null.cov.cutoff is validated", {
  # Unvalidated, an NA gave an internal "missing value where TRUE/FALSE needed",
  # and a character or length-2 value was accepted silently, the latter
  # recycling into two concatenated messages.
  fit <- fixture_cs1_gaussian()
  for (bad in list(NA, "0.8", c(0.5, 0.8), -1)) {
    expect_error(
      fixture_cs1_model_set(fit = fit, null.cov.cutoff = bad),
      "null.cov.cutoff must be a single non-negative number"
    )
  }
})

test_that("a random effect written bs=c('re') is exempt like bs='re'", {
  # The exemption is read from the parsed call, not by matching the text. A
  # regex missed bs=c('re'), which mgcv accepts and which fits the identical
  # null model, so a nested factor was dropped and the call stopped with a
  # message naming max.predictors -- none of them the cause
  # (FSSgam_package#23).
  fit <- fixture_cs1_gaussian()
  fit$use.dat$fine <- factor(paste(fit$use.dat$ZONE, fit$use.dat$site))

  for (nt in c("s(site,bs='re')", "s(site,bs=c('re'))", 's(site, bs = "re")')) {
    model.set <- suppress_nnet_nans(fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "depth", pred.vars.fact = "fine",
      null.terms = nt, max.predictors = 2
    ))
    expect_true("fine" %in% model.set$included.vars)
    expect_null(model.set$null.term.correlations)
  }
})

test_that("a null.terms variable that is not a column of use.dat is skipped", {
  # null.terms accepts any formula fragment, so it may name a function or a term
  # written over several columns. A correlation is not defined for those.
  fit <- fixture_cs1_gaussian()

  expect_no_warning(model.set <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "complexity", pred.vars.fact = "ZONE",
    null.terms = "s(log(depth+1),k=3,bs='cr')+s(depth,k=3,bs='cr')+s(site,bs='re')"
  ))
  expect_identical(rownames(model.set$null.term.correlations), "depth")
})

test_that("an empty null.terms leaves null.term.correlations absent", {
  model.set <- fixture_cs1_model_set(null.terms = "")
  expect_null(model.set$null.term.correlations)
})

test_that("dropping predictors below max.predictors is reported as such", {
  # enumerate_candidate_models() stops with "max.predictors is greater than the
  # number of predictors", which a user has no reason to connect to the drop
  # they were just warned about (FSSgam_package#23).
  fit <- fixture_cs1_gaussian()
  base.cor <- fixture_cs1_model_set(fit = fit)$predictor.correlations
  supplied <- base.cor
  supplied["depth", "complexity"] <- 0.95
  supplied["complexity", "depth"] <- 0.95

  expect_error(
    suppressWarnings(fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "complexity", pred.vars.fact = "ZONE",
      null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')",
      cor.matrix = supplied, max.predictors = 2
    )),
    "After dropping complexity for correlation with a null.terms variable"
  )
})

test_that("a bare mgcv::gamm() test.fit is rejected with a usable message", {
  # It records no call, so stats::update() has nothing to re-evaluate and every
  # candidate refit fails. The error update() gives -- "need an object with call
  # component" -- names neither the argument nor the remedy, and the message
  # that followed advised the user to stop using uGamm, which is what they would
  # have had to use to succeed (FSSgam_package#34).
  use.dat <- fixture_cs1_data()
  test.fit <- mgcv::gamm(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
    random = list(site = ~1), data = use.dat
  )
  expect_true(is.null(test.fit$call))

  expect_error(
    generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("complexity", "depth"), max.predictors = 2
    ),
    "mgcv::gamm\\(\\) cannot be used"
  )
})

test_that("a uGamm test.fit of the same shape is accepted", {
  # The remedy the message gives has to work, or it is not a remedy.
  skip_if_not_installed("MuMIn")
  use.dat <- fixture_cs1_data()
  test.fit <- MuMIn::uGamm(
    Herbivore.abundance ~ s(depth, k = 4, bs = "cr"),
    random = list(site = ~1), data = use.dat
  )

  expect_no_error(generate_model_set(
    use.dat = use.dat, test.fit = test.fit, k = 4,
    pred.vars.cont = c("complexity", "depth"), max.predictors = 2
  ))
})

test_that("a variable that is both a predictor and a forced term is not screened against itself", {
  # Its correlation with itself is 1, so without excluding the forced terms from
  # the screened set it would be dropped for being correlated with itself
  # (FSSgam_package#23).
  fit <- fixture_cs1_gaussian()

  model.set <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = c("complexity", "depth"), pred.vars.fact = "ZONE",
    null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')"
  )

  expect_true("depth" %in% model.set$included.vars)
  expect_false("depth" %in% colnames(model.set$null.term.correlations))
})

test_that("a strong negative correlation with a forced term is screened", {
  # The screen takes the absolute value. check_correlations() returns
  # non-negative estimates for the pairs it computes, so only a supplied matrix
  # reaches this, and a hand-built matrix may hold a signed correlation
  # (FSSgam_package#23).
  fit <- fixture_cs1_gaussian()
  base.cor <- fixture_cs1_model_set(fit = fit)$predictor.correlations
  supplied <- base.cor
  supplied["depth", "complexity"] <- -0.95
  supplied["complexity", "depth"] <- -0.95

  expect_warning(
    screened <- fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "complexity", pred.vars.fact = "ZONE",
      null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')",
      cor.matrix = supplied, max.predictors = 1
    ),
    "complexity \\(depth, max 0.95\\)"
  )
  expect_false("complexity" %in% screened$included.vars)
})

test_that("null_term_variables keeps only names that are columns of use.dat", {
  # null.terms accepts any formula fragment, so it may name something that is
  # not a column. Indexing use.dat with such a name would error, and a
  # correlation is not defined for it (FSSgam_package#23).
  use.dat <- fixture_cs1_data()

  expect_identical(
    FSSgam:::null_term_variables("s(depth,k=3,bs='cr')", use.dat), "depth"
  )
  expect_identical(
    FSSgam:::null_term_variables("s(nosuchcolumn,k=3,bs='cr')", use.dat),
    character(0)
  )
  expect_identical(
    FSSgam:::null_term_variables("s(depth,k=3,bs='cr')+s(nosuchcolumn,k=3,bs='cr')",
                                 use.dat),
    "depth"
  )
})

test_that("a supplied cor.matrix overrides the computed estimate for a forced term", {
  # The point of reading a supplied matrix at all. A test whose supplied cell is
  # above the cutoff passes whether or not the supplied value is used, since the
  # computed one is above it too; this one supplies a value BELOW the cutoff
  # where the computed estimate is above, so it fails if the supplied value is
  # ignored (FSSgam_package#23).
  set.seed(1)
  use.dat <- fixture_cs1_data()
  use.dat$depth.copy <- use.dat$depth + rnorm(nrow(use.dat), 0, 0.01)
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  args <- list(
    use.dat = use.dat, test.fit = test.fit, k = 3,
    pred.vars.cont = c("complexity", "depth.copy"), pred.vars.fact = "ZONE",
    null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')", max.predictors = 1
  )

  # computed, depth against depth.copy is 1.00 and the predictor is dropped
  expect_warning(do.call(generate_model_set, args), "depth.copy \\(depth")

  # supplied as 0.1, it is not
  nms <- c("depth", "complexity", "depth.copy", "ZONE")
  supplied <- matrix(0, 4, 4, dimnames = list(nms, nms))
  diag(supplied) <- 1
  supplied["depth", "depth.copy"] <- 0.1
  supplied["depth.copy", "depth"] <- 0.1

  expect_no_warning(
    kept <- do.call(generate_model_set, c(args, list(cor.matrix = supplied)))
  )
  expect_true("depth.copy" %in% kept$included.vars)
  expect_identical(unname(kept$null.term.correlations["depth", "depth.copy"]), 0.1)
})

test_that("a forced term whose correlations cannot be computed is skipped, not fatal", {
  # The FSSgam_package#13 argument, asserted rather than only stated. A
  # predictor of a class check_correlations() cannot classify makes the
  # null-term computation fail; the screen is skipped with a warning and the
  # model set is still built, so that caller keeps working and is told their
  # forced terms were not screened.
  use.dat <- fixture_cs1_data()
  use.dat$when <- as.Date("2020-01-01") + seq_len(nrow(use.dat))
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  nms <- c("complexity", "when")
  supplied <- matrix(0, 2, 2, dimnames = list(nms, nms))
  diag(supplied) <- 1

  expect_warning(
    model.set <- generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("complexity", "when"),
      null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')",
      max.predictors = 1, cor.matrix = supplied
    ),
    "could not be computed, so those terms were not screened"
  )
  expect_true(model.set$n.mods > 0)
  expect_true(all(c("complexity", "when") %in% model.set$included.vars))
  # NULL rather than a row of zeros. A zero here would state a correlation for
  # the pairs the warning has just said were not screened (FSSgam_package#41).
  expect_null(model.set$null.term.correlations)
})

test_that("an interaction term is dropped when any of its parts is", {
  # keep.terms() splits on every separator, not only .by. -- a te() term built
  # on a dropped predictor must not re-enter the model set
  # (FSSgam_package#23).
  set.seed(1)
  n <- 200
  use.dat <- data.frame(
    y = rnorm(n), forced = rnorm(n), cB = rnorm(n), cC = rnorm(n)
  )
  use.dat$cA <- use.dat$forced + rnorm(n, 0, 0.01)
  test.fit <- mgcv::gam(y ~ s(forced, k = 3, bs = "cr"), data = use.dat)

  expect_warning(
    model.set <- generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("cA", "cB", "cC"),
      smooth.smooth.interactions = TRUE,
      null.terms = "s(forced,k=3,bs='cr')", max.predictors = 2
    ),
    "cA \\(forced"
  )
  expect_false(any(grepl("cA", names(model.set$mod.formula), fixed = TRUE)))
})

test_that("a failed computation keeps the forced terms the supplied matrix did name", {
  # Two forced terms, one named in the supplied matrix and one not, with a
  # predictor of a class check_correlations() cannot classify so the computation
  # for the unnamed one fails. The supplied row must survive: discarding it
  # would lose a screen the user asked for by supplying it
  # (FSSgam_package#23).
  use.dat <- fixture_cs1_data()
  use.dat$when <- as.Date("2020-01-01") + seq_len(nrow(use.dat))
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  nms <- c("complexity", "when", "depth")
  supplied <- matrix(0, 3, 3, dimnames = list(nms, nms))
  diag(supplied) <- 1
  supplied["depth", "complexity"] <- 0.2
  supplied["complexity", "depth"] <- 0.2

  expect_warning(
    model.set <- generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("complexity", "when"),
      # SCORE1 is a second fixed forced term, absent from the supplied matrix,
      # so its row must be computed -- which fails on the Date predictor
      null.terms = "s(depth,k=3,bs='cr')+s(SCORE1,k=3,bs='cr')+s(site,bs='re')",
      max.predictors = 1, cor.matrix = supplied
    ),
    "SCORE1 and the predictors could not be computed"
  )

  expect_false(is.null(model.set$null.term.correlations))
  expect_true("depth" %in% rownames(model.set$null.term.correlations))
  expect_identical(unname(model.set$null.term.correlations["depth", "complexity"]), 0.2)
  # and SCORE1, which was neither supplied nor computed, has no row at all
  # rather than a row of zeros (FSSgam_package#41)
  expect_false("SCORE1" %in% rownames(model.set$null.term.correlations))
})

test_that("a predictor correlated with any one of several forced terms is dropped", {
  # The screen takes any() across the forced terms, not all(). With all(), a
  # predictor at r = 1.00 with one forced term and uncorrelated with another is
  # retained in every candidate and nothing is reported -- FSSgam_package#23
  # reinstated. Every other drop test here uses a single forced term, so none
  # distinguishes the two.
  set.seed(1)
  n <- 200
  use.dat <- data.frame(
    y = rnorm(n), forcedA = rnorm(n), forcedB = rnorm(n), other = rnorm(n)
  )
  use.dat$copy <- use.dat$forcedA + rnorm(n, 0, 0.01)
  test.fit <- mgcv::gam(y ~ s(forcedA, k = 3, bs = "cr"), data = use.dat)

  expect_warning(
    model.set <- generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("copy", "other"),
      null.terms = "s(forcedA,k=3,bs='cr')+s(forcedB,k=3,bs='cr')",
      max.predictors = 1
    ),
    "copy \\(forcedA"
  )
  expect_false("copy" %in% model.set$included.vars)
  expect_true("other" %in% model.set$included.vars)
})

test_that("a factor predictor correlated with a forced term is dropped", {
  # The drop applies to pred.vars.fact as well as pred.vars.cont. Unfiltered, a
  # dropped factor re-entered the candidates while included.vars said it had
  # been dropped (FSSgam_package#23).
  set.seed(1)
  n <- 200
  use.dat <- data.frame(y = rnorm(n), other = rnorm(n))
  use.dat$forced <- factor(sample(c("a", "b"), n, TRUE))
  use.dat$echo <- use.dat$forced
  test.fit <- mgcv::gam(y ~ s(other, k = 3, bs = "cr"), data = use.dat)

  expect_warning(
    model.set <- suppress_nnet_nans(generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3, pred.vars.cont = "other",
      pred.vars.fact = "echo", null.terms = "forced", max.predictors = 1
    )),
    "echo \\(forced"
  )
  expect_false("echo" %in% model.set$included.vars)
  expect_false(any(grepl("echo", names(model.set$mod.formula), fixed = TRUE)))
})

test_that("a linear interaction term built on a dropped predictor is dropped", {
  # keep.terms() is applied to linear.interaction.terms as well as to
  # interaction.terms. Unfiltered, a .t. term re-entered built on a dropped
  # linear.vars entry (FSSgam_package#23).
  set.seed(1)
  n <- 200
  use.dat <- data.frame(y = rnorm(n), other = rnorm(n), forced = rnorm(n))
  use.dat$lin <- use.dat$forced + rnorm(n, 0, 0.01)
  use.dat$fac <- factor(sample(c("a", "b"), n, TRUE))
  test.fit <- mgcv::gam(y ~ s(other, k = 3, bs = "cr"), data = use.dat)

  expect_warning(
    model.set <- suppress_nnet_nans(generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3, pred.vars.cont = "other",
      pred.vars.fact = "fac", linear.vars = "lin",
      null.terms = "s(forced,k=3,bs='cr')", max.predictors = 2
    )),
    "lin \\(forced"
  )
  expect_false(any(grepl("lin", names(model.set$mod.formula), fixed = TRUE)))
})

test_that("an NA in one supplied direction does not stop the call", {
  # A hand-built matrix may leave a cell NA. Unfilled, it reaches the comparison
  # against null.cov.cutoff and stops the call with "missing value where
  # TRUE/FALSE needed", which is the failure FSSgam_package#27 removed from the
  # predictor screen (FSSgam_package#23). The cell reports the 0 the reverse
  # direction supplies, not the NA; with no value in either direction it would
  # be computed instead (FSSgam_package#41).
  fit <- fixture_cs1_gaussian()
  nms <- c("depth", "complexity", "ZONE")
  supplied <- matrix(0, 3, 3, dimnames = list(nms, nms))
  diag(supplied) <- 1
  supplied["depth", "complexity"] <- NA

  expect_no_error(
    model.set <- fixture_cs1_model_set(
      fit = fit, pred.vars.cont = "complexity", pred.vars.fact = "ZONE",
      null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')",
      cor.matrix = supplied, max.predictors = 1
    )
  )
  expect_identical(unname(model.set$null.term.correlations["depth", "complexity"]), 0)
})

test_that("a pair the supplied matrix gives no value for is computed", {
  # An NA in both directions was read as a correlation of zero: the pair was not
  # screened against null.cov.cutoff, the predictor appeared in every candidate
  # alongside the forced term, and nothing was reported -- decided by a cell the
  # user left empty. The same NA between two predictors stops the call
  # (FSSgam_package#27), and that inconsistency is FSSgam_package#41.
  #
  # Both shapes give no value for depth/a. The square matrix says NA in each
  # direction; the column-only matrix says NA in the one direction it has, and
  # is what exercises recording the supplied cells from the transposed block,
  # after the transposition rather than before it.
  use.dat <- fixture_cs1_data()
  use.dat$a <- use.dat$depth * 2
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  nms <- c("depth", "a", "complexity")
  square <- matrix(0, 3, 3, dimnames = list(nms, nms))
  diag(square) <- 1
  square["depth", "a"] <- NA
  square["a", "depth"] <- NA
  shapes <- list(
    both.directions = square,
    column.only = square[setdiff(nms, "depth"), , drop = FALSE]
  )

  # the supplied zero for depth/complexity is not the computed estimate, so the
  # assertions below distinguish a cell-wise merge from recomputing the row
  expect_gt(check_correlations(use.dat[, nms])["depth", "complexity"], 0.3)

  for (nm in names(shapes)) {
    expect_warning(
      model.set <- generate_model_set(
        use.dat = use.dat, test.fit = test.fit, k = 3,
        pred.vars.cont = c("a", "complexity"),
        null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')",
        max.predictors = 1, cor.matrix = shapes[[nm]]
      ),
      "a \\(depth", info = nm
    )
    expect_equal(
      unname(model.set$null.term.correlations["depth", "a"]), 1, info = nm
    )
    expect_identical(
      unname(model.set$null.term.correlations["depth", "complexity"]), 0, info = nm
    )
    expect_false("a" %in% model.set$included.vars, info = nm)
    # split on the separator rather than matching "a" anywhere in a candidate
    # name, which a one-letter fixed match does
    expect_false(
      "a" %in% unlist(strsplit(names(model.set$mod.formula), "+", fixed = TRUE)),
      info = nm
    )
  }
})

test_that("a failed computation leaves a pair with no supplied value at zero", {
  # The computation an unsupplied pair now triggers can fail, on the same
  # predictor class FSSgam_package#13 is about. The screen is skipped for that
  # pair rather than the call stopping, the pairs the matrix did give a value
  # for are kept, and the warning says which case it is (FSSgam_package#41).
  # The skipped pair reports zero rather than being absent, the row being kept
  # for the pair that was supplied; an NA there is what the screen stops on.
  use.dat <- fixture_cs1_data()
  use.dat$when <- as.Date("2020-01-01") + seq_len(nrow(use.dat))
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  nms <- c("depth", "complexity", "when")
  supplied <- matrix(0, 3, 3, dimnames = list(nms, nms))
  diag(supplied) <- 1
  supplied["depth", "complexity"] <- 0.2
  supplied["complexity", "depth"] <- 0.2
  supplied["depth", "when"] <- NA
  supplied["when", "depth"] <- NA

  expect_warning(
    model.set <- generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("complexity", "when"),
      null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')",
      max.predictors = 1, cor.matrix = supplied
    ),
    "the pairs the supplied cor.matrix gave no value for were not screened"
  )
  expect_identical(unname(model.set$null.term.correlations["depth", "complexity"]), 0.2)
  expect_identical(unname(model.set$null.term.correlations["depth", "when"]), 0)
  expect_true(all(c("complexity", "when") %in% model.set$included.vars))
})

test_that("an NA in the computed forced-term block is treated as zero", {
  # The computed block's own zero-fill, distinct from the supplied block's.
  # check_correlations() returns NA for a zero-variance predictor, and without
  # the fill that NA reaches the comparison against null.cov.cutoff and stops
  # the call with "missing value where TRUE/FALSE needed" -- the failure the
  # fill exists to prevent (FSSgam_package#23).
  set.seed(1)
  n <- 100
  use.dat <- data.frame(
    y = rnorm(n), forced = rnorm(n), flat = rep(1, n), other = rnorm(n)
  )
  expect_true(is.na(suppressWarnings(
    check_correlations(use.dat[, c("forced", "flat")])
  )["forced", "flat"]))

  test.fit <- mgcv::gam(y ~ s(other, k = 3, bs = "cr"), data = use.dat)

  expect_no_error(model.set <- suppressWarnings(generate_model_set(
    use.dat = use.dat, test.fit = test.fit, k = 3,
    pred.vars.cont = c("flat", "other"),
    null.terms = "s(forced,k=3,bs='cr')", max.predictors = 1
  )))
  expect_identical(unname(model.set$null.term.correlations["forced", "flat"]), 0)
})

test_that("null.term.correlations is NULL where every predictor is a forced term", {
  # One of the three documented cases for the element being NULL, and the only
  # one with no other test (FSSgam_package#23).
  fit <- fixture_cs1_gaussian()

  model.set <- fixture_cs1_model_set(
    fit = fit, pred.vars.cont = "depth", pred.vars.fact = NA,
    null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')", max.predictors = 1
  )

  expect_null(model.set$null.term.correlations)
})

test_that("the drop warning names the forced terms a predictor correlates with", {
  # The message is what a user acts on: which predictor was dropped, against
  # what, and how strongly. Naming the predictor alone would not say which
  # forced term to reconsider (FSSgam_package#23).
  set.seed(1)
  n <- 200
  use.dat <- data.frame(y = rnorm(n), forcedA = rnorm(n), other = rnorm(n))
  use.dat$copy <- use.dat$forcedA + rnorm(n, 0, 0.01)
  test.fit <- mgcv::gam(y ~ s(forcedA, k = 3, bs = "cr"), data = use.dat)

  expect_warning(
    generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("copy", "other"),
      null.terms = "s(forcedA,k=3,bs='cr')", max.predictors = 1
    ),
    "copy \\(forcedA, max 1\\)"
  )
})

test_that("the forced-term screen reads both directions of an asymmetric matrix", {
  # check_non_linear_correlations() returns a deliberately asymmetric matrix --
  # row is the response, column the predictor -- so the two directions answer
  # different questions, and only [predictor, forced] is large when a candidate
  # is a deterministic function of a forced term. Reading one direction admitted
  # such a candidate silently, which is the condition FSSgam_package#23 exists
  # to end. exceeds_cutoff() already reads both triangles for cov.cutoff.
  set.seed(1)
  n <- 300
  use.dat <- data.frame(y = rnorm(n), forcedA = rnorm(n), other = rnorm(n))
  use.dat$curve <- use.dat$forcedA^2 + rnorm(n, 0, 0.02)

  # the fixture must be asymmetric across the cutoff, or this pins nothing
  asym <- check_non_linear_correlations(use.dat[, c("forcedA", "curve")])
  expect_lt(asym["forcedA", "curve"], 0.8)
  expect_gt(asym["curve", "forcedA"], 0.8)

  test.fit <- mgcv::gam(y ~ s(other, k = 3, bs = "cr"), data = use.dat)

  expect_warning(
    model.set <- generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("curve", "other"),
      null.terms = "s(forcedA,k=3,bs='cr')", non.linear.correlations = TRUE,
      max.predictors = 1
    ),
    "curve \\(forcedA"
  )
  expect_false("curve" %in% model.set$included.vars)
  expect_gt(model.set$null.term.correlations["forcedA", "curve"], 0.8)
})

test_that("null.cov.cutoff = 1 admits an exact duplicate of a forced term", {
  # The documented way to turn the screen off. The comparison is strictly
  # greater than, so a correlation of exactly 1 does not exceed a cutoff of 1;
  # with >= it would, and the documented behaviour would fail on the one case a
  # user would test it with (FSSgam_package#23).
  set.seed(1)
  n <- 200
  use.dat <- data.frame(y = rnorm(n), forced = rnorm(n), other = rnorm(n))
  use.dat$exact <- use.dat$forced
  test.fit <- mgcv::gam(y ~ s(other, k = 3, bs = "cr"), data = use.dat)
  args <- list(
    use.dat = use.dat, test.fit = test.fit, k = 3,
    pred.vars.cont = c("exact", "other"),
    null.terms = "s(forced,k=3,bs='cr')", max.predictors = 1
  )

  expect_identical(
    unname(suppressWarnings(do.call(generate_model_set, args))$null.term.correlations["forced", "exact"]),
    1
  )
  expect_no_warning(kept <- do.call(generate_model_set, c(args, list(null.cov.cutoff = 1))))
  expect_true("exact" %in% kept$included.vars)
})

test_that("the drop warning names only the forced terms actually exceeded", {
  # The message subsets the forced terms by which of them exceeded the cutoff.
  # Without that subscript it names every forced term, sending the user to one
  # they need not reconsider (FSSgam_package#23).
  set.seed(1)
  n <- 200
  use.dat <- data.frame(
    y = rnorm(n), forcedA = rnorm(n), forcedB = rnorm(n), other = rnorm(n)
  )
  use.dat$copy <- use.dat$forcedA + rnorm(n, 0, 0.01)
  test.fit <- mgcv::gam(y ~ s(forcedA, k = 3, bs = "cr"), data = use.dat)

  w <- tryCatch(
    generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("copy", "other"),
      null.terms = "s(forcedA,k=3,bs='cr')+s(forcedB,k=3,bs='cr')",
      max.predictors = 1
    ),
    warning = function(x) conditionMessage(x)
  )
  expect_match(w, "copy (forcedA,", fixed = TRUE)
  expect_false(grepl("forcedB", w, fixed = TRUE))
})

test_that("the supplied forced-term block reads both directions", {
  # The supplied half of the both-directions fix. No test reached it: every
  # earlier block exercising asymmetry supplies no cor.matrix, so both
  # single-direction reductions of this path left the suite green
  # (FSSgam_package#23).
  use.dat <- fixture_cs1_data()
  use.dat$dep2 <- use.dat$depth * 2
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  nms <- c("depth", "dep2", "complexity")
  supplied <- matrix(0, 3, 3, dimnames = list(nms, nms))
  diag(supplied) <- 1
  # only the [predictor, forced] direction is large
  supplied["depth", "dep2"] <- 0.1
  supplied["dep2", "depth"] <- 0.99
  args <- list(
    use.dat = use.dat, test.fit = test.fit, k = 3,
    pred.vars.cont = c("dep2", "complexity"),
    null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')", max.predictors = 1
  )

  expect_warning(
    model.set <- do.call(generate_model_set, c(args, list(cor.matrix = supplied))),
    "dep2 \\(depth"
  )
  expect_identical(unname(model.set$null.term.correlations["depth", "dep2"]), 0.99)
  expect_false("dep2" %in% model.set$included.vars)
})

test_that("an NA in one supplied direction does not discard the other", {
  # pmax() defaults to na.rm = FALSE, and the supplied block was filled after
  # the two directions were combined rather than before. So an NA in one
  # direction discarded a real value in the other: the cell reported 0, the
  # predictor stayed in every candidate, and nothing was reported
  # (FSSgam_package#23).
  use.dat <- fixture_cs1_data()
  use.dat$dep2 <- use.dat$depth * 2
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  nms <- c("depth", "dep2", "complexity")
  supplied <- matrix(0, 3, 3, dimnames = list(nms, nms))
  diag(supplied) <- 1
  supplied["depth", "dep2"] <- NA
  supplied["dep2", "depth"] <- 0.99

  expect_warning(
    model.set <- generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("dep2", "complexity"),
      null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')",
      max.predictors = 1, cor.matrix = supplied
    ),
    "dep2 \\(depth"
  )
  expect_identical(unname(model.set$null.term.correlations["depth", "dep2"]), 0.99)
})

test_that("a supplied matrix naming a forced term in one dimension only is accepted", {
  # have comes from rownames() and cols from colnames(), so reading the reverse
  # direction as supplied[cols, have] assumed every such name appears in both.
  # A matrix naming a forced term in one dimension only -- a shape
  # build_predictor_correlation_matrix() documents itself as accepting --
  # aborted with a bare "subscript out of bounds" (FSSgam_package#23).
  use.dat <- fixture_cs1_data()
  use.dat$dep2 <- use.dat$depth * 2
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  nms <- c("depth", "dep2", "complexity")
  square <- matrix(0, 3, 3, dimnames = list(nms, nms))
  diag(square) <- 1
  square["depth", "dep2"] <- 0.99
  rows.only <- square[, setdiff(colnames(square), "depth"), drop = FALSE]

  expect_warning(
    model.set <- generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("dep2", "complexity"),
      null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')",
      max.predictors = 1, cor.matrix = rows.only
    ),
    "dep2 \\(depth"
  )
  expect_identical(unname(model.set$null.term.correlations["depth", "dep2"]), 0.99)
})

test_that("a supplied cor.matrix naming a forced term in either dimension is used", {
  # have was taken from rownames() alone, so a matrix naming a forced term only
  # as a column was discarded in full and a computed estimate substituted -- the
  # value the user supplied to control the screen ignored. The third
  # orientation of the same asymmetry: the value read, the row shape, and the
  # column shape (FSSgam_package#23).
  use.dat <- fixture_cs1_data()
  use.dat$a <- use.dat$depth * 2
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  args <- list(
    use.dat = use.dat, test.fit = test.fit, k = 3,
    pred.vars.cont = c("a", "complexity"),
    null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')", max.predictors = 1
  )
  nms <- c("depth", "a", "complexity")
  square <- matrix(0, 3, 3, dimnames = list(nms, nms))
  diag(square) <- 1
  # both directions, so that each shape below retains the value in whichever
  # orientation survives the subsetting
  square["a", "depth"] <- 0.99
  square["depth", "a"] <- 0.99

  # the same value, only the shape of the matrix differing
  shapes <- list(
    square = square,
    rows.only = square[, setdiff(colnames(square), "depth"), drop = FALSE],
    cols.only = square[setdiff(rownames(square), "depth"), , drop = FALSE]
  )
  for (nm in names(shapes)) {
    expect_warning(
      model.set <- do.call(generate_model_set, c(args, list(cor.matrix = shapes[[nm]]))),
      "a \\(depth", info = nm
    )
    expect_identical(
      unname(model.set$null.term.correlations["depth", "a"]), 0.99, info = nm
    )
  }
})

test_that("a bs given as a variable is exempt from the screen", {
  # mgcv accepts bs = mybs, and the exemption reads bs from the parsed call in
  # baseenv(), where such a name does not resolve. Screening a grouping factor
  # in error drops a legitimate predictor; not screening a fixed term in error
  # leaves the behaviour of earlier versions, so the unreadable case is exempt
  # (FSSgam_package#23).
  use.dat <- fixture_cs1_data()

  expect_length(FSSgam:::null_term_variables("s(site,bs=mybs)", use.dat), 0)
  # and a readable non-re bs is still screened
  expect_identical(
    FSSgam:::null_term_variables("s(depth,k=3,bs='cr')", use.dat), "depth"
  )
})

test_that("the larger of two differing supplied directions is used", {
  # fill() combines the two directions with pmax(). A plain assignment would
  # leave whichever direction is filled second, and every earlier test either
  # supplies one direction or supplies the same value both ways, so none
  # distinguishes them (FSSgam_package#23).
  use.dat <- fixture_cs1_data()
  use.dat$a <- use.dat$depth * 2
  test.fit <- mgcv::gam(
    log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
    data = use.dat
  )
  nms <- c("depth", "a", "complexity")
  supplied <- matrix(0, 3, 3, dimnames = list(nms, nms))
  diag(supplied) <- 1
  supplied["depth", "a"] <- 0.99   # forward
  supplied["a", "depth"] <- 0.10   # reverse, smaller

  expect_warning(
    model.set <- generate_model_set(
      use.dat = use.dat, test.fit = test.fit, k = 3,
      pred.vars.cont = c("a", "complexity"),
      null.terms = "s(depth,k=3,bs='cr')+s(site,bs='re')",
      max.predictors = 1, cor.matrix = supplied
    ),
    "a \\(depth"
  )
  expect_identical(unname(model.set$null.term.correlations["depth", "a"]), 0.99)
})

test_that("a forced term written without an explicit bs is screened", {
  # The exemption tests the term's bs argument, and a smooth with no bs at all
  # is a fixed term, not a random effect. No other test uses a null.terms smooth
  # written without one (FSSgam_package#23).
  use.dat <- fixture_cs1_data()

  expect_identical(
    FSSgam:::null_term_variables("s(depth,k=3)", use.dat), "depth"
  )
  expect_identical(
    FSSgam:::null_term_variables("s(depth,k=3)+s(site,bs='re')", use.dat), "depth"
  )
})
