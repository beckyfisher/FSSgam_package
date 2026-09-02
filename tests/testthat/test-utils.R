# Unexported helpers. testthat runs tests in an environment whose parent is the
# package namespace, so these are reachable by name.

# ---- classify_correlation_predictors ----------------------------------------

test_that("classify_correlation_predictors splits factor-like from continuous columns", {
  dat <- data.frame(
    f = factor(c("a", "b", "a")),
    ch = c("x", "y", "x"),
    i = 1:3,
    n = c(1.5, 2.5, 3.5),
    stringsAsFactors = FALSE
  )

  vars <- classify_correlation_predictors(dat)

  expect_equal(vars$fact.vars, c("f", "ch"))
  expect_equal(vars$cont.vars, c("i", "n"))
})

test_that("classify_correlation_predictors names the offending column and its class", {
  dat <- data.frame(n = 1:3)
  dat$when <- as.Date("2020-01-01") + 0:2

  expect_error(classify_correlation_predictors(dat), "when")
  expect_error(classify_correlation_predictors(dat), "Date")
  expect_error(classify_correlation_predictors(dat), "not supported")
})

test_that("classify_correlation_predictors accepts a single column of each kind", {
  expect_equal(
    classify_correlation_predictors(data.frame(n = 1:3))$cont.vars, "n"
  )
  expect_equal(
    classify_correlation_predictors(data.frame(f = factor(c("a", "b"))))$fact.vars, "f"
  )
})

# ---- resolve_candidate_family -----------------------------------------------

test_that("resolve_candidate_family returns NULL when test.fit has no explicit family", {
  # a gaussian gam carries no family= in its call, and update() then needs no
  # override at all
  fit <- fixture_cs1_gaussian()

  expect_null(fit$test.fit$call$family)
  expect_null(resolve_candidate_family(fit$test.fit))
})

test_that("resolve_candidate_family re-evaluates the unfitted family from test.fit's call", {
  fit <- fixture_cs1_tweedie()

  fam <- resolve_candidate_family(fit$test.fit)

  expect_s3_class(fam, "extended.family")
  expect_equal(fam$family, "Tweedie")
  # freshly constructed, so it must not carry test.fit's estimated parameter
  expect_false(identical(
    environment(fam$getTheta), environment(stats::family(fit$test.fit)$getTheta)
  ))
})

test_that("resolve_candidate_family returns an independent object on each call", {
  fit <- fixture_cs1_nb()

  a <- resolve_candidate_family(fit$test.fit)
  b <- resolve_candidate_family(fit$test.fit)

  expect_false(identical(environment(a$getTheta), environment(b$getTheta)))
  a$putTheta(2)
  expect_false(isTRUE(all.equal(a$getTheta(TRUE), b$getTheta(TRUE))))
})

test_that("resolve_candidate_family falls back to NULL when the family cannot be resolved", {
  # a family expression naming something that no longer exists is not an error:
  # returning NULL leaves update() to its own no-override behaviour, which is
  # never worse than the pre-fix state
  fit <- fixture_cs1_tweedie()
  fit$test.fit$call$family <- as.name("no_such_family_object")

  expect_null(resolve_candidate_family(fit$test.fit))
})

# ---- clone_independent_family -----------------------------------------------

test_that("clone_independent_family returns plain families unchanged", {
  # a plain glm family has no mutable per-instance state to clone
  fam <- stats::gaussian()
  expect_identical(clone_independent_family(fam), fam)

  fam <- stats::binomial()
  expect_identical(clone_independent_family(fam), fam)
})

test_that("clone_independent_family gives an extended family its own mutable state", {
  fam <- mgcv::nb()
  clone <- clone_independent_family(fam)

  expect_s3_class(clone, "extended.family")
  expect_equal(clone$family, fam$family)
  expect_false(identical(environment(clone$getTheta), environment(fam$getTheta)))

  # changing the clone must not reach back into the original
  original.theta <- fam$getTheta(TRUE)
  clone$putTheta(log(5))
  expect_equal(fam$getTheta(TRUE), original.theta)
})

test_that("clone_independent_family leaves namespaced environments shared", {
  # only anonymous (per-instance) closure environments are duplicated. Named
  # environments -- package namespaces -- must stay shared, or the family's own
  # functions lose access to the package internals they call.
  fam <- mgcv::nb()
  clone <- clone_independent_family(fam)

  is.fn <- vapply(clone, is.function, logical(1))
  envs <- lapply(clone[is.fn], environment)
  named <- vapply(envs, function(e) environmentName(e) != "", logical(1))
  for (nm in names(clone)[is.fn][named]) {
    expect_identical(environment(clone[[nm]]), environment(fam[[nm]]), info = nm)
  }
})

test_that("a cloned extended family still fits", {
  # the clone has to remain usable, not merely independent: nb()'s postproc
  # calls an unexported mgcv helper through its parent environment chain
  use.dat <- fixture_cs1_data()
  use.dat$Herbivore.abundance <- round(use.dat$Herbivore.abundance)

  fam <- clone_independent_family(mgcv::nb())
  fit <- mgcv::gam(
    Herbivore.abundance ~ s(depth, k = 3, bs = "cr"), family = fam, data = use.dat
  )

  expect_s3_class(fit, "gam")
  expect_true(is.finite(stats::AIC(fit)))
})

# ---- compute_model_weights --------------------------------------------------

test_that("compute_model_weights adds deltas and Akaike weights", {
  mod.data.out <- data.frame(
    modname = c("null", "a", "b"),
    AICc = c(100, 102, 110),
    BIC = c(105, 104, 120),
    r2.vals = c(0.1, 0.4, 0.2)
  )

  out <- compute_model_weights(mod.data.out, report.unique.r2 = FALSE)

  expect_equal(out$delta.AICc, c(0, 2, 10))
  expect_equal(out$delta.BIC, c(1, 0, 16))
  expect_equal(out$wi.AICc, round(wi(c(100, 102, 110)), 3))
  expect_equal(sum(out$wi.AICc), 1, tolerance = 1e-3)
  expect_false("r2.vals.unique" %in% colnames(out))
})

test_that("compute_model_weights subtracts the null r2 when asked", {
  mod.data.out <- data.frame(
    modname = c("null", "a", "b"),
    AICc = c(100, 102, 110),
    BIC = c(105, 104, 120),
    r2.vals = c(0.1, 0.4, 0.05)
  )

  out <- compute_model_weights(mod.data.out, report.unique.r2 = TRUE)

  expect_equal(out$r2.vals.unique, c(0, 0.3, -0.05))
})

test_that("compute_model_weights ignores NA when finding the minimum", {
  mod.data.out <- data.frame(
    modname = c("null", "a", "b"),
    AICc = c(100, NA, 110),
    BIC = c(105, NA, 120),
    r2.vals = c(0.1, NA, 0.2)
  )

  out <- compute_model_weights(mod.data.out, report.unique.r2 = FALSE)

  expect_equal(out$delta.AICc, c(0, NA, 10))
  expect_true(is.na(out$wi.AICc[2]))
  expect_equal(sum(out$wi.AICc, na.rm = TRUE), 1, tolerance = 1e-3)
})

# ---- compute_variable_importance --------------------------------------------

test_that("compute_variable_importance sums the best n models per predictor for 'min.n'", {
  # a is in three models, b in two, so min.n = 2 and a's score uses only its two
  # largest weights
  mod.data.out <- data.frame(
    wi.AICc = c(0.5, 0.3, 0.15, 0.05),
    wi.BIC = c(0.4, 0.4, 0.1, 0.1),
    a = c(1, 1, 1, 0),
    b = c(1, 0, 0, 1)
  )

  vi <- compute_variable_importance(mod.data.out, included.vars = c("a", "b"),
                                     VI.mods = "min.n")

  expect_equal(vi$aic$variable.weights.raw, c(a = 0.8, b = 0.55))
  expect_equal(vi$bic$variable.weights.raw, c(a = 0.8, b = 0.5))
})

test_that("compute_variable_importance sums every model for 'all'", {
  mod.data.out <- data.frame(
    wi.AICc = c(0.5, 0.3, 0.15, 0.05),
    wi.BIC = c(0.4, 0.4, 0.1, 0.1),
    a = c(1, 1, 1, 0),
    b = c(1, 0, 0, 1)
  )

  vi <- compute_variable_importance(mod.data.out, included.vars = c("a", "b"),
                                     VI.mods = "all")

  expect_equal(vi$aic$variable.weights.raw, c(a = 0.95, b = 0.55))
  expect_equal(vi$bic$variable.weights.raw, c(a = 0.9, b = 0.5))
})

# ---- worker_packages() ------------------------------------------------------

test_that("worker_packages carries what update() needs and not gamm4", {
  # FSSgam_package#14. mgcv is necessary: fit_mod_l() refits through
  # stats::update(), which for a gam() test.fit resolves the stored call's
  # `gam` symbol by lexical lookup ending at the worker's search path. MuMIn is
  # not necessary and is kept defensively, since a uGamm test.fit dispatches to
  # MuMIn's own update.gamm method, whose body runs inside that namespace.
  #
  # gamm4 is absent because nothing in this package calls it and naming it would
  # attach it on every worker for every model set.
  pkgs <- FSSgam:::worker_packages()

  expect_true(all(c("mgcv", "MuMIn", "FSSgam") %in% pkgs))
  expect_false("gamm4" %in% pkgs)
})

test_that("gamm4 is not a hard dependency of this package", {
  # The invariant that makes the above effective. An @importFrom directive makes
  # the FSSgam namespace load gamm4 wherever it is loaded, workers included, so
  # leaving gamm4 out of the .packages vector achieves nothing on its own. The
  # directive and the DESCRIPTION entry had to go together, and R CMD check
  # covers each direction separately: an Imports entry with no NAMESPACE import
  # raises the "All declared Imports should be used" NOTE, and a NAMESPACE
  # import with no DESCRIPTION entry raises "Namespace dependency missing from
  # DESCRIPTION Imports/Depends entries", which is an ERROR. This expectation
  # covers neither; it covers gamm4 being made a hard dependency again.
  desc <- utils::packageDescription("FSSgam")
  hard <- paste(desc$Depends, desc$Imports, desc$LinkingTo)
  expect_false(grepl("gamm4", hard, fixed = TRUE))
})
