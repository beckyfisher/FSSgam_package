# The model selection criterion: which log-likelihood AICc and BIC are built
# from. Covers FSSgam_package#42 (mgcv's censored families report a criterion
# not built from a censored log-likelihood) and FSSgam_package#44 (a
# quasi-likelihood has no criterion at all).

# ---- reading the family off a fit -------------------------------------------

test_that("fit_family_name reads the family of each supported fit class", {
  expect_equal(FSSgam:::fit_family_name(fixture_cs1_gaussian()$test.fit),
               "gaussian")
  # An extended family reports its estimated parameter as part of the name once
  # fitted, so this is not the bare "cnorm".
  expect_match(FSSgam:::fit_family_name(fixture_cnorm()$test.fit),
               "^cnorm\\(")
  expect_equal(FSSgam:::fit_family_name(fixture_cs1_quasipoisson()$test.fit),
               "quasipoisson")
})

test_that("fit_family_name reads a uGamm fit, which stats::family cannot", {
  # stats::family() re-evaluates the recorded call, whose first element is not a
  # function symbol for a uGamm(lme4 = FALSE) fit, so it errors rather than
  # returning a family. That is a property of MuMIn, not of this package, and is
  # not asserted here; what is asserted is that fit_family_name() reads the
  # family regardless.
  gamm.fit <- fixture_coral_gamm()$test.fit
  expect_equal(FSSgam:::fit_family_name(gamm.fit), "gaussian")

  ugamm4.fit <- fixture_coral_ugamm()$test.fit
  expect_equal(FSSgam:::fit_family_name(ugamm4.fit), "binomial")
})

test_that("fit_family_name falls back to stats::family for a plain glm", {
  glm.fit <- stats::glm(c(1, 3, 2, 5, 4) ~ c(1, 2, 3, 4, 5),
                        family = stats::poisson())
  expect_equal(FSSgam:::fit_family_name(glm.fit), "poisson")
})

test_that("fit_family_name returns NA where no family can be read", {
  expect_true(is.na(FSSgam:::fit_family_name(list(a = 1))))
  expect_true(is.na(FSSgam:::fit_family_name(structure(list(), class = "unknown"))))
})

test_that("is_censored_family matches a fitted censored family name", {
  expect_true(FSSgam:::is_censored_family("cnorm"))
  expect_true(FSSgam:::is_censored_family("clog"))
  expect_true(FSSgam:::is_censored_family("cnorm(0.524)"))
  expect_true(FSSgam:::is_censored_family("clog(0.275)"))
  expect_false(FSSgam:::is_censored_family("gaussian"))
  expect_false(FSSgam:::is_censored_family("cloglog"))
  expect_false(FSSgam:::is_censored_family("cnormal"))
  expect_false(FSSgam:::is_censored_family(NA_character_))
})

# ---- the censored log-likelihood --------------------------------------------

test_that("censored_loglik agrees with an independent computation from mgcv's own components", {
  # censored_loglik_oracle() reads the family's ls and dev.resids slots, which
  # share no code with censored_loglik()'s densities and are not among the parts
  # FSSgam_package#42 reports as defective. See helper-fixtures.R.
  for (fixture in list(fixture_cnorm, fixture_clog)) {
    for (data in list(fixture_censored_data(seed = 3, n = 200),
                      fixture_mixed_censored_data(seed = 4, n = 200))) {
      fit <- fixture(data = data)$test.fit
      expect_equal(FSSgam:::censored_loglik(fit),
                   censored_loglik_oracle(fit),
                   tolerance = 1e-8,
                   info = fit$family$family)
    }
  }
})

test_that("censored_loglik handles a response with nothing censored", {
  d <- fixture_censored_data(seed = 5, n = 200, cens.frac = 0)
  fit <- fixture_cnorm(data = d)$test.fit
  expect_true(all(attr(fit$y, "censor") == fit$y))
  expect_equal(FSSgam:::censored_loglik(fit), censored_loglik_oracle(fit),
               tolerance = 1e-8)
})

test_that("censored_loglik accounts for prior weights", {
  d <- fixture_mixed_censored_data(seed = 6, n = 200)
  set.seed(6)
  d$wt <- sample(c(1, 2, 3), nrow(d), replace = TRUE)
  fit <- mgcv::gam(ycens ~ s(x, k = 5, bs = "cr"), data = d,
                   family = mgcv::cnorm(), method = "REML", weights = d$wt)
  expect_false(all(fit$prior.weights == 1))
  expect_equal(FSSgam:::censored_loglik(fit), censored_loglik_oracle(fit),
               tolerance = 1e-8)
})

test_that("getTheta(TRUE) returns the scale itself, not its logarithm", {
  # The condition for shipping censored_loglik() at all: it reads the scale
  # through getTheta(TRUE), and theta is documented as a log scale parameter.
  # What this test would catch is an mgcv change to getTheta()'s semantics,
  # which would otherwise give a wrong scale silently, rather than a repair of
  # the aic slots FSSgam_package#42 reports.
  set.seed(11)
  n <- 4000
  d <- data.frame(x = stats::runif(n))
  d$y <- stats::rnorm(n, 1 + sin(2 * pi * d$x), 0.5)
  d$ycens <- cbind(d$y, d$y)
  fit <- mgcv::gam(ycens ~ s(x, k = 10, bs = "cr"), data = d,
                   family = mgcv::cnorm(), method = "REML")
  expect_equal(fit$family$getTheta(TRUE), 0.5, tolerance = 0.05)
  expect_equal(fit$family$getTheta(TRUE), exp(fit$family$getTheta()),
               tolerance = 1e-8)
})

# ---- building the criterion from a log-likelihood ---------------------------

test_that("criterion_from_loglik reproduces MuMIn::AICc and stats::BIC exactly", {
  # Only the log-likelihood is replaced: handed the log-likelihood the default
  # route uses, it must return the default route's numbers.
  for (fit in list(fixture_cs1_gaussian()$test.fit,
                   fixture_cs1_tweedie()$test.fit,
                   fixture_cs1_nb()$test.fit)) {
    crit <- FSSgam:::criterion_from_loglik(fit, mgcv_reported_loglik)
    expect_equal(crit$AICc, MuMIn::AICc(fit), tolerance = 1e-10)
    expect_equal(crit$BIC, stats::BIC(fit), tolerance = 1e-10)
  }
})

test_that("criterion_from_loglik returns NA where the log-likelihood cannot be computed", {
  fit <- fixture_cs1_gaussian()$test.fit
  expect_equal(FSSgam:::criterion_from_loglik(fit, function(f) stop("no")),
               list(AICc = NA, BIC = NA))
  expect_equal(FSSgam:::criterion_from_loglik(fit, function(f) NA_real_),
               list(AICc = NA, BIC = NA))
  expect_equal(FSSgam:::criterion_from_loglik(fit, function(f) c(1, 2)),
               list(AICc = NA, BIC = NA))
})

# ---- resolving which log-likelihood the set is ranked on --------------------

test_that("resolve_criterion leaves an ordinary family on the default route", {
  expect_null(FSSgam:::resolve_criterion(fixture_cs1_gaussian()$test.fit, NULL))
  expect_null(FSSgam:::resolve_criterion(fixture_cs1_tweedie()$test.fit, NULL))
})

test_that("resolve_criterion supplies the censored log-likelihood, with a message", {
  fit <- fixture_cnorm()$test.fit
  expect_message(resolved <- FSSgam:::resolve_criterion(fit, NULL),
                 "censored log-likelihood")
  expect_identical(resolved, FSSgam:::censored_loglik)
})

test_that("resolve_criterion stops on a test.fit whose criterion is undefined", {
  fit <- fixture_cs1_quasipoisson()$test.fit
  expect_true(is.na(MuMIn::AICc(fit)))
  expect_error(FSSgam:::resolve_criterion(fit, NULL), "quasipoisson")
  expect_error(FSSgam:::resolve_criterion(fit, NULL), "cannot be ranked")
})

test_that("resolve_criterion validates a supplied logLik.fn against the test.fit", {
  fit <- fixture_cs1_gaussian()$test.fit
  expect_error(FSSgam:::resolve_criterion(fit, "not a function"),
               "must be a function")
  expect_error(FSSgam:::resolve_criterion(fit, function(f) stop("broken")),
               "failed on the test.fit")
  expect_error(FSSgam:::resolve_criterion(fit, function(f) NA_real_),
               "single finite numeric")
  expect_error(FSSgam:::resolve_criterion(fit, function(f) c(1, 2)),
               "single finite numeric")
  expect_identical(FSSgam:::resolve_criterion(fit, mgcv_reported_loglik),
                   mgcv_reported_loglik)
})

test_that("a supplied logLik.fn overrides the censored default without a message", {
  fit <- fixture_cnorm()$test.fit
  expect_no_message(
    resolved <- FSSgam:::resolve_criterion(fit, mgcv_reported_loglik))
  expect_identical(resolved, mgcv_reported_loglik)
})

# ---- extract_mod_dat --------------------------------------------------------

test_that("extract_mod_dat is unchanged when logLik.fn is NULL", {
  fit <- fixture_cs1_gaussian()$test.fit
  expect_equal(extract_mod_dat(fit), extract_mod_dat(fit, logLik.fn = NULL))
  expect_equal(extract_mod_dat(fit)$AICc, MuMIn::AICc(fit))
  expect_equal(extract_mod_dat(fit)$BIC, stats::BIC(fit))
})

test_that("extract_mod_dat builds the criterion from logLik.fn when supplied", {
  fit <- fixture_cs1_gaussian()$test.fit
  # Handed the default route's own log-likelihood, nothing changes.
  same <- extract_mod_dat(fit, logLik.fn = mgcv_reported_loglik)
  expect_equal(same$AICc, MuMIn::AICc(fit), tolerance = 1e-10)
  expect_equal(same$BIC, stats::BIC(fit), tolerance = 1e-10)
  # Only the criterion columns respond to it; the rest of the row does not.
  shifted <- extract_mod_dat(fit, logLik.fn = function(f) mgcv_reported_loglik(f) - 10)
  expect_equal(shifted$AICc, same$AICc + 20, tolerance = 1e-10)
  expect_equal(shifted$BIC, same$BIC + 20, tolerance = 1e-10)
  expect_equal(shifted[c("r2.vals", "edf", "edf.less.1")],
               same[c("r2.vals", "edf", "edf.less.1")])
})

# ---- fit_model_set ----------------------------------------------------------

test_that("fit_model_set is unchanged on an ordinary family", {
  fit <- fixture_cs1_gaussian()
  model.set <- fixture_cs1_model_set(fit)
  default <- fit_quietly(model.set, parallel = FALSE, progress = FALSE)
  supplied <- fit_quietly(model.set, parallel = FALSE, progress = FALSE,
                          logLik.fn = mgcv_reported_loglik)
  expect_equal(default$mod.data.out$AICc, supplied$mod.data.out$AICc,
               tolerance = 1e-8)
  expect_equal(default$mod.data.out$BIC, supplied$mod.data.out$BIC,
               tolerance = 1e-8)
  expect_equal(default$variable.importance, supplied$variable.importance)
})

test_that("fit_model_set stops before fitting on a quasi-likelihood test.fit", {
  fit <- fixture_cs1_quasipoisson()
  model.set <- fixture_cs1_model_set(fit)
  expect_error(fit_quietly(model.set, parallel = FALSE, progress = FALSE),
               "quasipoisson")
})

test_that("a quasi-likelihood model set is fitted and ranked when logLik.fn is supplied", {
  fit <- fixture_cs1_quasipoisson()
  model.set <- fixture_cs1_model_set(fit)
  poisson.loglik <- function(f) {
    sum(stats::dpois(f$y, f$fitted.values, log = TRUE))
  }
  out <- fit_quietly(model.set, parallel = FALSE, progress = FALSE,
                     logLik.fn = poisson.loglik)
  expect_false(anyNA(out$mod.data.out$AICc))
  expect_false(anyNA(out$mod.data.out$BIC))
  expect_false(anyNA(out$mod.data.out$wi.AICc))
  expect_equal(sum(out$mod.data.out$wi.AICc), 1, tolerance = 1e-3)
})

test_that("a censored model set is ranked on the censored log-likelihood", {
  # The defect FSSgam_package#42 reports, on its own seed-1 model set: the
  # criterion mgcv reports selects x alone, and the censored log-likelihood
  # selects x+z, which is the generating model.
  fit <- fixture_cnorm()
  model.set <- generate_model_set(use.dat = fit$use.dat, test.fit = fit$test.fit,
                                  pred.vars.cont = c("x", "z", "w", "v"),
                                  max.predictors = 2, k = 5)
  expect_message(corrected <- fit_quietly(model.set, parallel = FALSE,
                                          progress = FALSE),
                 "censored log-likelihood")

  best <- function(out) out$mod.data.out$modname[which.min(out$mod.data.out$AICc)]
  expect_equal(best(corrected), "x+z")

  # Nothing here asserts what mgcv's own reported criterion gives. Measured on
  # 2026-09-05, mgcv 1.9-4, that route selected "x", put the spread of
  # delta.AICc at 35.8 against 310.8, admitted four models within
  # delta.AICc <= 5 against two, and ranked the noise predictor v at 0.253 above
  # z at 0.182. Asserting any of that would fail if mgcv repaired its aic slots,
  # which is a repair and not a regression.

  # The criterion column is the censored log-likelihood's, model by model.
  expected <- vapply(corrected$success.models, FSSgam:::censored_loglik, numeric(1))
  k <- vapply(corrected$success.models,
              function(f) attr(stats::logLik(f), "df"), numeric(1))
  n <- nrow(fit$use.dat)
  expect_equal(unname(corrected$mod.data.out$AICc),
               unname(-2 * expected + 2 * k + 2 * k * (k + 1) / (n - k - 1)),
               tolerance = 1e-8)

  # Variable importance follows it: z is in the generating model and w and v are
  # noise, and only the corrected ranking puts z above both.
  vi <- corrected$variable.importance$aic$variable.weights.raw
  expect_true(vi[["z"]] > vi[["w"]])
  expect_true(vi[["z"]] > vi[["v"]])
  expect_equal(unname(vi[["x"]]), 1, tolerance = 1e-3)
})

test_that("both fitting paths give the same censored criterion", {
  fit <- fixture_cnorm(data = fixture_censored_data(seed = 8, n = 150))
  model.set <- generate_model_set(use.dat = fit$use.dat, test.fit = fit$test.fit,
                                  pred.vars.cont = c("x", "z"),
                                  max.predictors = 2, k = 5)
  saved <- suppressMessages(fit_quietly(model.set, parallel = FALSE,
                                        progress = FALSE))
  unsaved <- suppressMessages(fit_quietly(model.set, parallel = FALSE,
                                          progress = FALSE,
                                          save.model.fits = FALSE))
  expect_equal(saved$mod.data.out$AICc, unsaved$mod.data.out$AICc,
               tolerance = 1e-8)
  expect_equal(saved$mod.data.out$BIC, unsaved$mod.data.out$BIC,
               tolerance = 1e-8)
})

test_that("a failed candidate still yields an all-NA row under a censored family", {
  fit <- fixture_cnorm(data = fixture_censored_data(seed = 9, n = 150))
  model.set <- generate_model_set(use.dat = fit$use.dat, test.fit = fit$test.fit,
                                  pred.vars.cont = c("x", "z"),
                                  max.predictors = 2, k = 5)
  model.set$mod.formula[["x+z"]] <-
    stats::as.formula("~ s(not_a_column, k = 5, bs = 'cr')")
  unsaved <- suppressMessages(fit_quietly(model.set, parallel = FALSE,
                                          progress = FALSE,
                                          save.model.fits = FALSE))
  expect_true("x+z" %in% names(unsaved$failed.models))
  broken <- unsaved$mod.data.out[unsaved$mod.data.out$modname == "x+z", ]
  expect_true(is.na(broken$AICc))
  expect_true(is.na(broken$BIC))
  expect_false(anyNA(unsaved$mod.data.out$AICc[unsaved$mod.data.out$modname != "x+z"]))
})

test_that("censored_loglik reads an interval pair in either column order", {
  # mgcv orders the pair itself -- its dev.resids slot takes pmin and pmax of
  # the two columns for a finite interval -- so a reversed pair fits identically
  # and must give the same log-likelihood. Read as lower and upper as given it
  # returned NaN, and with it an NA criterion for every candidate in the set.
  d <- fixture_mixed_censored_data(seed = 12, n = 200)
  reversed <- d
  reversed$ycens <- d$ycens[, c(2, 1)]
  # Only the finite intervals are reversed; an infinite bound in the first
  # column is a different coding and is not what this covers.
  finite.interval <- is.finite(d$ycens[, 2]) & d$ycens[, 2] != d$ycens[, 1]
  reversed$ycens[!finite.interval, ] <- d$ycens[!finite.interval, ]
  expect_true(any(finite.interval))

  fit <- fixture_cnorm(data = d)$test.fit
  fit.reversed <- fixture_cnorm(data = reversed)$test.fit
  expect_equal(stats::coef(fit.reversed), stats::coef(fit), tolerance = 1e-6)
  expect_false(is.na(FSSgam:::censored_loglik(fit.reversed)))
  expect_equal(FSSgam:::censored_loglik(fit.reversed),
               censored_loglik_oracle(fit.reversed), tolerance = 1e-8)
  expect_equal(FSSgam:::censored_loglik(fit.reversed),
               FSSgam:::censored_loglik(fit), tolerance = 1e-6)
})

test_that("censored_loglik refuses a coding whose rows mgcv drops", {
  # mgcv selects its censored cases on the second column alone -- yat == -Inf
  # and yat == Inf -- so a row with an infinity in the first column matches
  # none of clog's four index sets and is dropped from the fit without a
  # message. Scoring it here would count observations the fit ignores, and the
  # offset depends on the fitted mean, so it is not constant across candidates.
  d <- fixture_mixed_censored_data(seed = 31, n = 200)
  infinite <- is.infinite(d$ycens[, 2])
  expect_true(any(infinite))
  d$ycens[infinite, ] <- d$ycens[infinite, c(2, 1)]
  fit <- fixture_clog(data = d)$test.fit
  expect_true(is.na(FSSgam:::censored_loglik(fit)))
  # And is reported before anything is fitted rather than as an NA column.
  expect_error(FSSgam:::resolve_criterion(fit, NULL), "could not be evaluated")
  expect_error(FSSgam:::resolve_criterion(fit, NULL), "second column")
})

test_that("censored_loglik does not underflow on an interval far into a tail", {
  # Computed as log(F(upper) - F(lower)) the difference underflows to zero far
  # from the mean and the whole sum becomes -Inf, which gives an NA criterion.
  d <- fixture_censored_data(seed = 13, n = 200)
  d$ycens[1, ] <- c(60, 60 + 1e-9)
  fit <- fixture_cnorm(data = d)$test.fit
  expect_true(is.finite(FSSgam:::censored_loglik(fit)))
  expect_equal(FSSgam:::censored_loglik(fit), censored_loglik_oracle(fit),
               tolerance = 1e-6)
})

test_that("censored_loglik agrees with the oracle on a heavily censored response", {
  # 85 per cent left censored. Censoring every row leaves nothing for mgcv to
  # estimate the smooth from and the fit itself fails, so this is as far as the
  # comparison can be taken.
  d <- fixture_censored_data(seed = 14, n = 300, cens.frac = 0.85)
  fit <- fixture_cnorm(data = d)$test.fit
  expect_gt(mean(attr(fit$y, "censor") == -Inf), 0.8)
  expect_equal(FSSgam:::censored_loglik(fit), censored_loglik_oracle(fit),
               tolerance = 1e-8)
})

test_that("resolve_criterion refuses a censored fit the built-in cannot evaluate", {
  fit <- fixture_cnorm()$test.fit
  # A censored coding the built-in returns no number for must be reported here,
  # before the fitting loop, not left to produce an NA criterion column.
  local_mocked_bindings(censored_loglik = function(fit) NaN, .package = "FSSgam")
  expect_error(FSSgam:::resolve_criterion(fit, NULL), "could not be evaluated")
})

test_that("fit_model_set stops where every fitted candidate has an NA criterion", {
  fit <- fixture_cs1_gaussian()
  model.set <- fixture_cs1_model_set(fit)
  # A logLik.fn that passes the test.fit check and fails on every candidate:
  # nothing before the fitting loop can detect it.
  # A constructor, not one closure: the counter must start again for the second
  # call, whose first evaluation is resolve_criterion()'s check on the test.fit.
  make_fussy <- function() {
    seen <- 0L
    function(f) {
      seen <<- seen + 1L
      if (seen == 1L) as.numeric(stats::logLik(f)) else NA_real_
    }
  }
  expect_error(fit_quietly(model.set, parallel = FALSE, progress = FALSE,
                           logLik.fn = make_fussy()),
               "cannot be ranked")
  # Both fitting paths reach it: the unsaved path records a failed model from
  # the fit rather than from an NA criterion, so a candidate that fitted and
  # was given none is not misreported as one that did not fit.
  expect_error(fit_quietly(model.set, parallel = FALSE, progress = FALSE,
                           save.model.fits = FALSE, logLik.fn = make_fussy()),
               "cannot be ranked")
})

test_that("a candidate fitted and given no criterion is reported, on both paths", {
  fit <- fixture_cs1_gaussian()
  model.set <- fixture_cs1_model_set(fit)
  # NA for the candidates containing ZONE, a value for every other fit.
  partial <- function(f) {
    terms <- attr(stats::terms(stats::formula(f)), "term.labels")
    if (any(grepl("ZONE", terms))) NA_real_ else as.numeric(stats::logLik(f))
  }
  saved <- expect_warning(
    fit_quietly(model.set, parallel = FALSE, progress = FALSE,
                logLik.fn = partial),
    "given no criterion")
  unsaved <- expect_warning(
    fit_quietly(model.set, parallel = FALSE, progress = FALSE,
                save.model.fits = FALSE, logLik.fn = partial),
    "given no criterion")
  # A candidate given no criterion fitted, so it is not a failed model on
  # either path, and the two paths agree on which candidates succeeded.
  expect_equal(length(saved$failed.models), 0L)
  expect_equal(length(unsaved$failed.models), 0L)
  expect_equal(names(unsaved$success.models), names(saved$success.models))
})

# ---- full_subsets_gam -------------------------------------------------------

test_that("full_subsets_gam forwards logLik.fn", {
  expect_true("logLik.fn" %in% names(formals(full_subsets_gam)))
  fit <- fixture_cs1_gaussian()
  default <- full_subsets_quietly(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"), null.terms = "s(site,bs='re')",
    max.predictors = 2, k = 3, parallel = FALSE, progress = FALSE)
  shifted <- full_subsets_quietly(
    use.dat = fit$use.dat, test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"), null.terms = "s(site,bs='re')",
    max.predictors = 2, k = 3, parallel = FALSE, progress = FALSE,
    logLik.fn = function(f) mgcv_reported_loglik(f) - 10)
  expect_equal(shifted$mod.data.out$AICc, default$mod.data.out$AICc + 20,
               tolerance = 1e-8)
})

test_that("full_subsets_gam stops on a quasi-likelihood test.fit", {
  fit <- fixture_cs1_quasipoisson()
  expect_error(
    full_subsets_quietly(
      use.dat = fit$use.dat, test.fit = fit$test.fit,
      pred.vars.cont = c("complexity", "depth"), null.terms = "s(site,bs='re')",
      max.predictors = 2, k = 3, parallel = FALSE, progress = FALSE),
    "quasipoisson")
})

# ---- the parallel path ------------------------------------------------------

test_that("logLik.fn reaches the doSNOW workers", {
  skip_unless_parallel_opt_in()
  fit <- fixture_cnorm(data = fixture_censored_data(seed = 10, n = 150))
  model.set <- generate_model_set(use.dat = fit$use.dat, test.fit = fit$test.fit,
                                  pred.vars.cont = c("x", "z"),
                                  max.predictors = 2, k = 5)
  sequential <- suppressMessages(fit_quietly(model.set, parallel = FALSE,
                                             progress = FALSE,
                                             save.model.fits = FALSE))
  parallel.out <- suppressMessages(fit_quietly(model.set, parallel = TRUE,
                                               n.cores = 2, progress = FALSE,
                                               save.model.fits = FALSE))
  expect_equal(parallel.out$mod.data.out$AICc, sequential$mod.data.out$AICc,
               tolerance = 1e-8)
})
