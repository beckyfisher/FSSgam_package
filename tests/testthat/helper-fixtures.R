# Shared fixtures for the FSSgam test suite.
#
# testthat sources helper-*.R before every test file, so everything defined here
# is available to all of them. These constructors exist so that the eight-line
# use.dat / test.fit / model.set preamble is written once rather than repeated
# in every test.
#
# Each constructor returns a named list rather than just the fitted object, so a
# test can reach both the prepared data and the test.fit without rebuilding
# either. Fixtures are deliberately small (k = 3, two or three predictors) to
# keep the whole suite fast.
#
# Datasets are referenced as FSSgam::case_study1 rather than loaded with data(),
# because data()'s default envir is .GlobalEnv and CRAN checks require tests to
# leave nothing behind there.

# mgcv is attached, not just imported. Its extended-family constructors (tw(),
# nb(), ...) build their component closures in an environment whose parent is
# .GlobalEnv rather than the mgcv namespace, so functions such as ldTweedie()
# are only resolvable when mgcv is on the search path. Any fixture using tw() or
# nb() fails with "could not find function ldTweedie" otherwise. The package's
# own @examples attach mgcv for the same reason.
library(mgcv)

# ---- data preparation -------------------------------------------------------

# case_study1 with the columns the fixtures model as factors.
fixture_cs1_data <- function() {
  use.dat <- FSSgam::case_study1
  use.dat$site <- as.factor(use.dat$site)
  use.dat
}

# case_study2 as prepared in the companion repository's case-study-2 vignette:
# one taxon's rows, the grouping columns as factors, Abundance as the response.
fixture_cs2_data <- function(taxon = "BDS") {
  use.dat <- FSSgam::case_study2
  use.dat <- use.dat[use.dat$Taxa == taxon, ]
  use.dat$response <- use.dat$Abundance
  use.dat$Location <- factor(use.dat$Location)
  use.dat$Site <- factor(use.dat$Site)
  use.dat$Status <- factor(use.dat$Status)
  stats::na.omit(use.dat)
}

# case_study3 as prepared in the companion repository's case-study-3 vignette.
# month and lunar.date are the cyclic predictors.
fixture_cs3_data <- function() {
  use.dat <- FSSgam::case_study3
  use.dat$year <- factor(use.dat$year)
  use.dat$Sex <- factor(use.dat$Sex)
  use.dat$Species <- factor(use.dat$Species)
  use.dat
}

# coral_data as prepared in the companion repository's extra-examples vignette:
# a cbind(successes, failures) binomial response with Site as a random effect.
fixture_coral_data <- function() {
  use.dat <- stats::na.omit(FSSgam::coral_data[, c(
    "Site", "Survey", "bleach.pres", "dredge.pres", "av.wave", "Depth",
    "allcoral", "totalpoints"
  )])
  use.dat$successes <- use.dat$allcoral
  use.dat$failures <- use.dat$totalpoints - use.dat$allcoral
  use.dat$Site <- factor(use.dat$Site)
  use.dat
}

# ---- test.fit constructors --------------------------------------------------

# The canonical Gaussian gam fixture: the same test.fit the pre-existing tests
# built inline.
fixture_cs1_gaussian <- function() {
  use.dat <- fixture_cs1_data()
  list(
    use.dat = use.dat,
    test.fit = mgcv::gam(
      log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
      data = use.dat
    )
  )
}

# Tweedie: an extended family that estimates an extra parameter, so it also
# exercises resolve_candidate_family()/clone_independent_family().
fixture_cs1_tweedie <- function() {
  use.dat <- fixture_cs1_data()
  list(
    use.dat = use.dat,
    test.fit = mgcv::gam(
      Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
      family = mgcv::tw(), data = use.dat
    )
  )
}

# Negative binomial, on a rounded count response. This is the family from the
# beckyfisher/FSSgam#12 regression.
fixture_cs1_nb <- function() {
  use.dat <- fixture_cs1_data()
  use.dat$Herbivore.abundance <- round(use.dat$Herbivore.abundance)
  list(
    use.dat = use.dat,
    test.fit = mgcv::gam(
      Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
      family = mgcv::nb(), data = use.dat
    )
  )
}

# A Tweedie model set on case_study2, with two nested random effects and a
# factor predictor -- the factor-heavy fixture.
fixture_cs2_tweedie <- function() {
  use.dat <- fixture_cs2_data()
  list(
    use.dat = use.dat,
    test.fit = mgcv::gam(
      response ~ s(lobster, k = 3, bs = "cr") + s(Location, Site, bs = "re"),
      family = mgcv::tw(), data = use.dat
    )
  )
}

# The cyclic fixture. Gamma family, as in the case-study-3 vignette. k = 5
# because bs = 'cc' needs a larger basis than the default fixture k = 3, which
# mgcv silently raises with a warning.
fixture_cs3_cyclic <- function() {
  use.dat <- fixture_cs3_data()
  list(
    use.dat = use.dat,
    test.fit = mgcv::gam(
      GSI ~ s(lunar.date, k = 5, bs = "cc"),
      family = stats::Gamma(), data = use.dat
    )
  )
}

# The gamm4/uGamm fixture: a binomial mixed model fitted through MuMIn::uGamm
# with lme4 = TRUE, giving an object of class c("gamm4", "gamm").
#
# k = 4, not the k = 3 used everywhere else: gamm4 represents each smooth as a
# random effect whose grouping factor has (number of basis columns) levels, and
# for k = 3 that is a single level, which lme4's checkNlevels() rejects with
# "grouping factors must have > 1 sampled level". Any gamm4 fixture therefore
# needs k >= 4.
fixture_coral_ugamm <- function() {
  use.dat <- fixture_coral_data()
  list(
    use.dat = use.dat,
    test.fit = MuMIn::uGamm(
      cbind(successes, failures) ~ s(Depth, k = 4, bs = "cr"),
      family = stats::binomial(), random = ~ (1 | Site),
      data = use.dat, lme4 = TRUE
    )
  )
}

# The gamm (nlme) counterpart: uGamm with lme4 = FALSE, giving class "gamm".
fixture_coral_gamm <- function() {
  use.dat <- fixture_coral_data()
  use.dat$prop <- use.dat$allcoral / use.dat$totalpoints
  list(
    use.dat = use.dat,
    test.fit = MuMIn::uGamm(
      prop ~ s(Depth, k = 4, bs = "cr"),
      random = list(Site = ~1), data = use.dat, lme4 = FALSE
    )
  )
}

# ---- model set constructors -------------------------------------------------

# The canonical 8-model candidate set used by most tests. Any generate_model_set
# argument can be overridden through ..., so a test that needs one non-default
# argument does not have to restate the other six.
fixture_cs1_model_set <- function(fit = fixture_cs1_gaussian(), ...) {
  args <- list(
    use.dat = fit$use.dat,
    test.fit = fit$test.fit,
    pred.vars.cont = c("complexity", "depth"),
    pred.vars.fact = "ZONE",
    null.terms = "s(site,bs='re')",
    max.predictors = 2,
    k = 3
  )
  overrides <- list(...)
  args[names(overrides)] <- overrides
  do.call(generate_model_set, args)
}

# ---- utilities --------------------------------------------------------------

# fit_model_set() writes a txtProgressBar to stdout on every call, which buries
# the testthat reporter's own output. Warnings and errors are unaffected --
# capture.output() only redirects stdout -- so expect_warning()/expect_error()
# still work through this wrapper.
#
# Four groups of tests call fit_model_set()/full_subsets_gam() directly rather
# than through these wrappers, deliberately: the two beckyfisher/FSSgam#10/#12
# family resolution tests in test-fit_model_set.R, the
# used.data/predictor.correlations regression at the top of
# test-full_subsets_gam.R, the size= regression at the bottom of it, and
# everything in test-deprecated.R. FSSgam_package#5 asks for those to be
# preserved as they stand, so they keep their original call form, and with it
# their progress output.
fit_quietly <- function(...) {
  out <- NULL
  utils::capture.output(out <- fit_model_set(...))
  out
}

# Same, for full_subsets_gam().
full_subsets_quietly <- function(...) {
  out <- NULL
  utils::capture.output(out <- full_subsets_gam(...))
  out
}

# Replaces one candidate's formula with one that cannot be fitted, keeping its
# name. Used to exercise the partial-failure paths: injecting the failure into
# the model set is more reliable than trying to find real data that makes some
# candidates fail and not others.
#
# The default modname is the C-collation spelling. Candidate names are built by
# sorting term names (see FSSgam_package#8), so the same candidate is called
# "depth+ZONE" in a UTF-8 locale. testthat forces LC_COLLATE=C, so the default
# holds throughout the suite, but pass modname explicitly if calling this from
# an interactive session.
break_one_candidate <- function(model.set, modname = "ZONE+depth") {
  stopifnot(modname %in% names(model.set$mod.formula))
  model.set$mod.formula[[modname]] <-
    stats::as.formula("~ s(not_a_column, k = 3, bs = 'cr') + s(site, bs = 're')")
  model.set
}

# One-line deparse of a formula or call, with the same joining behaviour as
# deparse1(): width.cutoff = 500L, and a space between any lines that still wrap.
#
# deparse1() would do this, but it was added in R 4.0.0 and DESCRIPTION declares
# R (>= 3.5). The tests must run on the declared minimum, so they cannot use it.
deparse_one <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

# Muffles only nnet's "NaNs produced" warning, leaving every other warning to
# reach the test.
#
# check_correlations() takes its factor-factor estimate from a summary() of a
# multinom() fit, and that summary computes standard errors it then discards. On
# a perfectly separated fit -- which factor.factor.interactions guarantees, since
# a pasted interaction column is by construction perfectly predicted by its own
# components -- the variance-covariance matrix has negative diagonal entries and
# sqrt() warns (FSSgam_package#10). Whether nnet's optimiser reaches that state
# is version dependent: it does on the CI runners and does not on nnet 7.3-20
# here.
#
# A blanket suppressWarnings() would also hide the "no non-missing arguments to
# max" warnings that this package's own regression tests rely on, so the filter
# is on the message.
suppress_nnet_nans <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("NaNs produced", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

# Guard for every test that starts a real doSNOW cluster.
#
# Two conditions must hold before such a test is meaningful, and one before it
# is safe:
#
#  - FSSgam must be installed, not loaded by pkgload. A doSNOW worker is a fresh
#    R process that loads FSSgam with library() from the installed library path;
#    it never sees pkgload's in-memory copy, so running these against a
#    devtools::load_all() copy tests nothing (CLAUDE.md Section 6, Phase 12).
#  - The run must have opted in explicitly. doSNOW cluster startup has been
#    observed to stall indefinitely whenever another process is working the
#    machine hard -- including for a trivial foreach() containing no FSSgam code
#    -- so an unattended stall consumes the entire runtime of an R CMD check,
#    covr or CI job.
#
# To run them:
#
#   R CMD INSTALL .
#   FSSGAM_TEST_PARALLEL=true NOT_CRAN=true Rscript -e \
#     'library(FSSgam); testthat::test_dir("tests/testthat", package = "FSSgam")'
skip_unless_parallel_opt_in <- function() {
  if (!identical(tolower(Sys.getenv("FSSGAM_TEST_PARALLEL")), "true")) {
    testthat::skip("set FSSGAM_TEST_PARALLEL=true to run the parallel tests")
  }
  if (exists(".__DEVTOOLS__", envir = asNamespace("FSSgam"), inherits = FALSE)) {
    testthat::skip(
      "FSSgam is loaded via pkgload; parallel tests need an installed package"
    )
  }
}
