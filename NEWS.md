# FSSgam (development version)

* Completed the snake_case rename across the public API.
  `full.subsets.gam()` is now a deprecated wrapper around the new
  `full_subsets_gam()` (same treatment as `generate.model.set()` and
  `fit.model.set()`, since it's used directly in the companion docs
  repository's published case studies). `check.correlations()`,
  `check.non.linear.correlations()`, `build.inclusion.mat()`,
  `extract.mod.dat()`, and `fit.mod.l()` are renamed outright to
  `check_correlations()`, `check_non_linear_correlations()`,
  `build_inclusion_mat()`, `extract_mod_dat()`, and `fit_mod_l()` with no
  backward-compatible alias -- a breaking change for any code calling these
  five directly, but none of them appear anywhere in the companion repo's
  published vignettes/FAQ.
* Bug fix: `fit_model_set()` (and `full.subsets.gam()`/`full_subsets_gam()`,
  which call it) failed to fit most candidate models when `test.fit` used
  an mgcv extended family that estimates an extra parameter, such as
  `nb()` (negative binomial). Each candidate refit reused the already-fitted
  family object from `test.fit`, which carries that estimated parameter in
  its own mutable environment; re-using it warm-started every refit from an
  unrelated formula's estimate and destabilised the fit for most (but not
  all) candidates. Candidate models now get a fresh, independent family
  object for every fit, the same way the null model already did -- and,
  unlike an earlier attempt at this same fix, this now holds regardless of
  how `family` was originally specified on `test.fit` (a literal call such
  as `family = nb()`, or a variable/list element such as
  `family = my.families[[2]]`), and no longer crashes every candidate under
  `parallel = TRUE` when family was supplied via a variable.
* Bug fix: `generate_model_set()` (and `full_subsets_gam()`, which calls it)
  did not validate the `null.terms` argument, so passing anything other than
  a single character string (e.g. `NA`, `NULL`, a number, a logical, a
  character vector of length > 1, or a factor) produced a confusing,
  type-dependent low-level error -- or, for a bare numeric/logical value,
  was silently accepted and spliced into a nonsense model formula instead of
  failing. `null.terms` is now validated up front with a clear error message
  (beckyfisher/FSSgam#7).
* Bug fix: `generate_model_set()` (and `full_subsets_gam()`, which calls it)
  failed on model sets with exactly one predictor. `use.dat[, all.predictors]`
  dropped to a bare vector before being passed to `check_correlations()` /
  `check_non_linear_correlations()`, which report on the columns of a
  data.frame, so the call failed with an unrelated "argument is of length
  zero" (or, for `non.linear.correlations = TRUE`, "invalid 'row.names'
  length") error instead of building a 1x1 correlation matrix.
  `check_non_linear_correlations()` now also returns the 1x1 unit matrix for a
  single-column input, matching `check_correlations()`.
* Bug fix: `generate_model_set()` silently dropped the interaction term from
  every `linear.vars` predictor after the first. The pattern used to identify
  `.t.` interaction terms was built as `paste(linear.vars, ".t.")`, a vector
  whenever more than one linear predictor was supplied, and `grep()` used only
  its first element (with a warning). Affected candidates were still named as
  interactions but fitted as the corresponding factor main-effect model, so
  the model table reported duplicate fits under distinct names.
* Bug fix: supplying factor predictors with `pred.vars.cont = NA` -- the
  documented way to run without smooth predictors -- made
  `generate_model_set()` construct a phantom `NA.by.<factor>` interaction term
  out of the `NA` itself. With `max.predictors = 1` the term was discarded
  again before the model set was returned and the only symptom was two spurious
  "no non-missing arguments to max; returning -Inf" warnings, but from
  `max.predictors = 2` it survived into the candidate set as a model whose
  formula smoothed the literal `NA`.
* Bug fix: `case_study1` was documented as having 28 variables; it has 27.
* Bug fix: `generate_model_set()`'s `@return` documented a `generated.models`
  element that it has never returned, and omitted `n.mods`, `mod.formula`,
  `test.fit` and `included.vars`, which it does. `fit_model_set()`'s `@return`
  documented a `used.data` element that it has never returned -- `full_subsets_gam()` is
  what carries that.
* Bug fix: a candidate model named as an interaction between a factor and a
  linear predictor was fitted as the factor main effect alone, whenever the
  `.t.` term came from the list form of `factor.smooth.interactions` naming a
  linear predictor absent from `linear.vars`. The model table then reported two
  identical fits under different names.
* Substantially expanded the `testthat` suite, from 105 passing expectations to
  447, and from 71.33% to 94.27% line coverage (FSSgam_package#5). New coverage
  includes: every non-default argument of `generate_model_set()`,
  `fit_model_set()` and `full_subsets_gam()`; the `save.model.fits = FALSE`
  fitting path; `gamm4`, `uGamm` and `gamm` test.fits; the `cyclic.vars`,
  `linear.vars` and `bs.arg` formula construction; the bundled datasets; every
  documented error and warning message; and numerical snapshot tests pinning
  fitted values for five end-to-end scenarios.
* Added a test-coverage GitHub Actions workflow
  (`.github/workflows/test-coverage.yaml`, using `covr` + Codecov) and a
  coverage badge to the README (FSSgam_package#3).

* Bug fix: `fit_model_set()` (and `full_subsets_gam()`, which calls it)
  failed on a model set with exactly one predictor, with `'x' must be an
  array of at least two dimensions`. The variable importance calculation
  indexed the model table by column name without `drop = FALSE`, so a single
  predictor collapsed it to a vector before `colSums()`. `generate_model_set()`
  was fixed for the single-predictor case earlier in this development cycle;
  `fit_model_set()` was not.
* Bug fix: with `VI.mods = 'min.n'`, a predictor present in no successfully
  fitted model was given the weight of the single best model containing any
  predictor, rather than none. The number of models to sum over is zero in
  that case, and `1:0` counts backwards.
* Bug fix: with `save.model.fits = FALSE` and `parallel = TRUE`, a single
  candidate that could not be fitted aborted the whole run instead of being
  recorded in `failed.models`. The `foreach()` on that path had its
  `.errorhandling` disabled; the `save.model.fits = TRUE` path was unaffected.
* Bug fix: an unrecognised `r2.type` was accepted silently and produced a
  column of `NA` r2 values. `fit_model_set()` and `full_subsets_gam()` now
  reject anything other than `"r2.lm.est"`, `"r2"` or `"dev"`. `extract_mod_dat()`
  is unchanged, since it documents `r2.type.` as passed through.
* Bug fix: an unrecognised `VI.mods` failed with `object 'aic.var.weights'
  not found`. `fit_model_set()` and `full_subsets_gam()` now reject anything
  other than `"min.n"` or `"all"`.
* Bug fix: `full_subsets_gam()`'s deprecated `factor.interactions`,
  `smooth.interactions` and `size` arguments raised an error on values they
  are documented to accept. A character vector of length greater than one
  gave `the condition has length > 1`, and `NA` gave `missing value where
  TRUE/FALSE needed`. All three are now detected with `missing()`, which is
  correct for every type, length and `NA`. As a result
  `smooth.interactions = NA` now warns about the deprecation, which it did
  not do reliably before.
* Behaviour change: `full_subsets_gam()`'s `max.models` default changes from
  500 to 200, matching `fit_model_set()`. The two disagreed and neither
  documented its default, so a candidate set of between 201 and 500 models
  saved its model fits through `full_subsets_gam()` but not through
  `generate_model_set()` plus `fit_model_set()`. Anyone relying on the old
  wrapper default for a set in that range will now receive the
  "model fits will not be saved" warning and an empty `success.models`; pass
  `max.models = 500` explicitly to keep the previous behaviour. The default
  is now stated in the documentation of both functions.
* `full_subsets_gam()` gains `VI.mods`, which it previously did not forward
  to `fit_model_set()`, so `VI.mods = 'all'` was unreachable through the
  wrapper (FSSgam_package#7).
* `fit_model_set()` and `full_subsets_gam()` gain `progress`, defaulting to
  `interactive()`. The progress bar was previously written to stdout
  unconditionally, so suppressing it in a script or report required wrapping
  the call in `capture.output()` (FSSgam_package#9).

* `check_correlations()` and `check_non_linear_correlations()` no longer emit
  spurious "NaNs produced" warnings when a factor-factor pair is perfectly
  separated. The deviance is read directly off the `nnet::multinom()` fit
  rather than from its `summary()`, which computes standard errors that were
  discarded (FSSgam_package#10). Reported correlation values are unchanged.
* `check_correlations()`'s factor-factor diagonal is now exactly 1. It was
  previously a fitted multinomial pseudo-R2 of about 0.999999, because each
  factor was regressed on itself (FSSgam_package#12). Off-diagonal values are
  unchanged.

* Bug fix: `max.predictors = 1` with `factor.factor.interactions` built
  interaction columns it should not have. The enumeration was written as
  `for (i in 2:n)`, and at `n = 1` that sequence counts backwards to
  `c(2, 1)`, so a size-1 "combination" was enumerated as well. Its pasted name
  is just the original variable name, so `used.data` gained a duplicate column
  per factor, and four "no non-missing arguments to max; returning -Inf"
  warnings were emitted. No interaction terms are now generated at
  `max.predictors = 1`, with a warning that says so; on a two-factor set
  `used.data` loses three spurious columns.
* Bug fix: `factor.factor.interactions` given as a character vector failed with
  `arguments imply differing number of rows` when every named combination
  exceeded `cov.cutoff`. It now warns and continues, as the `TRUE` form always
  did. The warning text for both forms now names `cov.cutoff` and its value,
  replacing a message that referred to a non-existent `cor.cuttoff`.

* Behaviour change: the collinearity screen that decides which interaction
  terms are built now tests both halves of the correlation matrix in every
  case. Three of the four places that perform it tested only the upper
  triangle, so a pair could be admitted whose correlation exceeded
  `cov.cutoff` in the other direction. The matrices are not symmetric:
  `check_correlations()` estimates a factor-factor value by fitting
  `multinom()` separately in each direction, and
  `check_non_linear_correlations()` is asymmetric by construction. In
  practice this affects only a pair whose two estimates straddle
  `cov.cutoff`, which for factor pairs means a difference of under 0.001; a
  model set containing such a pair will now exclude the interaction it
  previously included.

* `factor.factor.interactions` given as a character vector no longer skips
  silently when every named pair exceeds `cov.cutoff`. An internal guard
  counted cells of the correlation matrix rather than pairs, so no interaction
  was built and nothing was reported; it now warns, as the `TRUE` form does.
  The model set produced is unchanged.

# FSSgam 1.0.0

Modernisation release ahead of CRAN submission.

* `DESCRIPTION` updated to use `Authors@R`, explicit `Imports` (`doSNOW`,
  `foreach`, `gamm4`, `mgcv`, `MuMIn`, `nnet`, `parallel`, `stats`, `utils`),
  and added the package repository to `URL`.
* Removed the unused `R.rsp` vignette builder declaration; this package does
  not ship vignettes (see the companion repository at
  <https://github.com/beckyfisher/FSSgam> for worked examples).
* Added runnable `@examples` to all exported functions.
* Added `data-raw/`, `build_package_FFSgam.R`, and `CLAUDE.md` to
  `.Rbuildignore`.
* Renamed the core implementations to `generate_model_set()` and
  `fit_model_set()`. `generate.model.set()` and `fit.model.set()` remain
  exported but are now deprecated wrappers (see `R/deprecated.R`) that warn
  via `.Deprecated()` and forward to the new names; `full.subsets.gam()`
  calls the new names internally so it does not itself emit deprecation
  warnings.
* Replaced `require()`/bare calls throughout `R/` with explicit
  `package::function()` calls, and `class(x) == "y"` checks with
  `inherits(x, "y")`. Note: smooth-term constructors (`s()`, `te()`, `t2()`)
  are kept unqualified inside model formulas, because mgcv resolves them by
  literal symbol name during formula parsing rather than through normal
  function dispatch -- `mgcv::s(...)` inside a formula does not work.
* Bug fix: `full.subsets.gam()`'s deprecated `size` argument (superseded by
  `max.predictors`) was a no-op -- it warned correctly but never actually
  fed its value into `max.predictors`. Calling `full.subsets.gam(size = n)`
  now actually constrains the model set to `n` predictors, as documented.
* Bug fix: `full.subsets.gam()`'s `used.data` and `predictor.correlations`
  output fields were always `NULL` (a field-name mismatch referenced
  `model.set$use.dat`/`$cor.matrix`, which don't exist -- the actual fields
  are `$used.data`/`$predictor.correlations`). Both are now populated as
  documented.
