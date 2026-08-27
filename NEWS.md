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
  (#7).
* Bug fix: `generate_model_set()` (and `full_subsets_gam()`, which calls it)
  failed on model sets with exactly one predictor. `use.dat[, all.predictors]`
  dropped to a bare vector before being passed to `check_correlations()` /
  `check_non_linear_correlations()`, which report on the columns of a
  data.frame, so the call failed with an unrelated "argument is of length
  zero" (or, for `non.linear.correlations = TRUE`, "invalid 'row.names'
  length") error instead of building a 1x1 correlation matrix.
  `check_non_linear_correlations()` now also returns the 1x1 unit matrix for a
  single-column input, matching `check_correlations()`.
* Bug fix: `generate_model_set()` emitted two spurious
  "no non-missing arguments to max; returning -Inf" warnings whenever factor
  predictors were supplied with `pred.vars.cont = NA` (the documented way to
  run without smooth predictors). A phantom `NA.by.<factor>` interaction term
  was being constructed from the `NA` itself; it was discarded again before
  the model set was returned, so only the warnings were visible.

* Added a test-coverage GitHub Actions workflow
  (`.github/workflows/test-coverage.yaml`, using `covr` + Codecov) and a
  coverage badge to the README (#3).

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
