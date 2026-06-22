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
