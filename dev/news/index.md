# Changelog

## FSSgam (development version)

## FSSgam 1.0.0

Modernisation release ahead of CRAN submission.

- `DESCRIPTION` updated to use `Authors@R`, explicit `Imports`
  (`doSNOW`, `foreach`, `gamm4`, `mgcv`, `MuMIn`, `nnet`, `parallel`,
  `stats`, `utils`), and added the package repository to `URL`.
- Removed the unused `R.rsp` vignette builder declaration; this package
  does not ship vignettes (see the companion repository at
  <https://github.com/beckyfisher/FSSgam> for worked examples).
- Added runnable `@examples` to all exported functions.
- Added `data-raw/`, `build_package_FFSgam.R`, and `CLAUDE.md` to
  `.Rbuildignore`.
- Renamed the core implementations to
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  and
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md).
  [`generate.model.set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate.model.set.md)
  and
  [`fit.model.set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit.model.set.md)
  remain exported but are now deprecated wrappers (see `R/deprecated.R`)
  that warn via
  [`.Deprecated()`](https://rdrr.io/r/base/Deprecated.html) and forward
  to the new names;
  [`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full.subsets.gam.md)
  calls the new names internally so it does not itself emit deprecation
  warnings.
- Replaced [`require()`](https://rdrr.io/r/base/library.html)/bare calls
  throughout `R/` with explicit `package::function()` calls, and
  `class(x) == "y"` checks with `inherits(x, "y")`. Note: smooth-term
  constructors ([`s()`](https://rdrr.io/pkg/mgcv/man/s.html),
  [`te()`](https://rdrr.io/pkg/mgcv/man/te.html),
  [`t2()`](https://rdrr.io/pkg/mgcv/man/t2.html)) are kept unqualified
  inside model formulas, because mgcv resolves them by literal symbol
  name during formula parsing rather than through normal function
  dispatch – `mgcv::s(...)` inside a formula does not work.
- Bug fix:
  [`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full.subsets.gam.md)’s
  deprecated `size` argument (superseded by `max.predictors`) was a
  no-op – it warned correctly but never actually fed its value into
  `max.predictors`. Calling `full.subsets.gam(size = n)` now actually
  constrains the model set to `n` predictors, as documented.
- Bug fix:
  [`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full.subsets.gam.md)’s
  `used.data` and `predictor.correlations` output fields were always
  `NULL` (a field-name mismatch referenced
  `model.set$use.dat`/`$cor.matrix`, which don’t exist – the actual
  fields are `$used.data`/`$predictor.correlations`). Both are now
  populated as documented.
