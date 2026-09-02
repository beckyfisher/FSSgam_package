# Changelog

## FSSgam (development version)

- `gamm4` moves from `Imports` to `Suggests`. **A fresh install of this
  package no longer installs `gamm4` or `lme4`. If you fit through
  `MuMIn::uGamm(lme4 = TRUE)`, as the companion repository’s case study
  2 does, install `gamm4` yourself.** Nothing in this package calls it:
  the only uses are `inherits(x, "gamm4")` class tests, which need no
  package, and an `@importFrom gamm4 gamm4` directive for a function
  never called. That directive is what made every parallel worker load
  `gamm4` and `lme4`, through the `FSSgam` namespace, for every model
  set regardless of the family fitted.

- `fit_model_set(parallel = TRUE)` no longer names `gamm4` in the
  packages it attaches on the `doSNOW` workers. Together with the change
  above this is the only direction recorded on FSSgam_package#14 that
  measured as reducing the dispatch stall, in which the call hangs
  indefinitely with no output and no error. Running
  `fit_model_set(parallel = TRUE)` on `case_study1` from installed
  copies of each version, alternating in blocks of ten so that both saw
  the same machine load, one cluster per fresh R process, gave 44 stalls
  in 140 runs before and 25 in 140 after (chi-squared p = 0.013, this
  host, 2026-09-01).

  The narrower comparison agrees. A `foreach()` whose body is `i * 2`,
  so that no package code runs at all, was run 50 times with `gamm4` in
  the packages the workers attach and 50 times without, alternating run
  by run: 34 stalls against 19 (chi-squared p = 0.005, same host and
  date).

  Three qualifications. That is a reduction and not a removal, and
  FSSgam_package#14 stays open. Both figures are one session’s
  measurement of a rate on one host, and an earlier run of the narrower
  comparison with its two arms taken one after the other rather than
  alternated gave a result an independent replication could not
  reproduce; only alternating designs are quoted here. And neither
  extends to a `uGamm(lme4 = TRUE)` test.fit: unserialising such an
  object loads `gamm4` and `lme4` on the worker by itself, before the
  packages are attached, so for those users the worker side is
  substantially unchanged.

  `check_correlations(parallel = TRUE)` benefits as well, without any
  change of its own: its workers load the `FSSgam` namespace and so used
  to load `gamm4` through it.

- Bug fix:
  [`check_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_correlations.md)
  estimated the correlation between two factor predictors from deviances
  computed on different sets of rows whenever either factor contained
  missing values. The fitted model was fitted on the pair’s complete
  cases and the intercept-only model on the whole column, so the ratio
  the estimate is built from was not comparable and the reported
  correlation was wrong by an amount that depended on how much of the
  column was missing. Both models are now fitted on the pair’s complete
  cases. Values are unchanged for data with no missing predictor values,
  which is every call made by
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md),
  since it rejects predictors containing `NA` (FSSgam_package#16).

- Bug fix:
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  screened the wrong pairs for collinearity when a supplied `cor.matrix`
  carried the same names in a different order in its rows and its
  columns. The correlation sub-matrix for a candidate was built from two
  positional indices, one in rowname order and one in colname order, so
  where the two orders differed the diagonal of that sub-matrix held
  cross correlations and the triangles held the ones from the original
  diagonal, and the candidate was dropped whatever its predictors’
  correlation. Both sites now index by name in one order. Measured on
  `master`, with a supplied matrix over `depth`, `complexity` and
  `ZONE`, `pred.vars.fact = "ZONE"` and everything else at its default:
  the matrix with its columns permuted returned 5 candidates where the
  same matrix unpermuted returned 9 at `max.predictors = 2` and 13 at
  `max.predictors = 3`. On this branch the two orders agree, at 9
  and 13. The same change means a factor named twice in `pred.vars.fact`
  no longer yields an interaction of that factor with itself
  (FSSgam_package#15).

- [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  no longer requires a user-supplied `cor.matrix` to name the hard coded
  factor interaction columns that `factor.factor.interactions` causes to
  be created. Which of those columns exist depends on the supplied
  matrix itself, through the collinearity screen that decides which
  factor combinations survive, so a user could not know in advance which
  names to supply; the call stopped partway through model set
  construction demanding them. Rows and columns for any interaction
  column the matrix does not carry are now computed from `use.dat` and
  appended, and every value the user supplied is kept as supplied. A
  predictor the user named and did not supply is still an error, so a
  misspelled name is still reported rather than quietly computed
  (FSSgam_package#15).

## FSSgam 1.1.0

- Completed the snake_case rename across the public API.
  [`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full.subsets.gam.md)
  is now a deprecated wrapper around the new
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)
  (same treatment as
  [`generate.model.set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate.model.set.md)
  and
  [`fit.model.set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit.model.set.md),
  since it’s used directly in the companion docs repository’s published
  case studies). `check.correlations()`,
  `check.non.linear.correlations()`, `build.inclusion.mat()`,
  `extract.mod.dat()`, and `fit.mod.l()` are renamed outright to
  [`check_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_correlations.md),
  [`check_non_linear_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_non_linear_correlations.md),
  [`build_inclusion_mat()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/build_inclusion_mat.md),
  [`extract_mod_dat()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/extract_mod_dat.md),
  and
  [`fit_mod_l()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_mod_l.md)
  with no backward-compatible alias – a breaking change for any code
  calling these five directly, but none of them appear anywhere in the
  companion repo’s published vignettes/FAQ.

- Bug fix:
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md)
  (and
  [`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full.subsets.gam.md)/[`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md),
  which call it) failed to fit most candidate models when `test.fit`
  used an mgcv extended family that estimates an extra parameter, such
  as [`nb()`](https://rdrr.io/pkg/mgcv/man/negbin.html) (negative
  binomial). Each candidate refit reused the already-fitted family
  object from `test.fit`, which carries that estimated parameter in its
  own mutable environment; re-using it warm-started every refit from an
  unrelated formula’s estimate and destabilised the fit for most (but
  not all) candidates. Candidate models now get a fresh, independent
  family object for every fit, the same way the null model already did –
  and, unlike an earlier attempt at this same fix, this now holds
  regardless of how `family` was originally specified on `test.fit` (a
  literal call such as `family = nb()`, or a variable/list element such
  as `family = my.families[[2]]`), and no longer crashes every candidate
  under `parallel = TRUE` when family was supplied via a variable.

- Bug fix:
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  (and
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md),
  which calls it) did not validate the `null.terms` argument, so passing
  anything other than a single character string (e.g. `NA`, `NULL`, a
  number, a logical, a character vector of length \> 1, or a factor)
  produced a confusing, type-dependent low-level error – or, for a bare
  numeric/logical value, was silently accepted and spliced into a
  nonsense model formula instead of failing. `null.terms` is now
  validated up front with a clear error message (beckyfisher/FSSgam#7).

- Bug fix:
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  (and
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md),
  which calls it) failed on model sets with exactly one predictor.
  `use.dat[, all.predictors]` dropped to a bare vector before being
  passed to
  [`check_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_correlations.md)
  /
  [`check_non_linear_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_non_linear_correlations.md),
  which report on the columns of a data.frame, so the call failed with
  an unrelated “argument is of length zero” (or, for
  `non.linear.correlations = TRUE`, “invalid ‘row.names’ length”) error
  instead of building a 1x1 correlation matrix.
  [`check_non_linear_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_non_linear_correlations.md)
  now also returns the 1x1 unit matrix for a single-column input,
  matching
  [`check_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_correlations.md).

- Bug fix:
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  silently dropped the interaction term from every `linear.vars`
  predictor after the first. The pattern used to identify `.t.`
  interaction terms was built as `paste(linear.vars, ".t.")`, a vector
  whenever more than one linear predictor was supplied, and
  [`grep()`](https://rdrr.io/r/base/grep.html) used only its first
  element (with a warning). Affected candidates were still named as
  interactions but fitted as the corresponding factor main-effect model,
  so the model table reported duplicate fits under distinct names.

- Bug fix: supplying factor predictors with `pred.vars.cont = NA` – the
  documented way to run without smooth predictors – made
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  construct a phantom `NA.by.<factor>` interaction term out of the `NA`
  itself. With `max.predictors = 1` the term was discarded again before
  the model set was returned and the only symptom was two spurious “no
  non-missing arguments to max; returning -Inf” warnings, but from
  `max.predictors = 2` it survived into the candidate set as a model
  whose formula smoothed the literal `NA`.

- Bug fix: `case_study1` was documented as having 28 variables; it has
  27.

- Bug fix:
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)’s
  `@return` documented a `generated.models` element that it has never
  returned, and omitted `n.mods`, `mod.formula`, `test.fit` and
  `included.vars`, which it does.
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md)’s
  `@return` documented a `used.data` element that it has never returned
  –
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)
  is what carries that.

- Bug fix: a candidate model named as an interaction between a factor
  and a linear predictor was fitted as the factor main effect alone,
  whenever the `.t.` term came from the list form of
  `factor.smooth.interactions` naming a linear predictor absent from
  `linear.vars`. The model table then reported two identical fits under
  different names.

- Substantially expanded the `testthat` suite, from 105 passing
  expectations at 1.0.0 to 568 at this release, and from 71.33% to 94.4%
  line coverage (FSSgam_package#5 covered the first 447 of those; the
  rest accompany the fixes below). New coverage includes: every
  non-default argument of
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md),
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md)
  and
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md);
  the `save.model.fits = FALSE` fitting path; `gamm4`, `uGamm` and
  `gamm` test.fits; the `cyclic.vars`, `linear.vars` and `bs.arg`
  formula construction; the bundled datasets; every documented error and
  warning message; and numerical snapshot tests pinning fitted values
  for five end-to-end scenarios.

- Added a test-coverage GitHub Actions workflow
  (`.github/workflows/test-coverage.yaml`, using `covr` + Codecov) and a
  coverage badge to the README (FSSgam_package#3).

- Bug fix:
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md)
  (and
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md),
  which calls it) failed on a model set with exactly one predictor, with
  `'x' must be an array of at least two dimensions`. The variable
  importance calculation indexed the model table by column name without
  `drop = FALSE`, so a single predictor collapsed it to a vector before
  [`colSums()`](https://rdrr.io/r/base/colSums.html).
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  was fixed for the single-predictor case earlier in this development
  cycle;
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md)
  was not.

- Bug fix: with `VI.mods = 'min.n'`, a predictor present in no
  successfully fitted model was given the weight of the single best
  model containing any predictor, rather than none. The number of models
  to sum over is zero in that case, and `1:0` counts backwards.

- Bug fix: with `save.model.fits = FALSE` and `parallel = TRUE`, a
  single candidate that could not be fitted aborted the whole run
  instead of being recorded in `failed.models`. The `foreach()` on that
  path had its `.errorhandling` disabled; the `save.model.fits = TRUE`
  path was unaffected.

- Bug fix: an unrecognised `r2.type` was accepted silently and produced
  a column of `NA` r2 values.
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md)
  and
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)
  now reject anything other than `"r2.lm.est"`, `"r2"` or `"dev"`.
  [`extract_mod_dat()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/extract_mod_dat.md)
  is unchanged, since it documents `r2.type.` as passed through.

- Bug fix: an unrecognised `VI.mods` failed with
  `object 'aic.var.weights' not found`.
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md)
  and
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)
  now reject anything other than `"min.n"` or `"all"`.

- Bug fix:
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)’s
  deprecated `factor.interactions`, `smooth.interactions` and `size`
  arguments raised an error on values they are documented to accept. A
  character vector of length greater than one gave
  `the condition has length > 1`, and `NA` gave
  `missing value where TRUE/FALSE needed`. All three are now detected
  with [`missing()`](https://rdrr.io/r/base/missing.html), which is
  correct for every type, length and `NA`. As a result
  `smooth.interactions = NA` now warns about the deprecation, which it
  did not do reliably before.

- Behaviour change:
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)’s
  `max.models` default changes from 500 to 200, matching
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md).
  The two disagreed and neither documented its default, so a candidate
  set of between 201 and 500 models saved its model fits through
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)
  but not through
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  plus
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md).
  Anyone relying on the old wrapper default for a set in that range will
  now receive the “model fits will not be saved” warning and an empty
  `success.models`; pass `max.models = 500` explicitly to keep the
  previous behaviour. The default is now stated in the documentation of
  both functions.

- [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)
  gains `VI.mods`, which it previously did not forward to
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md),
  so `VI.mods = 'all'` was unreachable through the wrapper
  (FSSgam_package#7).

- [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md)
  and
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)
  gain `progress`, defaulting to
  [`interactive()`](https://rdrr.io/r/base/interactive.html). The
  progress bar was previously written to stdout unconditionally, so
  suppressing it in a script or report required wrapping the call in
  [`capture.output()`](https://rdrr.io/r/utils/capture.output.html)
  (FSSgam_package#9).

- [`check_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_correlations.md)
  and
  [`check_non_linear_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_non_linear_correlations.md)
  no longer emit spurious “NaNs produced” warnings when a factor-factor
  pair is perfectly separated. The deviance is read directly off the
  [`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) fit
  rather than from its
  [`summary()`](https://rdrr.io/r/base/summary.html), which computes
  standard errors that were discarded (FSSgam_package#10). Reported
  correlation values are unchanged.

- [`check_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_correlations.md)’s
  factor-factor diagonal is now exactly 1. It was previously a fitted
  multinomial pseudo-R2 of about 0.999999, because each factor was
  regressed on itself (FSSgam_package#12). Off-diagonal values are
  unchanged.

- Bug fix: `max.predictors = 1` with `factor.factor.interactions` built
  interaction columns it should not have. The enumeration was written as
  `for (i in 2:n)`, and at `n = 1` that sequence counts backwards to
  `c(2, 1)`, so a size-1 “combination” was enumerated as well. Its
  pasted name is just the original variable name, so `used.data` gained
  a duplicate column per factor, and four “no non-missing arguments to
  max; returning -Inf” warnings were emitted. No interaction terms are
  now generated at `max.predictors = 1`, with a warning that says so; on
  a two-factor set `used.data` loses three spurious columns.

- Bug fix: `factor.factor.interactions` given as a character vector
  failed with `arguments imply differing number of rows` when every
  named combination exceeded `cov.cutoff`. It now warns and continues,
  as the `TRUE` form always did. The warning text for both forms now
  names `cov.cutoff` and its value, replacing a message that referred to
  a non-existent `cor.cuttoff`.

- Behaviour change: the collinearity screen that decides which
  interaction terms are built now tests both halves of the correlation
  matrix in every case. Three of the four places that perform it tested
  only the upper triangle, so a pair could be admitted whose correlation
  exceeded `cov.cutoff` in the other direction. The matrices are not
  symmetric:
  [`check_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_correlations.md)
  estimates a factor-factor value by fitting
  [`multinom()`](https://rdrr.io/pkg/mgcv/man/multinom.html) separately
  in each direction, and
  [`check_non_linear_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_non_linear_correlations.md)
  is asymmetric by construction. In practice this affects only a pair
  whose two estimates straddle `cov.cutoff`, which for factor pairs
  means a difference of under 0.001; a model set containing such a pair
  will now exclude the interaction it previously included.

- `factor.factor.interactions` given as a character vector no longer
  skips silently when every named pair exceeds `cov.cutoff`. An internal
  guard counted cells of the correlation matrix rather than pairs, so no
  interaction was built and nothing was reported; it now warns, as the
  `TRUE` form does. The model set produced is unchanged.

- Bug fix: a name in `cyclic.vars` was matched against the assembled
  model terms with an unanchored
  [`grep()`](https://rdrr.io/r/base/grep.html), so it captured every
  predictor whose name contained it, and the name was used as a regular
  expression. Declaring `depth` cyclic also fitted `depthx` with
  `bs = "cc"`, and a predictor name containing a full stop matched any
  character in that position. The affected models were fitted with the
  wrong smoothing basis, without error, under a candidate name that
  looked correct. The basis is now chosen from the variable names as
  each term is built.

- Bug fix: a [`te()`](https://rdrr.io/pkg/mgcv/man/te.html) term over
  three or more predictors was given only two `bs` values, because the
  code that assigned them took the first two variables. mgcv warns “bs
  wrong length and ignored” and substitutes its own default, so both
  `bs.arg` and any `cyclic.vars` were silently discarded for that term.
  Each marginal now carries its own basis. Reachable through
  `smooth.smooth.interactions` given as a character vector with
  `max.predictors >= 3`.

- Behaviour change: candidate model names are now sorted in byte order,
  via `sort(method = "radix")`, and no longer depend on the session’s
  collation locale (FSSgam_package#8). A model named `complexity+ZONE`
  in an `en_US.UTF-8` session was named `ZONE+complexity` in a C-locale
  one, so a saved analysis was not reproducible across machines with
  different locales. Names are now the C-locale form everywhere. The
  fitted results are identical, but `modname` values and the row order
  of `mod.data.out` change for anyone working in a non-C locale; a saved
  model table matched on `modname` needs regenerating.

- Behaviour change: a user-supplied `cor.matrix` now governs every stage
  that screens on correlation, not only the exclusion of assembled
  candidate models (FSSgam_package#13). There are three such stages:
  which factor-factor interaction columns are built, which
  [`te()`](https://rdrr.io/pkg/mgcv/man/te.html) smooth-smooth
  interaction terms are built, and which candidates survive. The first
  two previously recomputed correlations from `use.dat` and ignored the
  supplied matrix, so a user who supplied a matrix saying two predictors
  were uncorrelated still found their
  [`te()`](https://rdrr.io/pkg/mgcv/man/te.html) term absent. Users who
  do not supply a `cor.matrix` see no change: the continuous-continuous
  block does not depend on which other predictors are present, verified
  against both
  [`check_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_correlations.md)
  and
  [`check_non_linear_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_non_linear_correlations.md).
  A supplied matrix must now include any hard coded factor interaction
  it causes to be created, which
  [`?generate_model_set`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  already required; missing predictors are reported by name.

- Bug fix: a user-supplied `cor.matrix` did not replace the automatic
  estimate, it only overrode it.
  [`check_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_correlations.md)
  (or
  [`check_non_linear_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_non_linear_correlations.md))
  was called over every predictor and its result discarded whenever a
  matrix was supplied, so supplying one neither avoided the
  [`multinom()`](https://rdrr.io/pkg/mgcv/man/multinom.html) and
  [`gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) fits nor allowed a
  predictor of a class those functions reject – a `Date` column, say –
  to be used at all. The estimate is now computed only when no matrix is
  supplied (FSSgam_package#13).

- Bug fix: with `save.model.fits = FALSE`, a run in which no candidate
  fitted returned a model table of `NA` with `delta.AICc` and `wi.AICc`
  filled with `NaN`, rather than the “None of your models fitted
  successfully” error raised on the `save.model.fits = TRUE` path. Which
  of the two happened depended only on a memory setting. Both paths now
  raise the error.

- [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md)
  and
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)
  reject a `progress` value that is not a single `TRUE` or `FALSE`, and
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)
  now validates `VI.mods` at entry as it already did `r2.type`, so an
  unrecognised value is reported before the candidate set is built
  rather than after.

- `smooth.smooth.interactions` naming a column of `use.dat` that is not
  among the predictors is now reported as such. It was previously
  screened against an empty sub-matrix of the correlation matrix, so
  `combine_uncorrelated()` took
  [`max()`](https://rdrr.io/r/base/Extremes.html) of nothing, warned “no
  non-missing arguments to max”, and built the
  [`te()`](https://rdrr.io/pkg/mgcv/man/te.html) term regardless of its
  correlations. The error message names the variable and the requirement
  (it previously said the variable was not supplied in `use.dat`, which
  was not what was being checked).

- [`fit_mod_l()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_mod_l.md)
  is no longer exported. It is an internal helper documented as not
  called directly, and its arguments changed as recently as the family
  resolution fix. It remains available as `FSSgam:::fit_mod_l()`.
  [`wi()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/wi.md),
  [`extract_mod_dat()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/extract_mod_dat.md)
  and
  [`build_inclusion_mat()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/build_inclusion_mat.md)
  stay exported, being more plausible to call directly.

- Documentation corrections. Twenty-six spelling errors across the
  reference pages, of which three affected meaning: the correlation
  values are the square *root* of an R squared, not the “square-route”;
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)’s
  description of what it sums referred to “the ?i values” where a
  character had been lost from “the wi values”; and
  [`check_non_linear_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_non_linear_correlations.md)’s
  note told users to increase `cor.cutoff`, which is not an argument –
  `cov.cutoff` is. A non-breaking hyphen in the `case_study2`
  documentation is replaced with an ASCII one.

- The sixteen argument descriptions shared between
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  and
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)
  were duplicated verbatim and had already drifted; they are now
  inherited, so they cannot drift again.
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)’s
  `null.terms` entry gains the sentence about fitting a correlation
  structure that only
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  carried.

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
