# FSSgam (development version)

* Bug fix: a user-supplied `cor.matrix` with exactly one cell was silently
  discarded and recomputed from `use.dat`. A model set with a single predictor
  is the only case in which a supplied matrix is 1x1, and
  `length(cor.matrix) == 1` was the test for "nothing supplied", which a 1x1
  matrix also satisfies. The supplied value governed no screen and did not
  appear in `predictor.correlations`, which is the opposite of the documented
  contract that a supplied matrix replaces the automatic estimate
  (FSSgam_package#26).

  The test is now on the default value rather than on length, in both places it
  was made -- `build_predictor_correlation_matrix()` and the `factor_correlations()`
  helper inside `resolve_factor_interactions()`, which had its own copy. It tests
  dimensionality rather than a list of accepted classes, so the S4 matrix
  `Matrix::Matrix()` returns is accepted, as it was before this argument was
  validated at all. Note that `Matrix::nearPD()` is not such a case: it returns
  a list whose `mat` component is the matrix, and only that component was ever
  accepted.

  **Behaviour change:** `cor.matrix` is now validated, and a value that is
  neither a matrix, a data.frame, nor `NA` is rejected naming its class. What
  happened to such a value before depended on its length, and neither outcome
  was right. A length-one value of any class -- a string, or a list holding a
  matrix -- satisfied the sentinel and was silently treated as nothing supplied,
  so the model set was built from a computed matrix as though the argument had
  not been given. A value of any other length fell through to the supplied
  branch and failed against the missing-predictor check, reporting predictors
  rather than the argument. Calls of the first kind used to run and now error.

* **Behaviour change:** an `NA` between two terms of a supplied `cor.matrix` is
  now reported by `generate_model_set()`, naming the pairs, where it previously
  stopped the call with `missing value where TRUE/FALSE needed` — naming neither
  the matrix, the argument, nor the pair (FSSgam_package#27).

  The report is made where each sub-matrix is formed, immediately before the
  screen takes `max(abs(.))` of it, rather than at a single earlier point. Two
  earlier placements were tried and both were wrong, in opposite directions.
  Checking only in `build_predictor_correlation_matrix()` missed every pair
  screened before it, `resolve_factor_interactions()` running first. Checking a
  union of every name that might be screened then reported pairs that are not:
  a name given in the character form of `factor.factor.interactions` is screened
  against the others in that argument and never against `pred.vars.cont`.
  Reporting at the screening site covers exactly the pairs that are screened.

  **A consequence, and a departure from what the issue asked for.** The report
  covers exactly the pairs some screen compares, so it depends on
  `max.predictors` and on the interaction arguments. The issue asked for the
  report not to depend on `max.predictors`; that would mean rejecting a matrix
  over a cell nothing reads, which is a false failure. An `NA` is reported
  whenever it is read and not otherwise.

  Note that "read" is not the same as "in a candidate together". A `.by.` term
  is one term of a candidate but splits into two for the screen, so at
  `max.predictors = 1` a continuous predictor and a factor are still compared;
  two continuous predictors are not.

* Bug fix: the correlation matrix `resolve_factor_interactions()` computes when
  no `cor.matrix` is supplied is now zero-filled, as the one
  `build_predictor_correlation_matrix()` computes already was. Without it a
  computed matrix reached the factor-factor screen with `NA` in it --
  `check_correlations()` returns `NA` for a pair involving a single-level factor
  (FSSgam_package#33) -- and the call stopped with `missing value where
  TRUE/FALSE needed` (FSSgam_package#27).

  **Behaviour change, and a poor outcome that is better than the error it
  replaces.** An `NA` becoming zero means "uncorrelated", so a factor-factor pair
  that previously stopped the call is now admitted. Where the `NA` came from a
  single-level factor the interaction column is a copy of the other factor:
  `fa.I.const` is `paste(fa, "a")`, so `fa` and `fa.I.const` are the same model
  under two names, and both enter the AICc table and the variable importance
  scores. That is not a good model set, and the real defect is that a
  single-level factor is accepted as a predictor at all, which is
  FSSgam_package#33.

  Measured over a grid of 432 scenarios -- six factor sets, three continuous
  sets, both forms of `factor.factor.interactions`, `max.predictors` 1 to 3,
  both settings of `non.linear.correlations`, two forms of
  `factor.smooth.interactions` -- comparing this version against the previous
  one:

  | outcome | cells |
  |---|---|
  | identical model set | 104 |
  | identical error | 208 |
  | error before, model set now | 40 |
  | error before, a different error now | 80 |
  | model set before, a different model set now | 0 |
  | model set before, an error now | 0 |

  No call that produced a model set before produces a different one, or fails,
  now.

  Two things bound the change. Every one of the 40 cells that now build has
  `pred.vars.cont = NA`: name a continuous predictor and the call still stops,
  in `check_correlations()`, on the same single-level factor
  (FSSgam_package#33). And all 80 differing errors are single-level-factor
  scenarios reaching a later failure than before; the scenarios whose predictors
  contain `NA` are 144 of the 208 that error identically, rejected before any of
  this code runs.

  Two earlier counts of this grid were wrong and are recorded so the figures are
  not taken from them. The first reported 392 identical, of which 288 cells
  passed an invalid `null.terms` and were identical only in the error string
  `validate_null_terms()` returns, so no model set was compared in them at all.
  The second reported 144 identical, which double-counted the 40: 144 is the
  number of cells that build under this version, not the number that are
  unchanged.

  `combine_uncorrelated()` and `enumerate_candidate_models()` are otherwise
  unchanged: stopping on an `NA` is the right behaviour at both, the alternative
  being to drop the combination silently, and this only makes the message name
  the pair.

* **`DESCRIPTION` now declares `Depends: R (>= 4.4.0)`, raised from
  `R (>= 3.5)`.** The old floor was not reachable: `MuMIn` and `mgcv`, both
  hard `Imports`, each declare `Depends: R (>= 4.4.0)` at the versions current
  when this was measured (MuMIn 1.48.19, mgcv 1.9-4), so the package could not
  be installed on anything lower whatever this field said. Nothing about what
  the package can do changes; the declaration now states what was already true.
  A user pinned to an older `MuMIn` and `mgcv` may find the new floor excludes
  a configuration that did work, which is the reason the change is noted here
  rather than treated as a correction (FSSgam_package#31).

  The floor is now checked rather than asserted. `tools/check-r-floor.R`
  compares the declared floor against the `Depends` of every hard dependency,
  transitively, and fails if the declared one is lower; it runs in CI as a new
  `r-floor-consistency` job. The closure matters: `Matrix` declares
  `R (>= 4.4)` and is not a direct dependency of this package, arriving through
  `mgcv` and `MuMIn`, so a check of the direct `Imports` alone would miss it. The
  script also stops rather than passing wherever it cannot count a constraint:
  a dependency that is not installed, an R constraint written in a form it does
  not read, and a `DESCRIPTION` declaring no floor at all. `.github/workflows/R-CMD-check.yaml` also gains an
  `ubuntu-latest`, `R 4.4.0` matrix entry, so the package is built, installed
  and tested on an R at the declared floor for the first time.

  The two check different things, and only the first detects what was wrong
  here. Installing a package rejects a floor set above the running R and
  accepts any floor below it, so no matrix of R versions can see a floor that
  is too low: reverting `DESCRIPTION` to `R (>= 3.5)` leaves every matrix job
  green.

  Consequently the test suite uses `deparse1()` directly and the local
  `deparse_one()` in `tests/testthat/helper-fixtures.R`, which existed only
  because `deparse1()` arrived in R 4.0.0, is removed.

* `fit_model_set()`'s `r2.vals.unique` column is now documented. It was returned
  in `mod.data.out` and described nowhere. The documentation states that the
  column is the model R2 minus the null model's R2, so it is the variance
  explained beyond whatever was supplied in `null.terms`; that it is a per-model
  quantity rather than a variance partition among terms, and so is a per-predictor
  contribution only at `max.predictors = 1`; that it is on whatever scale
  `r2.type` produced, and so can be negative; and that values are comparable only
  within a model set sharing the same `null.terms` and `r2.type`
  (FSSgam_package#24).

* `check_non_linear_correlations()` now guards the intercept-only `multinom()`
  fit against `try()` failure, as `check_correlations()` already did. Where only
  the fitted model was tested, a failed null was passed to `round()` as the
  character vector `try()` returns, raising "non-numeric argument to mathematical
  function", which propagated out of the function rather than leaving the cell
  `NA`, which is how a failed fit is already handled. The guard is defensive: both
  models are fitted on the same rows, so no case reaching it is known, and no test
  accompanies it (FSSgam_package#19). This does not make the function degrade to
  `NA` in every failure. The `lm()` branch of the same function, and the
  equivalent one in `check_correlations()`, are not guarded at all, and a
  single-level factor predictor still aborts both: FSSgam_package#33.

* The two `github.io` URLs in `DESCRIPTION`, and the one of them also written in
  `README.md`, are rewritten with a trailing slash. Without it both redirect with
  a 301 and `R CMD check --as-cran` reports them under checking CRAN incoming
  feasibility, which `devtools::check()` does not run (FSSgam_package#20).

* `gamm4` moves from `Imports` to `Suggests`. **A fresh install of this package
  no longer installs `gamm4` or `lme4`. If you fit through
  `MuMIn::uGamm(lme4 = TRUE)`, as the companion repository's case study 2 does,
  install `gamm4` yourself.** Nothing in this package calls it: the only uses are
  `inherits(x, "gamm4")` class tests, which need no package, and an
  `@importFrom gamm4 gamm4` directive for a function never called. That directive
  is what made every parallel worker load `gamm4` and `lme4`, through the
  `FSSgam` namespace, for every model set regardless of the family fitted.

* `fit_model_set(parallel = TRUE)` no longer names `gamm4` in the packages it
  attaches on the `doSNOW` workers. Together with the change above this is the
  only direction recorded on FSSgam_package#14 that measured as reducing the
  dispatch stall, in which the call hangs indefinitely with no output and no
  error. Running `fit_model_set(parallel = TRUE)` on `case_study1` from installed
  copies of each version, alternating in blocks of ten so that both saw the same
  machine load, one cluster per fresh R process, gave 44 stalls in 140 runs
  before and 25 in 140 after (chi-squared p = 0.013, this host, 2026-09-01).

  The narrower comparison agrees. A `foreach()` whose body is `i * 2`, so that no
  package code runs at all, was run 50 times with `gamm4` in the packages the
  workers attach and 50 times without, alternating run by run: 34 stalls against
  19 (chi-squared p = 0.005, same host and date).

  Three qualifications. That is a reduction and not a removal, and
  FSSgam_package#14 stays open. Both figures are one session's measurement of a
  rate on one host, and an earlier run of the narrower comparison with its two
  arms taken one after the other rather than alternated gave a result an
  independent replication could not reproduce; only alternating designs are
  quoted here. And neither extends to a `uGamm(lme4 = TRUE)` test.fit:
  unserialising such an object loads `gamm4` and `lme4` on the worker by itself,
  before the packages are attached, so for those users the worker side is
  substantially unchanged.

  `check_correlations(parallel = TRUE)` benefits as well, without any change of
  its own: its workers load the `FSSgam` namespace and so used to load `gamm4`
  through it.

* Bug fix: `check_correlations()` estimated the correlation between two factor
  predictors from deviances computed on different sets of rows whenever either
  factor contained missing values. The fitted model was fitted on the pair's
  complete cases and the intercept-only model on the whole column, so the ratio
  the estimate is built from was not comparable and the reported correlation was
  wrong by an amount that depended on how much of the column was missing. Both
  models are now fitted on the pair's complete cases. Values are unchanged for
  data with no missing predictor values, which is every call made by
  `generate_model_set()`, since it rejects predictors containing `NA`
  (FSSgam_package#16).

* Bug fix: `generate_model_set()` screened the wrong pairs for collinearity when
  a supplied `cor.matrix` carried the same names in a different order in its rows
  and its columns. The correlation sub-matrix for a candidate was built from two
  positional indices, one in rowname order and one in colname order, so where the
  two orders differed the diagonal of that sub-matrix held cross correlations and
  the triangles held the ones from the original diagonal, and the candidate was
  dropped whatever its predictors' correlation. Both sites now index by name in
  one order. Measured on `master`, with a supplied matrix over `depth`,
  `complexity` and `ZONE`, `pred.vars.fact = "ZONE"` and everything else at its
  default: the matrix with its columns permuted returned 5 candidates where the
  same matrix unpermuted returned 9 at `max.predictors = 2` and 13 at
  `max.predictors = 3`. On this branch the two orders agree, at 9 and 13. The
  same change means a factor named twice in `pred.vars.fact` no longer yields an
  interaction of that factor with itself (FSSgam_package#15).

* `generate_model_set()` no longer requires a user-supplied `cor.matrix` to name
  the hard coded factor interaction columns that `factor.factor.interactions`
  causes to be created. Which of those columns exist depends on the supplied
  matrix itself, through the collinearity screen that decides which factor
  combinations survive, so a user could not know in advance which names to
  supply; the call stopped partway through model set construction demanding
  them. Rows and columns for any interaction column the matrix does not carry
  are now computed from `use.dat` and appended, and every value the user
  supplied is kept as supplied. A predictor the user named and did not supply is
  still an error, so a misspelled name is still reported rather than quietly
  computed (FSSgam_package#15).

# FSSgam 1.1.0

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
* Substantially expanded the `testthat` suite, from 105 passing expectations at
  1.0.0 to 568 at this release, and from 71.33% to 94.4% line coverage
  (FSSgam_package#5 covered the first 447 of those; the rest accompany the fixes
  below). New coverage includes: every non-default argument of `generate_model_set()`,
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

* Bug fix: a name in `cyclic.vars` was matched against the assembled model
  terms with an unanchored `grep()`, so it captured every predictor whose name
  contained it, and the name was used as a regular expression. Declaring
  `depth` cyclic also fitted `depthx` with `bs = "cc"`, and a predictor name
  containing a full stop matched any character in that position. The affected
  models were fitted with the wrong smoothing basis, without error, under a
  candidate name that looked correct. The basis is now chosen from the
  variable names as each term is built.
* Bug fix: a `te()` term over three or more predictors was given only two
  `bs` values, because the code that assigned them took the first two
  variables. mgcv warns "bs wrong length and ignored" and substitutes its own
  default, so both `bs.arg` and any `cyclic.vars` were silently discarded for
  that term. Each marginal now carries its own basis. Reachable through
  `smooth.smooth.interactions` given as a character vector with
  `max.predictors >= 3`.

* Behaviour change: candidate model names are now sorted in byte order, via
  `sort(method = "radix")`, and no longer depend on the session's collation
  locale (FSSgam_package#8). A model named `complexity+ZONE` in an
  `en_US.UTF-8` session was named `ZONE+complexity` in a C-locale one, so a
  saved analysis was not reproducible across machines with different locales.
  Names are now the C-locale form everywhere. The fitted results are
  identical, but `modname` values and the row order of `mod.data.out` change
  for anyone working in a non-C locale; a saved model table matched on
  `modname` needs regenerating.

* Behaviour change: a user-supplied `cor.matrix` now governs every stage that
  screens on correlation, not only the exclusion of assembled candidate models
  (FSSgam_package#13). There are three such stages: which factor-factor
  interaction columns are built, which `te()` smooth-smooth interaction terms
  are built, and which candidates survive. The first two previously recomputed
  correlations from `use.dat` and ignored the supplied matrix, so a user who
  supplied a matrix saying two predictors were uncorrelated still found their
  `te()` term absent. Users who do not supply a `cor.matrix` see no change: the
  continuous-continuous block does not depend on which other predictors are
  present, verified against both `check_correlations()` and
  `check_non_linear_correlations()`. A supplied matrix must now include any
  hard coded factor interaction it causes to be created, which
  `?generate_model_set` already required; missing predictors are reported by
  name.

* Bug fix: a user-supplied `cor.matrix` did not replace the automatic estimate,
  it only overrode it. `check_correlations()` (or
  `check_non_linear_correlations()`) was called over every predictor and its
  result discarded whenever a matrix was supplied, so supplying one neither
  avoided the `multinom()` and `gam()` fits nor allowed a predictor of a class
  those functions reject -- a `Date` column, say -- to be used at all. The
  estimate is now computed only when no matrix is supplied (FSSgam_package#13).
* Bug fix: with `save.model.fits = FALSE`, a run in which no candidate fitted
  returned a model table of `NA` with `delta.AICc` and `wi.AICc` filled with
  `NaN`, rather than the "None of your models fitted successfully" error raised
  on the `save.model.fits = TRUE` path. Which of the two happened depended only
  on a memory setting. Both paths now raise the error.
* `fit_model_set()` and `full_subsets_gam()` reject a `progress` value that is
  not a single `TRUE` or `FALSE`, and `full_subsets_gam()` now validates
  `VI.mods` at entry as it already did `r2.type`, so an unrecognised value is
  reported before the candidate set is built rather than after.
* `smooth.smooth.interactions` naming a column of `use.dat` that is not among
  the predictors is now reported as such. It was previously screened against an
  empty sub-matrix of the correlation matrix, so `combine_uncorrelated()` took
  `max()` of nothing, warned "no non-missing arguments to max", and built the
  `te()` term regardless of its correlations. The error message names the
  variable and the requirement (it previously said the variable was not supplied
  in `use.dat`, which was not what was being checked).

* `fit_mod_l()` is no longer exported. It is an internal helper documented as
  not called directly, and its arguments changed as recently as the family
  resolution fix. It remains available as `FSSgam:::fit_mod_l()`. `wi()`,
  `extract_mod_dat()` and `build_inclusion_mat()` stay exported, being more
  plausible to call directly.

* Documentation corrections. Twenty-six spelling errors across the reference
  pages, of which three affected meaning: the correlation values are the square
  *root* of an R squared, not the "square-route"; `full_subsets_gam()`'s
  description of what it sums referred to "the ?i values" where a character had
  been lost from "the wi values"; and `check_non_linear_correlations()`'s note
  told users to increase `cor.cutoff`, which is not an argument -- `cov.cutoff`
  is. A non-breaking hyphen in the `case_study2` documentation is replaced with
  an ASCII one.
* The sixteen argument descriptions shared between `generate_model_set()` and
  `full_subsets_gam()` were duplicated verbatim and had already drifted; they
  are now inherited, so they cannot drift again. `full_subsets_gam()`'s
  `null.terms` entry gains the sentence about fitting a correlation structure
  that only `generate_model_set()` carried.

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
