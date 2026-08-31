## Session: pre-CRAN refactor
Date: 2026-08-31
Model: Claude Opus 5 (claude-opus-5)

### Prompts and Responses

**User:** I want to refactor the package before getting it ready to submit to
cran. We have built a comprehensive test suite, during which several issues
were identified. Should we tackle the issues first? Or first build a plan for
the refactor?

**Claude:** Recommended planning first, on the grounds that the open issues are
not independent of each other: FSSgam_package#8 (locale-dependent model names)
and #13 (a supplied `cor.matrix` not governing `te()` construction) both change
outputs the numerical snapshots pin, and fixing them in isolation would mean
updating those snapshots twice.

---

**User:** can you review the full code base and let me know if there are
improvements you would suggest as part of any pre-cran re-factor. I want to
sort out as much as before before we submit to cran

**Claude:** Read all 1,947 lines in `R/` at commit `a8235e2` and reported
twenty-two items in six groups. Four defects that silently change or abort
results (A1 unanchored `grep()` matching of `cyclic.vars`; A2 `fit_model_set()`
failing on a single-predictor set; A3 `for(i in 2:n)` counting backwards at
n = 1; A4 `.errorhandling` commented out on the unsaved parallel path), four
arguments accepting invalid values without complaint (`r2.type`, `VI.mods`, the
three deprecated `full_subsets_gam()` sentinels, and disagreeing `max.models`
defaults), duplicated interaction-resolution code, unsafe idioms, and
documentation errors.

Each of the four defects and each validation gap was reproduced on this host
before being reported.

Recorded a correction in the process: the documentation items were first framed
as CRAN concerns. `devtools::check()` had been run with
`_R_CHECK_CRAN_INCOMING_ = FALSE`, which is where such checks surface, so
`R CMD check --as-cran` was re-run on a built tarball with the local incoming
checks enabled. It also returned Status: OK, so nothing in that group blocks
submission; it is fixed because it is published on the pkgdown reference site.

---

**User:** The refactor should implement fixes for all the items identified and
also restructure the interactio-resolution code.

**Claude:** Wrote the two planning documents required by `CLAUDE.md` Section 13:
`planning/pre-cran-refactor-human.md` (the collaborator-facing plan) and
`planning/pre-cran-refactor-claude.md` (the implementation specification).
Added `^planning$` to `.Rbuildignore`.

Two invariants underpinning the design were verified empirically before being
written down, rather than asserted from reading the code:

- `cor.matrix` contains no `.te.` name (row names `depth, SCORE2, complexity`
  against a candidate set containing `depth.te.SCORE2`), so
  `build_predictor_correlation_matrix()` can be moved ahead of
  `resolve_smooth_smooth_interactions()`.
- the continuous-continuous block of the correlation matrix is identical
  whether computed on the continuous subset or on the full predictor set, for
  both `non.linear.correlations` settings, maximum absolute difference 0. So
  subsetting the resolved matrix is behaviour-preserving for anyone not
  supplying their own.

The order of work puts the restructure at phase 3 and the behaviour changes at
phase 4, so that the restructure is verifiable as a pure no-op by golden master
and the behaviour changes then land on code that has already been simplified.

---

**User:** cor.matrix should also govern factor-factor column construction and
smooth.smooth.interactions should retain the current behaviour.

**Claude:** Recorded both decisions in the planning documents.

The first extends FSSgam_package#13 to the third stage that uses correlations:
`build_factor_interaction_columns()` will screen against a supplied
`cor.matrix` when one is given and compute its own only when none is. It must
use the *supplied* matrix rather than the resolved one, because the resolved
matrix is built over `all.predictors`, which includes the interaction columns
this stage is deciding whether to create.

The second rejects the specification's own recommendation to fix both
`smooth.smooth.interactions` branches at pairwise. `TRUE` builds bivariate
`te()` terms and a character vector builds terms of size 2 up to
`max.predictors`; the latter is a working capability nobody has asked to lose.
The documentation is corrected to describe both forms instead. This moves the
item out of phase 4, where it was a behaviour change, into phase 7, where it is
a documentation correction.

---

### Phase 1 — defects and validation outside the interaction code

**A2 — `fit_model_set()` fails on a single-predictor model set.**

Reproduced first:

```r
tf <- gam(log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr"), data = d)
ms <- generate_model_set(use.dat = d, test.fit = tf, pred.vars.cont = "depth",
                         max.predictors = 1, k = 3)
fit_model_set(ms)
#> Error in colSums(mod.data.out[, included.vars]) :
#>   'x' must be an array of at least two dimensions
```

`compute_variable_importance()` indexed `mod.data.out` by a character vector of
column names without `drop = FALSE`, so a single predictor collapsed the result
to a vector. Fixed at all three `colSums()` sites.

A second defect was found in the same function while fixing it and corrected in
the same commit. `min.mods` is 0 when some predictor appears in no surviving
model, and the `VI.mods = 'min.n'` branch indexed the sorted weights with
`[1:min.mods]`. `1:0` is `c(1, 0)`, so it selected the single largest weight
rather than none. Replaced with `seq_len(min.mods)`, which is empty at 0.

Both are pinned by new tests in `test-fit_model_set_options.R`. Confirmed they
fail against the unmodified function — the first with the error above, the
second reporting weights of `1.0 0.0 NA` where `0 0 0` is expected — and pass
after.

---

**A4 — the unsaved parallel fitting path lost its error handling.**

`.errorhandling='pass'` was commented out on the `foreach()` in
`fit_and_summarise_unsaved_models()`, while the saved path
(`fit_and_summarise_saved_models()`) has always had it. One unfittable
candidate therefore aborted the whole run instead of being recorded in
`failed.models`.

Re-enabling it is not a one-character change. With `.errorhandling='pass'` a
failing element comes back as the condition object rather than as the named
vector `unlist(extract_mod_dat(...))` returns, and `do.call("rbind", .)` over
the mixed list produces a malformed table. Added two unexported helpers,
`na_mod_dat_row()` and `normalise_mod_dat_rows()`, which replace any element
that is not one of those vectors with `extract_mod_dat()`'s all-NA row. The
candidate then reaches `failed.models` through the `is.na(AICc)` test the
caller already applies, which is how a fit returning a `try-error` is already
handled.

Deviated from the specification on one detail. It proposed testing elements
with `is.numeric(x) && length(x) == 6`; `is.atomic()` is used instead, because
the all-NA row is *logical*, not numeric, so `is.numeric()` would reject rows
that are already correct and replace them with themselves. A condition object
is a list, so `is.atomic()` still separates the two.

The replacement logic was factored into a named helper rather than left inline
specifically so that a sequential test can reach it: the sequential path cannot
produce a malformed element, so inline code would be exercised only under
`parallel = TRUE`, which `covr` never runs (Phase 13).

`normalise_mod_dat_rows()` is tested directly, and an end-to-end assertion was
added to `test-fit_model_set_parallel.R` behind `skip_unless_parallel_opt_in()`.

**Verification note.** That end-to-end assertion has not yet been observed to
pass. Two attempts against an installed copy stalled: the doSNOW workers
started and then sat idle, and the run was killed at its 400 s timeout. Another
R process belonging to an unrelated session was saturating a core throughout,
which is the loaded-machine condition recorded in the Phase 13 caution. To be
re-run in isolation before the pull request. The sequential suite is green at
467 passing expectations (baseline 438).

---

**B1 — `r2.type` accepted any value.**

`fit_model_set(model.set, r2.type = "rsq")` returned an all-NA `r2.vals`
column with no message: `extract_mod_dat()` matches `r2.type.` against three
literals in sequence and leaves the value at `NA` when none matches. Added
`match.arg(r2.type, c("r2.lm.est", "r2", "dev"))` to `fit_model_set()` and to
`full_subsets_gam()`.

Validated in the wrapper as well as in `fit_model_set()`, although the wrapper
forwards the argument, because `generate_model_set()` runs first: without the
earlier check the whole candidate set is built before the error is raised.

Deliberately *not* added to `extract_mod_dat()`. That function is exported and
documents `r2.type.` as being passed straight through, so validating there
would make an exported function stricter than it is today. The same reasoning
was later applied to `check_correlations()`.

`"r2"` is an exact match against the second element, so `match.arg()`'s partial
matching does not make it ambiguous with `"r2.lm.est"`; this is asserted rather
than assumed, by a test that runs all three documented values through.

---

**B2 — `VI.mods` accepted any value.**

`fit_model_set(model.set, VI.mods = "alll")` failed with `object
'aic.var.weights' not found`: `compute_variable_importance()` had two `if`
blocks and no `else`, so a value matching neither assigned nothing and the
failure surfaced further down as a missing object rather than as a bad
argument.

Added `match.arg(VI.mods, c("min.n", "all"))` to `fit_model_set()` and changed
the second `if` to `else`. The helper is unexported and only reached after that
`match.arg()`, so the two branches are exhaustive.

The corresponding assertion for `full_subsets_gam()` was written and then moved
to the FSSgam_package#7 commit below, because the wrapper has no `VI.mods`
argument at this point in the sequence — the test failed with `unused argument
(VI.mods = "alll")`, which is the defect #7 reports rather than the one under
test here.

---

**B3 — the deprecated `full_subsets_gam()` arguments failed on valid input.**

`factor.interactions`, `smooth.interactions` and `size` were declared with a
`"previous.arg"` sentinel default and tested with `!=` or `is.na()`. Both
comparisons fail on input the arguments are documented to accept:

```
full_subsets_gam(..., smooth.interactions = c("ZONE", "ZONE2"))
#> Error in if (is.na(smooth.interactions) == TRUE) { :
#>   the condition has length > 1
full_subsets_gam(..., factor.interactions = NA)
#> Error in if (factor.interactions != "previous.arg") { :
#>   missing value where TRUE/FALSE needed
```

Replaced the sentinel with `missing()`, which is correct for every type,
length and `NA`. The formals lose their defaults as a result.

This removes the special branch that handled `smooth.interactions = NA`. Under
`missing()`, `NA` is simply assigned, which is its documented meaning for
`factor.smooth.interactions` ("If specified as NA no factor-continuous
predictor interactions will be included"), and it now warns like any other use
of a deprecated argument, which it did not do reliably before.

Four tests added. Two pin the two errors above; both were confirmed to fail
against the unmodified function with exactly those messages. One pins that `NA`
still warns. The fourth pins that `missing()` keeps reporting correctly through
`full.subsets.gam()`, which forwards everything through `...` rather than
naming the arguments — it asserts that the alias emits its own `.Deprecated()`
warning but none of the three legacy-argument warnings, and that a legacy
argument passed through the alias is still detected.

A stale comment on a neighbouring test was corrected in the same commit: it
described the `is.na()` branch that no longer exists.

---
