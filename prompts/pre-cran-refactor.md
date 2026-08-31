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

**B4 — `max.models` had different defaults in the wrapper and the pair.**

`fit_model_set()` defaulted to 200 and `full_subsets_gam()` to 500,
undocumented in either, so a candidate set of 300 saved its model fits through
the wrapper and did not through `generate_model_set()` + `fit_model_set()`.
Set both to 200 and documented the default in both `@param` entries.

200 rather than 500 because it is the value the two primary user-facing
functions already used; raising `fit_model_set()` to 500 would have increased
the memory a default call can consume, which is what the argument exists to
bound. This is a behaviour change for anyone with a set between 201 and 500
calling `full_subsets_gam()`: they now get the warning and lose
`success.models`. Recorded in `NEWS.md`.

The test asserts the two defaults are equal rather than only that each is 200,
so the two cannot drift apart again without failing.

`NEWS.md` entries for all of phase 1 were written in one commit at the end of
the phase rather than one per fix, since the file is a per-release changelog.

---

**FSSgam_package#7 — `full_subsets_gam()` did not forward `VI.mods`.**

`fit_model_set()` takes `VI.mods` but the wrapper never passed it on, so the
`'all'` setting was unreachable through `full_subsets_gam()`. Added the
argument with the same `'min.n'` default, forwarded it, and copied the
`@param` entry.

The existing test asserting the absence (`expect_false("VI.mods" %in%
names(formals(full_subsets_gam)))`, written during Phase 13 to pin the
behaviour until the issue was addressed) was replaced.

Writing the replacement surfaced something worth recording. The obvious test —
run the wrapper with each setting and compare against `fit_model_set()` called
directly — passes on the standard fixture *even if nothing is forwarded*,
because the two settings produce identical variable importance there. Measured
on that model set: `min.mods` is 4, and every model beyond the fourth for each
predictor carries a `wi.AICc` of 0.000 once `compute_model_weights()` rounds to
three decimal places, so summing the best four and summing all seven give the
same total. Predictor counts of 4, 4, 5 and 7 are not enough to separate them.

The separating case is a response with no signal, which spreads the weights so
that the models beyond the fourth carry non-zero weight. A second test uses a
`set.seed(42)` normal deviate as the response and asserts that the two settings
differ, and that `'all'` is never smaller than `'min.n'`.

---

**FSSgam_package#9 — the progress bar could not be suppressed.**

`fit_model_set()` wrote a `txtProgressBar` to stdout unconditionally, so the
only way to keep it out of a report, a script or a test reporter was to wrap
the call in `capture.output()` — which is what the suite's `fit_quietly()`
helper exists to do.

Added `progress = interactive()` to `fit_model_set()` and
`full_subsets_gam()`, forwarded into both fitting helpers. When `FALSE`, `pb`
is left `NULL` and `.options.snow` is passed an empty list, so no callback is
registered and nothing reaches stdout at all. The default means the bar
appears at the console but not in scripts, reports or `R CMD check`.

The progress callback was renamed from `progress` to `update_pb()`, since the
old name now shadows the argument.

This also removes the double `close(pb)` in the parallel branch of both
fitting helpers: the bar was closed inside the `if (parallel)` block and again
after it.

`fit_quietly()`/`full_subsets_quietly()` were left as they are. They remain
correct with the bar off, and Phase 13 records four groups of tests that
deliberately keep their direct call form.

**A note on `DESCRIPTION`.** `devtools::document()` was run once at the end of
the phase rather than after each commit, so `man/fit_model_set.Rd` and
`man/full_subsets_gam.Rd` are regenerated here for all of the phase's
documentation changes. The roxygen2 installed in this environment is 8.0.0,
which rewrote `RoxygenNote: 7.3.2` to `Config/roxygen2/version: 8.0.0`. That
change was reverted: it reflects the toolchain on this machine rather than a
decision about the package, and adopting a major roxygen2 version is a
separate matter. The generated `.Rd` content itself shows no structural
difference attributable to the version — the diff contains only the intended
argument and text changes.

---

### Phase 2 — the correlation functions

**FSSgam_package#10 — spurious "NaNs produced" warnings out of `nnet`.**

`check_correlations()` and `check_non_linear_correlations()` took the
factor-factor deviance from `summary(nnet::multinom(...))$deviance`. The
summary method computes a variance-covariance matrix in order to report
standard errors, and on a perfectly separated fit its diagonal is negative, so
`sqrt()` warns. Only `$deviance` is used, and the fit carries it directly.
Changed both calls in both functions, and both branches in
`check_correlations()`, to read it off the fit.

Verified rather than assumed, in two ways.

`identical(summary(f)$deviance, f$deviance)` is `TRUE` for the fits involved.
A golden-master comparison over 14 scenarios — continuous only, factors only,
mixed, one factor, a single factor column, a single continuous column, and a
`case_study2` mixed set, each through both `check_correlations()` and
`check_non_linear_correlations()` — gives bit-identical matrices before and
after, at `tolerance = 0`.

Recorded honestly: the warning itself could not be reproduced on this machine,
either before or after. nnet here is 7.3-20, and its optimiser does not reach
the state that warns, even on a deliberately perfectly separated fit
(`ZONE ~ paste(ZONE, ZONE2)`). This matches what Phase 13 recorded. The change
is therefore justified by removing a computation whose result is discarded,
which is where the warning came from on the platforms that did show it — not by
a reproduction here.

The added test asserts no warning is emitted on a perfectly separated pair. It
does not fail against the unfixed code on this nnet version, which is stated in
a comment on the test itself so it is not mistaken for a regression test. The
`suppress_nnet_nans()` helper is retained, and its comment corrected: it
described the `summary()` call that no longer exists.

---

**FSSgam_package#12 — the factor-factor diagonal was not exactly 1.**

`fill_factor_factor_correlations()` enumerated all ordered pairs of factors,
including each factor with itself. Regressing a factor on itself gives a
pseudo-R2 of about 0.999999 rather than exactly 1, and costs one `multinom()`
fit per factor for a value known in advance. Self-pairs are now dropped from
the grid and the factor block's diagonal assigned directly.

Only the factor block's diagonal is assigned, via
`out.cor.mat[cbind(fact.vars, fact.vars)] <- 1`, rather than
`diag(out.cor.mat) <- 1` over the whole matrix. The continuous block comes from
`cor()`, which already returns exactly 1 there, and the narrower assignment
cannot disturb it.

Dropping the self-pairs makes the grid empty when there is exactly one factor,
which `check_correlations()` does reach. `1:nrow(lm.grid)` would then count
backwards, so the two loop bounds became `seq_len()` and the fill loop
`seq_along()`. This was necessary for correctness here, not idiom tidying.

Verified with the same 14-scenario golden master. Every off-diagonal cell is
identical at `tolerance = 0`; the only difference is on the diagonal, at most
6e-07, and the diagonal is exactly 1 in every scenario afterwards.

Three existing assertions were tightened from `> 0.99` to exact equality, as
Phase 13 anticipated when it wrote them deliberately loose so they would keep
passing once the issue was fixed.

---

**C4 — the null multinomial was refitted once per pair.**

`null.fit` in `fill_factor_factor_correlations()` depends only on `var.1`, but
was refitted inside the loop for every ordered pair. Hoisted to one fit per
factor, computed before the loop and looked up by name. `foreach()` exports it
to the workers automatically, as it already does for `dat` and `lm.grid`.

**A correction to the specification.** It stated that hoisting "changes `r.est`
for factor pairs with missing data", because the null fit used `dat[, var.1]`
while the pair fit used `dat.r`, the pair's complete cases. Re-reading the
original line shows the null was *already* computed on the whole `dat` column
regardless of `var.2`, so the hoist recomputes the identical value fewer times
and changes nothing. Confirmed against the 14-scenario golden master: identical
at `tolerance = 0` in every case. The `NEWS.md` entry the specification asked
for was therefore not written, and the code comment drafted from that claim was
corrected before committing.

The `dat` / `dat.r` asymmetry itself is real and pre-existing, and was left
alone: the null deviance comes from the whole column while the pair deviance
comes from `na.omit()` of the pair, so for a factor carrying `NA`s the two are
computed on different row sets. Changing it would change reported correlations
for a case `generate_model_set()` already rejects. Recorded in a comment at the
hoisted fit, which is now the one place it is visible.

**Measurement.** `multinom()` calls for a four-factor `case_study1` set,
counted with `trace()`: 32 on `master`, 24 after FSSgam_package#12 removed the
self-pairs, 16 after this change. Elapsed time is not a useful figure here —
0.022 s against 0.020 s over five replicates — because each fit on this data is
sub-millisecond. The count is the quantity that scales.

---

### Verifying A4 under `parallel = TRUE`, and a pre-existing intermittent stall

The A4 regression test could not be observed to pass for some time, and the
investigation produced two corrections worth recording.

**First correction: A4 was not the cause of the stalls.** From single runs it
appeared that re-enabling `.errorhandling='pass'` made the unsaved parallel path
hang: removing that one argument gave a clean 1.6 s run, leaving it gave a
timeout. That inference was wrong. Adding two `message()` calls to the same
unmodified function also produced a clean 1.6 s run, which no theory of
`.errorhandling` explains. The path stalls *intermittently*, so one observation
either side established nothing.

Measured properly, 12 repetitions of the same plain script against each build,
25 s timeout (a successful run takes 1.6 s):

| build | completed | stalled |
|---|---|---|
| this branch, A4 applied | 8/12 | 4/12 |
| `master` worktree | 6/12 | 6/12 |

`master` stalls at least as often, so A4 did not introduce this. A trivial
`foreach()` completes in 0.29 s, and `.errorhandling='pass'` with `.packages`
and `.options.snow = list()` completes in isolation, so it is not the cluster
machinery or the arguments themselves.

This extends the Phase 13 caution rather than repeating it. Phase 13 recorded
the stall for the *test file* and noted that the same contents "run from a plain
script completed six times out of six". That no longer holds: the plain script
stalls at roughly the same rate as the test file.

**Second correction: a defective probe, not a defective package.** The first
successful A4 run reported nine rows in a table for a model set of eight. The
probe injected the broken formula with
`ms$mod.formula[["ZONE+depth"]] <- ...`, and `[[<-` on a name that does not
exist *appends* rather than replaces. In a plain UTF-8 session the candidate is
named `depth+ZONE`; `ZONE+depth` is its C-collation spelling, which is what
testthat produces. This is FSSgam_package#8 seen from the other side, and it is
exactly what `break_one_candidate()`'s `stopifnot(modname %in% ...)` exists to
catch, so the committed test was never at risk.

**A4 verified.** With the probe corrected, under `parallel = TRUE` and
`save.model.fits = FALSE` with one unfittable candidate:

```
failed.models names : depth+ZONE
rows in table       : 8 of 8
broken row AICc isNA: TRUE
others all finite   : TRUE
matches sequential  : TRUE
```

Before the fix this configuration aborted the whole run.

---

**The stall is gamm4 loading on the workers, not this package.**

Instrumenting a stall showed the last marker reached is "entering foreach":
`makeCluster()` and `registerDoSNOW()` both complete, so it is not cluster
startup. A `foreach()` with a trivial body (`i * 2`), calling no FSSgam code at
all, reproduces it purely from `.packages`:

| `.packages` | completed | loads gamm4 |
|---|---|---|
| `mgcv` | 10/10 | no |
| `mgcv, MuMIn` | 10/10 | no |
| `gamm4` | 7/10 | yes |
| `mgcv, gamm4` | 3/10 | yes |
| `mgcv, FSSgam` | 6/10 | yes, via Imports |
| `mgcv, gamm4, MuMIn, FSSgam` | 8/10 | yes |

Every set that loads gamm4 stalls and no set without it does: 20/20 against
24/40. No ordering is read into the individual gamm4 rates; at n = 10 those
intervals overlap heavily, and the only supported statement is gamm4 present
against gamm4 absent.

Not established: whether this is specific to WSL2 or to this lme4/Matrix build,
and the mechanism. gamm4 loads lme4, Matrix and nlme, so a namespace or S4
method-table race on the workers is plausible but was not demonstrated.

This is outside the refactor's scope and predates all of it, so no code change
was made. A correction to `CLAUDE.md` is drafted, since three statements there
are now known to be wrong: that the stall is in cluster startup, that a plain
script completed six times out of six, and that it is contention from another
process.

---

### Phase 3 — restructuring interaction resolution

**The golden master.** Thirty-one `generate_model_set()` scenarios, captured
before any change and compared with `identical()` afterwards. Each records
`n.mods`, the candidate names, every deparsed formula, `predictor.correlations`,
the names and dimensions of `used.data`, `included.vars`, and every warning or
error raised. Errors and warnings are captured as values so that a scenario
which currently fails keeps failing in the same way.

Building it changed what I understood about the code in three ways, each found
by measuring rather than reading.

**The three-way `te()` is real, but the obvious scenario does not show it.**
With the default `cov.cutoff` of 0.28, `smooth.smooth.interactions = TRUE` and
the character form give byte-identical output at `max.predictors = 3`. The
reason is that `depth` and `complexity` correlate at 0.3336, so that pair is
screened out and the three-way combination containing it goes with it. At
`cov.cutoff = 0.4` all three pairs are admitted and the character form produces
`depth.te.SCORE2.te.complexity` while `TRUE` does not. Without that scenario I
would have "retained" a capability with nothing demonstrating it was still
there.

**A third divergence between the two factor branches, not in the review.** The
character branch is wrapped in `if(length(which(factor.correlations <
cov.cutoff)) > 1)`, which counts cells of the correlation matrix rather than
pairs. With two factors whose two off-diagonal cells straddle the cutoff, the
count is 1 and the entire block is skipped in silence. This defeated the first
attempt at a scenario separating the upper/lower triangle screening: both
branches produced nothing, for two unrelated reasons.

**A scenario that does separate them.** `check_correlations()` fits `multinom()`
separately in each direction, so its factor block is asymmetric — by at most
8.5e-4 over these columns. `ZONE3`/`ZONE5` measure 0.5425 one way and 0.5434 the
other, so `cov.cutoff = 0.543` falls between. With a third uncorrelated factor
present to clear the cell-counting guard, the logical branch drops
`ZONE3.I.ZONE5` and the character branch keeps it: 15 candidates against 16.

**The restructure.** `combination_sizes()` and `combine_uncorrelated()` replace
four copies of "enumerate every combination, drop the ones that are too
correlated". `resolve_factor_interactions()` is 140 lines, down from 161;
`resolve_smooth_smooth_interactions()` is 52, down from 73, with the `tt` data
frame that both copies built and never read now gone.

Two details preserved deliberately rather than tidied:

- `combination_sizes()` reproduces `2:n` counting backwards at `n = 1`, so
  size-1 "combinations" are still enumerated. That is the A3 defect, and it is
  corrected in its own commit so the restructure and the fix are reviewable
  apart.
- `combine_uncorrelated()` keeps the `if (max(...) > cov.cutoff) out <- NA`
  shape rather than reducing to a logical predicate. `max()` of an all-`NA`
  sub-matrix returns `NA`, and `if (NA)` raises an error; a `vapply()` over a
  logical would silently drop the combination instead, which is a different
  behaviour.

`screen.both` holds each caller at its current triangle-screening behaviour.
All 31 scenarios identical, suite green at 508.

---
