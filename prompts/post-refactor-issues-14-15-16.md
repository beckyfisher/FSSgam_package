## Session: post-refactor issues 14, 15 and 16
Date: 2026-09-01
Model: Claude Opus 5 (claude-opus-5)

### Prompts and Responses

**User:** We just finished a major refactor. This raised three new issues. Work
through each issue sequentially, in a loop: address the problem and post a PR,
then have an independent session draft a review comment, address what it raises,
draft a summary comment, and have a fresh session review that. Repeat for up to
10 cycles or until the review finds nothing substantive. Genuine issues outside
the scope of the one being worked on are raised on GitHub as new issues.

**Claude:** Identified the three issues as FSSgam_package#14 (`parallel = TRUE`
stalls when `gamm4` is loaded on the doSNOW workers), #15 (a supplied
`cor.matrix` must name the hard coded factor interaction columns), and #16
(`check_correlations()` compares factor-factor deviances computed on different
row sets). Each of #15 and #16 was raised with more than one candidate
resolution and no decision recorded, so the three decisions were put to the
user before any code was written:

- #16: fit the intercept-only model on the pair's complete cases, keeping the
  hoisted per-factor value for any pair that drops no extra rows.
- #15: accept a matrix covering only the predictors the user named, and compute
  rows and columns for the interaction columns the function itself creates.
- #14: measure the three untried directions recorded on the issue and adopt one
  only if it measurably lowers the stall rate.

---

### Issue #14 -- parallel = TRUE stalls when gamm4 is loaded on the workers

**Claude:** The issue records three untried directions. All three were measured
before anything was written, one cluster per fresh R process with a `foreach()`
body of `i * 2`, so no package code executes:

| variant | stalled |
|---|---|
| baseline, `.packages = c("mgcv", "gamm4", "MuMIn", "FSSgam")` | 17 of 70 |
| `makePSOCKcluster(setup_strategy = "sequential")` | 9 of 70 |
| `clusterCall()` preload, no `.packages` | 10 of 70 |
| both together | 7 of 50 |

None of those three separates from the baseline at that number of replicates, and
`setup_strategy` also postdates R 4.0 while `DESCRIPTION` declares R (>= 3.5).
The fourth direction, narrowing `.packages`, is the one the issue's own table
supports: `mgcv` alone completed 39 of 40 while `gamm4` present completed 29 of
50.

Narrowing it is not sufficient on its own. `gamm4` was under `Imports`, and
loading the `FSSgam` namespace on a worker loads its imports, so a worker given
`.packages = c("mgcv", "MuMIn", "FSSgam")` still loaded `gamm4`. Nothing in `R/`
calls `gamm4`; the only uses are `inherits(x, "gamm4")` class tests, and the
`@importFrom gamm4 gamm4` directive imported a function the package never calls.
It was moved to `Suggests` and the directive removed. `MuMIn` only suggests
`gamm4` as well, so nothing else pulls it in. Confirmed on the installed copy:
`library(FSSgam)` leaves `gamm4` and `lme4` unloaded, and a worker running a real
candidate refit reports neither in `loadedNamespaces()`.

`worker_packages()` (`R/utils.R`) then adds `gamm4` back for a `test.fit` that is
one. `uGamm(lme4 = TRUE)` gives class `c("gamm4", "gamm")`; `gam()` and
`uGamm(lme4 = FALSE)` do not, and neither reaches `gamm4` when refitted.

Two measurements support the change, both one cluster per fresh R process:

- with a `foreach()` body executing no package code, 23 stalls in 60 runs with
  `gamm4` on the workers and 0 in 60 without;
- through `fit_model_set(parallel = TRUE)` on `case_study1`, eight candidates,
  running the two versions alternately in blocks of ten from separate installed
  libraries so that both arms see the same machine load, 44 stalls in 140 before
  and 25 in 140 after, chi-squared p = 0.013.

The second figure is the one that matters and it is a reduction, not a removal.
The first 30 runs of each arm were indistinguishable and the difference accrued
over the remaining 110, which is why the measurement was extended rather than
stopped at the first result. Something other than `gamm4` accounts for the
residual rate; a probe of `loadedNamespaces()` on a worker during a real refit
rules out `gamm4` and `lme4` as that cause, and nothing further was established.

Three tests were added: that `worker_packages()` omits `gamm4` for a `gam()` and
for a `uGamm(lme4 = FALSE)` test.fit and includes it for a `gamm4` one, and that
`gamm4` is absent from `Imports`. The last is the invariant that makes the other
two effective, and it fails if anyone moves `gamm4` back.

---

#### Issue #14, review round 1 (independent session)

**Claude (reviewer):** One blocking finding, two substantive. Eight test blocks
build a `test.fit` with `MuMIn::uGamm(lme4 = TRUE)`, which fails outright when
`gamm4` is not installed; while `gamm4` was under `Imports` that could not
happen, and CRAN requires a package to check successfully with its suggested
packages absent. `R CMD check` on the branch tarball with `gamm4` hidden from
the library path reports 2 ERRORs, tests failing 8. The headline measurement did
not reproduce: run on the reviewing session's copy, the trivial-body reprex gave
5 stalls in 40 with `gamm4` against 11 in 40 without, in the opposite direction,
at a much lower overall rate than either the pull request or the issue records.
And the `gamm4` branch of `worker_packages()` is redundant, because
`MuMIn::uGamm()` calls `require("gamm4")` itself, so a worker given
`.packages = c("mgcv", "MuMIn", "FSSgam")` refits a `uGamm(lme4 = TRUE)` model
set with no failures and ends with `gamm4` and `lme4` loaded.

**Claude:** All three addressed.

`skip_if_not_installed("gamm4")` was added to `fixture_coral_ugamm()` and to the
one direct `uGamm(lme4 = TRUE)` call. Verified against a library view built by
symlinking every package in the site library except `gamm4`: 547 passing, 0
failing, 0 errors, 14 skipped, against 571 and 7 with `gamm4` present.

The `gamm4` branch was removed rather than defended. It buys no correctness, and
on this change's own account of the mechanism it restores the condition the
change identifies as harmful, for exactly the users who fit `gamm4` models.
`worker_packages()` now takes no argument and returns `c("mgcv", "MuMIn",
"FSSgam")`, with the reason for each recorded on the face of it: `update()`
resolves the stored call's function symbol, `gam` or `uGamm`, by lexical lookup
ending at the worker's search path, so neither `mgcv` nor `MuMIn` is reached
without being attached.

The measurement claim was withdrawn rather than defended. The trivial-body arms
were run one after the other rather than interleaved, so they are confounded with
time and with machine load, and a failure to reproduce them is the expected
result. The interleaved end-to-end comparison, 44 stalls in 140 against 25 in
140 with the two versions alternating in blocks of ten, is the measurement that
survives, and `NEWS.md` now carries that one alone, attributed to this host and
date. An attempt to re-run the trivial-body comparison interleaved was abandoned:
an unrelated job had taken the host's load average above 20 on 22 cores and every
run of both arms timed out, 10 of 10 each.

The `dsm` finding was raised as FSSgam_package#30.

---

#### Issue #14, review round 2 (fresh session)

**Claude (reviewer):** No blocking findings; every non-timing claim reproduces,
including `R CMD check` with `gamm4` hidden returning Status OK. Three
substantive points, all about the recorded account rather than the code. The
stated reason for keeping `MuMIn` in `worker_packages()` is wrong: `MuMIn`
registers `update.gamm`, whose body calls `update.default` from inside the
`MuMIn` namespace, so `uGamm` is never resolved from the search path and a
`uGamm` model set refits 8 of 8 with `.packages = "FSSgam"` alone; `mgcv` is
genuinely required, 0 of 4 without it. "No longer loaded onto the workers" holds
only for a non-`gamm4` test.fit: unserialising a `uGamm(lme4 = TRUE)` object
loads `gamm4` and `lme4` by itself, and `doSNOW` sends the export environment in
the same `clusterCall` that carries the package vector, before the `library()`
loop runs, so the comment claiming that naming `gamm4` "would only move that load
forward to the dispatch" is reversed. And the dependency change is user-visible
in a way `NEWS.md` did not state: a fresh install no longer brings `gamm4` or
`lme4`, so a case-study-2 style call errors on a clean machine.

The reviewer also recorded, unprompted, that the change is justified on
dependency-hygiene grounds alone and should stand even if the stall result never
replicates.

**Claude:** All three corrected. `worker_packages()` now records what is
necessary and what is defensive, with the measured figures. `NEWS.md` leads with
the install consequence in bold and states the `uGamm(lme4 = TRUE)` qualification
explicitly, so the 44-in-140 figure is not read as extending to those users.

The reviewer's finding that it was the `@importFrom` directive rather than the
`Imports` field that loaded `gamm4` was taken as well: the earlier wording
implied that `MuMIn`, `nnet` and `doSNOW` are loaded on a worker by the `FSSgam`
namespace, which they are not. The two had to change together regardless, since
an `Imports` entry with no `NAMESPACE` import raises "All declared Imports should
be used", and that is what covers the other half of the invariant the test
asserts.

One unclaimed beneficiary was added to `NEWS.md`: `check_correlations(parallel =
TRUE)` dispatches with `.packages = c('nnet')`, but its workers load the `FSSgam`
namespace and so used to load `gamm4` through it.

The mismatch between `Depends: R (>= 3.5)` and `MuMIn`'s `R (>= 4.4.0)` was
raised as FSSgam_package#31. It weakens the reason given for rejecting
`setup_strategy = "sequential"`, which is recorded there.

---

#### Issue #14, review round 3 (fresh session)

**Claude (reviewer):** No blocking or substantive finding. The round 2
corrections were reproduced from private installs of both refs: `library(FSSgam)`
loads `gamm4` and `lme4` on `master` and neither on the branch, and on the branch
the namespace loads only `mgcv`, so the correction that `MuMIn`, `nnet` and
`doSNOW` are not loaded that way is right. A doctored build with `gamm4` back
under `Imports` but no `NAMESPACE` import loads nothing, confirming it was the
directive. A worker simulated in a fresh R process, unserialising the payload and
then attaching, refits a `gam()` model set 8 of 8 with `mgcv` and 0 of 8 without,
and a `uGamm(lme4 = TRUE)` set with nothing attached at all. `doSNOW:::workerInit()`
takes `exportenv` as an argument and runs its `library()` loop in the body, so the
ordering holds by construction rather than by observation.

One nit: the test comment named the wrong check. Re-adding the `@importFrom`
while `gamm4` stays under `Suggests` raises "Namespace dependency missing from
DESCRIPTION Imports/Depends entries", an ERROR; the NOTE named in the comment
catches the reverse. Both were verified on doctored builds.

**Claude:** The comment now names both directions and records that this
expectation covers neither of them, only `gamm4` being made a hard dependency
again.

---

#### Issue #14, the withdrawn measurement remade

**Claude:** The host's load fell back to about 5 on 22 cores once the unrelated
job finished, and the narrower comparison was run again, this time alternating
run by run rather than taking the two arms one after the other: 34 stalls in 50
with `gamm4` in the packages the workers attach and 19 in 50 without, chi-squared
p = 0.005. It agrees with the end-to-end figure in direction and in rough
magnitude.

The lesson is about the design rather than the number. The original 23 in 60
against 0 in 60 was measured with one arm run to completion before the other
began, so it is confounded with anything that changed on the host in between, and
the reviewing session's failure to reproduce it is what that design predicts. The
figure now recorded was taken under the same design as the end-to-end one, which
is the only reason it can be quoted alongside it. Both remain one session's
measurement of a rate on one host.

---

### Issue #16 -- factor-factor deviances computed on different row sets

**Claude:** `fill_factor_factor_correlations()` estimates the correlation
between two factors as `sqrt(1 - fit/null.fit)`. `fit` is fitted on
`na.omit(dat[, c(var.1, var.2)])`; `null.fit` came from a per-factor
intercept-only model that `multinom()` had fitted on `var.1`'s own complete
cases. When `var.2` carries missing values the two are computed on different
numbers of rows and the ratio is not meaningful.

The correction chosen matches the convention the rest of the matrix already
follows -- `cor(use = "pairwise.complete.obs")` for the continuous block and a
per-pair `lm()` for the factor-continuous block are both pairwise complete-case.
Fitting the null per pair unconditionally would have undone the earlier
reduction from n(n-1) null fits to n, so `complete.counts` records how many rows
each factor's own null was fitted on and the null is refitted only for a pair
whose complete-case count is lower. The pair's rows are a subset of the
factor's, so equal counts mean an identical row set and the hoisted value is
exact. Data with no missing values therefore costs nothing extra and returns
identical values.

`null.fit` is now guarded against being a `try-error` alongside `fit`. Before
this the sequential loop would have raised an error out of `round()` if the null
fit failed while the pair fit succeeded; the cell is now left `NA`, which is how
a failed pair fit was already handled.

Two tests were added. The first constructs two factors where `f2` reproduces
`f1` on three rows in four and is missing for most of one level, so the row sets
differ by 25 of 90 rows; the estimate moves from 0.848 to 0.728. It asserts the
matched value and separately asserts the unmatched value is absent, so reverting
the fix fails on both. The second asserts the cached path: with `f2` as the
response the pair's rows are exactly `f2`'s own complete cases, and the hoisted
null is already correct. Confirmed the first test fails against `origin/master`
and the second passes, which is what distinguishes them.

Full suite after the change: 571 passing, 0 failing, 7 skipped, against 568 on
`master`.

---

#### Issue #16, review round 1 (independent session)

**Claude (reviewer):** No blocking defect. The `complete.counts` shortcut was
confirmed exactly equivalent to refitting the null for every pair, by running a
variant with the condition replaced by `if(TRUE)` over 150 generated data frames
and comparing with `identical()`; all 150 matched, and the 10 frames with no
missing values matched `origin/master` as well. `complete.counts` was confirmed
to reach the doSNOW workers. Two substantive points were raised: the null is
refitted before the `try-error` guard, so a pair whose fitted model has already
failed pays for a discarded fit and emits its warnings (measured: 5 warnings
against master's 3 on a pair with no overlapping complete cases); and the second
test's name claims it detects whether a refit occurred, which it cannot.

**Claude:** Both were addressed. The refit condition now also requires `fit` to
have succeeded, which restores the warning count to 3. The second test was
renamed to what it asserts, and a third block was added that does pin the
shortcut: six ordered pairs over three factors with missingness in different
rows, each compared against the null refitted for that pair, so any divergence
between the retained value and the refitted one fails. The defensive nature of
the `null.fit` guard is now stated on the face of the code, the `na.action`
comment was qualified, and `check_non_linear_correlations()` gained the matching
`@details` paragraph.

The sibling function's own missing `null.fit` guard was raised as
FSSgam_package#19 rather than fixed here.

`DESCRIPTION` was bumped to 1.1.0.9000 so that the development `NEWS.md` section
and the declared version agree.

---

#### Issue #16, review round 2 (fresh session)

**Claude (reviewer):** No blocking defect. The three claims singled out in the
round-1 response were re-measured and hold: the warning count is 3 and matches
`master`, the third test block fails under two independent breaks of the
shortcut, and the suite figures reproduce. Four points were raised: the new test
block's comment misstates which of its six pairs take the retained value (two,
not four); the pull request body still described the state at the first commit;
no test reaches the parallel copy of the changed lines, because the only
`parallel = TRUE` fixture is built from `case_study1`, which has no missing
values; and where every factor is missing in rows the others are not, none of the
hoisted nulls is used, so 15 `multinom()` calls are made where refitting per pair
with no hoist would be 12.

**Claude:** The comment was corrected against the block's own counts. The
`parallel = TRUE` fixture now gives `ZONE2` missing values, so both copies of the
refit condition are compared; verified against an installed copy with
`FSSGAM_TEST_PARALLEL=true`, two runs completing and agreeing and a third
stalling, which is FSSgam_package#14.

The hoisted nulls are now fitted only for the factors some pair retains, decided
from the missingness pattern before anything is fitted rather than lazily, which
a cache could not do across doSNOW workers. A factor's null is read only by a
pair that drops none of the rows on which that factor is present, so the
condition is exact. Measured with `trace(nnet::multinom)` on three three-level
factors over 90 rows: 9 calls with no missing values, matching `master`, and 12
with disjoint missingness, down from 15 and equal to the no-hoist minimum. The
unread fits' warnings go with them: the no-overlap frame now raises 2 warnings
where `master` raises 3.

The version bump was kept and `cran-comments.md` was corrected instead. It now
records that a development version is carried between releases, that this adds
the "Version contains large components" NOTE under `R CMD check --as-cran`, that
`devtools::check()` does not run that check at all, and that the version is set
to the release number before a submission is built.

The redirecting `DESCRIPTION` URLs were raised as FSSgam_package#20.

---

#### Issue #16, review round 3 (fresh session)

**Claude (reviewer):** No blocking or substantive defect. Every measurement in
the round-2 response reproduced. The `retains` precomputation was established
exact rather than accepted: the reading condition and the fitting condition are
complementary over the same quantity, since `nrow(na.omit(dat[, c(v, w)]))` is
`sum(present[[v]] & present[[w]])` for every column class
`classify_correlation_predictors()` admits. An instrumented build raising an
error on a `NULL` null with a successful pair fit was run over 400 generated
frames and never fired, recording 200 benign `NULL` reads, all with the pair fit
already failed. A variant with the refit deleted from the `%dopar%` body alone
was built, installed and run against the branch's tests, and reported an error,
which confirms the new fixture reaches the parallel copy.

Four nits: the guard is safe by `&&` short-circuit rather than by construction,
since `inherits(NULL, "try-error")` is `FALSE`; the comment claiming no case
reaching the guard is known is falsifiable; an entirely missing factor column is
an exception to fitting the null only where a pair retains, because zero equals
zero; and `cran-comments.md` describes the version bump as adding a NOTE where
it adds a line to one.

**Claude:** All four applied. The claim about the guard was the interesting one:
under `options(na.action = "na.fail")` the pair model is fitted on `dat.r`, whose
incomplete rows are already removed, while a retained null is fitted on the raw
column and fails, so the guard is load-bearing rather than defensive. Confirmed
directly: under `na.fail` this branch returns 0.7281 and an `NA` for the other
direction where `master` raises "non-numeric argument to mathematical function".
The comment now records that rather than the opposite.

`!is.null(null.fit)` was added to both copies of the guard so it does not depend
on the order of its operands, and the hoisted null is now skipped for a column
with no observed values, which takes the all-missing case from 3 fits to 2
against `master`'s 4.

---

#### Issue #16, review round 4 (fresh session, confirmation)

**Claude (reviewer):** No blocking and no substantive defect. The `na.action`
table reproduced exactly, and `na.pass` was found to behave as `na.fail` does.
Both clauses added in the previous round were established unreachable with a
successful pair fit, by argument and then by fuzzing 400 generated frames
against a build with the shortcut deleted and the null refitted per pair
unconditionally: 400 of 400 matrices `identical()`, no error on any build, and
no cell `NA` on the branch where `master` returned a number. Three nits: the
round-3 comment's claim that the warning count is the same is wrong, the branch
raises 2 where `master` raises 3; neither new clause is covered by a test, both
having been deleted in variant builds with the suite still passing; and a
misaligned continuation line.

**Claude:** The warning claim was wrong and is corrected on the pull request. A
test now pins the zero-count clause by counting `multinom()` calls with
`trace()`, since the returned matrix is the same either way: two pair fits, both
failing on a zero-row frame, and no null at all. It fails with 3 calls when the
clause is removed. The `!is.null()` clause is left untested, because it is
unreachable under the default `na.action` and a test for it would have to set a
global option to reach a defensive branch. The continuation line was realigned.

Rounds 3 and 4 both found nothing substantive, so the cycle on this pull request
is closed.

---
