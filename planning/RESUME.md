# Pre-CRAN refactor — complete

Branch `pre-cran-refactor`, off `master` at `a8235e2`. All seven phases of
`pre-cran-refactor-human.md` are done. Nothing is pushed; nothing is merged.

## State

| | baseline | now |
|---|---|---|
| `testthat` | 438 passing, 11 skipped | 546 passing, 12 skipped |
| `devtools::check()` | 0 errors, 0 warnings, 0 notes | unchanged |
| `R CMD check --as-cran`, incoming checks on | Status OK | Status OK |
| `resolve_factor_interactions()` | 161 lines | 109 |
| `resolve_smooth_smooth_interactions()` | 73 lines | 55 |

Thirty-three commits, one per item, each behaviour fix carrying its regression
test in the same commit. Two exceptions, both recorded in their commit
messages: A3 is committed with the second branch divergence because fixing one
made the other reachable, and `NEWS.md` is written per phase rather than per
commit, being a per-release document.

## What was verified, and how

Phases 3 to 5 were checked against a golden master of 31 `generate_model_set()`
scenarios, compared with `identical()`. The scripts are in
`planning/golden-master/`, with a README explaining why four of the scenarios
exist — the obvious ones do not discriminate between the duplicated branches,
and each of the four was built against a measured correlation value.

Phase 2 used a separate 14-scenario master over `check_correlations()` and
`check_non_linear_correlations()`.

Every behaviour change was confirmed to fail before the fix, except
FSSgam_package#10, where the warning being removed is not reproducible on
nnet 7.3-20 — stated on the face of the test so it is not mistaken for a
regression test.

## Open issues raised during the work

- **FSSgam_package#14** — `parallel = TRUE` stalls on roughly half of runs when
  `gamm4` is loaded on the doSNOW workers. Not a defect in this package: a
  `foreach()` with a body of `i * 2` reproduces it. Pre-existing; `master` stalls
  at the same rate. Includes a reprex and the measurements.
- **FSSgam_package#15** — a supplied `cor.matrix` must now name the hard coded
  factor interaction columns, which a user cannot know in advance. Follows from
  #13 and is a usability problem rather than a defect; three options given.
- **FSSgam_package#16** — `check_correlations()` compares a pair deviance
  computed on complete cases against a null deviance computed on the whole
  column. Long-standing, reachable only by calling the function directly, since
  `generate_model_set()` rejects predictors with `NA`.

Still open from Phase 13 and unaffected by this work: #6, #9 (now fixed by the
`progress` argument — close it), #12 (fixed here — close it).

## The parallel path

Verified against the final code, one cluster per process:

```
path=saved    AICc matches sequential=TRUE  failed=0  rows=8
path=unsaved  AICc matches sequential=TRUE  failed=1  rows=9
```

The test file itself could not be run to completion: it creates about seven
clusters, and at the roughly one-in-two per-cluster stall rate of
FSSgam_package#14 the chance of all seven completing is about one per cent.
Measured, 10 attempts at the whole file, 10 stalls. Check the parallel path with
a single-cluster script rather than the test file, and do not read a whole-file
stall as a failure of the code under test. Recorded on the issue.

## Before submitting to CRAN

1. Decide whether the behaviour changes warrant `1.1.0` rather than `1.0.0`.
   Four are user-visible and two change output: candidate model names in a
   non-C locale, and which interaction terms are built when a `cor.matrix` is
   supplied. `NEWS.md` describes each.
3. Update `cran-comments.md`.
4. Open a pull request against `master`. Every measurement quoted in the prompt
   log was taken on this WSL host under R 4.6.1; re-measure rather than citing
   those numbers.
