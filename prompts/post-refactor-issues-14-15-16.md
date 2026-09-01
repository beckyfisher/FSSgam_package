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

**Claude:** Identified the three issues as FSSgam_package#14, #15 and #16, and
put the behaviour decisions in #15 and #16 to the user before writing any code,
since each was raised with more than one candidate resolution and none recorded.
This branch carries the work on #14 alone.

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
