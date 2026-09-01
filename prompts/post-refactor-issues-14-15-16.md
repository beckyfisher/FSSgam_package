
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
