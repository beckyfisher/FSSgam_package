## Session: fix-family-parallel-issue-10
Date: 2026-06-22
Model: Claude Sonnet 4.6

### Prompts and Responses

**User:** can you write a reprex for this issue?
https://github.com/beckyfisher/FSSgam/issues/10.

**Claude:** Fetched issue #10 (filed by the package author, 2022-06-29):
"When family is passed as an element of a vector, all models fail to fit
when parallel is TRUE." Her own comment on the issue: "the failure message
relates to being unable to find the object being used to pass the family
to the test model."

Read `fit_mod_l()` (`R/functions_supporting.R`) and `fit_model_set()`'s
two per-candidate loops (`R/fit-model-set.R`): candidate refits call
`update(test.fit., formula=formula.l, data=use.dat)` with no explicit
`family=`, so `update()` re-evaluates the original, unevaluated `family`
expression captured in `test.fit$call` — wherever that refit happens to
run. Under `parallel = TRUE`, that's a doSNOW worker process that never
had the calling session's variables (e.g. `family.vec`) exported to it.

Wrote a reprex (`issue10_reprex.R`, package root) using `family =
family.vec[[2]]` on `case_study1`. First attempt against `library(FSSgam)`
showed no failure in either mode — turned out the *installed* package
predated the Phase 9 snake_case rename entirely (still had `fit.mod.l()`,
not `fit_mod_l()`). Switching to `devtools::load_all(".")` reproduced it
exactly: 8/8 candidates fit with `parallel = FALSE`, 0/8 with `parallel =
TRUE`, erroring `object 'family.vec' not found`.

**User:** (interrupting a tangent into `doSNOW`/MPI cluster-type
debugging) Apologies, I should have mentioned that a previous bug (issue
12) that we fixed recently may have fixed this also.

**Claude:** Checked the opposite direction first: `git log --follow` on
`R/functions_supporting.R` found `a0a901d` ("fix for #10", same day as the
issue), which added `family = family(test.fit.)` to the `update()` call —
exactly the line Phase 11 (today, earlier this session block) removed to
fix #12. The two fixes directly conflict; Phase 11 had reopened #10, not
fixed it. Confirmed by re-running the reprex against current `master`
(post-#12-fix): same 8/8 vs 0/8 split as above.

**User:** Yes, implement the fix for the #10/#12 conflict.

**Claude:** Implemented `resolve_candidate_family()` and
`clone_independent_family()` (`R/utils.R`). The former evaluates
`test.fit$call$family` once, in `environment(stats::formula(test.fit))`
rather than wherever the caller happens to be; the latter deep-copies
every *anonymous* closure environment reachable from the family object's
function components (skipping shared/named ones like the `stats`/`mgcv`
namespaces, identified via `environmentName(e) == ""`), so refits never
share mutable extended-family state (theta etc.) regardless of how
`family` was specified. Discovered along the way that `match()` does not
reliably compare environments by reference identity in this R version
(silently returned the same index for every element of a list of distinct
environments) — had to write the environment-to-clone lookup by hand using
`identical()` instead.

`fit_model_set()`'s two per-candidate loops now resolve one independent
family per candidate *before* `parallel::makeCluster()` runs, passing the
result into `fit_mod_l()`'s new `family.` parameter explicitly, so nothing
crossing to a worker is ever an unresolved expression.

Hit two more snags while verifying: (1) testing the actual `parallel =
TRUE` path requires the package to be genuinely *installed*
(`R CMD INSTALL .`) — a doSNOW worker loads `FSSgam` via `library()` from
the installed path, never from `devtools::load_all()`'s in-memory version,
so the stale install silently made the parallel path untestable until
reinstalled. (2) A new regression test (family supplied via a list
element, called from inside a `testthat::test_that()` block) failed for
an unrelated reason: `mgcv::gam()` normalises away the environment a
formula was written in, so `environment(stats::formula(test.fit))`
resolves to `globalenv()` even when `test.fit` was built inside another
function's local scope — true of the pre-fix mechanism too, not a
regression, but it meant the test's own local variable was invisible to
either mechanism. Fixed by assigning that variable to `globalenv()` via
`<<-` with `on.exit()` cleanup, matching how a real top-level user script
would have it in scope.

Added two regression tests (`test-fit_model_set.R`,
`test-functions_supporting.R`); verified the second one fails clearly
against the pre-fix (Phase-11-only) code via `git stash` before confirming
it passes against the fix. Full `testthat` suite: 99 tests (up from 93), 0
failures. `devtools::check()`: 0 errors/warnings, 1 unrelated NOTE
("unable to verify current time" — sandbox/network artifact). Updated
`NEWS.md` (amended the existing unreleased #12 bullet rather than adding a
confusing third entry, since both fixes land in the same unreleased
version), `CLAUDE.md` (Section 5 caution bullet + Section 6 Phase 12), and
`.Rbuildignore` (added `issue10_reprex.R`, following the existing
`build_package_FFSgam.R` precedent for ungit-ignored scratch scripts).

---
