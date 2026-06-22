## Session: fix-nb-family-issue-12
Date: 2026-06-22
Model: Claude Sonnet 4.6

### Prompts and Responses

**User:** Can you write a plan to create a reprex and address this issue?
https://github.com/beckyfisher/FSSgam/issues/12

**Claude:** Fetched issue #12 (filed by the package author, Dec 2023, no
comments/reprex attached): "There is an issue with fitting the full model
set when the family `nb()` is used. Most models (but not all) appear to
fail. Note the same issues occurs in the use of `lapply`."

Reproduced the bug directly against the installed package: fitting
`test.fit` with `family = mgcv::nb()` on `case_study1`, then running
`generate_model_set()` + `fit_model_set()`, failed 12 of 13 candidate
models (only the null model survived), under both `parallel = FALSE` and
`parallel = TRUE`.

Root-caused it to `fit_mod_l()` in `R/functions_supporting.R`: every
candidate refit was passed `family = stats::family(test.fit.)` — the
*already-fitted* family object. mgcv's extended families that estimate an
extra parameter (`nb()`, `tw()`, and anything using mgcv's generic
`getTheta()`/`putTheta()` convention) store that estimate (theta) inside the
family object's own mutable environment. Reusing that object across every
candidate formula warm-started each refit's theta search from `test.fit.`'s
unrelated estimate, destabilising mgcv's internal IRLS loop for most (but
not all) formulas. `build_null_model()` already omitted the `family=`
override for the null model, which is why it alone always succeeded — the
giveaway for the root cause.

Verified the fix (dropping the `family=` override so `update()` re-evaluates
a fresh family from `test.fit.`'s own call) on a properly *installed* copy
of the package (not just `devtools::load_all()`, to rule out namespace-
isolation artifacts), comparing original vs. patched side by side, under
both bare and namespace-qualified `nb()`, sequential and parallel (`doSNOW`)
execution. An initial apparent total failure (0/13) under `parallel = TRUE`
for *both* old and new code turned out to be a test-harness artifact (a
stale, previously-installed `FSSgam` 1.0.0 on the cluster workers' default
library path) — not a real bug; resolved by pinning worker `.libPaths()`.

Wrote a plan (`Fix GitHub issue #12: models fail to fit with family = nb()`)
covering: the `fit_mod_l()` fix with a decision-point comment, a matching
cross-reference comment in `build_null_model()`, a new regression test in
`tests/testthat/test-fit_model_set.R` mirroring the existing Tweedie test,
a `NEWS.md` bug-fix entry, and porting the fix to `master` individually
(per this repo's branch-porting convention from Phase 8).

Asked the user whether to also post the reprex/explanation as a GitHub
comment and close the issue, or keep that step manual. User chose: fix +
test only, with a draft comment handed back for manual review/posting.

Implemented the plan exactly as approved:
- `R/functions_supporting.R`: removed `family=stats::family(test.fit.)`
  from `fit_mod_l()`'s non-`dsm` branch; added a comment explaining the
  `.Theta`/extended-family mutable-state mechanism and referencing #12.
- `R/generate-model-set.R`: added a short cross-reference comment in
  `build_null_model()` so the two code paths don't drift out of sync again.
- `tests/testthat/test-fit_model_set.R`: added
  `"fit_model_set works with a negative binomial family (extended family
  with estimated theta)"`, mirroring the existing Tweedie test structure.
- `NEWS.md`: added a "Bug fix: ..." bullet under the development-version
  heading, matching the style of the existing Phase 7 bug-fix entries.

Verified on `dev`: full testthat suite (33 tests across 13 files, 0
failures), `devtools::check()` clean (0 errors/warnings/notes), and the
original 13-model reprex now at 0/13 failed. Committed on `dev`
(`29daf70`), then cherry-picked cleanly onto `master` (`9cea95f`) per the
plan's branch-porting step — re-verified the full suite and `devtools::check()`
there too (0 errors/0 warnings; one unrelated "unable to verify current
time" NOTE caused by the sandbox having no network access for R's clock
check, not by the code change). Neither branch has been pushed.

Drafted a GitHub comment for issue #12 (root cause, fix, and the reprex)
for the user to review and post/close herself, per her chosen scope.

---
