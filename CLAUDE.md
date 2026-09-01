# Project: FSSgam

This file provides **project-specific context only** — repo structure, dependencies,
analysis decisions, and prompt logging. It does not impose personal style preferences
on collaborators.

**For Claude:** If a `CLAUDE.md` exists in the parent directory of this repo,
read it before proceeding — it may contain user-specific environment and style
preferences. This repo-level file takes precedence over any parent-level file
where they conflict.

---

## 1. Repo Type

**R package** — standard usethis/devtools structure. CRAN-ready as of the
v1.0.0 modernisation (Section 6): roxygen2 docs, explicit `Imports`,
testthat suite, GitHub Actions R CMD check, `NEWS.md`, `cran-comments.md`.
Has a pkgdown reference site (function docs only, no Articles/vignettes —
see Section 5).

---

## 2. What This Repo Does

FSSgam implements full-subsets multiple regression for ecological data using
GAMs and GAMMs (via mgcv and gamm4). The two primary user-facing functions are
`generate_model_set()` (builds a candidate model set from a test fit and
predictor list) and `fit_model_set()` (fits and compares that model set by
AICc). `full_subsets_gam()` is a one-shot wrapper around both, recommended
only for small candidate sets since it does not let you inspect the model
set before fitting.

The full public API is now snake_case. `generate.model.set()`,
`fit.model.set()`, and `full.subsets.gam()` (dot-case) are deprecated
aliases retained for backward compatibility because they're used by
existing downstream code (the first two cited in the published paper, the
third used directly in the companion repo's case studies) — see Section 5.
`check.correlations()`, `check.non.linear.correlations()`,
`build.inclusion.mat()`, `extract.mod.dat()`, and `fit.mod.l()` were also
renamed to snake_case, but with no alias — see Section 5 for why that's a
deliberate, asymmetric decision rather than an oversight. `fit_mod_l()` was
then unexported altogether in 1.1.0 (Section 6 Phase 14) and is reachable
only as `FSSgam:::fit_mod_l()`; `wi()`, `extract_mod_dat()` and
`build_inclusion_mat()` are still exported.

`generate_model_set()` (in `R/generate-model-set.R`) and `fit_model_set()`
(in `R/fit-model-set.R`) are each decomposed into named, unexported helpers
— see Section 6, Phase 6 for what they are and how they were verified.

The companion publication repo with all case studies, vignettes, and full
worked examples lives at https://github.com/beckyfisher/FSSgam. That repo is
the permanent, citable reference and must not be modified.

---

## 3. Package Dependencies in Scope

Dependencies are defined in `DESCRIPTION` (Imports + Suggests). Read that file
to determine what is available — don't assume this list stays in sync.

- **Imports**: `doSNOW`, `foreach`, `gamm4`, `mgcv`, `MuMIn`, `nnet`,
  `parallel`, `stats`, `utils`
- **Suggests**: `covr`, `testthat (>= 3.1.5)`
- **Depends**: `R (>= 3.5)` only — no other package may go in `Depends`

`doSNOW` is not installed by default in fresh environments (e.g. this WSL
container) even though it's a long-standing declared dependency — ask before
installing it, but it genuinely is needed for `devtools::check()`/`load_all()`
to succeed.

---

## 4. Key Files and Structure

```
R/
  FSSgam.R                              # package-level doc ("_PACKAGE"), @importFrom directives
  data.R                                 # roxygen2 @docType data entries for bundled datasets
  deprecated.R                          # generate.model.set()/fit.model.set()/full.subsets.gam() —
                                         # .Deprecated() wrappers (Section 6 Phase 9)
  generate-model-set.R                  # generate_model_set() + 8 unexported helpers (Section 6 Phase 6)
  fit-model-set.R                       # fit_model_set() + 4 unexported helpers (Section 6 Phase 6)
  check-correlations.R                  # check_correlations() + 3 unexported helpers (Section 6 Phase 6b)
  check-non-linear-correlations.R       # check_non_linear_correlations() + 3 helpers (Section 6 Phase 6b)
  utils.R                                # classify_correlation_predictors() — shared by both check-*.R files;
                                         # resolve_candidate_family()/clone_independent_family() — shared by
                                         # fit_mod_l() and fit-model-set.R's two per-candidate loops (Section 6 Phase 12)
  full-subsets-gam.R                    # full_subsets_gam(), calls the snake_case names internally
                                         # (renamed from function_full_subsets_gam_v1.11.R, Section 6 Phase 9)
  functions_supporting.R                # wi(), extract_mod_dat(), build_inclusion_mat() (exported but
                                         # documented "not called directly"), and fit_mod_l(), which is
                                         # unexported as of 1.1.0 (Section 6 Phase 14)
tests/testthat/                         # test-<name>.R, testthat 3rd edition
man/                                    # auto-generated by roxygen2, do not edit manually
data-raw/
  case_study_datasets.R   # reads CSVs and calls usethis::use_data() for each
  case_study1_dataset.csv, case_study2_dataset.csv, case_study3_dataset.csv,
  extra_examples_coral_data.csv
                          # .Rbuildignore'd — not part of the built package
data/                     # .rda objects produced by data-raw/case_study_datasets.R
                          # datasets: case_study1, case_study2, case_study3, coral_data
.github/workflows/R-CMD-check.yaml  # r-lib check-standard template; runs on push to
                          # main/master/dev and on pull_request
.github/workflows/pkgdown.yaml      # r-lib pkgdown template; builds + deploys to gh-pages on
                          # push to main/master/dev. Sets the PKGDOWN_MODE env var from the
                          # branch (dev -> devel, else release) and passes it as a pkgdown
                          # build_site() override (Section 6 Phase 11)
_pkgdown.yml              # reference-only site config: no Articles. development.mode is a
                          # static "release" default for local builds; CI overrides it per
                          # branch (Section 6 Phase 11)
DESCRIPTION               # package metadata; Version is kept in sync between master and dev
                          # (Section 6 Phase 11) -- safe to merge/PR normally between branches
NAMESPACE                 # auto-generated by roxygen2, do not edit manually
NEWS.md                   # changelog (one entry per version)
cran-comments.md          # CRAN submission notes; .Rbuildignore'd
ignore/                   # gitignored scratch folder (unlike data-raw/ or prompts/, nothing
                          # here is tracked in git at all). Holds build_package_FFSgam.R
                          # (legacy dev script with a hardcoded Windows path — do not rely on
                          # it) and issue10_reprex.R (manual reproduction script for the
                          # #10/#12 family-resolution conflict, see Section 5 — re-run it by
                          # hand if you ever touch fit_mod_l()'s family resolution again, since
                          # it's the only place that actually exercises the parallel = TRUE
                          # doSNOW path for that regression). Both files used to live at the
                          # repo root (tracked + .Rbuildignore'd) and were moved here, dropping
                          # out of git tracking, in commit 32037a1 (2026-06-23).
prompts/                  # session prompt logs (Section 7); .Rbuildignore'd
planning/                 # the two planning documents for the 1.1.0 pre-CRAN refactor
                          # (pre-cran-refactor-human.md and -claude.md, per the parent
                          # CLAUDE.md Section 13) plus golden-master/, the before/after
                          # comparison scripts that phase was verified with. Tracked in git
                          # and .Rbuildignore'd, so none of it reaches the built package.
```

---

## 5. Known Constraints or Decisions

- **`generate.model.set()`, `fit.model.set()`, and `full.subsets.gam()` are
  permanent, deprecated aliases**, not just historical names to preserve.
  They live in `R/deprecated.R`, call `.Deprecated()`, and forward `...` to
  `generate_model_set()` / `fit_model_set()` / `full_subsets_gam()`. The
  first two are cited in the published paper (Fisher et al. 2018, Ecology
  and Evolution); `full.subsets.gam()` is used directly in the companion
  repo's case-study-2.Rmd (verified by reading that repo's vignettes during
  the Phase 9 rename — see Section 6). Do not remove any of the three, and
  do not let `full_subsets_gam()`, `generate_model_set()`, or
  `fit_model_set()` call through their own deprecated wrappers internally
  (they must call the snake_case names directly, or every call to them
  emits spurious deprecation warnings).

- **`check.correlations()`, `check.non.linear.correlations()`,
  `build.inclusion.mat()`, `extract.mod.dat()`, and `fit.mod.l()` were
  renamed to snake_case with no backward-compatible alias** (Section 6
  Phase 9) — a deliberate, asymmetric decision relative to the three
  aliases above, not an inconsistency. Before renaming, every vignette in
  the companion repo (`case-study-1/2/3.Rmd`, `extra-examples.Rmd`,
  `faq.Rmd`) was checked for direct calls to these five names; none were
  found, unlike `full.subsets.gam()`. If you ever find one of these five
  names called directly in the companion repo (or anywhere else
  load-bearing), that's new information that would justify adding a
  deprecated alias retroactively — don't assume the original check still
  holds without re-verifying.

- **`fit_mod_l()`'s family resolution (`R/functions_supporting.R`,
  `resolve_candidate_family()`/`clone_independent_family()` in
  `R/utils.R`) embodies a fix for two GitHub issues that directly conflict
  (#10 and #12) — see Section 6 Phase 12 for the full mechanism. Do not
  "simplify" this by passing `family = stats::family(test.fit.)` (breaks
  #12: shares one mutable, already-fitted family object across every
  candidate refit) or by dropping back to letting `update()` re-evaluate
  `test.fit.$call$family` directly wherever the refit happens to run
  (breaks #10 under `parallel = TRUE`, and silently fails to give
  independent state per refit whenever `family` was supplied as a
  variable/list element rather than a literal call, regardless of
  parallel). Both were tried, in that order, across this package's
  history, and both regressed the other issue. `ignore/issue10_reprex.R`
  (see Section 4) reproduces the #10 half of the conflict end to end under
  `parallel = TRUE` -- re-run it by hand if you touch this code, since the
  committed test suite deliberately doesn't spin up a real doSNOW cluster.

- **Do not change the GitHub remote URL** for this package repo. Users install
  via `devtools::install_github("beckyfisher/FSSgam_package")` and this must
  continue to work.

- **The publication repo (`beckyfisher/FSSgam`) is strictly read-only** for
  the purposes of this work. Do not commit to it or suggest changes to it.
  Its URLs are permanent references in the published paper (Fisher et al. 2018).

- **Vignettes live in the publication repo, not here.** Do not create a
  `vignettes/` folder in this package repo. Full worked examples and case
  studies are permanently hosted at https://github.com/beckyfisher/FSSgam.
  The `@examples` blocks in roxygen2 should be minimal, self-contained, and
  fast — use the bundled datasets (`case_study1`, `case_study2`, etc.). In
  practice every current example runs in well under 1 second, so none need
  `\donttest{}`; re-time new examples before assuming that's still true.

- **The pkgdown site here is Reference-only — no Articles/vignettes section.**
  Full worked examples and case studies stay at the companion publication
  repo; `_pkgdown.yml`'s `navbar.structure.left` deliberately omits
  `articles`. The `reference:` index is curated into named groups (Core
  functions, All-in-one wrapper, Correlation diagnostics, Supporting
  functions, Deprecated aliases, Datasets) — when adding a new exported
  function, add it to an existing group or a new one, don't leave it
  unlisted (`pkgdown::check_pkgdown()` will fail with "must be a known
  topic name" style errors if `_pkgdown.yml` references a function that no
  longer exists, or "topics not in any index" if a new exported function is
  missing from every group).

- **`master` and `dev` share the same `DESCRIPTION` `Version`.** This was
  not always true. Originally `dev` ran `1.0.0.9000` so pkgdown's
  `development: mode: auto` could resolve `master` to `release` and `dev`
  to `devel` purely from the Version string (Phase 8) — but that meant
  every ordinary merge/PR between the branches fought the deliberate
  divergence. It broke for real once (PR #2, 2026-06-22: merging `dev` into
  `master` dragged `1.0.0.9000` across, fixed by hand afterward) before
  being designed out altogether. Phase 11 moved the release/devel decision
  into `.github/workflows/pkgdown.yaml` (keyed off the branch name via the
  `PKGDOWN_MODE` env var, passed into a manual `as_pkgdown()`/`build_site()`
  call rather than `build_site_github_pages()`'s own override forwarding —
  see the Phase 11 caution note for why that distinction matters), so
  `_pkgdown.yml`'s `development.mode` is now a plain static `release`
  default and `DESCRIPTION`'s `Version` no longer needs to differ between
  branches. Normal merges/PRs between `master` and `dev` are safe again —
  both deploy to the same `gh-pages` branch non-destructively (the deploy
  action uses `clean: false`).

- **Bundled datasets are managed via `data-raw/case_study_datasets.R`.**
  Do not regenerate or modify the `.rda` files in `data/` manually. If new
  data is needed, update `data-raw/case_study_datasets.R` and re-run it.
  Each dataset must have a `@format` roxygen2 entry in `R/data.R` for CRAN
  compliance.

- **`s()`, `te()`, and `t2()` must stay unqualified inside any `gam()`/`gamm()`
  formula** — `mgcv::s(...)` inside a formula breaks `model.frame()`
  construction silently (caught by `try()`, surfaces as an unexplained `NA`).
  mgcv resolves these by literal symbol name during formula parsing, not via
  normal namespaced function dispatch. This bit us once during the 1.0.0
  refactor (see the comment in `function_check_non_linear_correlations_v1.00.R`)
  — when adding `package::function()` namespacing elsewhere, never apply it
  to a smooth constructor written inside a formula, even though it's safe
  (and required for `R CMD check`) everywhere else, including the outer
  `gam()`/`gamm()` call itself.

- **CRAN submission is the end goal.** Avoid patterns that generate NOTES in
  `R CMD check`: use `package::function()` for all non-base calls in `R/`
  (except smooth constructors inside formulas, above), replace `T`/`F` with
  `TRUE`/`FALSE`, prefer `inherits(x, "y")` over `class(x) == "y"`, and
  ensure all exported functions have complete roxygen2 documentation with
  runnable `@examples`.

- **Apache License 2.0** — retain existing licence; do not change it. The
  DESCRIPTION string `Apache License (== 2.0)` is already CRAN-canonical
  (verified via `tools:::analyze_license()`) and needs no accompanying
  `LICENSE` file.

- **The FAQ for this package lives at**
  `https://github.com/beckyfisher/FSSgam/blob/master/vignettes/faq.Rmd`
  (lowercase `faq.Rmd`, inside `vignettes/`, not at the repo root) and is
  the best source of truth for understanding function behaviour and edge
  cases. Read it before writing tests or `@examples`. The case study
  vignettes in the same `vignettes/` folder (e.g. `case-study-1.Rmd`) are
  the best source of real, working `generate.model.set()`/`fit.model.set()`
  call patterns to adapt for examples and tests.

---

## 6. Refactor Phases

### Phases 1–5 — Completed (v1.0.0 modernisation)

All of the following landed in the v1.0.0 modernisation and are now the
baseline state of the repo, not aspirational:

1. **DESCRIPTION and roxygen2** — `Authors@R`, explicit `Imports`, version
   1.0.0, `URL`/`BugReports`, runnable `@examples` on every exported
   function, `NEWS.md`.
2. **snake_case + deprecation wrappers** — `generate_model_set()` /
   `fit_model_set()` are the real implementations; dot-case names are
   `.Deprecated()` wrappers in `R/deprecated.R`. `require()`/bare calls
   replaced with explicit `package::function()` (see the `s()`/`te()`/`t2()`
   exception in Section 5); `T`/`F` → `TRUE`/`FALSE`;
   `class(x) == "y"` → `inherits(x, "y")`.
3. **Tests** — `tests/testthat/`, covering Gaussian and Tweedie model sets,
   the non-data.frame and NA-predictor error paths, and deprecation
   behaviour (including that `full.subsets.gam()` itself stays
   warning-free).
4. **GitHub Actions + badges** — `.github/workflows/R-CMD-check.yaml`
   (Ubuntu/macOS/Windows release + Ubuntu devel), README badges.
5. **CRAN prep** — `cran-comments.md`, confirmed `devtools::check()` is
   clean (0 errors/warnings/notes) and all examples run fast.

Do not re-do this work; build on top of it. If you find a regression in any
of the above, fix it in place rather than reverting to the pre-1.0 style.

### Phase 6 — Decompose the monolithic function bodies — Completed

Done. `generate_model_set()` now lives in `R/generate-model-set.R`, split
into eight unexported helpers (`validate_use_dat()`, `build_null_model()`,
`check_predictor_missingness()`, `resolve_factor_interactions()`,
`resolve_smooth_smooth_interactions()`, `build_predictor_correlation_matrix()`,
`enumerate_candidate_models()`, `build_model_formulas()`). `fit_model_set()`
now lives in `R/fit-model-set.R`, split into four (`fit_and_summarise_saved_models()`
and `fit_and_summarise_unsaved_models()` are kept separate rather than
unified — `save.model.fits=FALSE` exists specifically so fitted model
objects are never all held in memory at once, and a shared helper would
have undone that — plus `compute_model_weights()` and
`compute_variable_importance()`). `function_generate_model_set.R` and
`function_fit_model_set.R` no longer exist.

Public arguments (`use.dat`, `test.fit`, `pred.vars.cont`, etc.) were not
touched — only internal structure moved. No `R/utils.R` was created: after
decomposing, nothing unexported turned out to be genuinely shared across
files — everything outside these two files is either a deprecated wrapper
or an already-exported, independently-documented function
(`check.correlations()`, `wi()`, etc.) that already has a sensible home.

Verified with a golden-master snapshot (14 `generate_model_set()` scenarios
covering every branch, including ones not in the committed test suite —
list-form `factor.smooth.interactions`, `cyclic.vars`, `linear.vars`, a
user-supplied `cor.matrix`; plus 8 `fit_model_set()`-specific scenarios:
`save.model.fits=FALSE`, `report.unique.r2`, `VI.mods='all'`, the
`max.models` override, `parallel=TRUE`) compared before/after each
extraction, on top of the full testthat suite and `devtools::check()`.
If you decompose anything further, use the same approach: a scratch
before/after comparison script for branches the committed tests don't
reach, not just the committed suite alone.

### Phase 7 — Broader test coverage — Completed

Done. `tests/testthat/` now also covers: `check.correlations()` and
`check.non.linear.correlations()` (basic matrices + invalid-column-class
errors); `generate_model_set()`'s factor-factor-interaction,
smooth-smooth-interaction (both the included-`te()` and
excluded-by-`cov.cutoff` cases), and non-linear-correlation code paths
(previously only checked ad hoc, never committed as tests); `full.subsets.gam()`
including its `factor.interactions`/`smooth.interactions`/`size` legacy
arguments; and the `functions_supporting.R` helpers (`wi()`,
`extract.mod.dat()`, `build.inclusion.mat()`, `fit.mod.l()`).

Writing these tests surfaced two genuine pre-existing bugs in
`full.subsets.gam()` (both predated the 1.0.0 rename, now fixed with
regression tests in `test-full_subsets_gam.R`, see `NEWS.md`):
- the deprecated `size` argument was a no-op (it warned but never actually
  fed its value into `max.predictors`);
- `used.data` and `predictor.correlations` in the return value were always
  `NULL` (referenced `model.set$use.dat`/`$cor.matrix`, which don't exist —
  the real fields are `$used.data`/`$predictor.correlations`).

If you find another function whose documented behaviour doesn't match what
it actually does while working in this codebase, fix it in its own commit
with a regression test, the same way — don't fold it into an unrelated
change.

### Phase 6b — Decompose check.correlations()/check.non.linear.correlations() — Completed

Done. `check.correlations()` moved to `R/check-correlations.R`, split into
`build_continuous_correlation_matrix()`, `build_factor_continuous_skeleton()`,
and `fill_factor_factor_correlations()`. `check.non.linear.correlations()`
moved to `R/check-non-linear-correlations.R`, split into
`build_correlation_pair_grid()`, `estimate_non_linear_correlation()`, and
`assemble_non_linear_correlation_matrix()`.

Both functions had byte-identical column-validation/classification logic
(14 lines). That's the genuinely shared utility Phase 6 didn't find —
extracted into `classify_correlation_predictors()` in `R/utils.R`, which
now exists.

Verified the same way as Phase 6: a golden-master snapshot (7
`check.correlations()` + 5 `check.non.linear.correlations()` scenarios
spanning continuous-only/factor-only/mixed/single-column/two-dataset
cases, plus the invalid-column-class error path) compared before/after,
all exact matches. The `parallel=TRUE` path was checked separately in
isolation (see below) rather than inside the snapshot script.

**Caution for next time:** the `parallel=TRUE` scenario hung for several
minutes the first time it ran *while other `devtools::` commands were
running concurrently against the same package directory* (`document()`/
`test()` in another shell while the capture script's own
`devtools::load_all()` was still in flight) — likely resource contention,
not a real bug. Re-run in isolation with `timeout`, it completed in under
2 seconds with correct values. If a parallel/cluster test seems to hang,
check whether something else is concurrently touching the package before
assuming the code is broken.

### Phase 8 — Merge dev into master; pkgdown reference site with dev/release split — Completed

Done. `master` was 18 commits behind `dev` and still had the pre-1.0.0
codebase (`Version: 1.11`, no tests, no GitHub Actions) — none of Phases
1–7/6b had been merged there. Verified the merge was conflict-free first
(master's only unique commit, a `rework` branch PR merge, had an empty
diff against the common ancestor) before merging `dev` into `master` and
pushing.

Added a pkgdown site (Reference only, no Articles — see Section 5) with a
GitHub Actions workflow (`.github/workflows/pkgdown.yaml`, the standard
r-lib/actions template, deploying to `gh-pages`). Set up a real dev/release
docs split using pkgdown's `development: mode: auto`: bumped `dev`'s
`DESCRIPTION` `Version` to `1.0.0.9000` via `usethis::use_dev_version()`
(also added the NEWS.md dev heading), left `master` at the release `1.0.0`,
and added the same `development: mode: auto` line to `_pkgdown.yml` on
both branches individually (not via merge — see Section 5 for why).
Verified the mode resolution directly via `pkgdown:::dev_mode()`: `master`
resolves to `release`, `dev` to `devel`.

**Discovery during this phase:** the live site at
`beckyfisher.github.io/FSSgam_package` was GitHub's default Jekyll
rendering of `master`'s README, not pkgdown's output — confirmed by
diffing the rendered HTML against the raw README and checking the GitHub
Pages API. The `gh-pages` branch already had the correct pkgdown build
sitting in it; GitHub Pages just wasn't pointed at it. If pkgdown docs
ever look wrong again, check what branch GitHub Pages settings actually
point at before assuming the build is broken.

### Phase 9 — Complete the snake_case rename — Completed

Done. `full.subsets.gam()` renamed to `full_subsets_gam()` (file renamed
`function_full_subsets_gam_v1.11.R` → `full-subsets-gam.R` to match the
kebab-case convention from Phase 6), with `full.subsets.gam()` added to
`R/deprecated.R` as a third permanent alias. `check.correlations()`,
`check.non.linear.correlations()`, `build.inclusion.mat()`,
`extract.mod.dat()`, and `fit.mod.l()` renamed outright with no alias —
see Section 5 for the verification behind that asymmetry and the decision
to treat it as deliberate rather than inconsistent.

All internal call-sites updated (`R/generate-model-set.R`,
`R/fit-model-set.R`, `R/utils.R`, `R/FSSgam.R` package doc), all affected
tests updated to call the new names, and `test-deprecated.R` extended with
a `full.subsets.gam()` deprecation test matching the existing
`generate.model.set()`/`fit.model.set()` ones. `_pkgdown.yml`'s reference
index updated to match.

Verified with the full `testthat` suite (32 files, 0 failures) and a clean
`devtools::check()` (0 errors/warnings/notes) — no golden-master snapshot
needed, since this was a pure rename with no logic changes, unlike the
Phase 6/6b decompositions.

### Phase 10 — Suggested next priorities

With Phases 1–9 and 6b complete, candidates for what comes next (none started):
- Tag a release and actually submit to CRAN (or decide what's still
  blocking that).
- The companion docs repo (`beckyfisher/FSSgam`) still calls the deprecated
  dot-case names in its vignettes — see `FSSgam-docs-CLAUDE.md` drafted for
  that repo (a copy may already be in place there as `CLAUDE.md`). Now that
  `full.subsets.gam()` is also a deprecated alias rather than the real
  implementation, the same applies to it as to `generate.model.set()`/
  `fit.model.set()`.

### Phase 11 — Fix `family = nb()` candidate-refit bug (#12); replace the
Version-based pkgdown split with a CI branch override — Completed

Done. `fit_mod_l()` (`R/functions_supporting.R`) was passing
`family = stats::family(test.fit.)` into `update()` for every candidate
refit — the *already-fitted* family object. mgcv's extended families that
estimate an extra parameter (`nb()`, `tw()`, ...) store that estimate in the
family object's own mutable environment, so every refit was warm-started
from `test.fit`'s unrelated estimate, destabilising most (but not all)
fits. Fix: drop the explicit `family=` override so `update()` re-evaluates
a fresh family from `test.fit`'s own call, the same approach
`build_null_model()` (`R/generate-model-set.R`) already used for the null
model — which is why the null model alone always survived. Added a
regression test mirroring the existing Tweedie test in
`test-fit_model_set.R`, using `family = nb()` on `case_study1`. Verified
12/13 candidate models failing before the fix, 0/13 after, sequentially and
under `parallel = TRUE`.

Also replaced the Phase 8 Version-based pkgdown dev/release split (see
Section 5) with a CI branch override, after it caused a second incident in
one day (see the PR #2 note this section used to carry). `_pkgdown.yml`'s
`development.mode` is now a static `release` default; `.github/workflows/pkgdown.yaml`
computes a `PKGDOWN_MODE` env var from `github.ref_name`/`github.base_ref`
(`dev` → `devel`, else `release`).

**Caution for next time:** the obvious implementation —
`pkgdown::build_site_github_pages(override = list(development = list(mode = ...)))`
— silently does *not* work. That function calls `as_pkgdown(pkg, override =
list(destination = dest_dir))` on a plain path first (computing `dst_path`/
`development$in_dev` from whatever mode was active at that point), then
forwards any further `override=` into `build_site()`'s *own*
`as_pkgdown(pkg, override = override)` call -- but by then `pkg` is already
a `pkgdown` object, and `as_pkgdown()` short-circuits for that case
(`if (is_pkgdown(pkg)) { pkg$meta <- modify_list(pkg$meta, override); return(pkg) }`),
patching `$meta$development$mode` without recomputing `dst_path`/`in_dev`.
The override looked like it worked (`pkg$development$mode` showed the new
value) while devel docs silently kept building to `docs/` instead of
`docs/dev/`. The workflow instead calls
`pkgdown::as_pkgdown(".", override = list(destination = "docs", development = list(mode = ...)))`
directly -- passing both overrides in the *same* call, while `pkg` is still
a plain path string -- then `pkgdown::clean_site()`, `pkgdown::build_site()`,
and the unexported `pkgdown:::build_github_pages()` (the same three calls
`build_site_github_pages()` makes internally). Verified end-to-end locally
for both modes: `release` built to `docs/index.html`, `devel` to
`docs/dev/index.html`, non-destructively (matching the deploy action's
`clean: false`). If you ever touch this step again, re-verify `dst_path`
and `development$in_dev` directly rather than trusting that
`development$mode` alone reflects what got built.

`dev`'s `DESCRIPTION` `Version` was reset from `1.0.0.9000` back to `1.0.0`
to match `master`; going forward the two branches keep the same `Version`
and merge/PR normally.

### Phase 12 — Fix the conflict between issues #10 and #12 in
`fit_mod_l()`'s family resolution — Completed

Done. Asked for a reprex of issue #10 ("Models fail to fit when
`parallel = TRUE`"), reported 2022-06-29 by the repo owner: when `family`
is supplied to `test.fit` as a variable or list element (e.g.
`family = family.vec[[2]]`) rather than written out literally, every
candidate model failed to fit once `parallel = TRUE`. `git log --follow`
on `R/functions_supporting.R` found the original same-day fix
(`a0a901d`, "fix for #10"): it added
`family = family(test.fit.)` to the `update()` call in `fit_mod_l()` --
exactly the line Phase 11 removed four years later to fix #12. The two
fixes directly conflict: #10 needs `update()` to *not* re-evaluate
`test.fit.$call$family` on whatever process actually runs the refit (a
doSNOW worker never has the calling session's variables in scope); #12
needs every refit to get its *own* independent, unfitted family object
(reusing `test.fit.`'s already-fitted one, or naively re-evaluating an
expression that's really just a list lookup, both warm-start refits from
an unrelated estimate, or from each other's). Toggling Phase 11's one line
on or off can only ever fix one at the other's expense.

Verified the conflict empirically before fixing anything: confirmed via
`devtools::load_all()` (the installed copy in this environment had gone
stale, predating the Phase 9 rename -- see the caution below) that current
`master` reproduces #10 exactly (8/8 candidates fit with `parallel =
FALSE`, 0/8 with `parallel = TRUE`, erroring on `object 'family.vec' not
found`), and confirmed empirically that `R`'s `match()` does not reliably
compare environments by reference identity (it silently returned the same
index for every element of a list of distinct environments in this R
version) -- relevant because the fix below clones family objects by
walking their component environments.

Fix: added `resolve_candidate_family()` and `clone_independent_family()`
(`R/utils.R`, shared the same way `classify_correlation_predictors()` is).
`resolve_candidate_family(test.fit)` evaluates `test.fit$call$family`
once -- in `environment(stats::formula(test.fit))`, not wherever the
caller happens to be -- then clones any mutable state via
`clone_independent_family()`, which deep-copies every *anonymous*
closure environment (`environmentName(e) == ""`) reachable from the
family object's function components, leaving shared/named environments
(base/stats namespaces, the enclosing package namespace itself) untouched
so internal helper lookups (e.g. `nb()`'s `postproc` calling mgcv's
unexported `find.null.dev()`) keep working. `fit_model_set()`'s two
per-candidate loops (`fit_and_summarise_saved_models()`/
`fit_and_summarise_unsaved_models()`, `R/fit-model-set.R`) now call this
once per candidate *before* `parallel::makeCluster()` is ever invoked,
building a `family.list` that gets passed into `fit_mod_l()`'s new
`family.` parameter explicitly -- so by the time anything crosses to a
worker, `family.` is already a resolved, independent value, never an
expression needing further lookup. `fit_mod_l()`'s own default for
`family.` (`resolve_candidate_family(test.fit.)`) preserves its existing
`@examples` direct-call behaviour.

`build_null_model()` (`R/generate-model-set.R`) was deliberately *not*
touched -- it runs once, sequentially, always on the same process that
called `generate_model_set()`, so it was never exposed to either bug;
its existing comment already warns not to add `family=` back in.

Caution for next time: don't try to recover "the environment `test.fit`
was created in" from `test.fit` itself beyond what's done here. Checked
`environment(stats::formula(test.fit))`, `attr(test.fit$terms,
".Environment")`, and `environment(test.fit$formula)` -- for a plain
`gam()` fit, mgcv normalises *all three* away to `globalenv()`, even when
`test.fit` was actually constructed inside another function's local
scope (e.g. inside a `testthat::test_that()` block). This is true of the
pre-Phase-12 `update()`-based mechanism too (it reaches `globalenv()` via
namespace lookup, same destination, different route) -- it is not a
regression this fix introduced, and resolving it properly would mean
threading the real calling environment through `generate_model_set()`/
`fit_model_set()`'s public APIs, which is out of scope for what #10
actually reported (top-level script usage). The two new regression tests
that exercise `family` as a list element (`test-fit_model_set.R`,
`test-functions_supporting.R`) assign that list to `globalenv()` via
`<<-` with `on.exit()` cleanup specifically to route around this, with a
comment explaining why -- don't "simplify" those back to a local variable.

Also reinstalled the package from source (`R CMD INSTALL .`) partway
through this phase after discovering the installed copy was stale
(predated the Phase 9 rename -- still had `fit.mod.l()`, not
`fit_mod_l()`). A `doSNOW` worker always loads `FSSgam` via `library()`
from the installed library path, never from `devtools::load_all()`'s
in-memory version, so testing the `parallel = TRUE` path specifically
requires the installed copy to actually be current -- `load_all()` alone
is not sufficient and will misleadingly appear to reproduce or fix
nothing changed about the parallel path.

Verified with the full `testthat` suite (now 99 tests, up from 93, 0
failures) and a clean `devtools::check()` (0 errors/warnings, 1 NOTE
about verifying the current time -- a sandbox/network artifact, not
related to this change). The parallel = TRUE crash itself was verified
manually outside the committed suite (matching the precedent set in
Phase 6b for the same reason: CRAN check overhead/flakiness from spinning
up a cluster in tests) -- 8/8 candidates now fit under `parallel = TRUE`
in the `family.vec[[2]]` reprex, versus 0/8 before this fix.

### Phase 13 — Comprehensive test suite (FSSgam_package#5) — Completed

Done. The suite went from 105 passing expectations / 71.33% line coverage
to 447 / 94.27%, with `tests/testthat/helper-fixtures.R` added so the
eight-line `use.dat`/`test.fit`/`model.set` preamble is written once.
Runtime went from 4.8 s to about 15 s (this WSL host, R 4.6.1; it is
machine dependent, so re-measure before and after rather than comparing
against these).

**Counting expectations.** The figure appears in three places -- `NEWS.md`,
this file and the pull request body -- and went out of step between them in
five consecutive review rounds, every time because one was updated and the
others were not. Update all three together or not at all. Quote testthat's
own `PASS` figure, which is `sum(as.data.frame(res)$passed)`. `sum(...$nb)` is *not* that: it includes
skipped expectations, so with the six skipped parallel tests it reads six
high, and every count recorded during Phase 13 was wrong by exactly that
margin until the sixth review round caught it. `master` and `fix_issues`
both report 105 across 36 `test_that()` blocks and have no skips, so the
baseline is unaffected. The "93 tests" figure in the Phase 12 entry does not
reproduce at all.

New files: `helper-fixtures.R`, `test-generate_model_set_formulas.R`
(cyclic.vars, linear.vars, bs.arg, the gamm4 `t2()` branch),
`test-fit_model_set_options.R` (save.model.fits, max.models, r2.type,
report.unique.r2, VI.mods), `test-fit_model_set_parallel.R`,
`test-utils.R`, `test-data.R`, and `test-numerical-regression.R`
(`expect_snapshot_value()` snapshots of the model table, variable
importance and candidate formulas for five end-to-end scenarios).

Points to note before extending the suite again:

- **testthat forces `LC_COLLATE=C`.** Candidate model names are built by
  `sort()`ing term names, so under testthat a factor named `ZONE` sorts
  ahead of a lowercase continuous predictor -- `"ZONE+complexity"`, not
  the `"complexity+ZONE"` you get interactively in `en_US.UTF-8`. Any
  test asserting on a `modname` must be written against C collation.
  Raised as FSSgam_package#8, since it also means a saved analysis is not
  reproducible across machines with different locales.
- **`library(mgcv)` is called in `helper-fixtures.R`, deliberately.**
  mgcv's extended-family constructors (`tw()`, `nb()`) build their
  component closures in an environment whose parent is `.GlobalEnv`, not
  the mgcv namespace, so `ldTweedie()` and friends only resolve when mgcv
  is attached. Any fixture using `tw()` fails with "could not find
  function ldTweedie" otherwise.
- **gamm4 fixtures need `k >= 4`.** gamm4 represents each smooth as a
  random effect whose grouping factor has (number of basis columns)
  levels; at `k = 3` that is one level, which lme4's `checkNlevels()`
  rejects outright ("grouping factors must have > 1 sampled level"). This
  originates in gamm4 rather than in this package, but it does mean the
  `k = 3` used everywhere else in the suite is unusable for a
  `uGamm(lme4 = TRUE)` test.fit.
- **`fit_quietly()`/`full_subsets_quietly()`** wrap the fitting calls in
  `capture.output()` purely to keep `fit_model_set()`'s unconditional
  `txtProgressBar` out of the reporter's output (FSSgam_package#9).
  Warnings and errors pass through, so `expect_warning()`/`expect_error()`
  still work. Four groups of tests deliberately keep their original direct
  call form, because issue FSSgam_package#5 asks for them to be preserved
  as they stand: the two beckyfisher/FSSgam#10/#12 family resolution tests,
  the two `full_subsets_gam()` Phase 7 regression tests, and everything in
  `test-deprecated.R`.
- **The tests may not use `deparse1()`.** It arrived in R 4.0.0 and
  `DESCRIPTION` declares `R (>= 3.5)`, so the suite has to run without it.
  `deparse_one()` in `helper-fixtures.R` is the replacement, and matches
  `deparse1()`'s `width.cutoff = 500L`. `deparse()`'s own default of 60
  wraps 8 of the 39 formulas across three representative candidate sets,
  and a wrapped formula deparses to a different string, so the cutoff
  matters. `testthat (>= 3.1.5)` is declared for the same class of reason:
  `expect_no_warning()` arrived in that version.
- **The numerical snapshots run everywhere, with one exclusion.** The
  binomial `gamm4`/`uGamm` scenario omits `edf` from its numeric snapshot,
  because gamm4 reports a smooth's edf through lme4's random-effect
  machinery and that column therefore tracks the lme4 version: measured
  4.000/5.000/8.000 here against 4.260/5.050/8.350 on `ubuntu-latest`. It
  is asserted structurally in that scenario instead. Everything else is
  compared wherever the suite runs. The comparison is split into a group
  compared exactly (the columns `compute_model_weights()` rounds to three
  decimal places) and a group compared at 1e-6 relative (the unrounded fit
  statistics), because `expect_snapshot_value()` applies its tolerance
  element-wise and relatively, and one setting cannot serve both.
- **`break_one_candidate()`** injects an unfittable formula into a model
  set to exercise the partial-failure paths. This is more reliable than
  locating real data on which some candidates fail and others do not.
- **Every `parallel = TRUE` test is opt-in.** `skip_unless_parallel_opt_in()`
  requires `FSSGAM_TEST_PARALLEL=true` *and* an installed (not pkgload)
  copy -- a doSNOW worker always loads the installed package (Phase 12),
  so running these against a `load_all()` copy tests nothing. The opt-in
  exists because doSNOW cluster startup stalls indefinitely on a loaded
  machine (see the caution below); an unattended stall consumes the entire
  runtime of a `devtools::check()`, `covr` or CI job. Note that
  `devtools::check()` sets `NOT_CRAN=true`, so `skip_on_cran()` alone does
  not protect it. Consequence: `covr` never executes those branches, which
  is why `check-correlations.R` and `fit-model-set.R` sit below the rest on
  line coverage -- every uncovered line in both is inside a
  `parallel == TRUE` branch. `.github/workflows/parallel-tests.yaml` is the
  one place they run automatically; it installs the package, sets the
  variable, and carries its own `timeout-minutes` so a stalled cluster
  fails that job alone. Locally:

  ```
  R CMD INSTALL .
  FSSGAM_TEST_PARALLEL=true NOT_CRAN=true Rscript -e \
    'library(FSSgam); testthat::test_dir("tests/testthat", package = "FSSgam")'
  ```

**Issue numbering.** The `#10` and `#12` referenced throughout Phases 11-12,
in `R/utils.R`, `R/fit-model-set.R`, `R/functions_supporting.R` and their
tests, are issues in the *publication* repository `beckyfisher/FSSgam`, not
in this one. This package repository now has its own `#10` *and* its own
`#12`, both filed during Phase 13, and both about something else entirely --
so a bare `#10` here resolves to the wrong issue. Every reference in `R/`
and `tests/` was qualified during Phase 13 for that reason. Write
`beckyfisher/FSSgam#10` or `FSSgam_package#6` anywhere the reference has to
survive on its own -- in `R/`, in `tests/`, in `NEWS.md`. Bare numbers are
fine only in prose that has already named the repository, as the paragraph
below does.

Six bugs were fixed: single-predictor model sets failing in
`check_correlations()`; phantom `NA.by.<factor>` terms when
`pred.vars.cont = NA`; interaction terms silently dropped for every
`linear.vars` entry after the first; the `case_study1` `@format` variable
count; two `@return` blocks that did not describe what their functions
return; and a second route to the dropped-interaction defect, via the list
form of `factor.smooth.interactions`.

The Phase 7 precedent of one commit per bug was followed only partly, which
is a deviation worth recording rather than repeating: `f599a2b` carries two
fixes and `a47d114` carries two more, folded in with a round of review
responses. Each pair was found together and each is two instances of one
theme, which is the reason but not a justification -- separate commits would
have been easier to revert or cite individually.

Six further findings were raised as issues in this repository rather than
resolved here, being behaviour decisions rather than defects:
FSSgam_package#6 (`extract_mod_dat()` returns NA r2 for `gamm` fits under
the default `r2.type`), #7 (`full_subsets_gam()` does not
forward `VI.mods`), #8 (locale-dependent model names), #9 (progress bar
cannot be suppressed), #10 (spurious "NaNs produced" warnings out of
`nnet`), and #12 (the factor-factor correlation diagonal is not exactly 1).
Two of those are pinned by an expectation that must change when the issue is
fixed: #6 (`test-functions_supporting.R` asserts the `NA` r2) and #7
(`test-full_subsets_gam.R` asserts `VI.mods` is absent from `formals()`).

#8 is a partial case. Candidate names are asserted throughout the suite under
C collation, which testthat forces, so any change to how they are built shows
up -- but the natural fix for #8 is to sort by byte order, under which the
asserted names are exactly what they already are, and no test changes.

#9, #10 and #12 are recorded in comments or in this file, and fixing any of
them would change no test: #9 has no expectation depending on the progress
bar and is named here rather than in a comment; #10's
`suppress_nnet_nans()` tolerates the warning being absent, which it already
is on nnet 7.3-20; and #12's diagonal is asserted as `> 0.99`, which keeps
passing once it is exactly 1.

That asymmetry is deliberate: an expectation that fails when a defect is
corrected is worse than a comment. Do not read the comments as coverage.

**Caution for next time:** doSNOW stalled repeatedly while this work was
done, including for a trivial `foreach()` containing no FSSgam code.
Measured on this WSL host: `test-fit_model_set_parallel.R` completed
cleanly on two of four consecutive attempts (14 expectations, 0 failures)
and stalled past 180 seconds on the other two. Before concluding that a
`parallel = TRUE` path is broken, re-run in isolation and re-run more than
once.

**Three statements this caution originally carried were wrong, and were
corrected on 2026-09-01 after the stall was measured properly. Do not
re-derive them.** Each one sent a later session down a blind alley. See
FSSgam_package#14, which holds the measurements and a reprex.

- *"cluster startup stalled"*. Startup completes. Instrumenting
  `fit_and_summarise_unsaved_models()` and running until a stall showed
  `parallel::makeCluster()` and `doSNOW::registerDoSNOW()` both returning;
  the last marker reached is immediately before the `%dopar%`. The stall is
  in the dispatch, and during it both workers are alive at 1-2% CPU.
- *"the same file's contents run from a plain script completed six times
  out of six"*. A plain script stalls at the same rate as the test file.
  The distinction that matters is not script against test file but **one
  cluster per fresh R process**: six consecutive clusters inside a single
  session complete 6/6, because after the first the packages are already
  loaded on the workers. Any measurement of this must use a fresh process
  per repetition.
- *"resource contention, not a real bug"* (the Phase 6b caution). It
  reproduces on an idle machine.

What it actually is: loading `gamm4` onto the workers. With a loop body of
`i * 2` and no FSSgam code executed, `.packages` alone moves the rate from
1 stall in 50 runs to 21 in 50. A stall was seen once without `gamm4`, so
it is a strong predictor rather than a precondition. Not established:
whether this is specific to WSL2 or to this lme4/Matrix build, and the
mechanism.

Consequence for anything touching the parallel path: a stall is not
evidence that the code under test is broken, and a single clean run is not
evidence that it works. Both were concluded during the pre-CRAN refactor
from one observation each, and both were wrong. Measure a rate across fresh
processes.

### Phase 14 — Pre-CRAN refactor: defects, argument validation, interaction
restructure; release 1.1.0 (PR #17) — Completed

Done. Scoped from a full read of the 1,947 lines in `R/` against `a8235e2`,
planned in `planning/pre-cran-refactor-human.md` and its `-claude.md`
companion, and carried out in seven phases: defects and validation, the two
correlation functions, a no-behaviour-change restructure, the behaviour
changes, idioms, unexporting `fit_mod_l()`, documentation. `NEWS.md` describes
every user-visible change; do not restate it here.

The four things worth carrying forward:

- **Interaction resolution is no longer four copies of one screen.**
  `combine_uncorrelated()` and `combination_sizes()` (`R/generate-model-set.R`)
  are shared by both `factor.factor.interactions` branches and both
  `smooth.smooth.interactions` branches. The two `factor.smooth.interactions`
  forms are normalised to one triple by
  `normalise_factor_smooth_interactions()` and consumed by
  `build_factor_smooth_terms()`. Four divergences between the duplicated
  branches were held at their existing behaviour behind explicit arguments
  through the restructure, then resolved separately, so the restructure itself
  is verifiable as a no-op. Three were fixed; one is retained by decision and
  documented on the face of the code: `smooth.smooth.interactions = TRUE`
  builds bivariate `te()` terms while a character vector builds terms up to
  `max.predictors`.

- **The correlation matrix is resolved once, before the `te()` terms are
  chosen.** One matrix now governs every stage that screens on correlation:
  which factor interaction columns are built, which `te()` terms are built,
  and which candidates survive. A supplied `cor.matrix` replaces the automatic
  estimate outright — `check_correlations()` is not called at all — so a
  predictor of a class it rejects can now be used by supplying one. Do not
  reintroduce a per-stage recomputation from `use.dat`.

- **The smoothing basis is chosen from the variable name as each term is
  built** (`bs_for()` in `build_model_formulas()`), never by rewriting the
  assembled term strings afterwards. The old form used `grep(cyclic.vars[r],
  term)`, unanchored and as a regular expression, and needed a second block to
  repair the `te()` terms the first had damaged. A name in `cyclic.vars` that
  matches no continuous predictor is now silently ignored rather than
  accidentally matching another predictor; that is not validated.

- **`sort(method = "radix")`** is what makes candidate names locale-independent
  (FSSgam_package#8). Any new code that names a candidate by sorting its terms
  must pass it.

Verified with a golden master of 31 `generate_model_set()` scenarios and 14
over the two correlation functions, compared with `identical()`; the scripts
are kept in `planning/golden-master/`, with a README explaining why four of
the 31 exist (the obvious scenarios do not discriminate between the duplicated
branches). Every behaviour change was confirmed to fail before its fix, with
one stated exception: FSSgam_package#10, where the warning being removed is
not reproducible on nnet 7.3-20 — written on the face of that test so it is
not mistaken for a regression test.

**Measure the suite, do not quote it.** The Phase 13 note about the three
places the expectation count is written applies, and was broken again here:
the figures in the first PR body (546 passing, 12 skipped; baseline 438, 11)
reproduced nowhere. Measured at the branch head on this host (Debian WSL2,
R 4.6.1, testthat 3.3.2) with `testthat::test_local()` and
`sum(as.data.frame(res)$passed)`: 567 passing, 0 failing, 7 skipped blocks,
against 447 and 6 on `master`. The most likely cause of the earlier figures is
a run that was regenerating `tests/testthat/_snaps/numerical-regression.md`,
which skips five blocks and the nine expectations inside them. Re-measure on a
clean checkout, with the snapshot file present.

**Check the parallel path with a single-cluster script, not with
`test-fit_model_set_parallel.R`.** That file creates about seven clusters, and
at the per-cluster stall rate of FSSgam_package#14 the chance of all seven
completing is about one per cent: 10 attempts at the whole file on this host
gave 10 stalls. Both fitting paths were verified against the final 1.1.0 code
one cluster per process, and both reproduce the sequential `AICc`, the unsaved
path recording the injected failing candidate rather than aborting the run.

**Four files in `R/` use CRLF line endings**: `data.R`, `FSSgam.R`,
`full-subsets-gam.R` and `functions_supporting.R`. An editor that rewrites
them as LF turns a three-line change into a whole-file diff. Check with
`file R/*.R` before editing one.

---

## 7. Prompt Log

Session logs for this project are in `prompts/`. Use a short kebab-case
descriptor as the filename for each session (e.g. `cran-modernisation.md`).
If a file with that name already exists, append to it.

### What to log

**Log a session only if it changes the package itself** — code in `R/`, tests,
`man/`, `DESCRIPTION`/`NAMESPACE`, `NEWS.md`, `_pkgdown.yml`, or workflow
files. The publication-disclosure requirement covers work that could bear on a
publication; a session that leaves the package byte-identical cannot.

Changes confined to `CLAUDE.md` or to `prompts/` itself do not count — those
are session housekeeping, not package work.

This is deliberately stricter than the parent `CLAUDE.md` Section 10, and
overrides it. In particular, do **not** open or append to a log for:

- git housekeeping — fetching, merging, rebasing, fast-forwarding a branch,
  pushing, reporting what is behind or untracked, deleting merged branches;
- reading, listing, searching or summarising files with no change made;
- environment checks (which R is active, whether a package is installed);
- answering a question about how the package behaves;
- writing or commenting on a GitHub issue or pull request without a code change;
- editing this `CLAUDE.md`, or tidying a prompt log.

If such a session goes on to change the package, log it from that point —
the log records the work that produced the change, not the preamble.

### Log format

```
## Session: <descriptor>
Date: <YYYY-MM-DD>
Model: <model name and version>

### Prompts and Responses

**User:** <prompt text>

**Claude:** <summary of response — full code blocks where relevant, prose summarised>

---
```
