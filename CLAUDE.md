# Project: FSSgam

This file provides **project-specific context only** — repo structure,
dependencies, analysis decisions, and prompt logging. It does not impose
personal style preferences on collaborators.

**For Claude:** If a `CLAUDE.md` exists in the parent directory of this
repo, read it before proceeding — it may contain user-specific
environment and style preferences. This repo-level file takes precedence
over any parent-level file where they conflict.

------------------------------------------------------------------------

## 1. Repo Type

**R package** — standard usethis/devtools structure. CRAN-ready as of
the v1.0.0 modernisation (Section 6): roxygen2 docs, explicit `Imports`,
testthat suite, GitHub Actions R CMD check, `NEWS.md`,
`cran-comments.md`. Has a pkgdown reference site (function docs only, no
Articles/vignettes — see Section 5).

------------------------------------------------------------------------

## 2. What This Repo Does

FSSgam implements full-subsets multiple regression for ecological data
using GAMs and GAMMs (via mgcv and gamm4). The two primary user-facing
functions are
[`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate_model_set.md)
(builds a candidate model set from a test fit and predictor list) and
[`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.md)
(fits and compares that model set by AICc).
[`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full_subsets_gam.md)
is a one-shot wrapper around both, recommended only for small candidate
sets since it does not let you inspect the model set before fitting.

The full public API is now snake_case.
[`generate.model.set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate.model.set.md),
[`fit.model.set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit.model.set.md),
and
[`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.md)
(dot-case) are deprecated aliases retained for backward compatibility
because they’re used by existing downstream code (the first two cited in
the published paper, the third used directly in the companion repo’s
case studies) — see Section 5. `check.correlations()`,
`check.non.linear.correlations()`, `build.inclusion.mat()`,
`extract.mod.dat()`, and `fit.mod.l()` were also renamed to snake_case,
but with no alias — see Section 5 for why that’s a deliberate,
asymmetric decision rather than an oversight.

[`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate_model_set.md)
(in `R/generate-model-set.R`) and
[`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.md)
(in `R/fit-model-set.R`) are each decomposed into named, unexported
helpers — see Section 6, Phase 6 for what they are and how they were
verified.

The companion publication repo with all case studies, vignettes, and
full worked examples lives at <https://github.com/beckyfisher/FSSgam>.
That repo is the permanent, citable reference and must not be modified.

------------------------------------------------------------------------

## 3. Package Dependencies in Scope

Dependencies are defined in `DESCRIPTION` (Imports + Suggests). Read
that file to determine what is available — don’t assume this list stays
in sync.

- **Imports**: `doSNOW`, `foreach`, `gamm4`, `mgcv`, `MuMIn`, `nnet`,
  `parallel`, `stats`, `utils`
- **Suggests**: `testthat (>= 3.0.0)`
- **Depends**: `R (>= 3.5)` only — no other package may go in `Depends`

`doSNOW` is not installed by default in fresh environments (e.g. this
WSL container) even though it’s a long-standing declared dependency —
ask before installing it, but it genuinely is needed for
`devtools::check()`/`load_all()` to succeed.

------------------------------------------------------------------------

## 4. Key Files and Structure

    R/
      FSSgam.R                              # package-level doc ("_PACKAGE"), @importFrom directives
      data.R                                 # roxygen2 @docType data entries for bundled datasets
      deprecated.R                          # generate.model.set()/fit.model.set()/full.subsets.gam() —
                                             # .Deprecated() wrappers (Section 6 Phase 9)
      generate-model-set.R                  # generate_model_set() + 8 unexported helpers (Section 6 Phase 6)
      fit-model-set.R                       # fit_model_set() + 4 unexported helpers (Section 6 Phase 6)
      check-correlations.R                  # check_correlations() + 3 unexported helpers (Section 6 Phase 6b)
      check-non-linear-correlations.R       # check_non_linear_correlations() + 3 helpers (Section 6 Phase 6b)
      utils.R                                # classify_correlation_predictors() — shared by both check-*.R files
      full-subsets-gam.R                    # full_subsets_gam(), calls the snake_case names internally
                                             # (renamed from function_full_subsets_gam_v1.11.R, Section 6 Phase 9)
      functions_supporting.R                # wi(), extract_mod_dat(), build_inclusion_mat(), fit_mod_l()
                                             # (all exported but documented "not called directly")
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
    build_package_FFSgam.R    # legacy dev script with a hardcoded Windows path; .Rbuildignore'd,
                              # not part of the package — do not rely on it
    prompts/                  # session prompt logs (Section 7); .Rbuildignore'd

------------------------------------------------------------------------

## 5. Known Constraints or Decisions

- **[`generate.model.set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate.model.set.md),
  [`fit.model.set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit.model.set.md),
  and
  [`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.md)
  are permanent, deprecated aliases**, not just historical names to
  preserve. They live in `R/deprecated.R`, call
  [`.Deprecated()`](https://rdrr.io/r/base/Deprecated.html), and forward
  `...` to
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate_model_set.md)
  /
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.md)
  /
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full_subsets_gam.md).
  The first two are cited in the published paper (Fisher et al. 2018,
  Ecology and Evolution);
  [`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.md)
  is used directly in the companion repo’s case-study-2.Rmd (verified by
  reading that repo’s vignettes during the Phase 9 rename — see Section
  6). Do not remove any of the three, and do not let
  [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full_subsets_gam.md),
  [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate_model_set.md),
  or
  [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.md)
  call through their own deprecated wrappers internally (they must call
  the snake_case names directly, or every call to them emits spurious
  deprecation warnings).

- **`check.correlations()`, `check.non.linear.correlations()`,
  `build.inclusion.mat()`, `extract.mod.dat()`, and `fit.mod.l()` were
  renamed to snake_case with no backward-compatible alias** (Section 6
  Phase 9) — a deliberate, asymmetric decision relative to the three
  aliases above, not an inconsistency. Before renaming, every vignette
  in the companion repo (`case-study-1/2/3.Rmd`, `extra-examples.Rmd`,
  `faq.Rmd`) was checked for direct calls to these five names; none were
  found, unlike
  [`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.md).
  If you ever find one of these five names called directly in the
  companion repo (or anywhere else load-bearing), that’s new information
  that would justify adding a deprecated alias retroactively — don’t
  assume the original check still holds without re-verifying.

- **Do not change the GitHub remote URL** for this package repo. Users
  install via `devtools::install_github("beckyfisher/FSSgam_package")`
  and this must continue to work.

- **The publication repo (`beckyfisher/FSSgam`) is strictly read-only**
  for the purposes of this work. Do not commit to it or suggest changes
  to it. Its URLs are permanent references in the published paper
  (Fisher et al. 2018).

- **Vignettes live in the publication repo, not here.** Do not create a
  `vignettes/` folder in this package repo. Full worked examples and
  case studies are permanently hosted at
  <https://github.com/beckyfisher/FSSgam>. The `@examples` blocks in
  roxygen2 should be minimal, self-contained, and fast — use the bundled
  datasets (`case_study1`, `case_study2`, etc.). In practice every
  current example runs in well under 1 second, so none need
  `\donttest{}`; re-time new examples before assuming that’s still true.

- **The pkgdown site here is Reference-only — no Articles/vignettes
  section.** Full worked examples and case studies stay at the companion
  publication repo; `_pkgdown.yml`’s `navbar.structure.left`
  deliberately omits `articles`. The `reference:` index is curated into
  named groups (Core functions, All-in-one wrapper, Correlation
  diagnostics, Supporting functions, Deprecated aliases, Datasets) —
  when adding a new exported function, add it to an existing group or a
  new one, don’t leave it unlisted
  ([`pkgdown::check_pkgdown()`](https://pkgdown.r-lib.org/reference/check_pkgdown.html)
  will fail with “must be a known topic name” style errors if
  `_pkgdown.yml` references a function that no longer exists, or “topics
  not in any index” if a new exported function is missing from every
  group).

- **`master` and `dev` share the same `DESCRIPTION` `Version`.** This
  was not always true. Originally `dev` ran `1.0.0.9000` so pkgdown’s
  `development: mode: auto` could resolve `master` to `release` and
  `dev` to `devel` purely from the Version string (Phase 8) — but that
  meant every ordinary merge/PR between the branches fought the
  deliberate divergence. It broke for real once (PR \#2, 2026-06-22:
  merging `dev` into `master` dragged `1.0.0.9000` across, fixed by hand
  afterward) before being designed out altogether. Phase 11 moved the
  release/devel decision into `.github/workflows/pkgdown.yaml` (keyed
  off the branch name via the `PKGDOWN_MODE` env var, passed into a
  manual `as_pkgdown()`/`build_site()` call rather than
  `build_site_github_pages()`’s own override forwarding — see the Phase
  11 caution note for why that distinction matters), so `_pkgdown.yml`’s
  `development.mode` is now a plain static `release` default and
  `DESCRIPTION`’s `Version` no longer needs to differ between branches.
  Normal merges/PRs between `master` and `dev` are safe again — both
  deploy to the same `gh-pages` branch non-destructively (the deploy
  action uses `clean: false`).

- **Bundled datasets are managed via `data-raw/case_study_datasets.R`.**
  Do not regenerate or modify the `.rda` files in `data/` manually. If
  new data is needed, update `data-raw/case_study_datasets.R` and re-run
  it. Each dataset must have a `@format` roxygen2 entry in `R/data.R`
  for CRAN compliance.

- **`s()`, `te()`, and `t2()` must stay unqualified inside any
  `gam()`/`gamm()` formula** — `mgcv::s(...)` inside a formula breaks
  [`model.frame()`](https://rdrr.io/r/stats/model.frame.html)
  construction silently (caught by
  [`try()`](https://rdrr.io/r/base/try.html), surfaces as an unexplained
  `NA`). mgcv resolves these by literal symbol name during formula
  parsing, not via normal namespaced function dispatch. This bit us once
  during the 1.0.0 refactor (see the comment in
  `function_check_non_linear_correlations_v1.00.R`) — when adding
  `package::function()` namespacing elsewhere, never apply it to a
  smooth constructor written inside a formula, even though it’s safe
  (and required for `R CMD check`) everywhere else, including the outer
  `gam()`/`gamm()` call itself.

- **CRAN submission is the end goal.** Avoid patterns that generate
  NOTES in `R CMD check`: use `package::function()` for all non-base
  calls in `R/` (except smooth constructors inside formulas, above),
  replace `T`/`F` with `TRUE`/`FALSE`, prefer `inherits(x, "y")` over
  `class(x) == "y"`, and ensure all exported functions have complete
  roxygen2 documentation with runnable `@examples`.

- **Apache License 2.0** — retain existing licence; do not change it.
  The DESCRIPTION string `Apache License (== 2.0)` is already
  CRAN-canonical (verified via `tools:::analyze_license()`) and needs no
  accompanying `LICENSE` file.

- **The FAQ for this package lives at**
  `https://github.com/beckyfisher/FSSgam/blob/master/vignettes/faq.Rmd`
  (lowercase `faq.Rmd`, inside `vignettes/`, not at the repo root) and
  is the best source of truth for understanding function behaviour and
  edge cases. Read it before writing tests or `@examples`. The case
  study vignettes in the same `vignettes/` folder
  (e.g. `case-study-1.Rmd`) are the best source of real, working
  [`generate.model.set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate.model.set.md)/[`fit.model.set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit.model.set.md)
  call patterns to adapt for examples and tests.

------------------------------------------------------------------------

## 6. Refactor Phases

### Phases 1–5 — Completed (v1.0.0 modernisation)

All of the following landed in the v1.0.0 modernisation and are now the
baseline state of the repo, not aspirational:

1.  **DESCRIPTION and roxygen2** — `Authors@R`, explicit `Imports`,
    version 1.0.0, `URL`/`BugReports`, runnable `@examples` on every
    exported function, `NEWS.md`.
2.  **snake_case + deprecation wrappers** —
    [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate_model_set.md)
    /
    [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.md)
    are the real implementations; dot-case names are
    [`.Deprecated()`](https://rdrr.io/r/base/Deprecated.html) wrappers
    in `R/deprecated.R`.
    [`require()`](https://rdrr.io/r/base/library.html)/bare calls
    replaced with explicit `package::function()` (see the
    `s()`/`te()`/`t2()` exception in Section 5); `T`/`F` →
    `TRUE`/`FALSE`; `class(x) == "y"` → `inherits(x, "y")`.
3.  **Tests** — `tests/testthat/`, covering Gaussian and Tweedie model
    sets, the non-data.frame and NA-predictor error paths, and
    deprecation behaviour (including that
    [`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.md)
    itself stays warning-free).
4.  **GitHub Actions + badges** — `.github/workflows/R-CMD-check.yaml`
    (Ubuntu/macOS/Windows release + Ubuntu devel), README badges.
5.  **CRAN prep** — `cran-comments.md`, confirmed `devtools::check()` is
    clean (0 errors/warnings/notes) and all examples run fast.

Do not re-do this work; build on top of it. If you find a regression in
any of the above, fix it in place rather than reverting to the pre-1.0
style.

### Phase 6 — Decompose the monolithic function bodies — Completed

Done.
[`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate_model_set.md)
now lives in `R/generate-model-set.R`, split into eight unexported
helpers (`validate_use_dat()`, `build_null_model()`,
`check_predictor_missingness()`, `resolve_factor_interactions()`,
`resolve_smooth_smooth_interactions()`,
`build_predictor_correlation_matrix()`, `enumerate_candidate_models()`,
`build_model_formulas()`).
[`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.md)
now lives in `R/fit-model-set.R`, split into four
(`fit_and_summarise_saved_models()` and
`fit_and_summarise_unsaved_models()` are kept separate rather than
unified — `save.model.fits=FALSE` exists specifically so fitted model
objects are never all held in memory at once, and a shared helper would
have undone that — plus `compute_model_weights()` and
`compute_variable_importance()`). `function_generate_model_set.R` and
`function_fit_model_set.R` no longer exist.

Public arguments (`use.dat`, `test.fit`, `pred.vars.cont`, etc.) were
not touched — only internal structure moved. No `R/utils.R` was created:
after decomposing, nothing unexported turned out to be genuinely shared
across files — everything outside these two files is either a deprecated
wrapper or an already-exported, independently-documented function
(`check.correlations()`,
[`wi()`](https://beckyfisher.github.io/FSSgam_package/reference/wi.md),
etc.) that already has a sensible home.

Verified with a golden-master snapshot (14
[`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate_model_set.md)
scenarios covering every branch, including ones not in the committed
test suite — list-form `factor.smooth.interactions`, `cyclic.vars`,
`linear.vars`, a user-supplied `cor.matrix`; plus 8
[`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.md)-specific
scenarios: `save.model.fits=FALSE`, `report.unique.r2`, `VI.mods='all'`,
the `max.models` override, `parallel=TRUE`) compared before/after each
extraction, on top of the full testthat suite and `devtools::check()`.
If you decompose anything further, use the same approach: a scratch
before/after comparison script for branches the committed tests don’t
reach, not just the committed suite alone.

### Phase 7 — Broader test coverage — Completed

Done. `tests/testthat/` now also covers: `check.correlations()` and
`check.non.linear.correlations()` (basic matrices + invalid-column-class
errors);
[`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate_model_set.md)’s
factor-factor-interaction, smooth-smooth-interaction (both the
included-`te()` and excluded-by-`cov.cutoff` cases), and
non-linear-correlation code paths (previously only checked ad hoc, never
committed as tests);
[`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.md)
including its `factor.interactions`/`smooth.interactions`/`size` legacy
arguments; and the `functions_supporting.R` helpers
([`wi()`](https://beckyfisher.github.io/FSSgam_package/reference/wi.md),
`extract.mod.dat()`, `build.inclusion.mat()`, `fit.mod.l()`).

Writing these tests surfaced two genuine pre-existing bugs in
[`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.md)
(both predated the 1.0.0 rename, now fixed with regression tests in
`test-full_subsets_gam.R`, see `NEWS.md`): - the deprecated `size`
argument was a no-op (it warned but never actually fed its value into
`max.predictors`); - `used.data` and `predictor.correlations` in the
return value were always `NULL` (referenced
`model.set$use.dat`/`$cor.matrix`, which don’t exist — the real fields
are `$used.data`/`$predictor.correlations`).

If you find another function whose documented behaviour doesn’t match
what it actually does while working in this codebase, fix it in its own
commit with a regression test, the same way — don’t fold it into an
unrelated change.

### Phase 6b — Decompose check.correlations()/check.non.linear.correlations() — Completed

Done. `check.correlations()` moved to `R/check-correlations.R`, split
into `build_continuous_correlation_matrix()`,
`build_factor_continuous_skeleton()`, and
`fill_factor_factor_correlations()`. `check.non.linear.correlations()`
moved to `R/check-non-linear-correlations.R`, split into
`build_correlation_pair_grid()`, `estimate_non_linear_correlation()`,
and `assemble_non_linear_correlation_matrix()`.

Both functions had byte-identical column-validation/classification logic
(14 lines). That’s the genuinely shared utility Phase 6 didn’t find —
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
`test()` in another shell while the capture script’s own
`devtools::load_all()` was still in flight) — likely resource
contention, not a real bug. Re-run in isolation with `timeout`, it
completed in under 2 seconds with correct values. If a parallel/cluster
test seems to hang, check whether something else is concurrently
touching the package before assuming the code is broken.

### Phase 8 — Merge dev into master; pkgdown reference site with dev/release split — Completed

Done. `master` was 18 commits behind `dev` and still had the pre-1.0.0
codebase (`Version: 1.11`, no tests, no GitHub Actions) — none of Phases
1–7/6b had been merged there. Verified the merge was conflict-free first
(master’s only unique commit, a `rework` branch PR merge, had an empty
diff against the common ancestor) before merging `dev` into `master` and
pushing.

Added a pkgdown site (Reference only, no Articles — see Section 5) with
a GitHub Actions workflow (`.github/workflows/pkgdown.yaml`, the
standard r-lib/actions template, deploying to `gh-pages`). Set up a real
dev/release docs split using pkgdown’s `development: mode: auto`: bumped
`dev`’s `DESCRIPTION` `Version` to `1.0.0.9000` via
`usethis::use_dev_version()` (also added the NEWS.md dev heading), left
`master` at the release `1.0.0`, and added the same
`development: mode: auto` line to `_pkgdown.yml` on both branches
individually (not via merge — see Section 5 for why). Verified the mode
resolution directly via `pkgdown:::dev_mode()`: `master` resolves to
`release`, `dev` to `devel`.

**Discovery during this phase:** the live site at
`beckyfisher.github.io/FSSgam_package` was GitHub’s default Jekyll
rendering of `master`’s README, not pkgdown’s output — confirmed by
diffing the rendered HTML against the raw README and checking the GitHub
Pages API. The `gh-pages` branch already had the correct pkgdown build
sitting in it; GitHub Pages just wasn’t pointed at it. If pkgdown docs
ever look wrong again, check what branch GitHub Pages settings actually
point at before assuming the build is broken.

### Phase 9 — Complete the snake_case rename — Completed

Done.
[`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.md)
renamed to
[`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full_subsets_gam.md)
(file renamed `function_full_subsets_gam_v1.11.R` → `full-subsets-gam.R`
to match the kebab-case convention from Phase 6), with
[`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.md)
added to `R/deprecated.R` as a third permanent alias.
`check.correlations()`, `check.non.linear.correlations()`,
`build.inclusion.mat()`, `extract.mod.dat()`, and `fit.mod.l()` renamed
outright with no alias — see Section 5 for the verification behind that
asymmetry and the decision to treat it as deliberate rather than
inconsistent.

All internal call-sites updated (`R/generate-model-set.R`,
`R/fit-model-set.R`, `R/utils.R`, `R/FSSgam.R` package doc), all
affected tests updated to call the new names, and `test-deprecated.R`
extended with a
[`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.md)
deprecation test matching the existing
[`generate.model.set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate.model.set.md)/[`fit.model.set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit.model.set.md)
ones. `_pkgdown.yml`’s reference index updated to match.

Verified with the full `testthat` suite (32 files, 0 failures) and a
clean `devtools::check()` (0 errors/warnings/notes) — no golden-master
snapshot needed, since this was a pure rename with no logic changes,
unlike the Phase 6/6b decompositions.

### Phase 10 — Suggested next priorities

With Phases 1–9 and 6b complete, candidates for what comes next (none
started): - Tag a release and actually submit to CRAN (or decide what’s
still blocking that). - The companion docs repo (`beckyfisher/FSSgam`)
still calls the deprecated dot-case names in its vignettes — see
`FSSgam-docs-CLAUDE.md` drafted for that repo (a copy may already be in
place there as `CLAUDE.md`). Now that
[`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.md)
is also a deprecated alias rather than the real implementation, the same
applies to it as to
[`generate.model.set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate.model.set.md)/
[`fit.model.set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit.model.set.md).

### Phase 11 — Fix `family = nb()` candidate-refit bug (#12); replace the

Version-based pkgdown split with a CI branch override — Completed

Done.
[`fit_mod_l()`](https://beckyfisher.github.io/FSSgam_package/reference/fit_mod_l.md)
(`R/functions_supporting.R`) was passing
`family = stats::family(test.fit.)` into
[`update()`](https://rdrr.io/r/stats/update.html) for every candidate
refit — the *already-fitted* family object. mgcv’s extended families
that estimate an extra parameter (`nb()`, `tw()`, …) store that estimate
in the family object’s own mutable environment, so every refit was
warm-started from `test.fit`’s unrelated estimate, destabilising most
(but not all) fits. Fix: drop the explicit `family=` override so
[`update()`](https://rdrr.io/r/stats/update.html) re-evaluates a fresh
family from `test.fit`’s own call, the same approach
`build_null_model()` (`R/generate-model-set.R`) already used for the
null model — which is why the null model alone always survived. Added a
regression test mirroring the existing Tweedie test in
`test-fit_model_set.R`, using `family = nb()` on `case_study1`. Verified
12/13 candidate models failing before the fix, 0/13 after, sequentially
and under `parallel = TRUE`.

Also replaced the Phase 8 Version-based pkgdown dev/release split (see
Section 5) with a CI branch override, after it caused a second incident
in one day (see the PR \#2 note this section used to carry).
`_pkgdown.yml`’s `development.mode` is now a static `release` default;
`.github/workflows/pkgdown.yaml` computes a `PKGDOWN_MODE` env var from
`github.ref_name`/`github.base_ref` (`dev` → `devel`, else `release`).

**Caution for next time:** the obvious implementation —
`pkgdown::build_site_github_pages(override = list(development = list(mode = ...)))`
— silently does *not* work. That function calls
`as_pkgdown(pkg, override = list(destination = dest_dir))` on a plain
path first (computing `dst_path`/ `development$in_dev` from whatever
mode was active at that point), then forwards any further `override=`
into `build_site()`’s *own* `as_pkgdown(pkg, override = override)` call
– but by then `pkg` is already a `pkgdown` object, and `as_pkgdown()`
short-circuits for that case
(`if (is_pkgdown(pkg)) { pkg$meta <- modify_list(pkg$meta, override); return(pkg) }`),
patching `$meta$development$mode` without recomputing
`dst_path`/`in_dev`. The override looked like it worked
(`pkg$development$mode` showed the new value) while devel docs silently
kept building to `docs/` instead of `docs/dev/`. The workflow instead
calls
`pkgdown::as_pkgdown(".", override = list(destination = "docs", development = list(mode = ...)))`
directly – passing both overrides in the *same* call, while `pkg` is
still a plain path string – then
[`pkgdown::clean_site()`](https://pkgdown.r-lib.org/reference/clean.html),
[`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html),
and the unexported `pkgdown:::build_github_pages()` (the same three
calls `build_site_github_pages()` makes internally). Verified end-to-end
locally for both modes: `release` built to `docs/index.html`, `devel` to
`docs/dev/index.html`, non-destructively (matching the deploy action’s
`clean: false`). If you ever touch this step again, re-verify `dst_path`
and `development$in_dev` directly rather than trusting that
`development$mode` alone reflects what got built.

`dev`’s `DESCRIPTION` `Version` was reset from `1.0.0.9000` back to
`1.0.0` to match `master`; going forward the two branches keep the same
`Version` and merge/PR normally.

------------------------------------------------------------------------

## 7. Prompt Log

Session logs for this project are in `prompts/`. Create the folder if it
does not exist. Use a short kebab-case descriptor as the filename for
each session (e.g. `cran-modernisation.md`). If a file with that name
already exists, append to it.

Log format for each session entry:

    ## Session: <descriptor>
    Date: <YYYY-MM-DD>
    Model: <model name and version>

    ### Prompts and Responses

    **User:** <prompt text>

    **Claude:** <summary of response — full code blocks where relevant, prose summarised>

    ---

**Note:** the `prompts/` folder does not exist yet — the v1.0.0
modernisation session was not logged here. Start the folder and log file
on the next session that touches this repo.
