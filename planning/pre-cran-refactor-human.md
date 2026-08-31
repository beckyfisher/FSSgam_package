# Pre-CRAN refactor — plan

This document is for collaborators. It states what the refactor does, why, in
what order, and what changes for someone using the package. The companion
document `pre-cran-refactor-claude.md` holds the implementation
specification: file and line references, replacement code, edge cases, tests,
and the alternatives that were rejected. Every decision is stated in full in
the specification; this document gives the decision and a one-sentence reason.

---

## 1. What this refactor decides

It decides that the package's known defects, its unvalidated arguments, its
duplicated interaction-resolution code, and its six open issues are all fixed
before the first CRAN submission rather than after it.

The work was scoped from a full read of the 1,947 lines in `R/`, carried out on
2026-08-31 against commit `a8235e2`. Twenty-two items were identified. Four are
defects that silently change or abort results; four are arguments that accept
invalid values without complaint; the remainder are duplicated code, unsafe
idioms, and documentation errors. Six issues were already open on the
repository. All of them are in scope.

## 2. Why now rather than after release

`R CMD check` is clean. It was run twice on 2026-08-31 against commit `a8235e2`
on this WSL host under R 4.6.1: once through `devtools::check()`, and once as
`R CMD check --as-cran` with the local CRAN incoming checks enabled
(`_R_CHECK_CRAN_INCOMING_=TRUE`, remote checks disabled as they require network
access). Both returned 0 errors, 0 warnings and 0 notes. The `testthat` suite
is green at 438 passing expectations, 0 failures, 11 skipped.

Nothing in this plan is therefore required to pass the check. That is the
argument for doing it, not against it. A clean check confirms the package
builds, installs, documents itself, and runs its examples. It does not call
`generate_model_set()` with two predictors whose names overlap, or with a
single predictor, or with `max.predictors = 1`, or with an invalid `r2.type`.
Each of those is broken today and each passed the check.

The second reason is that release freezes the public interface. Two of the
changes below alter what users see in their output, and one removes an exported
function. After a version is on CRAN these become breaking changes to a
published package, which requires a deprecation cycle. Before it, they are
simply how the package works.

## 3. What changes for someone using the package

Four changes are visible to users. The rest are internal.

**Candidate model names change for anyone not working in the C locale.** Model
names are built by sorting the term names, and sorting currently follows the
session's collation setting. A model that is named `complexity+ZONE` in an
`en_US.UTF-8` session is named `ZONE+complexity` in a C-locale session. After
this change it is `ZONE+complexity` everywhere. The fitted results are
identical; the name and the row order of `mod.data.out` are not. Anyone holding
a saved model table and matching on `modname` will need to regenerate it.

**A supplied `cor.matrix` now governs every stage that uses correlations.**
There are three such stages: which factor-factor interaction columns are built,
which `te()` smooth-smooth interaction terms are built, and which assembled
candidates survive the collinearity screen. At present a supplied matrix governs
only the third; the other two recompute correlations from the data and ignore
it. A user who supplies a matrix saying two predictors are uncorrelated
currently still finds their `te()` term absent. After this change the supplied
matrix decides all three. Users who do not supply a `cor.matrix` see no change,
and that is verified rather than assumed (specification sections 5.3 and 11).

A supplied `cor.matrix` that omits one of the named predictors will now be
rejected with a message naming what is missing, rather than failing later with a
subscript error.

**Invalid argument values now produce errors.** `r2.type = "rsq"` currently
returns a column of `NA` with no message. `VI.mods = "alll"` currently fails
with `object 'aic.var.weights' not found`. Both will name the argument and list
the values it accepts.

**`fit_mod_l()` is no longer exported.** It is an internal helper documented as
"not called directly", and its arguments changed as recently as Phase 12. It
remains available as `FSSgam:::fit_mod_l()`. `wi()`, `extract_mod_dat()` and
`build_inclusion_mat()` stay exported, being more plausible to call directly.

One new argument is added: `progress`, on `fit_model_set()` and
`full_subsets_gam()`, defaulting to `interactive()`. This turns off the
progress bar in scripts, reports and checks while keeping it at the console.
`VI.mods` is added to `full_subsets_gam()`, which currently cannot reach it.

## 4. The four defects

Each was reproduced on this host before being written down. The reproductions
are in specification section 3.

**A1 — a cyclic variable captures every predictor whose name contains it.**
Cyclic variables are matched into formula terms with an unanchored `grep()`, so
declaring `depth` cyclic also makes `depthx` cyclic. The name is used as a
regular expression, so a predictor containing a full stop matches any
character. The consequence is a model fitted with the wrong smoothing basis,
with no error and a candidate name that looks correct. This is the most serious
item in the plan.

**A2 — `fit_model_set()` fails on a single-predictor model set.** The variable
importance calculation indexes a data frame in a way that collapses it to a
vector when one predictor is present, and `colSums()` then errors.
`generate_model_set()` was fixed for this case during Phase 13;
`fit_model_set()` was not.

**A3 — `max.predictors = 1` builds interactions it should not.** A loop written
as `for(i in 2:n)` counts backwards when `n` is 1, so it also runs at `i = 1`.
With two factors and `max.predictors = 1`, `used.data` gains duplicate columns
and four spurious warnings are emitted.

**A4 — the large-model-set path loses its error handling in parallel.** When
model fits are not being saved and `parallel = TRUE`, the line that tells
`foreach()` to record a failure rather than abort is commented out. One
unfittable candidate therefore aborts the whole run instead of appearing in
`failed.models`. The other fitting path does this correctly.

## 5. The restructure

Its purpose is to remove the duplication that produced A3, and to make the two
functions that build interaction terms readable enough that the next defect is
visible on inspection. It changes no behaviour.

`resolve_factor_interactions()` is 161 lines and does four jobs. Two of them are
implemented twice, once for when the user passes `TRUE` and once for when they
pass a character vector, and the two copies have drifted: one screens both
halves of the correlation matrix and the other only the upper half; one warns
when nothing survives and the other is silent.
`resolve_smooth_smooth_interactions()` is 73 lines with the same two-copy
structure, and both copies build a data frame of pasted columns that is then
discarded without being read.

The restructure normalises each argument to one internal form first, then runs
one implementation. The repeated "enumerate every combination of these
variables and drop the ones that are too correlated" block, which currently
appears three times, becomes a single helper. Full design in specification
section 5.

It is verified the way Phases 6 and 6b were: a snapshot of outputs across every
argument combination, taken before and compared after, including combinations
the committed test suite does not reach.

## 6. Order of work

Seven phases. The order is chosen so that each phase is verifiable on its own,
and so the restructure lands on behaviour that has stopped moving.

| Phase | Purpose | Done when |
|---|---|---|
| 1 | Fix the defects and validation outside the interaction code: A2, A4, B1–B4, and issues #7 and #9 | Each has a regression test; suite green |
| 2 | Fix the correlation functions: issues #10 and #12, and the redundant refits | The two tests written to tolerate the defects become exact |
| 3 | Restructure interaction resolution, changing no behaviour: the duplicated branches and the dead code | Golden-master snapshot identical before and after |
| 4 | Apply the behaviour changes on the restructured code: A1, A3, and issues #8 and #13 | Numerical snapshots updated deliberately, with the change recorded in `NEWS.md` |
| 5 | Replace the unsafe idioms | Suite green; no output change |
| 6 | Unexport `fit_mod_l()` | `NAMESPACE`, `_pkgdown.yml` and tests updated |
| 7 | Documentation, spelling and DESCRIPTION | `R CMD check --as-cran` clean |

Phases 1 and 2 are independent of each other and of everything after them.
Phase 3 must precede phase 4, so that the behaviour changes are made in code
that has already been simplified and snapshot-verified. Phases 5 to 7 can
follow in any order.

## 7. Decisions taken

Five decisions were put to the maintainer on 2026-08-31 and answered. Each is
stated in full in the specification section named.

- **A supplied `cor.matrix` governs both the `te()` terms and the candidate
  screen**, because that is what `?generate_model_set` already documents.
  Specification 6.4.
- **It governs the factor-factor interaction columns as well**, so that one
  argument governs all three stages rather than two of three. Specification 11.
- **Model names are sorted in byte order** via `sort(method = "radix")`, because
  that is the smallest change that makes a saved analysis reproducible across
  machines. Specification 6.2.
- **`smooth.smooth.interactions` keeps its current behaviour** — `TRUE` builds
  pairwise `te()` terms, a character vector builds terms up to `max.predictors`
  — and the documentation is corrected to describe that, because the
  alternative removed a working capability nobody had asked to lose.
  Specification 6.3.
- **`fit_mod_l()` is unexported and the other three helpers stay**, because it
  is the one whose signature is most likely to move again. Specification 8.

The fourth of these moves an item out of phase 4 and into phase 7: it is now a
documentation correction, not a behaviour change. No other phase boundary moves.

## 8. Nothing outstanding

The two questions this document previously carried have been answered and appear
in section 7. Implementation can proceed through all seven phases without a
further decision.

Two points will nonetheless need a judgement while phase 3 is underway, and are
called out here so they are not settled silently. The specification records the
intended answer for each, but each rests on reading code rather than on a
measurement, so either could turn out otherwise once the restructure exposes it.

The first is the pair of divergences between the duplicated interaction-building
branches, described in section 5: one screens both halves of the correlation
matrix and the other only the upper half, and one warns when nothing survives
while the other is silent. The specification holds both at their current
behaviour through phase 3, behind an argument, and resolves them in phase 4 in
favour of the stricter and noisier form. If the golden master shows that
resolution changes more than the two scenarios expected, it should be raised
rather than absorbed.

The second is whether the `cor.matrix` completeness check described in section 3
belongs in `generate_model_set()` alone or also in `check_correlations()`. The
specification places it in `generate_model_set()` only, because
`check_correlations()` is exported and adding validation to an exported function
makes it stricter than it is today. That reasoning is the same as the one
applied to `extract_mod_dat()` in specification 3.3, and would only need
revisiting if a caller is found that relies on the looser behaviour.
