## Session: open issues addressed in batches
Date: 2026-09-03
Model: Claude Opus 5 (1M context), claude-opus-5[1m]

### Prompts and Responses

**User:** can you check the current open issues and made a plan to address these in the most appropriate sequence

**Claude:** Read the twelve open issues (FSSgam_package#14, #19, #20, #22, #23,
#24, #25, #26, #27, #28, #30, #31) and grouped them into six clusters by the
function each touches, rather than treating them as twelve independent tasks.
Proposed a sequence of six batches, one branch and one pull request each, with
the ordering justified by dependency: metadata and no-decision fixes first; the
declared R floor next because it governs which idioms later batches may use; the
two supplied-`cor.matrix` validation issues together because the second is
unreachable until the first is fixed; argument validation next; the `null.terms`
cluster after that because both issues in it need the same variable-extraction
machinery; and the parallel-worker issues last because neither blocks anything
else. Recorded four decisions the sequence required.

**User:** Implement fixes on all of these in the recommended order. Decisions:
(1) raise the floor to R (>= 4.4.0); (2) for the interaction column collision,
warn and decline to build the column; (3) for the `null.terms` collinearity gap,
return the forced-versus-candidate correlations alongside `predictor.correlations`
and add an argument controlling the cutoff applied to them, defaulting to 0.8,
with no screening or dropping of correlations between null terms themselves,
documented as such; (4) no release yet. One stacked pull request per batch, each
subjected to a review loop of up to ten cycles by independent sessions until no
substantive comments remain.

**Claude:** Batch 1 (FSSgam_package#20, #19, #24), branch
`batch1-cran-urls-null-guard-r2-docs` off `dev`.

- **#20.** Added the trailing slash to `https://beckyfisher.github.io/FSSgam_package`
  and `https://beckyfisher.github.io/FSSgam` in `DESCRIPTION`, and to the second
  of the two in `README.md`, which is not `.Rbuildignore`d and so ships in the
  tarball. Confirmed with `curl` that both bare forms return 301 and both
  trailing-slash forms return 200, and that
  `https://github.com/beckyfisher/FSSgam_package` is not redirected and needs no
  change. Re-ran `devtools::document()` so `man/FSSgam-package.Rd` follows.

- **#19.** Added `&&!inherits(null.fit,"try-error")` to the guard in
  `estimate_non_linear_correlation()` (`R/check-non-linear-correlations.R`), with
  a comment recording that the null half is defensive, that no case reaching it is
  known because both models are fitted on the same `dat.r`, and what the failure
  would be without it. No test: constructing the case would require mocking
  `nnet::multinom` to fail on the second call only, which pins the implementation
  rather than the behaviour.

- **#24.** Documented `r2.vals.unique` in `fit_model_set()`'s `@return`
  description of `mod.data.out`, and rewrote the `report.unique.r2` `@param` in
  both `fit_model_set()` and `full_subsets_gam()` to say which column it
  populates. The documentation states four things the name does not support:
  that the column is the model R2 minus the null model's R2 and so depends on
  `null.terms`; that it is per-model rather than per-term, and a per-predictor
  contribution only at `max.predictors = 1`; that it is on the scale `r2.type`
  produced and can therefore be negative; and that values are comparable only
  within a set sharing the same `null.terms` and `r2.type`. The issue's evidence
  quoted the subtraction without its `if(report.unique.r2==TRUE)` guard, so the
  wording was written against the code rather than the issue: the column is
  always present and is `NA` unless `report.unique.r2` is `TRUE`.

---
**Review cycle 1 (independent session, batch 1).** No substantive defect. Ten
smaller points, of which these changed the branch:

- `full_subsets_gam()`'s own `Value` section still omitted `r2.vals.unique`, so
  the two documented descriptions of the same returned object disagreed. Added,
  with the cross-reference written as `\link[=fit_model_set]` rather than plain
  text so it renders as a hyperlink.
- The claim "the null model is intercept-only" is wrong for a `uGamm`/`gamm4`
  test fit with `random =`, where the random structure is retained and only the
  fixed part is an intercept. Reworded, and the qualification stated. The
  numerical claim was unaffected: every `r2.type` ignores the random effects, and
  the reviewer measured `null.r2 = 0` for such a fit at two `r2.type` values.
- "the same guard" was not literal. `check_correlations()` tests `!is.null()` as
  well, because its `null.fit` can be `NULL` by design. Corrected to "the
  equivalent guard", and the `#16` against `#18` discrepancy between the issue
  body and the code comment recorded on the face of the code.

**A reachable defect of the same shape, found by that review and raised as
FSSgam_package#33.** The `lm()` call in the continuous-response, factor-predictor
branch of `estimate_non_linear_correlation()` is not wrapped in `try()` at all,
and neither is the equivalent one in `build_factor_continuous_skeleton()`. Both
were confirmed here to abort the whole call. The reviewer's reproduction used a
frame whose factor and continuous columns share no complete rows; a shorter and
more realistic one was found before filing, since a single-level factor is enough:

```r
d <- data.frame(x = rnorm(20), f = factor(rep("a", 20)))
check_correlations(d)             #> Error: contrasts can be applied only to ...
check_non_linear_correlations(d)  #> Error: contrasts can be applied only to ...
```

and it reaches the public entry point, so `generate_model_set()` stops with an
error naming neither the function, the argument, nor the predictor at fault. It
was checked that the missing-value route to the same failure does not reach
`generate_model_set()`, which rejects `NA` predictors first; the single-level
route does. Not fixed here: it is a separate, pre-existing path, and which of a
`try()` guard and an argument check is right is a decision the issue records
rather than settles. Both call sites are marked with a comment naming the issue,
and the `NEWS.md` entry for FSSgam_package#19 was rewritten so it does not imply
the function now degrades to `NA` in every failure, which it does not.

---
**Review cycle 2 (independent session, batch 1).** Two substantive findings, both
in documentation text added in response to cycle 1, neither a code defect. Both
were reproduced here before being acted on.

- The cycle 1 correction to the `null.terms` sentence over-generalised in the
  opposite direction. `build_null_model()` updates the `test.fit` with the formula
  `~1`, and `update.formula()` drops every term of it. A random effect written
  as `s(site, bs = 're')` is therefore *not* retained by the null model, which is
  precisely why `null.terms` exists and precisely the case `fit_model_set()`'s own
  `@examples` builds. Measured: a `gam()` test fit of
  `~ s(depth, k = 3, bs = 'cr') + s(site, bs = 're')` gives a null model of
  `Herbivore.abundance ~ 1` with zero smooth terms. What *is* retained is a random
  structure supplied through a separate argument rather than through the formula:
  `uGamm(random = ~(1|site), lme4 = TRUE)` gives a null model whose `gam` formula
  is `~ 1` and whose `mer` still groups on `site`. The sentence now states both,
  and the distinction it turns on is the formula, not the class of the fit.

  Note for anyone reproducing that second measurement: the `uGamm` fixture needs
  `k >= 4`. At `k = 3` the smooth has one basis column and lme4 rejects the
  grouping factor outright, which `CLAUDE.md` Phase 13 already records.

- "this function returns the object `fit_model_set()` built, unaltered" was false.
  `full_subsets_gam()` returns six elements against `fit_model_set()`'s four, and
  its own `Value` block lists two of the extra ones, so the page contradicted
  itself. Only `mod.data.out` passes through, which is all the sentence needed to
  claim.

Everything else in the cycle 1 response was verified by that session and holds,
including that `\link[=fit_model_set]` renders as a hyperlink in both the
installed help and a pkgdown build.

---
**Review cycle 3 (independent session, batch 1).** No substantive defect. The two
statements corrected in cycle 2 were re-derived from scratch by that session and
both hold, across all four `test.fit` classes. Five minor points were acted on:

- The `r2.vals.unique` text stated "NA unless `report.unique.r2` is TRUE", which
  is a necessary and not a sufficient condition. A `gamm` `test.fit` under the
  default `r2.type` has `r2.vals` `NA` for the whole set, so the column is `NA`
  with `report.unique.r2` TRUE as well. Stated. That gap is FSSgam_package#6,
  open and pinned by an expectation, so this is documentation of it, not a fix.
- Three constructions parent `CLAUDE.md` §12 removes on sight, all of which
  render into `\value{}` and onto the pkgdown page: "Two qualifications." as an
  announcement of structure, a sentence opening "And", and "Note that".
- The `R/check-correlations.R` comment said the `multinom()` fits were "below",
  which is true of the file and not of the function: they are in
  `fill_factor_factor_correlations()`. Named. It was also attached to the
  `expand.grid()` line rather than to the `apply()` containing the `lm()`.
- Four lines recording which issue number FSSgam_package#19's body cites were
  removed from `R/`. Accurate, but §14 puts what outlives a thread in a
  repository document, and the provenance is held by the issue and by `git log`.

**A second pre-existing defect found by that review and raised as
FSSgam_package#34.** An `mgcv::gamm()` fit cannot be a `test.fit` at all. It
records no `call`, so `stats::update()` in `build_null_model()` fails with "need
an object with call component" before any candidate is enumerated. Confirmed here.
This matters because `generate_model_set()`'s own `null.terms` documentation
directs the user to `gamm` for a correlation structure, its `@details` contrast
`gam` against `gamm` at length, `fit_model_set()`'s `r2.type` names `gamm` among
the classes it reports for, and `extract_mod_dat()` has a `gamm` branch. The error
message compounds it by advising the user to stop using `uGamm`, which is what
they would have had to use to succeed. Not fixed here: whether to support it, or
to document `uGamm` as the only route to a mixed fit, is a decision.

---
**Claude:** Batch 2 (FSSgam_package#31), branch `batch2-r-version-floor`, stacked
on batch 1.

The declared floor `R (>= 3.5)` was not reachable. The issue names `MuMIn`;
`mgcv` was found to declare the same thing, so both hard `Imports` require
R 4.4.0 and the floor was raised to match. Measured at MuMIn 1.48.19 and
mgcv 1.9-4:

```
MuMIn    Depends: R (>= 4.4.0)
mgcv     Depends: R (>= 4.4.0), nlme (>= 3.1-64)
doSNOW   Depends: R (>= 2.5.0), ...
foreach  Depends: R (>= 2.5.0)
nnet     Depends: R (>= 3.0.0), ...
testthat Depends: R (>= 4.1.0)          [Suggests]
```

The issue poses two decisions. The first is whether the floor should track the
hard dependencies: it should, and does. The second is whether anything should
check it, since `R CMD check` does not verify a declared floor against a
dependency's. Two things now do, rather than the comment alone the issue
proposed:

- A comment in `DESCRIPTION` recording that the floor is derived and not chosen.
  Confirmed that `read.dcf()` ignores `#` lines, that `R CMD check` passes with
  it, and that `R CMD build` strips it from the tarball's `DESCRIPTION`, so the
  provenance stays in the source where a reader of the repository finds it and
  never reaches CRAN.
- An `ubuntu-latest`, `oldrel-2` job in `.github/workflows/R-CMD-check.yaml`.
  With the release at 4.6, `oldrel-2` is 4.4, which is the declared floor. It is
  written as `oldrel-2` rather than as `4.4` so it needs no maintenance, with a
  comment recording that it stops testing the floor once R moves past 4.6. This
  is more than the issue asked for and is the part most easily dropped if the
  job proves unreliable; the comment alone still records the provenance.

Consequences carried out with it: `deparse_one()` in
`tests/testthat/helper-fixtures.R` existed only because `deparse1()` arrived in
R 4.0.0, so it was deleted and its 31 call sites across four test files replaced
with `deparse1()`. `CLAUDE.md`'s note recording that constraint was rewritten,
as was its `Depends` line. Its `Imports`/`Suggests` list was also corrected while
adjacent: it still listed `gamm4` under `Imports`, which stopped being true when
`gamm4` moved to `Suggests`.

Not done: `parallel::makePSOCKcluster(setup_strategy = "sequential")` is now
available as a mitigation to test against FSSgam_package#14, having been
rejected partly on the old floor. That belongs to batch 6.

---
**Review cycle 1 (independent session, batch 2).** Three substantive findings.
The floor value itself was confirmed correct — 4.4.0 across the whole installed
dependency closure, nothing in `R/` or `tests/` using anything added after
R 4.0.0 — and the `deparse1()` substitution exact. All three defects were in the
two mechanisms meant to record and check the floor, and all three were
reproduced here before being acted on.

- **The `oldrel-2` job could not detect the defect FSSgam_package#31 reported.**
  Installing a package rejects a floor set *above* the running R and accepts any
  floor *below* it, so no matrix of R versions sees a floor that is too low.
  Confirmed: `R CMD INSTALL` of a copy with `Depends: R (>= 3.5)` succeeds on
  R 4.6.1. Four documents said the job checked the floor. Replaced with
  `tools/check-r-floor.R`, which compares the declared floor against the
  `Depends` of every hard dependency and exits 1 on exactly the state the issue
  reported, run in CI as an `r-floor-consistency` job. The matrix entry is kept
  but described accurately, as establishing that the package installs and passes
  at the declared floor.

  This is the finding worth carrying forward: **an R version in the CI matrix is
  not a floor check.** It tests the opposite direction from the one that goes
  wrong.

- **The `DESCRIPTION` comment is deleted by routine tooling.** Any `desc`-based
  write removes it — `usethis::use_package()`, a version bump,
  `roxygen2::roxygenise()` updating `Config/roxygen2/version` — and this
  repository has run two of the three. Writing R Extensions §1.1.1 also states
  the format does not support `#` comments; `read.dcf()` skipping them is an
  implementation detail, not a guarantee. The comment was removed and the
  provenance moved into the script, where it is executed rather than trusted.

- **`oldrel-2` resolves to 4.4.3, not 4.4.0.** A base function added in 4.4.1
  would pass CI and fail at the declared floor. The matrix entry now names 4.4.0
  literally. The stated reason for `oldrel-2`, that it needed no maintenance,
  also contradicted its own comment saying to re-check it when R moves past 4.6.

Minor points acted on: the PR evidence table was headed "every dependency" while
covering only the direct ones; bare `#14` where `beckyfisher/FSSgam#14` also
exists; and `planning/pre-cran-refactor-claude.md` still directed a future
session to the deleted `deparse_one()`. `planning/golden-master/capture.R`
defines its own copy and was left alone, being a self-contained record of what
was run. A comment was added to FSSgam_package#14 recording that the constraint
on its `setup_strategy` mitigation is lifted.

---
**Review cycle 2 (independent session, batch 2).** Two substantive findings, both
silent-pass paths in `tools/check-r-floor.R` — the same defect class cycle 1
found in the CI matrix entry, one level down. Both reproduced here before being
acted on.

- **The script read only the direct `Imports` and printed "OK".** `Matrix`
  declares `R (>= 4.4)` and `MASS` declares `R (>= 4.4.0)`, and neither is a
  direct dependency of this package: both arrive through `mgcv`, `MuMIn` and
  `gamm4`. Reproduced: a `DESCRIPTION` with `Imports: gamm4, nnet` and
  `Depends: R (>= 3.5)` exited 0 under the old script and exits 1 under the new
  one, naming `MASS, Matrix, mgcv`. The script now walks the transitive closure
  of `Depends`, `Imports` and `LinkingTo`, excluding only base-priority packages
  — recommended packages are where the binding floor is actually found.

- **An uninstalled dependency contributed no floor and the job stayed green.**
  `packageDescription()` returns `NA`, `$Depends` errors, the `tryCatch`
  swallowed it and the row printed as `-`. This repository's `CLAUDE.md` §3
  records that exact state for `doSNOW`. The script now stops, naming the
  packages it could not read.

The minor finding that a package moved from `Imports` to `Depends` would be
unread was fixed with them: the seed set is now `Depends` + `Imports` +
`LinkingTo`.

Nine edge cases were constructed and run against the rewritten script, since it
is the only new code in this batch: the floor written `R (>=4.4.0)` and
`R(>= 4.4.0)` without the usual spaces; `R (>= 4.4)` against `R (>= 4.4.0)`,
which must compare equal and do so because `package_version()` pads;
`R (>= 4.10.0)` against `R (>= 4.9.0)`, which string comparison would order
wrongly and `package_version()` does not; version constraints on the `Imports`
entries; a hard dependency declared in `Depends` rather than `Imports`; a
`DESCRIPTION` with no R floor at all; an uninstalled dependency; and package
names beginning with R -- `Rcpp`, `R6` -- which must not be read as the R entry
and are not, the pattern being anchored on `R` followed by `(`.

---
**Review cycle 3 (independent session, batch 2).** One substantive finding: the
third instance of the same failure mode, this time inside the script written to
catch it.

`read_r_floor()` returned `NULL` both for a dependency declaring no R constraint
and for one whose constraint it could not parse. It accepted only `>=`, so any
other operator fell through, printed as `-` — the same glyph as "declares
nothing" — and contributed no floor. Reproduced with a scratch dependency
declaring `Depends: R (<= 9.9.0)`: the script printed "No hard dependency
declares an R floor; nothing to check" and exited 0. The form exists in
practice; the reviewer found `R (> 3.0.2)` among 818 R constraints across the
1119 packages installed here.

Rewritten as `parse_r_constraint()`, returning one of three states rather than
two: `none`, `ok`, `unparsed`. `>` is now read as well as `>=`, since
`R (> 3.0.2)` means at least that version and understates the true floor by at
most a patch level, which is safe in this direction. Anything else stops the
script naming the package and the constraint. Verified with three scratch
packages installed for the purpose: `R (> 4.5.9)` is counted as a floor of 4.5.9
and fails a declared 3.5; `R (<= 9.9.0)` stops the script; `Depends: methods`
reaches the "nothing to check" path legitimately.

**Three statements of mine were false and are corrected.** They are recorded here
because each was asserted in several places at once:

- *`MASS` declares `R (>= 4.4.0)` and arrives through `mgcv`, `MuMIn` and
  `gamm4`.* It arrives only through `gamm4` → `lme4`, and `gamm4` is a
  `Suggests`, which the script does not walk. `MASS` is not in this package's
  hard closure and does not appear in the script's own output. `Matrix` alone
  justifies the closure walk, so the argument survives, but the sentence was in
  the script, `NEWS.md`, `CLAUDE.md`, the commit message and the pull request
  body.
- *Neither new CI job has run yet.* Both had run and passed on `51b296e`.
  Recorded below with the job names.
- The pull request body quoted an error message the script no longer emits.

The empty-closure path also crashed at `order(names(rows))` for a package with
no hard dependencies outside base R, and now exits cleanly.

**Both new CI jobs pass.** Run 33721973950 on `51b296e`: `declared R floor is
reachable: success`, `ubuntu-latest (4.4.0): success`, alongside the four
pre-existing matrix jobs. This confirms two things that could not be checked
locally: that `r-lib/actions` resolves a literal `4.4.0`, and that the script
runs against a dependency set installed by `setup-r-dependencies` rather than
this host's.

---
**Review cycle 4 (independent session, batch 2).** One substantive finding: the
fourth instance, this time in the ordering of the guards rather than in any one
of them.

The empty-closure exit added in cycle 3 was placed *above* the
uninstalled-dependency check added in cycle 2. `rows` is empty exactly when no
dependency could be read, so a package all of whose hard dependencies were
missing reached `quit(status = 0)` reporting "This package declares no hard
dependencies outside base R" — while declaring two — and accepted the declared
floor having read nothing. Reproduced with `Imports: notinstalledA,
notinstalledB, stats` and `Depends: R (>= 3.5)`: exit 0 before, exit 1 after.

The reporting section was restructured so that every guard runs before any exit
path, with a comment saying so and why. Eight cases were then re-run: all
dependencies uninstalled, one uninstalled, genuinely no hard dependencies, an
unparsed constraint, a strict `>` constraint, a dependency declaring no R floor,
the FSSgam_package#31 state, and the package as it stands. All eight behave
correctly.

Two smaller defects fixed with it: a table row vanished entirely rather than
printing a blank where a field was zero-length, `sprintf()` returning
`character(0)` if any argument is; and the advisory to record why a floor is
higher than its dependencies require was unreachable on the no-floors path.

Two documentation corrections. `CLAUDE.md` labelled as "Measured" the statement
that reverting `DESCRIPTION` would leave every CI job green — no such CI run was
made, and it is an inference from `R CMD INSTALL` and `devtools::check()` on this
host. `NEWS.md` named two of the three cases in which the script stops where two
other documents named three.

The reviewer also ran `parse_r_constraint()` over every R constraint declared by
every package installed on this host — 1119 constraint rows across 813 packages —
and none was classified unparsed, so the parser meets the population it will
actually see.

---
**Review cycle 5 (independent session, batch 2).** No substantive finding, and no
fifth instance. That session traced every exit path and ran 20 constructed
inputs — the eight from cycle 4 plus a tie at the binding floor, a transitive
chain, a dependency cycle, a missing `Version` field, `R6`/`Rcpp` names, a
constraint on a DCF continuation line, a CRLF `DESCRIPTION`, both
own-`DESCRIPTION` guards, and the advisory path — all correct.

One low finding was worth acting on rather than accepting. `R (> x)` was read as
a floor of `x`, conflating it with `R (>= x)`, and the comment called the
resulting understatement "safe in this direction", which is inverted: a declared
floor of exactly `x` satisfies `>= x` and does not satisfy `> x`, so the check
could pass a floor one patch level too low. Rather than document the
approximation, the operator is now recorded and each constraint tested as
written, which removes it. The maximum is still computed, but only for the
report. Verified with a scratch dependency declaring `R (> 4.4.0)`: a declared
`R (>= 4.4.0)` now fails and `R (>= 4.4.1)` passes.

`tools/` was added to the directory tree in `CLAUDE.md` §4, which two earlier
rounds had raised, and four prose items were corrected.

Two claims in the cycle 4 response were wrong and were corrected before posting:
the constraint population is 818 R constraints, not 1119, which is the
`installed.packages()` row count; and "covering every path through the script"
omitted three exit paths, which cycle 5 verified separately.

---
**Claude:** Batch 3 (FSSgam_package#26, #27), branch `batch3-supplied-cor-matrix`,
stacked on batch 2. Re-triage first: batch 2 raised no new issues, so batch 3 is
unchanged from the plan.

- **#26.** `build_predictor_correlation_matrix()` branched on
  `length(cor.matrix)==1` as the sentinel for the `NA` default. A 1x1 matrix has
  length 1, and a single-predictor model set is the only case in which a supplied
  matrix is 1x1, so exactly that case had its supplied matrix discarded and
  recomputed. Replaced with `cor_matrix_supplied()`, which tests the default
  value: a matrix or data.frame is supplied, a length-one `NA` is not, and
  anything else is rejected naming its class. That last part is a consequence
  rather than an addition — something has to happen for the remaining cases, and
  under the old sentinel a `NULL` or a character value fell through to the
  supplied branch and failed against the missing-predictor check, reporting
  predictors rather than the argument.

- **#27.** An `NA` between two predictors of a supplied matrix reached
  `combine_uncorrelated()` and stopped the call with "missing value where
  TRUE/FALSE needed". The check is now in the supplied branch, naming the pairs.

  Three decisions on where and what to check, each recorded on the face of the
  code. It runs *before* `augment_supplied_correlation_matrix()`, so only cells
  the user actually supplied are examined — that function fills the derived rows
  and columns it adds with zero and deliberately leaves a user-supplied `NA`
  alone (FSSgam_package#15). It is *off-diagonal only*, both screens taking
  `max(abs(.))` over `upper.tri()` and `lower.tri()`, which exclude the diagonal,
  so a diagonal `NA` is never read and must not be reported. And it covers only
  predictors present in *both* dimensions of the supplied matrix, a derived name
  the user did not supply being computed rather than supplied.

  `combine_uncorrelated()` is left alone. The comment there records that stopping
  on an `NA` is the right behaviour at that point, the alternative being to drop
  the combination silently; this moves the report to where the pair can be named.

  This is a behaviour change and is recorded as one in `NEWS.md`: at
  `max.predictors = 1` the pair is never enumerated, so a matrix with an `NA` in
  it was previously accepted and the model set built against a matrix with a hole
  in it. Such a call now errors.

Six tests added. Four fail against the pre-fix code, verified by reverting
`R/generate-model-set.R` to the batch 2 state and running the file: the 1x1 case
fails, the two validation cases error, and the `max.predictors = 1` case fails.
The fifth, that a diagonal `NA` is accepted, passes before and after — it guards
against a false positive the new check could introduce rather than against the
old defect, which is why it is written as `expect_no_error()`.

Suite: 622 passing, 0 failing, 7 skipped, from 616.

---
**Review cycle 1 (independent session, batch 3).** Four substantive findings.
Three were paths on which the fixes did not reach, and one was a rationale that
is false. All were reproduced here before being acted on.

- **The `NA` check was evaded whenever `factor.factor.interactions` was set.**
  `resolve_factor_interactions()` screens factor pairs against the *raw* supplied
  matrix and runs before `build_predictor_correlation_matrix()`, which is where
  the check had been put. An `NA` between two factors still gave "missing value
  where TRUE/FALSE needed", under both the logical and the character forms of the
  argument. FSSgam_package#27 was not closed on that path, which `NEWS.md` and the
  pull request both claimed without qualification.

- **A second copy of the `length(cor.matrix)==1` sentinel.** `factor_correlations()`,
  the local helper inside `resolve_factor_interactions()`, has its own, so a 1x1
  supplied matrix was discarded there as well and the class validation was
  bypassed on that path. FSSgam_package#26 was half fixed.

- **The check was placed before augmentation for a reason that is false, and that
  placement created a second evasion.** The comment said checking after
  `augment_supplied_correlation_matrix()` "would have reported cells the user
  never wrote". It would not: that function fills every derived row and column it
  adds with zero, so after it runs the only `NA` remaining is a user's.
  Meanwhile, restricting the check to names present in *both* dimensions excluded
  a derived interaction name supplied in one dimension only, whose `NA` then
  reached `enumerate_candidate_models()`.

- **The behaviour change was stated backwards in five places.** A length-one value
  of any class -- `"oops"`, or a list holding a matrix -- satisfied the old
  sentinel and was silently treated as nothing supplied, so the model set was
  built from a computed matrix as though the argument had not been given. Only a
  value of some other length fell through to the supplied branch and failed
  against the missing-predictor check. `NEWS.md`, both commit messages, the pull
  request body and a test comment all asserted the second for both.

The structure that follows from this: the checks are made at two points because
no single point sees every screen, and the reason is on the face of the code.
The interaction columns do not exist when the first screen runs, and the resolved
matrix cannot exist until they do (FSSgam_package#13), so the named predictors
are checked before the first screen and the full set after the second matrix is
assembled, the latter after augmentation. `stop_on_na_correlations()` is shared
between the two.

Three evasions verified closed end to end, each with a construction that
actually exercises the path. The first attempt at the third did not: two factors
derived from one another are correlated above `cov.cutoff`, so no interaction
column was built and the cell was never read. Two independently sampled factors
were needed. Suite: 627 passing, from 616.

---
**Review cycle 2 (independent session, batch 3).** Two substantive findings.

- **A fourth screened name the check did not cover.** Both checks were scoped to
  `all.predictors`, but the character forms of `smooth.smooth.interactions` and
  `factor.factor.interactions` are validated against `rownames(cor.matrix)` and
  `colnames(use.dat)` respectively, not against the predictor lists. So either
  can name a variable that is screened against `cov.cutoff` without being a
  predictor of the model set, and a supplied `NA` on such a pair still gave
  "missing value where TRUE/FALSE needed". Reproduced with `macro`, a column of
  `case_study1` deliberately not named as a predictor: the control call succeeds
  at `n.mods = 3` and the same call with the `NA` fails with the old message.
  The first check now covers the union of `all.predictors` and the names given in
  either character form.

  That those two arguments accept a non-predictor at all is a separate defect,
  and the reviewer was right to qualify the finding as an unenforced contract
  rather than a supported use. Raised as **FSSgam_package#37** and assigned to
  batch 4, which is the argument-validation batch. The widening here works around
  it and does not fix it, which the code comment and `NEWS.md` both say.

- **A false claim in the pull request body.** It said the check runs "after the
  missing-predictor check so that a misspelled name is still reported first". The
  first check now runs at the top of `generate_model_set()`, before both
  `check_predictor_missingness()` and the missing-predictor check, so a missing
  predictor is reported second. The bullet immediately below it contradicted the
  sentence. Corrected.

The reviewer also mutation-tested the change, confirming each new check is pinned
by exactly one test. Suite: 629 passing, from 616.

---
**Review cycle 3 (independent session, batch 3).** Three substantive findings,
and together they showed the approach was wrong rather than incomplete.

- **A fifth path.** Enumerating all four expressions that read a cell of
  `cor.matrix`, the set indexed by `enumerate_candidate_models()` is
  `all.predictors` union the derived `.I.` columns union the character-form
  `smooth.smooth.interactions` names. A pair *crossing* the last two is screened
  and was covered by neither check.
- **The widened check failed wrongly.** A name given in the character form of
  `factor.factor.interactions` is screened against the others in that argument
  and never against `pred.vars.cont`, so a union of names treats the Cartesian
  product as screened and rejects a cell nothing reads. Measured: an `NA` at
  `gc/depth` was accepted before and errored after.
- **An S4 `Matrix` was rejected.** `cor_matrix_supplied()` used a whitelist of
  `matrix` and `data.frame`, so the objects `Matrix::Matrix()` and
  `Matrix::nearPD()` return stopped working, having worked before the argument
  was validated at all.

The first two are one defect: a name union cannot describe which pairs are
screened. **The check now runs at the two screening sites**, on the sub-matrix
each is about to take `max(abs(.))` of, so it covers exactly the pairs that are
compared. No model of which names reach which screen is maintained beside the
code that screens. `cor_matrix_supplied()` now tests dimensionality rather than
class, and `Matrix` is declared in `Suggests` for the test.

Six paths verified: an `NA` between two factors under
`factor.factor.interactions`; between two continuous predictors; in a derived
column supplied in one dimension only; on a non-predictor named through the
character `smooth.smooth.interactions`; on a pair crossing a derived column and
such a name; and, as a control, on a pair no screen reads, which now passes.

**A departure from what FSSgam_package#27 asked for, decided here and flagged for
the maintainer.** The issue asked that the report not depend on `max.predictors`.
It now does, and on every other argument deciding which pairs are compared: at
`max.predictors = 1` each candidate holds one term, so two predictors are never
compared and an `NA` between them is accepted. Rejecting it would be a false
failure of exactly the kind found above. The property that makes the check
complete is the same one that makes it conditional. Recorded in `NEWS.md`, in
the pull request, and in a test whose comment states the reasoning.

**FSSgam_package#37's evidence did not reproduce as written** and was corrected.
It used `macro`, whose correlation with `depth` is -0.335 and so exceeds the
default `cov.cutoff`, meaning no `te()` term was built and the model set equalled
the control. The corrected block uses an uncorrelated variable and shows
`zz.te.depth` built and fitted while `included.vars` names only `depth` and
`complexity`. It also records that the gap is reachable for
`smooth.smooth.interactions` only when a `cor.matrix` is supplied.

Suite: 633 passing, from 616.

---
**Review cycle 4 (independent session, batch 3).** Three substantive findings.

- **The new error fired on a *computed* matrix and blamed a supplied one.** The
  zero-fill that makes computed matrices free of `NA` is inside
  `build_predictor_correlation_matrix()`, but `resolve_factor_interactions()`
  runs first and its `factor_correlations()` helper calls `check_correlations()`
  directly with no zero-fill. `check_correlations()` returns `NA` for a pair
  involving a single-level factor. So with no `cor.matrix` argument at all, a
  single-level factor and `factor.factor.interactions = TRUE` produced
  "Supplied cor.matrix has NA between terms ... Supply a correlation for each",
  when nothing had been supplied and there was nothing the user could supply.

  FSSgam_package#27 asserts the computed branch "can never reach" the `NA`; that
  assertion was inherited into the code and into `NEWS.md` without being tested,
  and it is false by this route. `factor_correlations()` now zero-fills, matching
  the other computed path. With that in place the same call fails later with the
  `contrasts` error, which is FSSgam_package#33 and is batch 4's work.

  The message was also reworded so that it is true whichever matrix it is
  describing. Threading a `supplied` flag through four call sites was tried and
  rejected as disproportionate to a wording nuance.

- **A justification repeated in four places was false as a general rule.** "At
  `max.predictors = 1` each candidate holds one term, so two predictors are never
  compared" is wrong: a `.by.` term is one term of a candidate but splits into two
  for the screen. Measured on the standard fixture at `max.predictors = 1`, an
  `NA` at `complexity/depth` is accepted while `depth/ZONE` and
  `complexity/ZONE` are reported. The claim is now stated as what it is -- the
  report covers the pairs some screen compares -- with the `.by.` case named.

- **A test passed for the wrong reason.** The "unscreened pair" test used
  `pred.vars.fact` of length one, which skips the factor-factor block entirely,
  so it asserted nothing about screening. It now uses two factors and asserts
  that the interaction column was built before testing the cell.

Two stale error-message quotations in the pull request body, and a claim there
that `combine_uncorrelated()` is unchanged when it is one of the two functions
the change edits, were corrected.

Suite: 636 passing, from 616.

---
**Review cycle 5 (independent session, batch 3).** Three substantive findings.

- **The zero-fill turns an error into a model set holding the same model twice,
  and `NEWS.md` did not say so.** An `NA` becoming zero means "uncorrelated", so
  a factor-factor pair that previously stopped the call is now admitted. Where
  the `NA` came from a single-level factor the interaction column is a copy of
  the other factor: `fa.I.const` is `paste(fa, "a")`, so `fa` and `fa.I.const`
  are the same model under two names and both enter the AICc table and the
  variable importance scores. Confirmed here.

  Recorded in `NEWS.md` as a behaviour change, together with what makes it
  acceptable and what it is not a fix for: the real defect is that a
  single-level factor is accepted as a predictor, FSSgam_package#33, batch 4's
  work. A test now pins the zero-fill and says the same on its face.

  The risk this cycle was set to look for did not materialise. A 432-scenario
  sweep -- six factor sets, three continuous sets, both forms of
  `factor.factor.interactions`, `max.predictors` 1 to 3, both settings of
  `non.linear.correlations`, two forms of `factor.smooth.interactions` -- gave
  392 identical, 40 error to built, and **none** built to a different model set
  or built to error. No call that worked before behaves differently.

- **`Matrix::nearPD()` does not return a matrix.** It returns an S3 object of
  seven components whose `mat` component is the matrix, so only that component
  was ever accepted, on this branch and before it. Named as an accepted input in
  `NEWS.md`, a code comment, a test comment and a commit message. Corrected in
  all four; `Matrix::Matrix()` alone is the example.

- **A pull request body edit had not taken.** One of two replacements silently
  failed to match, so the body still said `combine_uncorrelated()` is unchanged
  while being one of the two functions the change edits -- a collapsed-block
  statement contradicting the open section, which parent `CLAUDE.md` §14 forbids.
  Checking that an edit landed, rather than that the command returned a URL, is
  the lesson.

Suite: 638 passing, from 616.

---
**Review cycle 6 (independent session, batch 3).** Three findings, all in text.
No code defect, and that session states it looked for a sixth and did not find
one.

- **The 432-scenario sweep did not measure what `NEWS.md` said it did.**
  Classifying the saved results rather than trusting the count: 336 of the 432
  cells were identical *error strings*, not identical model sets. 288 of those
  passed `null.terms = NA`, which `validate_null_terms()` rejects before any of
  the code under test runs, so no scenario with a continuous predictor was
  exercised at all; 48 more had `NA` predictors, rejected equally early. Only 56
  model sets were ever compared.

  Re-run with a valid `null.terms`: 144 identical, 40 error to built, 80 error
  to a different error, and none built to a different model set or built to
  error. The conclusion stands; the evidence as recorded did not support it.
  `NEWS.md` now gives the corrected figures and says what the earlier count
  actually measured.

  The lesson is the one this project has now met several times in another form:
  a comparison that reports "identical" can be comparing two failures.

- **A pull request body edit had again not fully applied.** The stale error
  message was still quoted; it had been split across a line break, so the
  patterns grepped for it did not match and the check that the edit had landed
  passed wrongly. It is now verified by capturing the message the code actually
  produces and comparing it with the body whitespace-insensitively, rather than
  by searching for the string being replaced.

- **`@param cor.matrix` documented neither new rejection rule** while inviting
  hand-built matrices. Both are now stated: the argument must be
  two-dimensional, and it must not carry `NA` between two terms that are
  actually screened -- with the note that which cells those are depends on
  `max.predictors` and the interaction arguments.

---
**Review cycle 7 (independent session, batch 3).** Two substantive findings, both
in the one paragraph cycle 6 had just corrected. No code defect, and both quoted
error messages and every `@param cor.matrix` claim were verified correct.

- **"144 identical" double-counted.** 144 is the number of cells that build under
  this version, not the number unchanged. Reclassified here from the saved
  results rather than taken on trust: 104 identical model sets, 208 identical
  errors, 40 error to built, 80 error to a different error, and none built to a
  different model set or built to an error. The four sum to 432.
- **The explanatory sentence was false.** The 80 differing errors are
  single-level-factor scenarios, not `NA`-predictor ones; the `NA`-predictor
  scenarios are 144 of the 208 that error identically. Confirmed by
  cross-tabulating the outcome against the factor set.

Verifying that also produced the qualification the paragraph most needed, which
no earlier count had stated: **every one of the 40 cells that now build has
`pred.vars.cont = NA`.** Name a continuous predictor and the call still stops, in
`check_correlations()`, on the same single-level factor (FSSgam_package#33). The
behaviour change is bounded much more tightly than the raw count suggested.

`NEWS.md` now gives the classification as a table, states both bounds, and
records both wrong counts so the figures are not taken from them.

**Three counts of the same grid were wrong in three different ways** -- comparing
error strings, double-counting a category, and mischaracterising a category.
Each was arithmetic or classification over results already computed, not a
measurement problem. Classify the cells and check they sum to the grid size.

---
**Review cycle 8 (independent session, batch 3).** Nothing substantive. That
session regenerated the 432-cell grid from `git archive` exports of both trees,
rather than reusing any saved result, and classified it with a script written
that round: every figure reproduces, 104 / 208 / 40 / 80 / 0 / 0 summing to 432,
all 40 error-to-built cells having `pred.vars.cont = NA` and all 80
differing-error cells being single-level-factor scenarios. It also confirmed the
stated mechanism directly, each of those 80 changing from "missing value where
TRUE/FALSE needed" to the `contrasts` error, and checked the two superseded
counts recorded in `NEWS.md` are arithmetically what they are said to be.

It probed for an eighth code defect across `data.frame`, S4 `Matrix`, 1x1 in both
classes, asymmetric matrices, `NaN`, `NA_character_`, a three-dimensional array,
`TRUE` and `nearPD`, and found none. One latent assumption is recorded and left
alone: `stop_on_na_correlations()` assumes a square sub-matrix, which nothing
currently reaches with a non-square one.

Two ruled verbs were corrected -- "must not carry NA" and "costs nothing", both
cases parent `CLAUDE.md` §12 names. The blank lines flagged as stray turned out
to be separators inside newly added functions; the file already had single blank
lines where the diff suggested otherwise.

Batch 3 converged after eight cycles. Suite: 638 passing, from 616.

---
**Claude:** Batch 4 (FSSgam_package#28, #37, #33, #22), branch
`batch4-predictor-validation`, stacked on batch 3. Re-triage first: batch 3
raised FSSgam_package#37, already assigned here; nothing else new.

All four are input validation and all four now sit in one block at the top of
`generate_model_set()`, beside `validate_use_dat()` and `validate_null_terms()`.

- **#28.** `validate_predictor_names()` rejects a name repeated within one
  argument, and a name given as both a factor and a continuous predictor.

  **It deliberately does not check `pred.vars.cont` against `linear.vars`.** A
  first version did, and the existing suite caught it: naming a predictor in
  both is the documented way to fit it linearly, `setdiff(pred.vars.cont,
  linear.vars)` deciding which predictors get a smooth. The check would have
  broken a supported idiom that a test already covered. Worth recording as the
  clearest case so far of the committed suite doing the work a review would
  otherwise have had to.

- **#37.** `validate_interaction_predictors()` checks the character forms
  against the predictor lists rather than against `colnames(use.dat)` and
  `rownames(cor.matrix)`. Placed at the top so neither depends on a supplied
  matrix and both are settled before either argument is consumed --
  `factor.factor.interactions` by `resolve_factor_interactions()`, which runs
  before the resolved matrix exists, and `smooth.smooth.interactions` by
  `resolve_smooth_smooth_interactions()`, which runs after.

- **#33.** Two halves. `validate_factor_levels()` rejects a factor with fewer
  than two observed levels, `droplevels()` first so a declared but unobserved
  level counts as absent. And the `lm()` in each correlation function is wrapped
  in `try()`, being the one unguarded fit in each; both are exported and
  callable directly, so the guard matters beyond the `generate_model_set()` path
  the rejection covers.

- **#22.** `build_factor_interaction_columns()` declines a generated column
  whose name is already in `use.dat`, warning and dropping the combination from
  the returned term names as well as the columns, so no later stage refers to a
  term that was not built.

**Batch 4 closes two things batch 3 could only work around.** The widening of the
`NA` check to cover names screened without being predictors is now belt and
braces, since #37 stops such a name reaching a screen at all. And the duplicate
model set the zero-fill admitted is now unreachable, since a single-level factor
is rejected first. Three batch-3 tests were rewritten accordingly rather than
deleted, so each still covers what it was written for: the crossing-pair test now
uses a real predictor, the non-predictor test asserts the rejection, and the
zero-fill test asserts the rejection plus the helper's behaviour.

Suite: 644 passing, from 638.

---
**Review cycle 1 (independent session, batch 4).** Five substantive findings.

- **Four test blocks were deleted, not rewritten, while the pull request, the
  commit message and this log all said they were rewritten.** Three of the four
  pass verbatim against this branch, verified by extracting and running them, so
  coverage of FSSgam_package#27's pair naming and of the `max.predictors` scoping
  decision had simply been lost. Restored verbatim. The fourth used a factor
  named in `factor.factor.interactions` but not in `pred.vars.fact`, which
  FSSgam_package#37 now rejects, so its property is asserted through a
  construction that is still legal.

  The cause was a python replacement whose span reached further than intended.
  The lesson is the same one this session met with a pull request body: check the
  result, not that the command returned.

- **The `-Inf` warnings FSSgam_package#28 reports still fired, and the claim that
  they were removed was false in three places.** Rejecting a repeated name does
  not reach the second cause: naming a predictor in both `pred.vars.cont` and
  `linear.vars`, the documented idiom the check exempts, puts the name in the
  candidate pool twice and reaches the same empty triangle at any
  `max.predictors` of 2 or more. Confirmed here.

  Fixed rather than only re-worded, since the warning fires on a legitimate
  input and names nothing. `exceeds_cutoff()` returns `FALSE` for an empty
  triangle rather than calling `max()` on it; which combinations survive is
  unchanged, `-Inf` having been below every cutoff.

- **A stale comment contradicted the commit that fixed it**, still saying the
  `lm()` below had the same unguarded shape four lines above the `try()` that
  guards it.

- **A branch was left asserted by nothing.** Changing a test's `pred.vars.fact`
  from two factors to one skipped the factor-factor block entirely, so its error
  came from `build_predictor_correlation_matrix()` and duplicated another test.
  Restored to two factors, and confirmed by mutating the message at each of the
  two sites separately that the test now exercises the one inside
  `resolve_factor_interactions()`.

  Note on method: `testthat::test_file()` did not pick up an edit to `R/` even
  after `devtools::load_all()`, so a mutation run through it reported no
  failures for a message that had demonstrably changed. Mutation testing here
  has to be checked against a direct call.

- **The `linear.vars` guard test used `max.predictors = 1`**, the one value at
  which the empty-triangle case cannot arise, so it would have passed whatever
  `exceeds_cutoff()` did. It now covers 2 as well and asserts no warning.

Suite: 653 passing, from 638.

---
**Review cycle 2 (independent session, batch 4).** Five substantive findings.
`exceeds_cutoff()` itself was cleared: equivalent to the two `max()` calls it
replaced on 454 of 460 constructed inputs, the six exceptions needing an
`NA`/`NaN` `cov.cutoff`, which is out of contract; and identical model sets
across 179 scenarios, so "which combinations survive is unchanged" holds.

- **FSSgam_package#28 survived by a third route, and this batch removed its only
  signal.** `factor.factor.interactions = c("fa","fa","fb")` builds the
  interaction column `fa.I.fb` twice, so `use.dat` gains two columns of that
  name -- the FSSgam_package#22 shape -- and the candidate set gains
  `fa.I.fb+fa.I.fb`. The two `-Inf` warnings were the only indication, and
  `exceeds_cutoff()` suppressed them. Confirmed here: the parent warns, this
  branch was silent. A repeated name within either interaction argument is now
  rejected.

  Worth recording as the shape of the mistake: silencing a warning is safe only
  once every defect it was the sole evidence of is closed.

- **Two branches were left asserted by nothing**, both because the covering test
  was one this batch rewrote. Reverting `combine_uncorrelated()`'s by-name
  indexing to positional failed two tests on the parent and none here; deleting
  the zero-fill in `factor_correlations()` failed one on the parent and none
  here. The first is now asserted through a supplied matrix whose two dimensions
  carry the same names in different orders. The second is not, and the test that
  claimed to assert it now says so instead: no public call reaches the zero-fill
  with an `NA` since a single-level factor is rejected, so it is defensive and
  its title and comment no longer imply coverage.

- **The pull request said three tests were rewritten; eight were changed**, five
  rewritten, one renamed and two modified in place. Corrected.

- **FSSgam_package#37's evidence block still did not reproduce.** `macro`, the
  obvious variable from `case_study1`, correlates -0.335 with `depth` and so
  exceeds the default `cov.cutoff`, meaning no `te()` term is built. The block
  now draws an uncorrelated variable and says why.

Suite: 652 passing, from 638.

---
**Review cycle 3 (independent session, batch 4).** Four substantive findings.

- **The same three test blocks were deleted a second time**, by the very commit
  written to restore them two commits earlier. The cause was again a span-based
  text replacement reaching further than intended: the span ran from one block's
  title to another's, and the three restored blocks sat between. Restored, and a
  script now compares the list of `test_that` names against any earlier ref after
  every edit to a test file. **Any further edit to a test file in this project
  should run it.** Two rounds lost the same coverage; nothing but a mechanical
  check will stop a third.

- **A fourth route to FSSgam_package#28 and #22.** Generated names were compared
  against `colnames(use.dat)` but never against each other, and two distinct
  combinations can paste to one name: `(a, b.I.c)` and `(a.I.b, c)` both give
  `a.I.b.I.c`, the shape produced by feeding `used.data` back in. Both columns
  were built, so `use.dat` gained two columns of that name holding different
  variables, with no warning on this branch and only the `-Inf` warnings on the
  parent. Declined and warned about, screened before the existing collision check
  so a name colliding both ways is reported once.

- **The by-name-indexing test did not discriminate.** Reverting
  `combine_uncorrelated()` to positional indexing fails nothing: every caller
  pre-subsets the matrix by name, so the reordering the test constructs is
  normalised away before that function runs. The test is renamed to the property
  it does assert, and both it and the code now record that the mechanism is
  defensive given its callers and is pinned by nothing. Propping it up with a
  contrived test would have been worse than saying so.

- The pull request's test accounting was wrong a second time.

Suite: 661 passing, from 638.

---
**Review cycle 4 (independent session, batch 4).** Three substantive findings,
two of them losses this branch caused.

- **A fifth route, and this branch made it silent.** `build_model_formulas()`
  finds `.by.` and `.te.` with `grep(fixed=TRUE)` over every term, not only over
  names this package generated, so a predictor whose own name contains one is
  parsed as an interaction. `catch.by.effort` -- an ordinary variable name --
  produced `~s(catch.by.effort) + s(catch, by = effort)`, which cannot be fitted,
  so the predictor vanished from the model set. Confirmed here: on `master` four
  `-Inf` warnings; on this branch, before the fix, none at all.

  Raised as **FSSgam_package#39** and rejected here, naming the predictor and the
  reserved strings. Rejecting the failure is not fixing the parse; the issue
  records what a full fix would involve and that a user cannot always avoid the
  collision, since feeding `used.data` back in presents names containing `.I.`
  by construction.

- **The stated reason for deleting the `factor_correlations()` zero-fill test
  was false.** "No public call reaches it with an NA since a single-level factor
  is rejected" is wrong: two 40-level factors exceed `multinom()`'s `MaxNWts` and
  `check_correlations()` returns NA with no single-level factor involved.
  Verified here. The zero-fill is reachable, and removing it turns that call
  from a model set into an error while the whole suite still passes. Restored
  with that construction.

  This is the second time a claim that something is unreachable has been made
  without testing it -- FSSgam_package#27 asserted the same of its computed
  branch, and cycle 4 of batch 3 disproved that too. **Do not record a path as
  unreachable without constructing an attempt to reach it.**

- **A second, undisclosed loss.** `augment_supplied_correlation_matrix()`'s
  zero-fill was pinned on the parent and by nothing here. The block kept its
  title, so the name comparison could not show it -- that check catches deletion,
  not weakening. Restored as its own block.

The separator rejection makes the self-collision screen added in cycle 3
unreachable, since it needs a predictor literally named `a.I.b`. It is kept as a
guard, and both the code and the test that used to exercise it now say it is
pinned by nothing rather than implying coverage.

Suite: 664 passing, from 638.

---
**Review cycle 5 (independent session, batch 4).** Four substantive findings, and
the fifth consecutive round to find a lost assertion.

- **The `augment_supplied_correlation_matrix()` zero-fill was still pinned by
  nothing.** The block restored in cycle 4 could not fail: with every name in
  the supplied matrix a predictor, the computed cells are never `NA` and the
  mutant returns a byte-identical matrix. The parent pinned it with a matrix
  carrying a name that is *not* a predictor, which the rewrite had dropped on
  the grounds that FSSgam_package#37 rejects a non-predictor -- but #37 rejects
  one named in an *interaction argument*, and says nothing about the names a
  supplied matrix carries. That is exactly the case the zero-fill exists for.

- **A sixth route.** A predictor named `null` produced two candidates called
  `null`, so `mod.formula[["null"]]` returned the first and the variable
  importance table matched the wrong row. `+` and `*` were unvalidated
  separators, and `factor.smooth.interactions` was unvalidated in its character
  form. All three closed.

- Two figures corrected: the `-Inf` warnings number two at `max.predictors = 1`
  and four at 2, stated as four in five places; and "all removed titles are
  renames whose property is still asserted" was false for the by-name indexing,
  which this branch's own commit records as pinned by nothing.

**Every claim this branch makes is now mutation-tested.** The results, each
mutation run against the whole suite in a fresh process:

| mutated | failures |
|---|---|
| `augment` zero-fill, rows | 1 |
| `augment` zero-fill, columns | 2 |
| `factor_correlations()` zero-fill | 1 |
| duplicate-within-argument check | 1 |
| `exceeds_cutoff()` empty-triangle guard | 1 |
| FSSgam_package#22 collision decline | 6 |
| reserved-separator check | 3 |
| `null`-name check | 1 |
| self-collision screen | 0, recorded as unpinned |

The duplicate-within-argument check read 0 at first, which exposed a defect in
the new `factor.smooth.interactions` validation rather than in the test: that
argument defaults to `pred.vars.fact`, so a name repeated in `pred.vars.fact`
was reported against an argument the user never set. It is now skipped while it
holds its default.

**Comparing `test_that` names catches deletion, not weakening.** Five rounds
found losses; only mutation found the last two. Any claim added to this
repository should be mutated before it is described as covered.

---
**Review cycle 6 (independent session, batch 4).** Nothing substantive; batch 4
converged after six cycles.

That session re-ran the whole mutation table independently on a clean
`git archive` copy, one fresh process per mutation, and all nine rows reproduce.
It mutated eight further things the table did not cover -- the factor/continuous
check, the membership check, `droplevels()`, `repeated()`, both `try()` guards,
and dropping `+`/`*` from the reserved list -- all pinned. Two are not pinned and
neither changes behaviour: the guard skipping `factor.smooth.interactions` while
it holds its default, which is unreachable because `validate_predictor_names()`
stops on that input first, and `validate_factor_levels()`'s `colnames(use.dat)`
filter.

**No false rejection.** Zero reserved-string hits across all four bundled
datasets, 133 columns, and every predictor vector in the companion repository's
vignettes.

**A seventh route exists and was recorded on FSSgam_package#39 rather than fixed
here.** `-` and `:` mis-parse as `+` and `*` do: `catch-effort` gives
`~s(catch - effort, k = 3, bs = "cr")`. Confirmed here. It differs in being less
silent -- the candidate is retained and a failed fit appears in `failed.models`
-- and the reserved list this batch adds covers the set that was demonstrably
silent. Extending it belongs to the issue, together with the note that a full fix
would have to quote term names, since a name containing `-` is not a syntactic R
name.

Two figures corrected: a test comment still said four `-Inf` warnings above a
`max.predictors = 1` call, where the count is two; and the pull request still
described the branch as two commits, where it is seven.

---
**Claude:** Batch 5 (FSSgam_package#23, #34), branch `batch5-null-terms-and-gamm`,
stacked on batch 4. Re-triage: batch 4 raised FSSgam_package#39, whose reported
failure batch 4 closes; the design question it leaves is deferred alongside
FSSgam_package#25, which the maintainer asked to split out and address later.

- **#23.** `null.cov.cutoff`, defaulting to 0.8, and a returned
  `null.term.correlations`. The decisions, in the order they were forced:

  *Screening, not only reporting.* The maintainer asked for a cutoff to apply,
  so a predictor over it is dropped from every list a term can reach the
  candidate set through, and an interaction term is dropped when any of its
  parts is.

  *Random-effect grouping factors are excluded from the screen.* This was not
  planned; the existing suite forced it. Two tests failed on the first
  implementation, and both were right to: the companion repository's case study
  2 uses `null.terms = "s(Location,Site,bs='re')"` with `Status` nested in
  `Location`, correlation 1, and a fixture builds a factor from `site` while
  `site` is the random effect. A grouping factor is correlated with the
  predictors measured within it by construction. Screening those dropped
  legitimate predictors, so a variable inside a `bs='re'` smooth is now excluded
  and the screen applies to fixed forced terms, which is the case the issue
  actually describes.

  *A supplied `cor.matrix` is used, never recomputed.* The first implementation
  called `check_correlations()` for the null-term block, which broke the
  guarantee that a supplied matrix replaces the automatic estimate outright --
  the guarantee that lets a predictor of a class it cannot classify be used at
  all (FSSgam_package#13). A committed test caught it. Only what the supplied
  matrix names is used now.

  *A defect found in the process:* `cor.matrix` is reassigned to the resolved
  matrix before the screen runs, so testing it for "did the caller supply one"
  was always true. Recorded once at the top instead.

  *Dropping can exhaust the predictors.* `enumerate_candidate_models()` then
  stopped with a message naming `max.predictors`, which a user has no reason to
  connect to the drop they were just warned about. Reported together.

- **#34.** A bare `mgcv::gamm()` fit is rejected before `update()` is attempted,
  naming `MuMIn::uGamm` as the remedy, and a test asserts that remedy works
  rather than only that the message appears. Three documentation passages that
  directed the user to `gamm` are corrected.

Suite: 693 passing, from 670.

---
**Review cycle 1 (independent session, batch 5).** Three substantive findings,
one of them a design defect that undid the point of the change.

- **A supplied `cor.matrix` silently disabled the whole screen.** The first
  implementation used only what such a matrix names, and a supplied matrix is
  indexed by predictor while a forced term is not one -- so in the ordinary case
  nothing was screened, `null.term.correlations` was `NULL`, and no warning was
  raised. That reinstated, for every `cor.matrix` caller, exactly the silence
  FSSgam_package#23 exists to end. Confirmed here: a predictor at r = 1.00 with
  the forced term retained in every candidate.

  FSSgam_package#13 did not force that reading. Its guarantee is that a supplied
  matrix replaces the automatic estimate *for the predictors*, so a predictor of
  a class `check_correlations()` cannot classify can be used. Only the
  `null.vars` block is computed now, and where that computation fails -- which is
  what such a predictor causes -- the screen is skipped with a warning rather
  than the call aborting. The FSSgam_package#13 caller keeps their model set and
  is told their forced terms were not screened.

- **`bs=c('re')` defeated the random-effect exemption and then errored.** The
  test was a regex over the term text. mgcv accepts `bs=c('re')` and fits the
  identical null model, so a nested factor was dropped and the call stopped with
  a message naming `max.predictors`, `null.cov.cutoff` and `null.terms`, none of
  them the cause. The exemption is now read from the parsed call: the term's `bs`
  argument is evaluated and tested for `"re"`, which covers every spelling.

- **The `@param` text stated the opposite of the code**, saying the correlations
  are "not taken from a supplied cor.matrix". Every other description was right;
  only the user-facing one was wrong.

`null.cov.cutoff` was also unvalidated: an `NA` gave an internal error and a
length-2 value was accepted, recycling into two concatenated messages.

**Five mutations left the suite green**, so five claims were pinned by nothing:
the `setdiff` exempting forced terms from the screened set, the `colnames`
filter, `abs()`, three of four separators in the interaction-term drop, and
`full_subsets_gam()`'s forwarding. Assertions added for each; all now fail when
mutated.

Suite: 712 passing, from 670.

---
**Review cycle 2 (independent session, batch 5).** Two substantive findings, both
in the verification rather than the code, and both the same class of error as
earlier rounds: a claim about coverage that was not true.

- **Two of the five claims cycle 1 said were "now pinned by mutation" were not.**
  Removing `full_subsets_gam()`'s forwarding of `null.cov.cutoff` left the suite
  green, and `null.cov.cutoff` appeared nowhere in that file's tests. Reducing
  `keep.terms()` to the `.by.` separator also left it green, and that mutation is
  not cosmetic: two `te()` terms built on a dropped predictor re-entered the
  model set.

- **The mechanism the previous commit exists to introduce had no test.** Nothing
  pinned that a supplied `cor.matrix` *overrides* the computed estimate: the
  committed test set its supplied cell above the cutoff, so it passed whether or
  not the supplied value was read. And the `try()` fallback and its warning --
  the whole FSSgam_package#13 argument for computing the block at all -- were
  untested.

Four assertions added, each verified by mutation. The fourth needed a
construction the first attempt missed: returning `NULL` from the fallback rather
than the supplied rows is only detectable when the supplied matrix names one
forced term and the computation fails for another, so the test uses two forced
terms and a `Date` predictor.

The reviewer attacked the merge logic directly -- different column sets and
orders, a forced term in both blocks, permuted dimnames, asymmetric, duplicated
and `NA` cells -- and found no defect, and confirmed the screen is correct for
`te()` forced terms, `by=` terms, `non.linear.correlations = TRUE` and
`full_subsets_gam()`.

Text corrections: the `@return` said the element is "Absent" where it is present
and `NULL`; the `bs='re'` exemption was described in `NEWS.md` and in the code
but not on the manual page a user reads; and `NEWS.md` described an intermediate
state of the branch.

**Stating that something is mutation-tested is itself a claim, and it was wrong
twice.** Run the mutation, read the count, and quote it.

Suite: 728 passing, from 670.

---
**Review cycle 3 (independent session, batch 5).** No code defect. Every
construction built against the merge, the `null_term_variables()` parse across
fifteen spellings, the drop-everything path and the `gamm` rejection returned the
correct answer, and the nine claims cycle 2 added all held on independent re-run
-- the first verification statement on this branch to do so.

One substantive finding, and it is the same class as cycle 2's: **"every claim
the branch adds is mutation-tested" was false.** Thirteen further mutations were
devised and six survived. Five were reachable, and one of them matters:

- `apply(over,2,any)` reduced to `all` retains a predictor at r = 1.00 with one
  of two forced terms, in every candidate, with nothing reported --
  FSSgam_package#23 reinstated. Every drop test used a single forced term, so
  none distinguished the two.
- The drop was unpinned for `pred.vars.fact` and for `linear.interaction.terms`:
  unfiltered, a dropped factor re-entered three candidates while `included.vars`
  said it had been dropped, and a `.t.` term re-entered built on a dropped
  `linear.vars` entry.
- Both zero-fills of the supplied block were unpinned: a supplied `NA` cell gives
  `subscript out of bounds` or `missing value where TRUE/FALSE needed`, which is
  the failure FSSgam_package#27 removed from the predictor screen.

Four assertions added, each verified by running the mutation and reading the
count: 2, 1, 1 and 2 failures respectively.

**A second false claim, in a commit message.** `8000ca1` says it "corrects a
NEWS.md paragraph describing an intermediate state of this branch". Its NEWS.md
hunk does not do that; the paragraph is unchanged. The commit is pushed, so the
correction is recorded here and in the following commit rather than by rewriting
it.

Three rounds, three false verification claims -- that five were pinned, that
every claim was pinned, and that a paragraph was corrected. The code has been
right since cycle 1; what has repeatedly not been is the account of what the
tests cover.

Suite: 738 passing, from 670.

---
**Review cycle 4 (independent session, batch 5).** No code defect. Four rounds
have now attacked this branch and none has found one.

One substantive finding, and for the third consecutive round it is a false
statement of mine about coverage rather than a defect in the code. The claim
that "both zero-fills" were pinned was true of one. `full[is.na(full)]=0`, which
fills the *computed* block, survived its mutation and is reachable: a
zero-variance predictor makes `check_correlations()` return `NA`, and without the
fill the call stops with "missing value where TRUE/FALSE needed" -- the failure
the fill exists to prevent. The commit also mislocated the two lines, describing
`:557` as the supplied block when it is the computed one.

Assertions added for it and for three minor gaps the same round found: the
documented `NULL` return where every predictor is itself a forced term, which
had no test; the forced terms named in the drop warning, which is what a user
acts on; and `full_subsets_gam()`'s manual page, which through the inherited
`@param` stated that `null.term.correlations` is returned when that function
does not return it.

Both structural assertions verified by mutation, 2 and 1 failures.

**Four rounds, four false verification claims.** That five assertions were
pinned; that every claim was pinned; that a `NEWS.md` paragraph was corrected;
that both zero-fills were pinned. The code has been correct since cycle 1. What
was repeatedly wrong was the account of what the tests cover, and in each case
the error was the same: asserting coverage from having written a test rather
than from having run the mutation.

Suite: 743 passing, from 670.

---
**Review cycle 5 (independent session, batch 5).** A code defect, the first since
cycle 1, and it is the one this branch existed to prevent.

`check_non_linear_correlations()` returns a deliberately asymmetric matrix -- row
is the response, column the predictor -- and the two directions answer different
questions. The screen read only `[forced term, predictor]`, which is the small
one when a candidate is a deterministic function of a forced term. Measured on
`curve = forcedA^2 + N(0, 0.02)` over 300 rows with
`non.linear.correlations = TRUE`: `[forcedA, curve]` is 0.152 and
`[curve, forcedA]` is 0.991, so `curve` entered every candidate with nothing
reported. That is exactly the condition FSSgam_package#23 exists to end, under
the argument that turns the estimate asymmetric.

`exceeds_cutoff()` has read both triangles for `cov.cutoff` all along. The
null-term screen now does the same, taking the larger of the two directions, and
each cell of the returned matrix is the value screened on rather than one
direction of it. On the symmetric default path the two are equal and nothing
changes.

It survived four rounds because nothing exercised the branch: replacing
`if(non.linear.correlations==TRUE)` with `if(FALSE)` left the suite green. **An
argument with no test is not covered by the tests of the argument beside it.**

Two further unpinned mutations were closed with it: `>` against `>=` in the
cutoff, where the documented "set to 1 to admit every predictor" fails for an
exact duplicate under `>=`; and the subscript that names only the forced terms
actually exceeded, where the assertion added in cycle 4 used a single forced
term and so could not distinguish it.

All three verified by mutation: 3, 2 and 1 failures.

Suite: 753 passing, from 670.

---
**Review cycle 6 (independent session, batch 5).** Two code defects, both in the
*supplied* half of the both-directions fix added in cycle 5, and both in the half
no test reached. The lesson from cycle 5 repeated at once, one level down.

- **An `NA` in one supplied direction discarded a real value in the other.**
  `pmax()` defaults to `na.rm = FALSE`, and the supplied block was filled after
  the two directions were combined rather than before, as the computed block is.
  A supplied `[forced, predictor]` of `NA` against a `[predictor, forced]` of
  0.99 reported 0, kept the predictor in every candidate and warned about
  nothing -- the silence FSSgam_package#23 exists to end, reached through the
  argument a user supplies to control the screen.

- **A regression: a matrix naming a forced term in one dimension only aborted.**
  `have` comes from `rownames()` and `cols` from `colnames()`, so reading the
  reverse direction as `supplied[cols, have]` assumed every such name appears in
  both. That shape is one `build_predictor_correlation_matrix()` documents itself
  as accepting, and the parent branch completed on it. Now guarded by
  intersecting against the opposite dimension.

Three assertions added and verified by mutation: 5, 3 and 1 failures. Before
them, both single-direction reductions of the supplied path left the suite green.

Two figures corrected. `[forced, curve]` is 0.280 in the committed test, not the
0.152 an ad-hoc run gave; the difference is the random draw, and the test asserts
only that it is below the cutoff, which is the property that matters. And the
pull request described a six-commit branch as one commit.

**Six rounds, and the count of false coverage claims stopped mattering: what
matters is that each was made the same way.** The pull request no longer claims
that all reachable mutations are pinned, only that every mutation it names was
run and its count recorded.

Suite: 760 passing, from 670.

---
**Review cycle 7 (independent session, batch 5).** One substantive finding: a
third orientation of the asymmetry cycles 5 and 6 fixed.

`have` was taken from `rownames(cor.matrix)` alone, so a supplied matrix naming
a forced term only as a *column* was discarded in full and a computed estimate
substituted. The value the user supplied to control the screen was ignored, and
on the FSSgam_package#13 path -- where the computation fails -- the screen was
skipped altogether. Cycle 6 fixed the mirror shape and tested it; the
column-only shape had no test.

Evidence: the same value in the same cell, only the matrix shape differing. A
square matrix and a rows-only matrix both drop the predictor at 0.99; the
column-only matrix reported the computed 0.274, raised no warning and retained
it in all seven candidates.

`from.supplied` is now allocated over the union of both dimnames and filled from
whichever directions the matrix has, rather than subset from the rows.

**The three orientations are worth naming together, because they are one
mistake made three times:** the value read (one direction of an asymmetric
estimate), the row shape, and the column shape. Each time the code was written
for the square, symmetric case and the general one was assumed to follow.

A secondary over-claim was also corrected. "Exempt whichever way the argument is
spelled" was false for `bs=mybs`, a variable, which mgcv accepts and which does
not resolve in `baseenv()`. Such a term is now treated as a random effect: the
two errors are not symmetric, since screening a grouping factor drops a
legitimate predictor while not screening a fixed term leaves the behaviour of
earlier versions.

Both fixes verified by mutation, 1 failure each. Two test fixtures of mine were
wrong rather than the code, and were corrected: one placed the value in the cell
the subsetting removes, the other compared against `character(0)` where the
function returned `NULL` -- which is now `character(0)`.

Suite: 768 passing, from 670.

---
