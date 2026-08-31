# Pre-CRAN refactor — where this stands

Branch `pre-cran-refactor`, off `master` at `a8235e2`. Last session
2026-08-31/09-01. Nothing is pushed; nothing is merged.

## Done

Phases 1, 2 and 3 of `pre-cran-refactor-human.md`. Eighteen commits, one per
item, each with its regression test in the same commit.

| Phase | Items | State |
|---|---|---|
| 1 | A2, A4, B1, B2, B3, B4, FSSgam_package#7, #9 | complete |
| 2 | FSSgam_package#10, #12, redundant null refits | complete |
| 3 | restructure interaction resolution, no behaviour change | complete |

Suite: 508 passing, 0 failing, 12 skipped (baseline was 438/0/11).
`devtools::check()` after phase 1: 0 errors, 0 warnings, 0 notes. Not re-run
after phases 2 and 3 — do that first on resuming.

## Next

Phase 4, then 5, 6 and 7. Specification sections 6, 7, 8 and 9 of
`pre-cran-refactor-claude.md`.

Phase 4 applies A1 (cyclic variables matched by unanchored `grep()`), A3
(`for(i in 2:n)` counting backwards at n = 1), FSSgam_package#8 (locale-dependent
candidate names) and #13 (a supplied `cor.matrix` governing all three stages).
It also resolves the two branch divergences that phase 3 deliberately held at
their current behaviour behind arguments:

- `screen.both` in `combine_uncorrelated()` — resolve in favour of screening
  both triangles, then remove the argument;
- `warn.when.empty` in `build_factor_interaction_columns()` — resolve in favour
  of warning, then remove the argument. The `FALSE` path currently reaches
  `cbind(use.dat, <0-column data.frame>)` and raises "arguments imply differing
  number of rows".

Two further divergences found during phase 3 and not yet decided:

- the character `factor.factor.interactions` branch has an outer guard,
  `length(which(factor.correlations < cov.cutoff)) > 1`, that counts matrix
  cells rather than pairs and silently skips the whole block. The logical
  branch has no counterpart. Decide whether to keep, fix or remove it.
- the `te()` arity difference is retained by decision (2026-08-31) and needs
  only the documentation rewrite in phase 7, specification 6.3.

**Before changing anything in phase 4**, re-capture the golden master on the
current branch tip and keep it as the reference:

```
Rscript planning/golden-master/capture.R /tmp/gm_before.rds
```

`planning/golden-master/README.md` lists which scenarios are expected to change
and why, and why four of them exist at all.

## Open, needing a decision from the maintainer

**The `gamm4` worker stall.** `fit_model_set(parallel = TRUE)` stalls on roughly
half of runs on this host, on `master` as well as on this branch. Measured, all
n = 10 or 12: a `foreach()` with a trivial body and no FSSgam code reproduces
it purely from `.packages` — `mgcv` 10/10, `mgcv, MuMIn` 10/10, `gamm4` 7/10,
`mgcv, gamm4` 3/10, `mgcv, FSSgam` 6/10, all four 8/10. Every set that loads
gamm4 stalls; no set without it does, 20/20 against 24/40. Instrumenting a
stall shows `makeCluster()` and `registerDoSNOW()` both complete and the stall
is inside the `%dopar%` dispatch.

This is not a defect in this package, it predates the refactor, and it is
outside the plan's scope. Nothing was changed. Two follow-ups were proposed and
are awaiting a decision:

1. correct `CLAUDE.md`, whose Phase 13 caution states that the stall is in
   cluster startup, that a plain script completed six times out of six, and
   that it is contention from another process — all three are now known to be
   wrong. A drafted correction with the measurements is in the last session's
   prompt log.
2. raise it as an issue in `FSSgam_package`.

Not established: whether this is specific to WSL2 or to this lme4/Matrix build,
and the mechanism.

**One environment note.** roxygen2 here is 8.0.0, and `devtools::document()`
rewrites `DESCRIPTION`'s `RoxygenNote: 7.3.2` to
`Config/roxygen2/version: 8.0.0`. That line was reverted each time, since it
reflects this machine rather than a decision about the package. The generated
`.Rd` files show no structural difference from the version change. Decide
whether to adopt roxygen2 8 deliberately.

## Working notes for whoever resumes

- Four files in `R/` use CRLF line endings — `data.R`, `FSSgam.R`,
  `full-subsets-gam.R`, `functions_supporting.R`. An editor that rewrites them
  as LF turns a three-line change into a whole-file diff. Check
  `git diff --numstat` after editing those four.
- The parallel tests need `R CMD INSTALL .` first; a `pkgload` copy tests
  nothing, because a doSNOW worker loads the installed package. Expect to
  retry: see the stall above.
- Every measurement quoted in the prompt log was taken on this WSL host under
  R 4.6.1. Re-measure rather than comparing against those numbers.
