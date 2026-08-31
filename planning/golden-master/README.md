# Golden-master scripts for the pre-CRAN refactor

These are the verification scripts for the phase 3 restructure and the phase 4
behaviour changes. They are kept here, rather than in a scratch directory, so
that a later session can reproduce the comparison. `planning/` is in
`.Rbuildignore`, so nothing here reaches the built package.

## What each script does

- **`capture.R`** — runs 31 `generate_model_set()` scenarios and saves the
  result to an `.rds`. Each scenario records `n.mods`, the candidate names,
  every deparsed formula, `predictor.correlations`, the names and dimensions of
  `used.data`, `included.vars`, and every warning or error raised. Errors and
  warnings are stored as values, so a scenario that currently fails must keep
  failing in the same way.
- **`compare.R`** — compares two such files with `identical()` and prints the
  differing field of any scenario that moved.
- **`capture-correlations.R`** — the equivalent for `check_correlations()` and
  `check_non_linear_correlations()`, 14 scenarios. Used for phase 2.

## Running them

From the repository root:

```
Rscript planning/golden-master/capture.R /tmp/gm_before.rds
# ... make the change ...
Rscript planning/golden-master/capture.R /tmp/gm_after.rds
Rscript planning/golden-master/compare.R /tmp/gm_before.rds /tmp/gm_after.rds
```

The reference used through phase 3 was captured on `master` at `a8235e2` and
matched at every step up to and including `0b5d991`. To regenerate it, check out
that commit, run `capture.R`, and check the branch back out.

## Scenarios that exist for a specific reason

Four were added after measuring the baseline, because the obvious scenarios did
not discriminate. Do not delete them without reading this.

- **`ssi_true_3way` / `ssi_char_3way`** use `cov.cutoff = 0.4` rather than the
  default 0.28. At the default, `smooth.smooth.interactions = TRUE` and the
  character form produce byte-identical output at `max.predictors = 3`, because
  `depth` and `complexity` correlate at 0.3336 and screening that pair out
  screens the three-way combination containing it out as well. At 0.4 all three
  pairs are admitted and the character form produces
  `depth.te.SCORE2.te.complexity` while `TRUE` does not. These are the only
  scenarios that demonstrate the arity difference between the two branches.

- **`ffi_asym_screen` / `ffi_char_asym_screen`** use `cov.cutoff = 0.543` and
  three factors. `check_correlations()` fits `multinom()` separately in each
  direction, so its factor block is asymmetric — by at most 8.5e-4 over these
  columns. `ZONE3`/`ZONE5` measure 0.5425 one way and 0.5434 the other, so the
  cutoff falls between them, and the logical branch (which screens both
  triangles) drops `ZONE3.I.ZONE5` while the character branch (upper only)
  keeps it: 15 candidates against 16. The third factor is required — with only
  two, the character branch is skipped by a separate guard and both branches
  produce nothing for unrelated reasons.

- **`ffi_char_cellcount_guard`** pins that separate guard,
  `length(which(factor.correlations < cov.cutoff)) > 1`, which counts matrix
  cells rather than pairs.

## Scenarios that will change in phase 4

These are expected to differ, and each difference must be inspected and recorded
rather than compared for identity:

- `max_pred_1` — A3. Currently produces duplicated `used.data` columns and four
  "no non-missing arguments to max" warnings.
- `max_pred_1_ssi`, `single_pred`, `single_fact` — same defect, other shapes.
- `supplied_cor_matrix`, `supplied_cor_matrix_ffi` — FSSgam_package#13.
  `supplied_cor_matrix_ffi` currently errors.
- `ffi_asym_screen` / `ffi_char_asym_screen` — when the triangle-screening
  divergence is resolved.
- every scenario, if the candidate-name sort order changes (FSSgam_package#8).
  Under `sort(method = "radix")` the names should not move, because testthat
  already forces `LC_COLLATE=C`; that is the thing to verify rather than assume.
