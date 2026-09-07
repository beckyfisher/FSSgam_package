# Pre-CRAN refactor — specification

Implementation specification for the plan in `pre-cran-refactor-human.md`.
Written for a Claude Code session to implement from. Line references are
against commit `a8235e2` on `master`.

---

## 1. Baseline

Measured 2026-08-31, this WSL host, R 4.6.1, commit `a8235e2`:

- `devtools::check(document = FALSE, args = c("--no-manual", "--no-build-vignettes"))`
  → 0 errors, 0 warnings, 0 notes, 51 s.
- `R CMD check --as-cran --no-manual` on the built tarball with
  `_R_CHECK_CRAN_INCOMING_=TRUE _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`
  → Status: OK.
- `testthat::test_local()` → 447 passed, 0 failed, 6 skipped blocks.

  Corrected on 2026-09-01. This was recorded as 438 passed, 11 skipped, which
  does not reproduce on a clean checkout. Regenerating
  `tests/testthat/_snaps/numerical-regression.md` skips five further blocks and
  the nine expectations inside them, which accounts for the difference exactly.
  Measure with the snapshot file present.

Re-measure before starting; these move.

**Counting expectations.** Phase 13's rule stands. Quote testthat's own `PASS`
figure, `sum(as.data.frame(res)$passed)`, not `sum(...$nb)`, which includes
skips. The figure appears in `NEWS.md`, `CLAUDE.md` and the pull request body:
update all three together or not at all.

## 2. Working arrangements

- Branch from `master`; do not commit to `master` directly.
- One commit per item. Phase 13 deviated from this by folding pairs of fixes
  into single commits and recorded it as a deviation not to repeat.
- Each behaviour fix carries a regression test in the same commit.
- `planning/` must be added to `.Rbuildignore` as `^planning$`, or
  `R CMD check` raises a "non-standard file/directory found at top level"
  note. Do this in the first commit.
- Prompt logging: repo `CLAUDE.md` Section 7 requires a log once the package
  itself changes. Open `prompts/pre-cran-refactor.md` at phase 1, not before.

---

## 3. Phase 1 — defects and validation outside the interaction code

### 3.1 A2 — `fit_model_set()` fails on a single-predictor model set

`R/fit-model-set.R:272`, `:291`, `:295`, inside `compute_variable_importance()`.

`mod.data.out[, included.vars]` drops to a vector when `included.vars` has
length 1, and `colSums()` errors.

Reproduced:

```r
tf <- gam(log.Herbivore.biomass ~ s(depth, k = 3, bs = "cr"), data = d)
ms <- generate_model_set(use.dat = d, test.fit = tf, pred.vars.cont = "depth",
                         max.predictors = 1, k = 3)
fit_model_set(ms)
#> Error in colSums(mod.data.out[, included.vars]) :
#>   'x' must be an array of at least two dimensions
```

Fix: `drop = FALSE` at all three sites.

```r
min.mods <- min(colSums(mod.data.out[, included.vars, drop = FALSE]))
variable.weights.raw <- colSums(mod.data.out[, included.vars, drop = FALSE] * mod.data.out$wi.AICc)
variable.weights.raw <- colSums(mod.data.out[, included.vars, drop = FALSE] * mod.data.out$wi.BIC)
```

Edge case to fix in the same commit: when a predictor appears in no surviving
model, `min.mods` is 0 and `sort(...)[1:0]` returns the first element rather
than none, so its weight is counted. Guard with `seq_len(min.mods)`, which is
empty at 0 and yields `sum(numeric(0)) == 0`.

Tests: single continuous predictor at `max.predictors = 1`, both `VI.mods`
values, asserting the returned importance vector is named and length 1.

### 3.2 A4 — `.errorhandling` disabled on the unsaved parallel path

`R/fit-model-set.R:218`. `.errorhandling='pass'` is commented out;
`R/fit-model-set.R:154` has it active on the saved path.

Fix: uncomment. Then handle the consequence, which is the reason it is not a
one-character change. With `.errorhandling='pass'`, a failing element of the
`foreach()` is the condition object, not the named numeric vector that
`unlist(extract_mod_dat(...))` returns, and `do.call("rbind", mod.dat)` on a
mixed list produces a malformed table.

After the `foreach()`, replace any element that is not a numeric vector of the
expected length with the all-`NA` template `extract_mod_dat()` returns for a
`try-error`:

```r
na.template <- unlist(list(AICc = NA, BIC = NA, r2.vals = NA,
                           r2.vals.unique = NA, edf = NA, edf.less.1 = NA))
bad <- !vapply(mod.dat, function(x) is.numeric(x) && length(x) == length(na.template),
               logical(1))
mod.dat[bad] <- list(na.template)
```

`failed.models` is then derived from `is.na(mod.data.out$AICc)` as it already
is, so a failed candidate lands in `failed.models` on both paths.

Note `fit_mod_l()` returns a `try-error` rather than throwing, and
`extract_mod_dat()` handles that class, so the `foreach()` body rarely throws.
`.errorhandling='pass'` is a backstop for anything that gets past both.

Testing: `break_one_candidate()` (`tests/testthat/helper-fixtures.R`) injects an
unfittable formula. The parallel assertion must sit behind
`skip_unless_parallel_opt_in()`, per Phase 13. Add a sequential test of the
post-processing in the same commit so `covr` reaches it.

### 3.3 B1 — `r2.type` accepts anything

Reproduced: `fit_model_set(ms, r2.type = "rsq")` returns `r2.vals` of `NA NA NA`
with no message, because `extract_mod_dat()` matches `r2.type.` against three
literals and leaves `tempOut` at `NA` when none match.

Fix: validate in `fit_model_set()` and `full_subsets_gam()`:

```r
r2.type <- match.arg(r2.type, c("r2.lm.est", "r2", "dev"))
```

Do **not** add this to `extract_mod_dat()`. It stays exported (section 8), its
`r2.type.` argument is documented as being passed through, and adding
validation there would make the exported function stricter than it is now.

Note `"r2"` is an exact match against the second element, so `match.arg()`'s
partial matching does not make it ambiguous with `"r2.lm.est"`.

### 3.4 B2 — `VI.mods` accepts anything

Reproduced: `fit_model_set(ms, VI.mods = "alll")` fails with
`object 'aic.var.weights' not found`, because `compute_variable_importance()`
has two `if` blocks and no `else`.

Fix: `VI.mods <- match.arg(VI.mods, c("min.n", "all"))` in `fit_model_set()`,
and change the second `if` in `compute_variable_importance()` to `else`.

### 3.5 B3 — the deprecated argument checks fail on valid input

`R/full-subsets-gam.R:137-161`. Three arguments use the sentinel default
`"previous.arg"` and are tested with `!=` or `is.na()`.

Reproduced:

```r
full_subsets_gam(..., smooth.interactions = c("ZONE", "Z2"))
#> Error in if (is.na(smooth.interactions) == TRUE) { :
#>   the condition has length > 1
full_subsets_gam(..., factor.interactions = NA)
#> Error in if (factor.interactions != "previous.arg") { :
#>   missing value where TRUE/FALSE needed
```

Fix: drop the sentinel and use `missing()`, which is correct for every type,
length and `NA`.

```r
full_subsets_gam <- function(..., factor.interactions, smooth.interactions, size) {
  if (!missing(factor.interactions)) {
    factor.factor.interactions <- factor.interactions
    warning(...)
  }
  if (!missing(smooth.interactions)) {
    factor.smooth.interactions <- smooth.interactions
    warning(...)
  }
  if (!missing(size)) {
    max.predictors <- size
    warning(...)
  }
```

The formals lose their defaults. `full.subsets.gam()` forwards through `...`,
so `missing()` still reports correctly through the deprecated alias; assert
that in a test.

This also removes the existing special case where `smooth.interactions = NA`
was handled by a separate branch. Under `missing()`, `NA` is simply assigned
to `factor.smooth.interactions`, which is its documented meaning ("If
specified as NA no factor-continuous predictor interactions will be
included"). Confirm the warning still fires for `NA`, which it did not
reliably before.

### 3.6 B4 — `max.models` defaults disagree

`fit_model_set()` defaults to 200 (`R/fit-model-set.R:75`);
`full_subsets_gam()` to 500 (`R/full-subsets-gam.R:130`). Undocumented in
either. A candidate set of 300 therefore saves its fits through the wrapper and
does not through the pair.

Fix: set both to 200. Record in `NEWS.md` as a behaviour change, since a user
with a set between 201 and 500 calling `full_subsets_gam()` will now get the
warning and lose `success.models`.

### 3.7 Issue #7 — `full_subsets_gam()` does not forward `VI.mods`

Add `VI.mods = 'min.n'` to the formals, forward it to `fit_model_set()`, and
copy the `@param` entry. Update the test in
`tests/testthat/test-full_subsets_gam.R` named "full_subsets_gam has no
VI.mods argument and always uses the min.n default" — it asserts the current
absence and must be replaced by a test that the argument is forwarded.

### 3.8 Issue #9 — the progress bar cannot be suppressed

Add `progress = interactive()` to `fit_model_set()` and `full_subsets_gam()`,
forwarded through. Guard the three `txtProgressBar`/`setTxtProgressBar`/`close`
call sites in each of `fit_and_summarise_saved_models()` and
`fit_and_summarise_unsaved_models()`, and the `.options.snow = opts` argument,
which must be `list()` when `progress` is `FALSE`.

Also fixes the double `close(pb)` in the parallel branch of
`fit_and_summarise_saved_models()` (`R/fit-model-set.R:159` and `:165`) — the
bar is closed inside the `if` and again after it.

`fit_quietly()`/`full_subsets_quietly()` in `helper-fixtures.R` can stay as
they are; they wrap `capture.output()` and remain correct with the bar off.
Phase 13 records that four groups of tests deliberately keep their direct call
form — do not change those.

---

## 4. Phase 2 — the correlation functions

### 4.1 Issue #10 — spurious "NaNs produced" warnings

`R/check-correlations.R:117` and `:131`. `summary.multinom()` computes the
variance-covariance matrix to report standard errors; on a perfectly separated
fit its diagonal is negative and `sqrt()` warns. Only `$deviance` is used.

Fix in both the parallel and sequential branches:

```r
fit <- try(nnet::multinom(dat.r[, var.1] ~ dat.r[, var.2], trace = FALSE)$deviance, silent = TRUE)
```

and the same for `null.fit`.

`suppress_nnet_nans()` in the test helpers tolerates the warning being absent,
so no test breaks. Leave the helper in place.

### 4.2 Issue #12 — the factor-factor diagonal is not exactly 1

Two changes in `fill_factor_factor_correlations()`:

1. Skip the self-pairs in the `lm.grid` loop
   (`lm.grid <- lm.grid[lm.grid$fact.var1 != lm.grid$fact.var2, ]`), which also
   removes `length(fact.vars)` redundant `multinom()` fits per call.
2. `diag(out.cor.mat) <- 1` after the block is filled, matching
   `assemble_non_linear_correlation_matrix()`.

Invariant: the loop currently writes all `n²` ordered pairs and therefore fills
every cell. Removing the `n` diagonal pairs leaves exactly the diagonal
unwritten, which step 2 then sets. No off-diagonal cell changes.

`tests/testthat/test-check_correlations.R` asserts `all(diag(cm) > 0.99)`;
change to `expect_equal(diag(cm), rep(1, n))`.

### 4.3 C4 — the null multinomial is refitted once per pair

`R/check-correlations.R:118`, `:132`. `null.fit` depends only on `var.1`, but is
refitted for every ordered pair. Hoist to one fit per factor, computed once
before the loop and looked up by name.

Note the existing `null.fit` uses `dat[, var.1]` while `fit` uses
`dat.r[, var.1]`, where `dat.r` has had incomplete cases removed for that pair.
The deviances are therefore computed on different row sets when the pair has
missing values. Hoisting makes the null fit unambiguously per-variable on the
full column. Record this in `NEWS.md`: it changes `r.est` for factor pairs with
missing data. `generate_model_set()` rejects predictors containing `NA`
(`check_predictor_missingness()`), so this is reachable only through a direct
call to `check_correlations()`.

---

## 5. Phase 3 — restructure interaction resolution (behaviour-preserving)

### 5.1 What is wrong now

`resolve_factor_interactions()` (`R/generate-model-set.R:243-403`, 161 lines)
does four jobs, two of them implemented twice:

| Job | Implementations | Divergence |
|---|---|---|
| Build factor-factor interaction columns | `logical` branch (255-291), `character` branch (293-320) | logical screens `upper.tri` and `lower.tri`; character screens only `upper.tri`. Logical warns when nothing survives; character is silent. |
| Coerce factor predictors to `factor` | one (323-324) | — |
| Resolve `.by.` factor-smooth terms | `character` branch (328-354), `list` branch (356-393) | list form validates against `all.predictors`; character form against `pred.vars.fact` |
| Resolve `.t.` factor-linear terms | inside both of the above | — |

`resolve_smooth_smooth_interactions()` (`:404-476`, 73 lines) has the same
two-copy structure, and both copies build `tt` (`:433`, `:465`) — a data frame
of every surviving pair of columns pasted together — and never read it.

The same "enumerate combinations, drop the too-correlated ones" block appears
three times, at `:262-274`, `:300-312` and `:420-431`.

### 5.2 Target structure

New unexported helpers in `R/generate-model-set.R`:

```r
# Every combination of `vars` of size in `sizes`, dropping any whose
# sub-matrix of `cor.matrix` exceeds `cov.cutoff` off the diagonal.
# Returns a list of character vectors, possibly empty.
combine_uncorrelated <- function(vars, sizes, cor.matrix, cov.cutoff)

# Normalises factor.factor.interactions (TRUE / character / FALSE) to the
# character vector of factors to combine, or character(0).
normalise_factor_factor_interactions <- function(factor.factor.interactions,
                                                 pred.vars.fact, use.dat)

# Builds the pasted interaction columns and appends them to use.dat.
# Returns list(use.dat =, interaction.terms =).
build_factor_interaction_columns <- function(use.dat, fact.vars, max.predictors,
                                             cov.cutoff, cor.matrix)

# Normalises factor.smooth.interactions (character / list / NA) to
# list(fact.vars =, cont.vars =, linear.vars =).
normalise_factor_smooth_interactions <- function(factor.smooth.interactions,
                                                 pred.vars.fact, pred.vars.cont,
                                                 linear.vars, all.predictors)

# Builds the ".by." and ".t." term names from that triple.
# Returns list(interaction.terms =, linear.interaction.terms =).
build_factor_smooth_terms <- function(spec, pred.vars.fact)
```

`resolve_factor_interactions()` becomes a short sequence of calls to these.
`resolve_smooth_smooth_interactions()` becomes one call to
`normalise_*` plus one to `combine_uncorrelated()`, with the dead `tt` removed.

Both divergences in the table above must be resolved deliberately, not
absorbed. Specify: screen **both** `upper.tri` and `lower.tri` (the stricter
behaviour, and correct for the asymmetric matrix
`check_non_linear_correlations()` returns), and **warn** when nothing survives.
The character branch therefore becomes stricter. That is a behaviour change,
so it belongs in phase 4, not here — implement the shared helper in phase 3
with a `screen.both` argument preserving each caller's current behaviour, and
remove the argument in phase 4.

### 5.3 Verification

Golden master, following Phases 6 and 6b. Capture `generate_model_set()`
output — `n.mods`, `names(mod.formula)`, every `deparse()`d formula,
`predictor.correlations`, `names(used.data)`, `included.vars` — for at least
these scenarios, before the restructure and after, and compare with
`identical()`:

1. continuous only; 2. factors only; 3. mixed; 4. `factor.factor.interactions
= TRUE`; 5. `factor.factor.interactions` as a character vector; 6.
`smooth.smooth.interactions = TRUE`; 7. `smooth.smooth.interactions` as a
character vector; 8. `factor.smooth.interactions` as a list; 9.
`factor.smooth.interactions = NA`; 10. `cyclic.vars`; 11. `linear.vars`; 12.
`linear.vars` with a list-form `factor.smooth.interactions`; 13. a supplied
`cor.matrix`; 14. `non.linear.correlations = TRUE`; 15. `max.predictors = 1`;
16. a single predictor.

Scenarios 15 and 16 will differ after phase 4 and must be captured again at
that point. Use `deparse1()`, not `deparse()` — the latter's default
`width.cutoff` of 60 wraps long formulas and a wrapped formula deparses to a
different string. This said `deparse_one()`, a local helper in
`helper-fixtures.R` that existed only because `deparse1()` postdated the
declared R floor. The floor was raised and the helper deleted
(FSSgam_package#31); `planning/golden-master/capture.R` still defines its own
copy and is left alone, being a self-contained record of what was run.

---

## 6. Phase 4 — behaviour changes

### 6.1 A1 — cyclic variables matched by unanchored `grep()`

`R/generate-model-set.R:620-634`. After all terms are assembled as strings, a
loop rewrites `bs=` for any term matching `grep(cyclic.vars[r], term)`.

Reproduced, `pred.vars.cont = c("depth", "depthx", "complexity")`,
`cyclic.vars = "depth"`:

```
depth        ~s(depth, k = 4, bs = "cc") + s(site, bs = "re")
depthx       ~s(depthx, k = 4, bs = "cc") + s(site, bs = "re")   # wrong
complexity   ~s(complexity, k = 4, bs = "cr") + s(site, bs = "re")
```

Two defects: the match is a substring match, and the pattern is a regular
expression, so a predictor name containing `.` matches any character.

Fix by choosing the basis when each term is constructed, from the variable
names in `mod.m`, which are known exactly — not by rewriting assembled strings
afterwards. Delete the whole post-hoc rewrite block, including the second block
at `:627-634` that repairs `te()` terms the first block damaged.

```r
bs_for <- function(v) if (v %in% stats::na.omit(cyclic.vars)) "'cc'" else bs.arg

# single smooths
paste0("s(", cont.smooths, ",k=", k, ",bs=", vapply(cont.smooths, bs_for, character(1)), ")")

# by smooths: the smooth variable is the part before ".by."
sm.var <- sub("\\.by\\..*$", "", by.smooths)
paste0("s(", gsub(".by.", ",by=", by.smooths, fixed = TRUE), ",k=", k,
       ",bs=", vapply(sm.var, bs_for, character(1)), ")")

# te/t2: one bs per variable
te.vars <- strsplit(te.smooths, ".te.", fixed = TRUE)
bs.vec  <- vapply(te.vars, function(v) paste0("bs=c(", paste(vapply(v, bs_for, character(1)),
                                                             collapse = ","), ")"),
                  character(1))
```

Invariant: for any model set where no predictor name is a substring of another
and no name contains a regular-expression metacharacter, the formulas are
byte-identical before and after. Assert that over the golden-master scenarios.

Tests: overlapping names (`depth`/`depthx`); a name containing `.`; a cyclic
variable inside a `te()` pair, which is the case the deleted second block
existed to repair.

### 6.2 A3 — `for(i in 2:n)` counts backwards at `n = 1`

`R/generate-model-set.R:262`, `:300`, `:420`, `:452`.

Reproduced, two factors, `factor.factor.interactions = TRUE`,
`max.predictors = 1`: `used.data` grows from 28 to 31 columns with `ZONE` and
`Z2` duplicated, and four warnings are emitted:

```
In max(abs(cor.mat.m[upper.tri(cor.mat.m)])) :
  no non-missing arguments to max; returning -Inf
```

The `i = 1` pass builds single-variable "interactions" whose pasted names are
just the original variable names, and `cbind()` appends them again.

Fix: `combine_uncorrelated()` takes `sizes` and computes it as
`if (max.size < 2) integer(0) else 2:max.size`, so the loop does not run.
The `-Inf` warnings disappear with it, since they came from taking `max()` of
the empty off-diagonal of a 1×1 sub-matrix.

Test: assert `ncol(used.data)` is unchanged and `anyDuplicated(names(used.data))`
is 0 at `max.predictors = 1`, and that no warning is emitted.

### 6.3 C2 — `te()` arity differs between the two branches

**Decision (maintainer, 2026-08-31): retain the current behaviour; correct the
documentation to describe it.**

`R/generate-model-set.R:418` reads `cont.cmbns.max.predictors=2#max.predictors`.
`smooth.smooth.interactions = TRUE` therefore builds only pairwise `te()`
terms, while a character vector builds combinations of size 2 up to
`max.predictors`. `?generate_model_set` describes both as "bivariate calls to
te", which is accurate for the `TRUE` form and wrong for the character form.

The recommendation to fix both branches at pairwise was rejected: it removes a
capability that is reachable today, and no user has asked for it to go. The
code is therefore unchanged and the documentation is corrected instead.

**This item is no longer a behaviour change, so it leaves phase 4.** Two pieces
of work remain, in two different phases:

1. **Phase 3.** The literal `2` survives the restructure as the `sizes`
   argument to `combine_uncorrelated()` for the `TRUE` branch, and
   `2:max.predictors` for the character branch. Replace the trailing
   `#max.predictors` comment, which reads as an unfinished edit, with a comment
   stating that the arity difference is intentional and pointing at the
   documented behaviour. Do not "restore" `max.predictors` here.
2. **Phase 7.** Rewrite the `smooth.smooth.interactions` `@param` entry to give
   both forms separately: `TRUE` builds bivariate `te()` terms over all
   continuous predictors; a character vector builds `te()` terms over every
   combination of the named predictors of size 2 up to `max.predictors`, so a
   three-way `te(a, b, c)` is produced when `max.predictors >= 3`.

No `NEWS.md` entry: nothing about the package's behaviour changes.

Note the interaction with A3 (section 6.2). The `TRUE` branch's hard-coded 2
means its loop was never the backwards `2:1` case, so only the character branch
at `R/generate-model-set.R:452` is affected there. Retaining the hard-coded 2
also means `smooth.smooth.interactions = TRUE` with `max.predictors = 1` still
constructs pairwise `te()` terms, which `enumerate_candidate_models()` then
excludes for exceeding `max.predictors`. That is current behaviour and stays.

### 6.4 Issue #13 — a supplied `cor.matrix` must govern `te()` construction

**Decision (maintainer, 2026-08-31): one matrix governs both stages.**

Reorder `generate_model_set()` so `build_predictor_correlation_matrix()` runs
before `resolve_smooth_smooth_interactions()`, and pass the resolved matrix in:

```r
factors.l   <- resolve_factor_interactions(...)          # unchanged position
all.predictors <- stats::na.omit(unique(c(all.predictors, pred.vars.fact)))
cor.matrix  <- build_predictor_correlation_matrix(...)   # moved up
smooth.smooth.interaction.terms <- resolve_smooth_smooth_interactions(
                                     cor.matrix = cor.matrix, ...)
use.mods    <- enumerate_candidate_models(cor.matrix = cor.matrix, ...)
```

The reorder is safe: `build_predictor_correlation_matrix()` needs
`all.predictors` including the factor-interaction columns, which
`resolve_factor_interactions()` has already created, and does **not** need the
`te()` term names.

**Verified 2026-08-31**, `case_study1`, `pred.vars.cont = c("depth", "SCORE2",
"complexity")` with `smooth.smooth.interactions = TRUE`: `cor.matrix` has
row names `depth, SCORE2, complexity` and contains no `.te.` name, while the
candidate set contains `depth.te.SCORE2`, `depth.te.complexity` and
`SCORE2.te.complexity`. `enumerate_candidate_models()` splits `.te.` terms into
their components before indexing the matrix.

`resolve_smooth_smooth_interactions()` then subsets `cor.matrix` to the
continuous predictors instead of calling `check_correlations()` itself.

**Invariant, verified 2026-08-31** on `case_study1` with
`pred.vars.cont = c("depth", "SCORE2", "complexity")` and one factor (`ZONE`)
added for the full set, for both `non.linear.correlations` settings:

```r
identical(f(use.dat[, pred.vars.cont]), f(use.dat[, all.predictors])[pred.vars.cont, pred.vars.cont])
#> TRUE   for f = check_correlations
#> TRUE   for f = check_non_linear_correlations
```

Maximum absolute difference 0 in both cases. The continuous-continuous block is
therefore unaffected by which other predictors are present, so subsetting the
full matrix is behaviour-preserving whenever no `cor.matrix` is supplied. This
was checked on one dataset with one factor; re-check if a scenario is added
with several factors or with a continuous predictor that has missing values,
since `cor()` uses `pairwise.complete.obs`.

**One known behaviour difference.** `build_predictor_correlation_matrix()`
replaces `NaN` and `NA` with 0 (`R/generate-model-set.R:489-491`); the local
call inside `resolve_smooth_smooth_interactions()` does not. A zero-variance
continuous predictor yields `NaN` from `cor()`, and
`max(abs(NaN)) > cov.cutoff` is `NA`, so the current code raises
`missing value where TRUE/FALSE needed`. After the change the cell is 0 and the
pair is admitted. This is an improvement, but it is a change: record it in
`NEWS.md` and add a test with a constant column.

### 6.5 Issue #8 — locale-dependent candidate model names

**Decision (maintainer, 2026-08-31): sort in byte order.**

`R/generate-model-set.R:570`:

```r
use.mods = unique(lapply(use.mods, FUN = sort, method = "radix"))
```

`method = "radix"` sorts in the C locale regardless of `LC_COLLATE`.

The committed tests already assert C-collation names, because testthat forces
`LC_COLLATE=C`, so no test changes. Verify explicitly by running the suite once
under `LC_COLLATE=en_US.UTF-8` with testthat's own setting overridden, and
confirming the names are unchanged — that is the only thing that actually
demonstrates the fix.

`NEWS.md` must state plainly that `modname` values and the row order of
`mod.data.out` change for users in non-C locales, and that saved model tables
matched on `modname` need regenerating.

Close issue #8 with the measurement, not just the commit reference.

---

## 7. Phase 5 — idioms

Behaviour-preserving. Suite must stay green with no snapshot change.

| ID | Pattern | Sites | Replacement |
|---|---|---|---|
| D1 | `if(length(cont.smooths>0))` | `R/generate-model-set.R:608,610,612,614,616,618,619` | `if(length(cont.smooths) > 0)` |
| D2 | `class(x)[1] == "y"` | 14 sites across `R/` | `inherits(x, "y")` |
| D3 | `cor.matrix[which(cor.matrix=="NaN")]=0` | `R/generate-model-set.R:489` | `cor.matrix[is.nan(cor.matrix)] <- 0` |
| D4 | `1:length(x)`, `1:nrow(x)` | 18 sites | `seq_along(x)`, `seq_len(nrow(x))` |
| D5 | `max(is.na(x))==1` | several | `anyNA(x)` |
| D6 | `wi()`'s element-by-element loop | `R/functions_supporting.R:28-31` | `exp(-0.5 * AICc.vals.adj)` |

D1 is safe because `length(x > 0)` equals `length(x)` for any atomic `x`, so
`if(length(x > 0))` and `if(length(x) > 0)` are already equivalent. Confirm by
running the golden master, not by argument alone.

D2 needs care at two sites. `validate_use_dat()` (`R/generate-model-set.R:189`)
uses `class(use.dat)[1] != "data.frame"` deliberately, to reject tibbles —
`inherits()` would accept them. Use
`!identical(class(use.dat), "data.frame")` and comment why. `build_model_formulas()`
tests `class(test.fit)[[1]] == "gamm4"` to choose between `te()` and `t2()`;
`gamm4` objects are plain lists with that class, so `inherits()` is equivalent,
but verify against a `uGamm(lme4 = TRUE)` fixture.

D6: `wi()` returns `NaN` when every input is `NA`, both before and after.
Leave that; changing it is a behaviour change with no reported need.

---

## 8. Phase 6 — unexport `fit_mod_l()`

**Decision (maintainer, 2026-08-31): unexport `fit_mod_l()`; keep `wi()`,
`extract_mod_dat()` and `build_inclusion_mat()` exported.**

Rationale: `fit_mod_l()` gained a `family.` argument in Phase 12 and its
internals are the most likely to move again; the other three are plausible to
call interactively.

Steps:

1. Replace `@export` with `@keywords internal` in `fit_mod_l()`'s roxygen block.
   Keep the documentation — it holds the beckyfisher/FSSgam#10/#12 reasoning.
   Do not use `@noRd`, which would discard it.
2. Re-run `devtools::document()`; confirm `export(fit_mod_l)` leaves `NAMESPACE`.
3. Remove `fit_mod_l` from `_pkgdown.yml`'s reference index. `pkgdown::check_pkgdown()`
   errors with "must be a known topic name" if it is left in, and with "topics
   not in any index" if an exported function is missing — run it.
4. The `@examples` block calls `fit_mod_l()` directly. Examples of an
   internal-keyword function still run under `R CMD check`, so qualify it as
   `FSSgam:::fit_mod_l()` or drop the block.
5. Tests in `tests/testthat/test-functions_supporting.R` call it directly.
   `test_dir(package = "FSSgam")` evaluates tests in the package namespace, so
   unqualified calls keep working — confirm rather than assume, and against an
   installed copy as well as `load_all()`.
6. `NEWS.md`: state that `fit_mod_l()` is no longer exported and remains
   available as `FSSgam:::fit_mod_l()`.

---

## 9. Phase 7 — documentation and DESCRIPTION

None of this is required by `R CMD check`, which passes today with all of it in
place, including with the local CRAN incoming checks enabled. It is published
on the pkgdown reference site, which is the reason to fix it.

1. `R/full-subsets-gam.R:97` reads "we summed the ?i values". A character was
   lost; it should be `wi`. This renders into the published HTML.
2. `R/data.R:20` contains a U+2011 non-breaking hyphen, which reaches
   `man/case_study2.Rd`. Replace with an ASCII hyphen.
3. Spelling, 23 occurrences: `calcualte`, `apporoximate` (×2), `multnomial`
   (×2), `assymetric`, `replationships`, `Pleasure` (for "Please"), `cuttoff`,
   `amonst` (×2), `acheived`, `charactervector` (×2), `square-route` (×4, for
   "square root"), `interrrogated`, `superceded`, `Institue`, `pacakge`,
   `fulls subsets`, `overriten`, `as all as` (for "as well as"), `TRUE of you`
   (×2, for "TRUE if you"). Two matter beyond tidiness: `square-route` sits
   inside the definition of what the correlation values are, and
   `check_non_linear_correlations()`'s `@note` tells users to increase
   `cor.cutoff`, which is not an argument — `cov.cutoff` is.
4. Fourteen `@param` blocks are duplicated verbatim between
   `generate_model_set()` and `full_subsets_gam()`. Replace the duplicates with
   `@inheritParams generate_model_set` so they cannot drift.
5. DESCRIPTION: quote software names ('mgcv', 'MuMIn', 'gamm4') and give the
   reference as `Fisher et al. (2018) <doi:10.1002/ece3.4134>`. Neither is
   flagged by the local check; both are conventions CRAN reviewers apply by
   hand.
6. Update `?generate_model_set`'s `cor.matrix` entry for the 6.4 decision, and
   `smooth.smooth.interactions` for 6.3.

---

## 10. Rejected alternatives

**10.1 Fixing A1 by anchoring the regular expression.** Replacing
`grep(cyclic.vars[r], term)` with an anchored, escaped pattern would fix the
substring and metacharacter defects but keeps the design of rewriting assembled
strings, which is what required the second repair block at `:627-634`.
Constructing the basis with the term removes both blocks.

**10.2 Fixing #13 by passing `use.dat` through unchanged and documenting the
limitation.** This was the second option in issue #13 and was rejected by the
maintainer: `?generate_model_set` already documents `cor.matrix` as governing
exclusion generally, so the documentation would have to be narrowed rather than
the code corrected.

**10.3 Fixing #8 by preserving the supplied predictor order.** More natural to
read, but changes names for every user in every locale rather than only in
non-C locales, and requires updating the numerical snapshots. Rejected by the
maintainer in favour of `method = "radix"`.

**10.4 Unexporting all four supporting functions.** Rejected: only
`fit_mod_l()` has an unstable signature, and `wi()` in particular is a
self-contained calculation a user could reasonably call.

**10.5 Simplifying `fit_mod_l()`'s family resolution.** Out of scope and
explicitly forbidden by repo `CLAUDE.md` Section 5. Both simplifications were
tried across this package's history and each regressed the issue the other
fixed. Do not touch `resolve_candidate_family()` or
`clone_independent_family()`. If phase 6 or any other phase appears to require
it, stop and ask.

---

## 11. Issue #13, third stage — factor-factor interaction columns

**Decision (maintainer, 2026-08-31): a supplied `cor.matrix` governs this stage
too. Option A below.**

Section 6.4 settles that a supplied `cor.matrix` governs `te()` term
construction as well as candidate exclusion. A third stage screens against
correlations it computes itself: `build_factor_interaction_columns()` calls
`check_correlations(use.dat[, pred.vars.fact])` to decide which factor-factor
interaction columns to build.

The documentation for `factor.factor.interactions` states that a supplied
`cor.matrix` "must include these hard coded interactions", which is now
consistent with the code rather than aspirational.

**Option A, adopted.** Use the supplied matrix when one is given; compute one
only when none is.

Implementation, in `build_factor_interaction_columns()`:

```r
# `cor.matrix` here is the user-supplied matrix or NA, not the resolved one:
# the resolved matrix cannot exist yet, because it is built over the
# interaction columns this function creates.
fact.cor <- if (!identical(cor.matrix, NA)) {
  cor.matrix[fact.vars, fact.vars, drop = FALSE]
} else {
  check_correlations(use.dat[, fact.vars, drop = FALSE])$r.est
}
```

Three points of care.

1. **Ordering.** The screen must use the *supplied* matrix directly, not the
   resolved one from `build_predictor_correlation_matrix()`. The resolved matrix
   is built over `all.predictors`, which includes the interaction columns this
   function is deciding whether to create, so it cannot exist yet. This is why
   the argument is named `cor.matrix` in section 5.2's signature and carries the
   raw user value.
2. **Missing factors.** A user-supplied matrix that omits one of
   `pred.vars.fact` currently reaches the collinearity screen at
   `enumerate_candidate_models()` and errors there on the subscript. Reaching it
   earlier changes only where the error is raised, not whether. Validate once,
   in `generate_model_set()`, that every name in `pred.vars.fact` and
   `pred.vars.cont` is present in a supplied `cor.matrix`, and raise a message
   naming the missing ones. Test it.
3. **`NEWS.md`.** This is a behaviour change for anyone supplying a
   `cor.matrix` together with `factor.factor.interactions`: the set of
   interaction columns built can differ from what the data would have given.
   Record it in the same entry as 6.4, since it is the same decision.

**Option B, rejected.** Leaving it to compute its own matrix, and documenting
that `cor.matrix` governs `te()` construction and candidate exclusion but not
factor-factor column construction. Rejected for the reason in 10.2: it narrows
the documentation to match the code rather than correcting the code, and it
leaves one argument governing two of three stages, which is harder to explain
than either extreme.

Golden-master consequence: scenario 13 of section 5.3 (a supplied `cor.matrix`)
must be extended with a case that also sets `factor.factor.interactions`, and
that scenario will differ after phase 4. Capture it again at that point,
alongside scenarios 15 and 16.
