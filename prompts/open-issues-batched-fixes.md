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
