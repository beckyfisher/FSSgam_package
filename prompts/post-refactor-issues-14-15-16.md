## Session: post-refactor issues 14, 15 and 16
Date: 2026-09-01
Model: Claude Opus 5 (claude-opus-5)

### Prompts and Responses

**User:** We just finished a major refactor. This raised three new issues. Work
through each issue sequentially, in a loop: address the problem and post a PR,
then have an independent session draft a review comment, address what it raises,
draft a summary comment, and have a fresh session review that. Repeat for up to
10 cycles or until the review finds nothing substantive. Genuine issues outside
the scope of the one being worked on are raised on GitHub as new issues.

**Claude:** Identified the three issues as FSSgam_package#14 (`parallel = TRUE`
stalls when `gamm4` is loaded on the doSNOW workers), #15 (a supplied
`cor.matrix` must name the hard coded factor interaction columns), and #16
(`check_correlations()` compares factor-factor deviances computed on different
row sets). Each of #15 and #16 was raised with more than one candidate
resolution and no decision recorded, so the three decisions were put to the
user before any code was written:

- #16: fit the intercept-only model on the pair's complete cases, keeping the
  hoisted per-factor value for any pair that drops no extra rows.
- #15: accept a matrix covering only the predictors the user named, and compute
  rows and columns for the interaction columns the function itself creates.
- #14: measure the three untried directions recorded on the issue and adopt one
  only if it measurably lowers the stall rate.

---

### Issue #16 -- factor-factor deviances computed on different row sets

**Claude:** `fill_factor_factor_correlations()` estimates the correlation
between two factors as `sqrt(1 - fit/null.fit)`. `fit` is fitted on
`na.omit(dat[, c(var.1, var.2)])`; `null.fit` came from a per-factor
intercept-only model that `multinom()` had fitted on `var.1`'s own complete
cases. When `var.2` carries missing values the two are computed on different
numbers of rows and the ratio is not meaningful.

The correction chosen matches the convention the rest of the matrix already
follows -- `cor(use = "pairwise.complete.obs")` for the continuous block and a
per-pair `lm()` for the factor-continuous block are both pairwise complete-case.
Fitting the null per pair unconditionally would have undone the earlier
reduction from n(n-1) null fits to n, so `complete.counts` records how many rows
each factor's own null was fitted on and the null is refitted only for a pair
whose complete-case count is lower. The pair's rows are a subset of the
factor's, so equal counts mean an identical row set and the hoisted value is
exact. Data with no missing values therefore costs nothing extra and returns
identical values.

`null.fit` is now guarded against being a `try-error` alongside `fit`. Before
this the sequential loop would have raised an error out of `round()` if the null
fit failed while the pair fit succeeded; the cell is now left `NA`, which is how
a failed pair fit was already handled.

Two tests were added. The first constructs two factors where `f2` reproduces
`f1` on three rows in four and is missing for most of one level, so the row sets
differ by 25 of 90 rows; the estimate moves from 0.848 to 0.728. It asserts the
matched value and separately asserts the unmatched value is absent, so reverting
the fix fails on both. The second asserts the cached path: with `f2` as the
response the pair's rows are exactly `f2`'s own complete cases, and the hoisted
null is already correct. Confirmed the first test fails against `origin/master`
and the second passes, which is what distinguishes them.

Full suite after the change: 571 passing, 0 failing, 7 skipped, against 568 on
`master`.

---

#### Review round 1 (independent session)

**Claude (reviewer):** No blocking defect. The `complete.counts` shortcut was
confirmed exactly equivalent to refitting the null for every pair, by running a
variant with the condition replaced by `if(TRUE)` over 150 generated data frames
and comparing with `identical()`; all 150 matched, and the 10 frames with no
missing values matched `origin/master` as well. `complete.counts` was confirmed
to reach the doSNOW workers. Two substantive points were raised: the null is
refitted before the `try-error` guard, so a pair whose fitted model has already
failed pays for a discarded fit and emits its warnings (measured: 5 warnings
against master's 3 on a pair with no overlapping complete cases); and the second
test's name claims it detects whether a refit occurred, which it cannot.

**Claude:** Both were addressed. The refit condition now also requires `fit` to
have succeeded, which restores the warning count to 3. The second test was
renamed to what it asserts, and a third block was added that does pin the
shortcut: six ordered pairs over three factors with missingness in different
rows, each compared against the null refitted for that pair, so any divergence
between the retained value and the refitted one fails. The defensive nature of
the `null.fit` guard is now stated on the face of the code, the `na.action`
comment was qualified, and `check_non_linear_correlations()` gained the matching
`@details` paragraph.

The sibling function's own missing `null.fit` guard was raised as
FSSgam_package#19 rather than fixed here.

`DESCRIPTION` was bumped to 1.1.0.9000 so that the development `NEWS.md` section
and the declared version agree.

---

#### Review round 2 (fresh session)

**Claude (reviewer):** No blocking defect. The three claims singled out in the
round-1 response were re-measured and hold: the warning count is 3 and matches
`master`, the third test block fails under two independent breaks of the
shortcut, and the suite figures reproduce. Four points were raised: the new test
block's comment misstates which of its six pairs take the retained value (two,
not four); the pull request body still described the state at the first commit;
no test reaches the parallel copy of the changed lines, because the only
`parallel = TRUE` fixture is built from `case_study1`, which has no missing
values; and where every factor is missing in rows the others are not, none of the
hoisted nulls is used, so 15 `multinom()` calls are made where refitting per pair
with no hoist would be 12.

**Claude:** The comment was corrected against the block's own counts. The
`parallel = TRUE` fixture now gives `ZONE2` missing values, so both copies of the
refit condition are compared; verified against an installed copy with
`FSSGAM_TEST_PARALLEL=true`, two runs completing and agreeing and a third
stalling, which is FSSgam_package#14.

The hoisted nulls are now fitted only for the factors some pair retains, decided
from the missingness pattern before anything is fitted rather than lazily, which
a cache could not do across doSNOW workers. A factor's null is read only by a
pair that drops none of the rows on which that factor is present, so the
condition is exact. Measured with `trace(nnet::multinom)` on three three-level
factors over 90 rows: 9 calls with no missing values, matching `master`, and 12
with disjoint missingness, down from 15 and equal to the no-hoist minimum. The
unread fits' warnings go with them: the no-overlap frame now raises 2 warnings
where `master` raises 3.

The version bump was kept and `cran-comments.md` was corrected instead. It now
records that a development version is carried between releases, that this adds
the "Version contains large components" NOTE under `R CMD check --as-cran`, that
`devtools::check()` does not run that check at all, and that the version is set
to the release number before a submission is built.

The redirecting `DESCRIPTION` URLs were raised as FSSgam_package#20.

---
