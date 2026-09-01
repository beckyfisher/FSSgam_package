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

**Claude:** Identified the three issues as FSSgam_package#14, #15 and #16, and
put the behaviour decisions in #15 and #16 to the user before writing any code,
since each was raised with more than one candidate resolution and none recorded.
This branch carries the work on #15 alone.

### Issue #15 -- a supplied cor.matrix had to name columns the user cannot know

**Claude:** With `factor.factor.interactions` set, `resolve_factor_interactions()`
appends a hard coded interaction column for every combination of factors that
survives the collinearity screen, and that screen reads the supplied
`cor.matrix`. `build_predictor_correlation_matrix()` then validated the same
matrix against `all.predictors`, which by that point includes those columns, and
stopped with them named. The names are deterministic, but which of them exist
depends on the matrix being validated, so a user cannot construct the matrix
without first knowing what the matrix will produce.

Option 2 of the three recorded on the issue was chosen: accept a matrix covering
only the predictors the user named, and compute rows and columns for the
interaction columns the function itself created. Option 1, falling back for any
missing predictor, was rejected because it would silently compute a value for a
predictor the user misspelled instead of reporting it, which is what the existing
check is for. Option 3, documenting the construction, leaves the error arriving
partway through construction.

`resolve_factor_interactions()` now returns the names of the columns it appended,
which is what lets `build_predictor_correlation_matrix()` tell a name the user
could have supplied from one the function invented. Missing predictors are split
on that: a named one is still an error, and a derived one is passed to
`augment_supplied_correlation_matrix()`, which computes a matrix over
`all.predictors` and splices in only the rows and columns involving a derived
column. Every supplied cell is kept as supplied, which was asserted by supplying
zeros for a pasted interaction column that the data say is perfectly collinear
with its components and checking that a candidate model holding both survives.

Computing those rows reads every predictor's data, so the FSSgam_package#13
allowance for a predictor of a class `check_correlations()` rejects does not
extend to this case. That is reported as such, naming both the interaction
columns and the underlying class error, rather than as the class error alone.

Five tests were added. Three fail against `origin/master` with the error the
issue reports. The other two pin behaviour the fix must not take away -- that a
supplied value for an interaction column is still used rather than recomputed,
and that a predictor the user named and omitted is still reported -- and pass
against the pre-fix code, which is recorded on the face of each so neither is
mistaken for a regression test.

Full suite: 583 passing, 0 failing, 7 skipped, against 568 on `master`.

---

#### Issue #15, review round 1 (independent session)

**Claude (reviewer):** One blocking defect and three substantive points. A
`cor.matrix` supplied as a `data.frame` crashes in the new path: assigning a
`data.frame` into a sub-block of a matrix replaces the matrix with a plain list,
dropping `dim` and `dimnames`, so the next assignment fails with "incorrect
number of subscripts on matrix". A `data.frame` works on `master` and works on
the branch whenever the matrix is complete. The comment claiming the `NA` cells
are never read is false: the character form of `smooth.smooth.interactions`
validates against the rownames of the resolved matrix rather than against the
predictors, so a supplied name that is not a predictor reaches
`combine_uncorrelated()`, where `if (NA)` raises an error naming nothing. No test
fails when the named/derived split is removed, because the blocks written for it
omit a factor, which `resolve_factor_interactions()` rejects several stages
earlier. And `@param cor.matrix` states the new exception and then retains,
unqualified, the FSSgam_package#13 sentence saying `check_correlations` is not
called at all.

Confirmed alongside: the golden master gives 30 of 31 scenarios `identical()`,
the one difference being the intended error message; the splice is correct with
reordered dimnames and with a partially supplied set of derived columns; and
`resolve_factor_interactions()` returns `character(0)` in every branch that
builds nothing.

**Claude:** All four fixed. The right hand side of the splice is coerced with
`as.matrix()`. The cells the function itself created are set to zero after the
splice, matching what the branch that computes the whole matrix already does with
the `NA` values `check_correlations()` returns, and leaving an `NA` the user
supplied among their own predictors alone. The supplied names are now kept in
their own order per dimension, with the derived ones appended to each, so a
matrix whose two dimnames differ is no longer reshaped into a square one with an
invented all-`NA` row. The `@param` sentence is qualified rather than deleted,
since it remains true for a complete matrix.

Four tests were added. Three fail against the previous commit on this branch, one
with the `data.frame` crash, one with "missing value where TRUE/FALSE needed",
and one on the reshaped dimnames. The fourth omits a continuous predictor
alongside a missing derived column, which is the case that reaches
`build_predictor_correlation_matrix()` with both kinds missing; it fails under
the reviewer's mutation folding named-missing into derived-missing.

The duplicate column a generated interaction name can create when it collides
with an existing column was raised as FSSgam_package#22. The redirecting URLs the
reviewer found in `cran-comments.md` are FSSgam_package#20, already open.

---

#### Issue #15, review round 2 (fresh session)

**Claude (reviewer):** All four behavioural claims from round 1 reproduce, and
the suite figures are exact. One substantive defect remains: `derived.missing`
was the union over the two dimensions, so a derived name the user supplied as a
column and not as a row was treated as missing from both, and the splice
overwrote the column they did supply. That contradicts the guarantee the change
is justified by. Two figures quoted in the text do not reproduce: the mutation
described in the round 1 comment fails a pre-existing block rather than passing
the suite, and the count of blocks failing against `master` is higher than
stated. Six nits and two out-of-scope findings.

**Claude:** The missing names are now tracked per dimension and each splice is
indexed with its own set. A test supplies `ZONE.I.ZONE.copy` as a column only,
with values 0.11, 0.22 and 0.33, and asserts they survive; it fails against the
previous commit.

One correction to the reviewer's account of the consequence. The values in
`predictor.correlations` were indeed overwritten, but the candidate set was not
changed in that configuration: `enumerate_candidate_models()` screens both
triangles, and the row was not supplied, so the computed 1 excludes
`ZONE+ZONE.I.ZONE.copy` whichever way the column is resolved. The test records
that rather than asserting a candidate-set difference that does not occur.

A duplicated name in a supplied matrix is now reported rather than collapsed by
`unique()`, which would have kept one of the two silently.

Both figures were re-measured rather than adjusted: 9 blocks of the branch's
interaction test file fail against `master`'s `R/`, not 3.

The two out-of-scope findings were raised as FSSgam_package#26, a 1x1 supplied
`cor.matrix` being silently discarded through the `length(cor.matrix)==1`
sentinel, and FSSgam_package#27, an `NA` between two predictors in a supplied
matrix failing with a message naming nothing. Both reproduce unchanged on
`master`.

---

#### Issue #15, review round 3 (fresh session)

**Claude (reviewer):** Every claim verified, including the round 2 correction
about the candidate set, which the previous review had asserted the opposite of.
Two substantive findings. The derived names are appended to the end of each
dimension, so a user who supplied one as a column anywhere other than last gets a
resolved matrix whose two dimnames hold the same set in different orders;
`combine_uncorrelated()` and `enumerate_candidate_models()` index the sub-matrix
by position in each dimension separately, so the cells screened are not the pairs
intended and the candidate set moves. Measured, 18 candidates with the derived
column last and 16 with it first, from matrices with identical numeric content.
And the duplicated-name error is broader than its justification: a duplicate
among names the model set never indexes stops a call that succeeds on `master`.

**Claude:** Both fixed. The two screens now index by name, both dimensions in one
order. The defect is pre-existing -- a complete supplied matrix with its columns
permuted returns 5 candidates on `master` where the unpermuted matrix returns 7 --
but the partial supply this change accepts is what produces the divergence from a
matrix the user wrote consistently, so it is fixed here rather than deferred. The
golden master is 31 of 31 identical across the change, which is what makes it
safe: every scenario there supplies dimnames in one order, so the fix is a no-op
for all of them. Three regression tests were added, and all three fail against
the previous commit.

The duplicated-name check now reports only names in `all.predictors`, since those
are the only ones indexed.

---

#### Issue #15, review round 4 (fresh session)

**Claude (reviewer):** Nothing blocking or substantive. The name-based indexing
was attacked directly: membership is identical between the two forms, the
ordering difference is a symmetric permutation wherever the dimnames agree, and
both triangles are screened, so the outcome is unchanged, which is why the golden
master holds at 31 of 31. Non-square selection is unreachable at both sites.
Four nits, including two figures that do not reproduce and one undocumented
behaviour change: a factor named twice in `pred.vars.fact` used to select a 1x1
sub-matrix whose empty triangles gave `max()` an empty vector, so the combination
survived and a `ZONE.I.ZONE` column was built; by name it is the 2x2 whose
off-diagonal is 1, and the combination is dropped.

**Claude:** The behaviour change is now tested and recorded in `NEWS.md`. The
residual half of it is not fixed here: `enumerate_candidate_models()`
deduplicates a candidate's terms before indexing, so it still reaches a 1x1
sub-matrix for a duplicated term and still admits a nonsense `ZONE+ZONE`
candidate with two `-Inf` warnings. That is pre-existing and unchanged, and was
raised as FSSgam_package#28, with the observation that rejecting a duplicated
entry in `pred.vars.cont`, `pred.vars.fact` or `linear.vars` up front is a
clearer error than anything downstream can give.

Both figures were re-measured: 10 of the 15 new blocks fail against `master`'s
`R/`, 8 of them with the missing-predictors error, and the permuted-matrix
candidate counts in `NEWS.md` now state the `max.predictors` they were measured
at, since 7 against 5 holds at 2 and 8 against 5 at the default of 3.

Rounds 3 and 4 both found nothing blocking, and round 4 nothing substantive, so
the cycle on this pull request is closed.

---
