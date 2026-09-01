
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
