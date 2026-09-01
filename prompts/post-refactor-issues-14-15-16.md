
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
