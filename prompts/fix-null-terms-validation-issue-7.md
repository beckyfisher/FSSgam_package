## Session: fix-null-terms-validation-issue-7
Date: 2026-06-22
Model: Claude Sonnet 4.6

### Prompts and Responses

**User:** can you help me resolve this issue? https://github.com/beckyfisher/FSSgam/issues/7

**Claude:** Fetched issue #7 (filed by the package author, Jan 2021, no
comments): "Warnings when null terms is not type character" — requesting
"Add some error catching around the null.term passed to generate.model.set".
The issue lives in the read-only companion docs repo, but the bug it
describes is in `generate_model_set()`, which lives in this package repo, so
the fix belongs here.

Traced `null.terms` through `R/generate-model-set.R`: it's checked via
`nchar(null.terms) > 0` in `build_null_model()` and twice more in
`build_model_formulas()`, with no upfront type/length validation. Reproduced
six distinct failure modes by passing non-character values for `null.terms`:
- `NA` → `nchar(NA)` errors "missing value where TRUE/FALSE needed"
- `NULL` → "argument is of length zero"
- a numeric or logical scalar (e.g. `123`, `TRUE`) → silently *accepted*,
  `nchar()` coerces it to a string and it gets spliced straight into
  `as.formula(paste("~ ", null.terms))`, producing a nonsense formula with no
  error at all
- a character vector of length > 1 → "the condition has length > 1"
- a factor → "'nchar()' requires a character vector"

None of these point the user at the actual problem. Added
`validate_null_terms()` (`R/generate-model-set.R`, same pattern as the
existing `validate_use_dat()`), requiring `null.terms` to be a single,
non-`NA` character string (`""` allowed, matching the default), called right
after `validate_use_dat(use.dat)` in `generate_model_set()`. Since
`full_subsets_gam()` calls `generate_model_set()` directly (not through the
deprecated wrapper), this covers both public entry points without touching
`full-subsets-gam.R`.

Added a regression test in `tests/testthat/test-generate_model_set.R`
looping over all six bad inputs above and asserting the new, clear error
message. Updated `NEWS.md`. Verified: full `testthat` suite green (105
tests, up from 99) and a clean `devtools::check()` (0 errors/warnings, the
known sandbox time-verification NOTE only).

---
