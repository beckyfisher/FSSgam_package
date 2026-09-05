## Session: criterion columns (FSSgam_package#42, #44)
Date: 2026-09-05
Model: Claude Opus 5 (1M context)

### Prompts and Responses

**User:** Address FSSgam_package#44 and #42 in a single PR, then review the PR
over up to five cycles, fixing what each review finds.

**Claude:** Both issues are failures of the same extractor. `extract_mod_dat()`
reads `AICc` from `MuMIn::AICc()` and `BIC` from `stats::BIC()`, and both
resolve to the fitted family's own `aic` slot. A quasi-likelihood has no
log-likelihood, so that slot returns `NA` and the whole model set is unrankable
(#44). `mgcv`'s `cnorm()` and `clog()` return a number that is not a censored
log-likelihood, so the set is ranked on a criterion that is wrong rather than
missing (#42).

Three decisions were put to the user before implementing, being the ones #42
lists as to be asked rather than guessed. The answers: the corrected values
overwrite `AICc` and `BIC` rather than arriving as new columns, so the delta
values, the weights and variable importance are corrected with them; the
censored log-likelihood is applied automatically with one message per model
set, so a user who never reads the issue gets a correct ranking and the change
in results between versions is visible; and the refusal for an undefined
criterion is raised in `fit_model_set()` rather than `generate_model_set()`,
because it is `logLik.fn` that decides whether such a `test.fit` can be used at
all and `generate_model_set()` never sees that argument.

Implemented as a new `R/criterion.R` holding `resolve_criterion()`, which
decides once per model set which log-likelihood the set is ranked on, plus
`censored_loglik()`, `criterion_from_loglik()` and the two small predicates they
need. `fit_model_set()`, `full_subsets_gam()` and `extract_mod_dat()` gain
`logLik.fn`.

**Reasoning recorded, since it is not recoverable from the diff.**

The censored log-likelihood is computed from the fitted mean, the scale and the
response's censoring coding rather than from `mgcv`'s components, so it is
independent of the three defects #42 documents and stays correct if they are
repaired. `mgcv`'s own `ls` and `dev.resids` slots would have given a shorter
implementation through the identity `deviance = 2 * (saturated - loglik)`, and
that route was rejected as the implementation and adopted as the *test oracle*
instead: the two share no code, so their agreement is evidence rather than a
tautology. They agree to 1e-8 or better for both families over left, right,
interval and no censoring, either column order, and non-unit prior weights.

Three properties of `mgcv` that the implementation depends on were measured
rather than assumed. `getTheta(TRUE)` returns the scale itself and not its
logarithm, despite `theta` being documented as a log scale parameter. A prior
weight divides the scale by its square root, matching `mgcv`'s own
`th <- theta - log(wt)/2`. And `attr(logLik(fit), "df")` is finite even for a
quasi-likelihood fit, whose `logLik()` value is `NA`, because `logLik.gam`
computes the degrees of freedom before it reads the `aic` slot; that is what
makes `logLik.fn` a usable route for a quasi `test.fit` rather than a dead end.

Two censoring codings were found during review to give no number and are
handled differently, for a reason worth keeping. A finite interval written with
the upper bound first fits identically under `mgcv`, which takes `pmin`/`pmax`
of the pair, so it is read the same way here. An *infinite* bound written first
is refused: `mgcv` selects a censored case on the second column alone, so such a
row matches none of `clog`'s index sets and is dropped from the fit without a
message, and scoring it would count observations the fit ignores. Measured on
`mgcv` 1.9-4, 60 of 200 rows written that way: -204.99 against the fit's own
-88.42.

The refusal for an unusable criterion is on the value `MuMIn::AICc()` returns
for the `test.fit`, not on a list of family names, so any family whose `aic`
slot is undefined is covered by the same route rather than only the two quasi
families the issue names.

**Two guards were added during review**, either side of the fitting loop.
`resolve_criterion()` puts the built-in through the same check a supplied
`logLik.fn` is put through, because a censored response it cannot evaluate
otherwise reached the fitting loop and produced an `NA` criterion for every
candidate under a message saying the ranking had been corrected -- which is
issue #44's output arrived at by another route. `fit_model_set()` then stops
where every candidate that fitted was given no criterion, and warns naming them
where only some were.

Reaching the second of those required separating "did not fit" from "fitted and
given no criterion" on the `save.model.fits = FALSE` path, which had defined a
failed model by an `NA` criterion. The two fitting paths disagreed in
consequence: the same partial failure returned eight successes with five `NA`
rows on one path and three successes with five `failed.models` on the other.
The unsaved loop now returns a `fit.ok` flag alongside each summary row, read
and dropped before `mod.data.out` is assembled, so both paths classify a failed
model by whether `fit_mod_l()` returned a fitted object.

**Rejected.** Quasi-AIC through `MuMIn::QAIC()`/`QAICc()` was not offered.
It needs a dispersion estimate and a decision about which fit supplies it, and
quasi-AIC comparisons are valid only within a common dispersion estimate, which
is a modelling decision rather than a defect fix. `logLik.fn` is what a user who
wants one would build it from.

**Not asserted in the tests.** No expectation states what `mgcv`'s own reported
criterion gives, so a repair upstream would fail none of them. The values that
route produced on 2026-09-05 are recorded in a comment beside the test instead.
The one deliberate exception is the scale-recovery test, which is aimed at a
change in `getTheta()`'s semantics rather than at a repair of the `aic` slots.

---
