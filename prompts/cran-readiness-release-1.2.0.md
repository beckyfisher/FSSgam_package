## Session: CRAN readiness audit and release 1.2.0
Date: 2026-09-07
Model: Claude Opus 5 (1M context), claude-opus-5[1m]

### Prompts and Responses

**User:** All issues are resolved. I am going to work on updating the
documentation package. While I do that, can you go through and check
everything in this part is ready for cran.

**Claude:** Audited the package against `R CMD check --as-cran` on a built
tarball, with `_R_CHECK_CRAN_INCOMING_REMOTE_=TRUE`, plus independent checks
that `R CMD check` does not run. Measured on Debian WSL2, R 4.6.1,
testthat 3.3.2, mgcv 1.9-4, MuMIn 1.48.19, roxygen2 8.0.0.

The check reported 1 ERROR, 1 WARNING and 3 NOTEs, of which only one NOTE was
the package: "Version contains large components (1.1.0.9000)". The ERROR and
WARNING were the PDF manual failing on a missing `inconsolata.sty`, and on the
`times` metrics once inconsolata was removed; this machine's TeX Live has five
texlive packages. `R CMD Rd2pdf` under `R_RD4PDF="hyper"` produced a 179 KB
manual, exit 0, no LaTeX errors, so the Rd sources are sound. That run is the
only coverage of the manual: `.github/workflows/R-CMD-check.yaml` passes
`--no-manual`. The remaining two NOTEs were the left-over `.tex` from that
failure and a missing `tidy`.

Reported as ready: tests (910 passing, 0 failing, 13 skipped under the check;
919 and 8 skipped blocks locally at `NOT_CRAN=true`), example timings (slowest
0.54 s), all eight declared Imports used, no `T`/`F`, `<<-`, `set.seed()`,
`options()`/`par()` mutation or `\dontrun`, `tools/check-r-floor.R` passing,
every export documented with `\value` and runnable `\examples` and present in a
`_pkgdown.yml` group.

Reported as outstanding: the development version and the 555-line
"(development version)" `NEWS.md` section; `master` 59 commits behind `dev`
with no tags and no releases anywhere in the repository; a stale
`cran-comments.md`; "in R" in the `Title`; no `cph` in `Authors@R` though every
file in `R/` states the AIMS copyright; two URLs `urlchecker::url_check()`
reported invalid (the README CRAN badge, 404, and the Wiley link, 403); the
typos `calaculation` and `FFSgam`; and a 33-line README with no installation
instruction or usage example.

**User:** go ahead and apply this fixes on dev, and then merge to master

**Claude:** Applied on a branch off `dev`, merged as pull request #51, then
`dev` to `master` as #52 and tagged `v1.2.0`.

Package metadata: `Title` drops "in R"; `Authors@R` gains
`person("Australian Institute of Marine Science", role = "cph")`; the
Fisher et al. (2018) reference becomes `\doi{10.1002/ece3.4134}` in
`R/FSSgam.R`, `R/fit-model-set.R`, `R/full-subsets-gam.R` and
`R/generate-model-set.R`; `calaculation` corrected in
`R/check-correlations.R`; the README rewritten with the method, installation, a
worked two-step call run against `case_study1` before committing, and the
citation, with the CRAN badge removed until acceptance.

**The audit recommended `Language: en-GB` and the measurement overturned it.**
The reasoning had been that the documentation uses Australian spelling, which
an en-US dictionary flags. Measured with `hunspell` 3.0.5 on the tree: `en_GB`
flags six words in the `Title` and `Description` against `en_US`'s four, adding
"Generalized" and "GAMs", and 57 against 56 across `man/`. `DESCRIPTION` names
"Generalized Additive Models", the standard term, so no dictionary is clean for
this package, and CRAN reports the `DESCRIPTION` words. Declaring `en-GB` would
report more of them. The field was added, then removed; the measurement is
recorded in `NEWS.md` and `cran-comments.md` so it is not revisited.

A second self-inflicted finding: the rewritten README first wrote the DOI as
`<https://doi.org/10.1002/ece3.4134>`, which markdown auto-links and which
redirects to the publisher that returns 403, so it reintroduced in the README
the finding the `\doi{}` conversion had removed from the Rd. Written as a
plain-text `doi:10.1002/ece3.4134` instead. README.md ships in the tarball, so
its URLs are checked.

Release: `Version: 1.2.0`, minor rather than patch because the closed section
records several behaviour changes, two new arguments (`null.cov.cutoff` and
`logLik.fn`), a raised R floor and `gamm4` moving to `Suggests`. `NEWS.md`
heading renumbered; `tools:::.build_news_db_from_package_NEWS_md()` now returns
1.2.0, 1.1.0 and 1.0.0, having previously returned only the latter two, so the
development section had been recorded against no version.

Verified against the built `FSSgam_1.2.0.tar.gz`: `R CMD check --as-cran` with
remote incoming checks reports 0 errors, 0 warnings and 1 note, that note being
"New submission" alone; 910 passing, 0 failing, 13 skipped;
`urlchecker::url_check()` reports all 11 URLs correct; `Rd2pdf` builds the
manual. All nine checks passed on #51 and #52 across macOS release, Windows
release and Ubuntu release, devel and 4.4.0, plus the R floor consistency job,
the opt-in parallel tests, coverage and pkgdown.

---

## Session: CRAN review response, ::: in documentation
Date: 2026-09-17
Model: Claude Opus 5 (1M context)

### Prompts and Responses

**User:** Supplied the CRAN reviewer's reply to the 1.2.0 submission. It quotes
the policy on `foo:::f`, reports one occurrence, `man/fit_mod_l.Rd` calling
`FSSgam:::fit_mod_l(...)` in its example, and asks for a fix and resubmission.

**Claude:** `fit_mod_l()` was unexported in 1.1.0 (Phase 14) and is documented
with `@keywords internal`, so no example can call it without `:::`. The two ways
out are to export it again or to drop the example; re-exporting reverses a
deliberate decision to keep an internal helper out of the public API, so the
example was removed. The `@examples` block was deleted from the roxygen block in
`R/functions_supporting.R` and `man/fit_mod_l.Rd` regenerated with
`roxygen2::roxygenise()`. That file uses CRLF endings, so the deletion was made
with `sed` rather than an editor that would have rewritten the whole file
(Phase 14 note). `grep` over `R/` and `man/` finds no remaining `:::`; the
occurrences in `tests/` are not documentation and were left.

The version number stays at 1.2.0, on the user's ruling, 1.2.0 never having been
published. `NEWS.md` is therefore unchanged, and `cran-comments.md` gains a
`## Resubmission` section quoting the reviewer's finding and stating what was
changed.

Two figures in `cran-comments.md` were re-measured rather than copied from the
previous submission, the example count having changed from thirteen to twelve. Measured on this host
(Debian WSL2, R 4.6.1) against the built `FSSgam_1.2.0.tar.gz`:
`R CMD check --as-cran` with both incoming-check variables set reports 0 errors,
0 warnings and 1 note, that note being "New submission" alone, with 910 passing,
0 failing and 13 skipped; the slowest of the twelve examples takes 0.52 seconds
elapsed, read from `FSSgam-Ex.timings`. The check was repeated with
`_R_CHECK_DEPENDS_ONLY_=true`, giving `Status: OK` with 887 passing, 0 failing
and 20 skipped, which reproduces the figures already recorded for 1.2.0.

---
