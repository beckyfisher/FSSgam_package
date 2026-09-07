## Test environments

* Local: Debian GNU/Linux (WSL2), R 4.6.1, `R CMD check --as-cran` on the built
  tarball with `_R_CHECK_CRAN_INCOMING_=TRUE` and
  `_R_CHECK_CRAN_INCOMING_REMOTE_=TRUE`
* GitHub Actions (`.github/workflows/R-CMD-check.yaml`): ubuntu-latest
  (release, devel and 4.4.0, the declared floor), macos-latest (release),
  windows-latest (release)
* R-hub was not used for this submission.

## R CMD check results

0 errors | 0 warnings | 1 note, for 1.2.0, measured on the local environment
above:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Rebecca Fisher <r.fisher@aims.gov.au>'

New submission
```

That note is the whole of it. The test suite reports 910 passing, 0 failing and
13 skipped under the check.

The check was run with `--no-manual`, this machine's TeX Live having neither
`inconsolata` nor the `times` metrics. The manual was built separately with
`R CMD Rd2pdf` under `R_RD4PDF="hyper"`, which completed with no LaTeX errors.

`urlchecker::url_check()` reports all 11 URLs in the package as correct.

The package was also checked with `_R_CHECK_DEPENDS_ONLY_=true`, so that the
`Suggests` are unavailable: `Status: OK`, with 887 tests passing, 0 failing and
20 skipped. `gamm4` is the only `Suggests` the tests use, and the fixture that
needs it calls `skip_if_not_installed()`.

## Notes for CRAN reviewers

* `full.subsets.gam()`, `generate.model.set()` and `fit.model.set()` are
  intentionally retained as exported, deprecated aliases (`R/deprecated.R`) for
  `full_subsets_gam()`, `generate_model_set()` and `fit_model_set()`. The dotted
  names are those cited in Fisher et al. (2018, Ecology and Evolution,
  <https://doi.org/10.1002/ece3.4134>) and are used by existing downstream code
  and teaching material; removing them would break published, citable
  workflows. Each wrapper calls `.Deprecated()` and forwards all arguments to
  its snake_case replacement.
* This package does not ship vignettes. Worked examples and case studies are
  maintained in the companion repository
  <https://github.com/beckyfisher/FSSgam>, which is cited in the package
  documentation (`@references`, `URL`).
* All `@examples` are runnable; none are wrapped in `\donttest{}`. The slowest
  of the thirteen takes 0.60 seconds elapsed, read from the 1.2.0 check's own
  `FSSgam-Ex.timings`. Three of them are the examples of the deprecated aliases
  named above, and each emits one deprecation warning, that being the behaviour
  the example documents.
* No `Language` field is declared. The documentation uses Australian spelling,
  but `DESCRIPTION` names "Generalized Additive Models", the standard term.
  Measured with `hunspell` 3.0.5: `en_GB` flags six words in the `Title` and
  `Description` against `en_US`'s four, so declaring `en-GB` would report more
  of them rather than fewer.
* A small number of tests are skipped on CRAN. Those covering `parallel = TRUE`
  start a real `doSNOW` cluster, and cluster dispatch has been observed to stall
  indefinitely on at least one platform when `gamm4` is loaded onto the workers.
  The stall reproduces with a `foreach()` loop containing no code from this
  package, so it is not a defect here, but an unattended stall would consume the
  entire runtime of a check. Those tests require an explicit opt-in environment
  variable and are run in a dedicated workflow with its own timeout. The code
  paths they cover are also exercised sequentially.

## Downstream dependencies

This is a new CRAN submission; there are no reverse dependencies to check.
