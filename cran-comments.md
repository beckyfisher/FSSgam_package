## Test environments

* Local: Debian GNU/Linux (WSL2), R 4.6.1, via `devtools::check()` and via
  `R CMD check --as-cran` on the built tarball with `_R_CHECK_CRAN_INCOMING_=TRUE`
  (remote incoming checks disabled, as they require network access)
* GitHub Actions (`.github/workflows/R-CMD-check.yaml`): ubuntu-latest
  (release and devel), macos-latest (release), windows-latest (release)
* R-hub was not used for this submission (not installed / not run).

## R CMD check results

0 errors | 0 warnings | 0 notes, under `devtools::check()`, for the 1.1.0
release.

`devtools::check()` sets `_R_CHECK_CRAN_INCOMING_=FALSE`, so it does not run
the CRAN incoming feasibility check and that figure does not cover it.
`R CMD check --as-cran` on the built tarball raises one incoming feasibility
NOTE, whose body lists "New submission" together with any other item that check
finds.

Between releases `DESCRIPTION` carries a development version (currently
1.1.0.9000) to match the development section of `NEWS.md`, which adds "Version
contains large components" to that NOTE's body. The version is set to the
release number before a submission is built, at which point that line does not
arise.

## Notes for CRAN reviewers

* `full.subsets.gam()`, `generate.model.set()` and `fit.model.set()` are intentionally retained as
  exported, deprecated aliases (`R/deprecated.R`) for `generate_model_set()`
  and `fit_model_set()`. The dotted names are the function names cited in
  Fisher et al. (2018, Ecology and Evolution,
  <https://doi.org/10.1002/ece3.4134>) and are used by existing downstream
  code and teaching material; removing them would break published, citable
  workflows. Each wrapper calls `.Deprecated()` and forwards all arguments
  to its snake_case replacement.
* This package does not ship vignettes. Worked examples and case studies
  are maintained in the companion repository
  <https://github.com/beckyfisher/FSSgam>, which is cited in the package
  documentation (`@references`, `URL`).
* All `@examples` are runnable (none are wrapped in `\donttest{}`); the
  slowest, `full_subsets_gam()`, completes in well under 1 second.
* A small number of tests are skipped on CRAN. Those covering
  `parallel = TRUE` start a real `doSNOW` cluster, and cluster dispatch has
  been observed to stall indefinitely on at least one platform when `gamm4`
  is loaded onto the workers. The stall is reproducible with a `foreach()`
  loop containing no code from this package, so it is not a defect here, but
  an unattended stall would consume the entire runtime of a check. Those
  tests require an explicit opt-in environment variable and are run in a
  dedicated workflow with its own timeout. The code paths they cover are
  also exercised sequentially.

## Downstream dependencies

This is a new CRAN submission; there are no reverse dependencies to check.
