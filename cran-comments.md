## Test environments

* Local: Debian GNU/Linux (WSL2), R 4.5.2, via `devtools::check()`
* GitHub Actions (`.github/workflows/R-CMD-check.yaml`): ubuntu-latest
  (release and devel), macos-latest (release), windows-latest (release)
* R-hub was not used for this submission (not installed / not run).

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes for CRAN reviewers

* `generate.model.set()` and `fit.model.set()` are intentionally retained as
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
  slowest, `full.subsets.gam()`, completes in well under 1 second.

## Downstream dependencies

This is a new CRAN submission; there are no reverse dependencies to check.
