# Package index

## Package overview

- [`FSSgam-package`](https://beckyfisher.github.io/FSSgam_package/reference/FSSgam-package.md)
  [`FSSgam`](https://beckyfisher.github.io/FSSgam_package/reference/FSSgam-package.md)
  : FSSgam: Full Subsets Multiple Regression in R Using GAMs

## Core functions

Build a candidate model set from a test fit and predictor list, then fit
and compare that model set by AICc.

- [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate_model_set.md)
  : generate_model_set
- [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.md)
  : fit_model_set

## All-in-one wrapper

One-shot wrapper around generate_model_set() and fit_model_set();
recommended only for small candidate sets since it does not let you
inspect the model set before fitting.

- [`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/reference/full.subsets.gam.md)
  : full.subsets.gam

## Correlation diagnostics

- [`check.correlations()`](https://beckyfisher.github.io/FSSgam_package/reference/check.correlations.md)
  : check.correlations
- [`check.non.linear.correlations()`](https://beckyfisher.github.io/FSSgam_package/reference/check.non.linear.correlations.md)
  : check.non.linear.correlations

## Supporting functions

Used internally by full.subsets.gam() and fit_model_set(). Exported but
not intended to be called directly.

- [`wi()`](https://beckyfisher.github.io/FSSgam_package/reference/wi.md)
  : wi
- [`extract.mod.dat()`](https://beckyfisher.github.io/FSSgam_package/reference/extract.mod.dat.md)
  : extract.mod.dat
- [`build.inclusion.mat()`](https://beckyfisher.github.io/FSSgam_package/reference/build.inclusion.mat.md)
  : build.inclusion.mat
- [`fit.mod.l()`](https://beckyfisher.github.io/FSSgam_package/reference/fit.mod.l.md)
  : fit.mod.l

## Deprecated aliases

Legacy dot.case names retained for backward compatibility. These are the
function names cited in Fisher et al. (2018, Ecology and Evolution) and
used by existing downstream code. New code should call the snake_case
functions above directly.

- [`generate.model.set()`](https://beckyfisher.github.io/FSSgam_package/reference/generate.model.set.md)
  : generate.model.set
- [`fit.model.set()`](https://beckyfisher.github.io/FSSgam_package/reference/fit.model.set.md)
  : fit.model.set

## Datasets

- [`case_study1`](https://beckyfisher.github.io/FSSgam_package/reference/case_study1.md)
  : Case study 1 dataset
- [`case_study2`](https://beckyfisher.github.io/FSSgam_package/reference/case_study2.md)
  : Case study 2 dataset
- [`case_study3`](https://beckyfisher.github.io/FSSgam_package/reference/case_study3.md)
  : Case study 3 dataset
- [`coral_data`](https://beckyfisher.github.io/FSSgam_package/reference/coral_data.md)
  : Coral_data
