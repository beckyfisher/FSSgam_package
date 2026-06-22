# Package index

## Package overview

- [`FSSgam-package`](https://beckyfisher.github.io/FSSgam_package/dev/reference/FSSgam-package.md)
  [`FSSgam`](https://beckyfisher.github.io/FSSgam_package/dev/reference/FSSgam-package.md)
  : FSSgam: Full Subsets Multiple Regression in R Using GAMs

## Core functions

Build a candidate model set from a test fit and predictor list, then fit
and compare that model set by AICc.

- [`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
  : generate_model_set
- [`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md)
  : fit_model_set

## All-in-one wrapper

One-shot wrapper around generate_model_set() and fit_model_set();
recommended only for small candidate sets since it does not let you
inspect the model set before fitting.

- [`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)
  : full_subsets_gam

## Correlation diagnostics

- [`check_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_correlations.md)
  : check_correlations
- [`check_non_linear_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_non_linear_correlations.md)
  : check_non_linear_correlations

## Supporting functions

Used internally by full_subsets_gam() and fit_model_set(). Exported but
not intended to be called directly.

- [`wi()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/wi.md)
  : wi
- [`extract_mod_dat()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/extract_mod_dat.md)
  : extract_mod_dat
- [`build_inclusion_mat()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/build_inclusion_mat.md)
  : build_inclusion_mat
- [`fit_mod_l()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_mod_l.md)
  : fit_mod_l

## Deprecated aliases

Legacy dot.case names retained for backward compatibility.
generate.model.set() and fit.model.set() are the function names cited in
Fisher et al. (2018, Ecology and Evolution); full.subsets.gam() is
retained because it’s used directly in the companion docs repository’s
published case studies. New code should call the snake_case functions
above directly.

- [`generate.model.set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate.model.set.md)
  : generate.model.set
- [`fit.model.set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit.model.set.md)
  : fit.model.set
- [`full.subsets.gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full.subsets.gam.md)
  : full.subsets.gam

## Datasets

- [`case_study1`](https://beckyfisher.github.io/FSSgam_package/dev/reference/case_study1.md)
  : Case study 1 dataset
- [`case_study2`](https://beckyfisher.github.io/FSSgam_package/dev/reference/case_study2.md)
  : Case study 2 dataset
- [`case_study3`](https://beckyfisher.github.io/FSSgam_package/dev/reference/case_study3.md)
  : Case study 3 dataset
- [`coral_data`](https://beckyfisher.github.io/FSSgam_package/dev/reference/coral_data.md)
  : Coral_data
