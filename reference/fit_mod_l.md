# fit_mod_l

Supporting function for functions full_subsets_gam and fit_model_set.
Not called directly.

## Usage

``` r
fit_mod_l(
  formula.l,
  test.fit.,
  use.dat,
  family. = resolve_candidate_family(test.fit.)
)
```

## Arguments

- formula.l:

  A model formula

- test.fit.:

  A dsm, gam or uGamm fitted model object

- use.dat:

  the data used to fit test.fit#

- family.:

  The family to refit formula.l with. Defaults to a fresh, independent
  re-evaluation of the family test.fit. itself used (see
  resolve_candidate_family in R/utils.R), so repeated calls never share
  mutable extended-family state (e.g. mgcv's nb()/tw() estimated theta).
  fit_model_set() resolves this once per candidate up front, on the
  calling process, and passes it in explicitly rather than relying on
  this default – see the comment below for why.

## Value

An updated dsm, gam or uGamm fitted model object

## Details

Generates an updated model fit based on the supplied formula. This
wrapper was required to allow full_subsets_gam and fit_model_set to be
applied to dsm models

## Examples

``` r
library(mgcv)
data(case_study1)
base.fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
                 family = tw(), data = case_study1)
fit_mod_l(formula.l = ~ s(complexity, k = 3, bs = "cr"),
          test.fit. = base.fit, use.dat = case_study1)
#> 
#> Family: Tweedie(p=1.284) 
#> Link function: log 
#> 
#> Formula:
#> Herbivore.abundance ~ s(complexity, k = 3, bs = "cr")
#> 
#> Estimated degrees of freedom:
#> 1.79  total = 2.79 
#> 
#> REML score: 262.1906     
```
