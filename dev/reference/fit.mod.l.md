# fit.mod.l

Supporting function for functions full.subsets.gam and fit_model_set.
Not called directly.

## Usage

``` r
fit.mod.l(formula.l, test.fit., use.dat)
```

## Arguments

- formula.l:

  A model formula

- test.fit.:

  A dsm, gam or uGamm fitted model object

- use.dat:

  the data used to fit test.fit#

## Value

An updated dsm, gam or uGamm fitted model object

## Details

Generates an updated model fit based on the supplied formula. This
wrapper was required to allow full.subsets.gam and fit_model_set to be
applied to dsm models

## Examples

``` r
library(mgcv)
data(case_study1)
base.fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
                 family = tw(), data = case_study1)
fit.mod.l(formula.l = ~ s(complexity, k = 3, bs = "cr"),
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
