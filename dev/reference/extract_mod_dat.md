# extract_mod_dat

Supporting function for functions full_subsets_gam and fit_model_set.
Not called directly.

## Usage

``` r
extract_mod_dat(mod.fit, r2.type. = "r2.lm.est")
```

## Arguments

- mod.fit:

  A dsm, gam or uGamm fitted model object

- r2.type.:

  The type of r2 to extract. Passed through arguments supplied to
  fit_model_set

## Value

A list of model fit parameters

## Details

Extracts model fit parameters from a dsm, gam or uGamm fitted model
object

## Examples

``` r
library(mgcv)
#> Loading required package: nlme
#> This is mgcv 1.9-4. For overview type '?mgcv'.
library(MuMIn)
data(case_study1)
fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
           family = tw(), data = case_study1)
extract_mod_dat(fit, r2.type. = "r2")
#> $AICc
#> [1] 602.3223
#> 
#> $BIC
#> [1] 611.5786
#> 
#> $r2.vals
#> [1] -0.00721
#> 
#> $r2.vals.unique
#> [1] NA
#> 
#> $edf
#> [1] 2.32
#> 
#> $edf.less.1
#> [1] 0
#> 
```
