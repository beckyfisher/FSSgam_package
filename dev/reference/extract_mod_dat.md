# extract_mod_dat

Supporting function for functions full_subsets_gam and fit_model_set.
Not called directly.

## Usage

``` r
extract_mod_dat(mod.fit, r2.type. = "r2.lm.est", logLik.fn = NULL)
```

## Arguments

- mod.fit:

  A dsm, gam or uGamm fitted model object

- r2.type.:

  The type of r2 to extract. Passed through arguments supplied to
  fit_model_set

- logLik.fn:

  A function of one argument, a fitted model, returning a single
  log-likelihood value, or NULL (the default) to read AICc from
  MuMIn::AICc and BIC from stats::BIC as before. When supplied, AICc and
  BIC are built from the value it returns, at the degrees of freedom and
  sample size the default route uses, so only the log-likelihood
  changes. fit_model_set passes this through, and supplies it itself for
  a test.fit fitted with one of mgcv's censored families.

## Value

A list of model fit parameters

## Details

Extracts model fit parameters from a dsm, gam or uGamm fitted model
object. Called directly, this function reads AICc and BIC from
MuMIn::AICc and stats::BIC whatever the fitted family is, so a censored
fit gives the value mgcv reports rather than one built from a censored
log-likelihood. It is fit_model_set that resolves which log-likelihood a
model set is ranked on and passes it here as logLik.fn.

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
