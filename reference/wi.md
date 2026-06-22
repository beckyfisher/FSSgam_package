# wi

Supporting function for functions full.subsets.gam and fit_model_set.
Not called directly.

## Usage

``` r
wi(AIC.vals)
```

## Arguments

- AIC.vals:

  vector of AICc, AIC or BIC values

## Value

A vector of Akaike weights

## Details

Calculates Akaike weight values from a vector of AICc, AIC or BIC values

## Examples

``` r
wi(c(100, 102, 105, 110))
#> [1] 0.686482022 0.252542623 0.056349876 0.004625479
```
