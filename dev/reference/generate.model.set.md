# generate.model.set

Deprecated alias for
[`generate_model_set`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md).
Retained because `generate.model.set` is the function name cited in
Fisher et al. (2018, Ecology and Evolution) and is used by existing
downstream code. New code should call
[`generate_model_set`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
directly.

## Usage

``` r
generate.model.set(...)
```

## Arguments

- ...:

  Arguments passed on to
  [`generate_model_set`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md).

## Value

See
[`generate_model_set`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md).

## Examples

``` r
library(mgcv)
data(case_study1)
use.dat <- case_study1
use.dat$site <- as.factor(use.dat$site)
test.fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
                 family = tw(), data = use.dat)
model.set <- generate.model.set(
  use.dat = use.dat,
  test.fit = test.fit,
  pred.vars.cont = c("complexity", "depth"),
  pred.vars.fact = "ZONE",
  null.terms = "s(site,bs='re')",
  max.predictors = 2,
  k = 3
)
#> Warning: 'generate.model.set' is deprecated.
#> Use 'generate_model_set' instead.
#> See help("Deprecated")
```
