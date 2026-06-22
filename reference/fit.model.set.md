# fit.model.set

Deprecated alias for
[`fit_model_set`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.md).
Retained because `fit.model.set` is the function name cited in Fisher et
al. (2018, Ecology and Evolution) and is used by existing downstream
code. New code should call
[`fit_model_set`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.md)
directly.

## Usage

``` r
fit.model.set(...)
```

## Arguments

- ...:

  Arguments passed on to
  [`fit_model_set`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.md).

## Value

See
[`fit_model_set`](https://beckyfisher.github.io/FSSgam_package/reference/fit_model_set.md).

## Examples

``` r
library(mgcv)
data(case_study1)
use.dat <- case_study1
use.dat$site <- as.factor(use.dat$site)
test.fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
                 family = tw(), data = use.dat)
model.set <- generate_model_set(
  use.dat = use.dat,
  test.fit = test.fit,
  pred.vars.cont = c("complexity", "depth"),
  pred.vars.fact = "ZONE",
  null.terms = "s(site,bs='re')",
  max.predictors = 2,
  k = 3
)
fit.model.set(model.set, parallel = FALSE)
#> Warning: 'fit.model.set' is deprecated.
#> Use 'fit_model_set' instead.
#> See help("Deprecated")
#>   |                                                                              |                                                                      |   0%  |                                                                              |=========                                                             |  12%  |                                                                              |==================                                                    |  25%  |                                                                              |==========================                                            |  38%  |                                                                              |===================================                                   |  50%  |                                                                              |============================================                          |  62%  |                                                                              |====================================================                  |  75%  |                                                                              |=============================================================         |  88%  |                                                                              |======================================================================| 100%
#> $mod.data.out
#>                                         modname
#> null                                       null
#> complexity                           complexity
#> depth                                     depth
#> ZONE                                       ZONE
#> ZONE+complexity                 ZONE+complexity
#> ZONE+depth                           ZONE+depth
#> ZONE+complexity.by.ZONE ZONE+complexity.by.ZONE
#> ZONE+depth.by.ZONE           ZONE+depth.by.ZONE
#>                                                                                        formula
#> null                                                                        s(site, bs = "re")
#> complexity                                s(complexity, k = 3, bs = "cr") + s(site, bs = "re")
#> depth                                          s(depth, k = 3, bs = "cr") + s(site, bs = "re")
#> ZONE                                                                 ZONE + s(site, bs = "re")
#> ZONE+complexity                    s(complexity, k = 3, bs = "cr") + ZONE + s(site, bs = "re")
#> ZONE+depth                              s(depth, k = 3, bs = "cr") + ZONE + s(site, bs = "re")
#> ZONE+complexity.by.ZONE s(complexity, by = ZONE, k = 3, bs = "cr") + ZONE + s(site, bs = "re")
#> ZONE+depth.by.ZONE           s(depth, by = ZONE, k = 3, bs = "cr") + ZONE + s(site, bs = "re")
#>                             AICc      BIC r2.vals r2.vals.unique  edf
#> null                    601.8165 613.9550 0.13574             NA 2.93
#> complexity              530.3521 546.5630 0.55313             NA 5.15
#> depth                   604.2478 620.2653 0.17579             NA 4.98
#> ZONE                    603.0187 616.7110 0.13529             NA 4.02
#> ZONE+complexity         531.1829 548.2208 0.55773             NA 5.94
#> ZONE+depth              604.7738 621.8040 0.17570             NA 5.94
#> ZONE+complexity.by.ZONE 531.6302 549.9272 0.53534             NA 6.83
#> ZONE+depth.by.ZONE      600.4762 619.5704 0.23486             NA 7.56
#>                         edf.less.1 delta.AICc delta.BIC wi.AICc wi.BIC
#> null                             0     71.464    67.392   0.000  0.000
#> complexity                       0      0.000     0.000   0.457  0.616
#> depth                            0     73.896    73.702   0.000  0.000
#> ZONE                             0     72.667    70.148   0.000  0.000
#> ZONE+complexity                  0      0.831     1.658   0.302  0.269
#> ZONE+depth                       0     74.422    75.241   0.000  0.000
#> ZONE+complexity.by.ZONE          0      1.278     3.364   0.241  0.115
#> ZONE+depth.by.ZONE               0     70.124    73.007   0.000  0.000
#>                         complexity depth ZONE
#> null                             0     0    0
#> complexity                       1     0    0
#> depth                            0     1    0
#> ZONE                             0     0    1
#> ZONE+complexity                  1     0    1
#> ZONE+depth                       0     1    1
#> ZONE+complexity.by.ZONE          1     0    1
#> ZONE+depth.by.ZONE               0     1    1
#> 
#> $failed.models
#> named list()
#> 
#> $success.models
#> $success.models$null
#> 
#> Family: Tweedie(p=1.469) 
#> Link function: log 
#> 
#> Formula:
#> Herbivore.abundance ~ s(site, bs = "re")
#> 
#> Estimated degrees of freedom:
#> 1.93  total = 2.93 
#> 
#> REML score: 297.6204     
#> 
#> $success.models$complexity
#> 
#> Family: Tweedie(p=1.265) 
#> Link function: log 
#> 
#> Formula:
#> Herbivore.abundance ~ s(complexity, k = 3, bs = "cr") + s(site, 
#>     bs = "re")
#> 
#> Estimated degrees of freedom:
#> 1.73 2.42  total = 5.15 
#> 
#> REML score: 261.5252     
#> 
#> $success.models$depth
#> 
#> Family: Tweedie(p=1.465) 
#> Link function: log 
#> 
#> Formula:
#> Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re")
#> 
#> Estimated degrees of freedom:
#> 1.42 2.57  total = 4.98 
#> 
#> REML score: 297.5266     
#> 
#> $success.models$ZONE
#> 
#> Family: Tweedie(p=1.469) 
#> Link function: log 
#> 
#> Formula:
#> Herbivore.abundance ~ ZONE + s(site, bs = "re")
#> 
#> Estimated degrees of freedom:
#> 2.02  total = 4.02 
#> 
#> REML score: 297.9798     
#> 
#> $success.models$`ZONE+complexity`
#> 
#> Family: Tweedie(p=1.267) 
#> Link function: log 
#> 
#> Formula:
#> Herbivore.abundance ~ s(complexity, k = 3, bs = "cr") + ZONE + 
#>     s(site, bs = "re")
#> 
#> Estimated degrees of freedom:
#> 1.71 2.23  total = 5.94 
#> 
#> REML score: 262.176     
#> 
#> $success.models$`ZONE+depth`
#> 
#> Family: Tweedie(p=1.465) 
#> Link function: log 
#> 
#> Formula:
#> Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + ZONE + s(site, 
#>     bs = "re")
#> 
#> Estimated degrees of freedom:
#> 1.47 2.47  total = 5.94 
#> 
#> REML score: 297.7074     
#> 
#> $success.models$`ZONE+complexity.by.ZONE`
#> 
#> Family: Tweedie(p=1.262) 
#> Link function: log 
#> 
#> Formula:
#> Herbivore.abundance ~ s(complexity, by = ZONE, k = 3, bs = "cr") + 
#>     ZONE + s(site, bs = "re")
#> 
#> Estimated degrees of freedom:
#> 1.78 1.00 2.05  total = 6.83 
#> 
#> REML score: 261.0315     
#> 
#> $success.models$`ZONE+depth.by.ZONE`
#> 
#> Family: Tweedie(p=1.45) 
#> Link function: log 
#> 
#> Formula:
#> Herbivore.abundance ~ s(depth, by = ZONE, k = 3, bs = "cr") + 
#>     ZONE + s(site, bs = "re")
#> 
#> Estimated degrees of freedom:
#> 1.71 1.00 2.85  total = 7.56 
#> 
#> REML score: 294.4957     
#> 
#> 
#> $variable.importance
#> $variable.importance$aic
#> $variable.importance$aic$variable.weights.raw
#> complexity      depth       ZONE 
#>      1.000      0.000      0.543 
#> 
#> 
#> $variable.importance$bic
#> $variable.importance$bic$variable.weights.raw
#> complexity      depth       ZONE 
#>      1.000      0.000      0.384 
#> 
#> 
#> 
```
