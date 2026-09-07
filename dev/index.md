# FSSgam

FSSgam constructs, fits and compares a complete model set of candidate
ecological or environmental predictors for a response variable of
interest. Models are generalized additive models, fitted with
[mgcv](https://cran.r-project.org/package=mgcv) or
[gamm4](https://cran.r-project.org/package=gamm4), and ranked by AICc
using [MuMIn](https://cran.r-project.org/package=MuMIn).

The approach admits more predictors than there are replicates, removes
models with correlated predictors automatically, and supports model sets
containing interactions between factors and smooth predictors as well as
smooth-by-smooth interactions through
[`te()`](https://rdrr.io/pkg/mgcv/man/te.html). The method is described
in Fisher et al. (2018), *Ecology and Evolution*,
<doi:10.1002/ece3.4134>.

## Installation

``` r

# install.packages("devtools")
devtools::install_github("beckyfisher/FSSgam_package")
```

## Usage

An analysis is built in two steps.
[`generate_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/generate_model_set.md)
builds the candidate set from a test fit and a list of predictors, and
[`fit_model_set()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/fit_model_set.md)
fits and compares that set. Separating them allows the candidate set to
be inspected before anything is fitted, and allows the fitted objects to
be discarded as they are summarised, which matters for large sets.

``` r

library(FSSgam)
library(mgcv)

data(case_study1)
use.dat <- case_study1
use.dat$site <- as.factor(use.dat$site)

test.fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
                family = tw(), data = use.dat)

model.set <- generate_model_set(
  use.dat        = use.dat,
  test.fit       = test.fit,
  pred.vars.cont = c("complexity", "depth"),
  pred.vars.fact = "ZONE",
  null.terms     = "s(site,bs='re')",
  max.predictors = 2,
  k              = 3
)

out <- fit_model_set(model.set)
out$mod.data.out          # the model table, ranked by AICc
out$variable.importance   # summed model weights per predictor, by criterion
```

[`full_subsets_gam()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/full_subsets_gam.md)
performs both steps in one call. It saves every fitted model, so it is
recommended only for small candidate sets.

[`check_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_correlations.md)
and
[`check_non_linear_correlations()`](https://beckyfisher.github.io/FSSgam_package/dev/reference/check_non_linear_correlations.md)
report the predictor correlations the model set is screened against, and
can be called before generating a set.

## Documentation

Function reference: <https://beckyfisher.github.io/FSSgam_package/>

Worked case studies, an FAQ and the material accompanying the
publication are maintained in the companion repository, which is the
citable reference for the method:
<https://github.com/beckyfisher/FSSgam> and
<https://beckyfisher.github.io/FSSgam/>

## Citation

Fisher R, Wilson SK, Sin TM, Lee AC, Langlois TJ (2018) A simple
function for full-subsets multiple regression in ecology with R.
*Ecology and Evolution* 8(12): 6104-6113. <doi:10.1002/ece3.4134>

# License

The code is released under the Apache License 2.0

``` R
Copyright 2020 Australian Institute of Marine Science

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at 

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
```
