#    Copyright 2020 Australian Institute of Marine Science
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

#' \tabular{ll}{
#' Package: \tab FSSgam\cr
#' Type: \tab Package\cr
#' Title: \tab FUll subsets multiple regression in R with gam(m4)\cr
#' Version: \tab 1.11\cr
#' Date: \tab 2018-09-14\cr
#' Author: \tab Rebecca Fisher\cr
#' Maintainer: \tab Rebecca Fisher\cr
#' License: \tab Apache 2\cr
#' LazyLoad: \tab yes\cr
#' Depends: \tab doSNOW, MuMIn, gamm4, mgcv, nnet\cr
#' }
#'
#' @details
#' Full subsets information theoretic approaches are becoming an increasingly popular tool for exploring predictive power and variable importance where a wide range of candidate predictors are being considered.
#'
#' This package provides simple function(s) that can be used to construct, fit and compare a complete model set of possible ecological or environmental predictors, given a response variable of interest. The function(s) are based on Generalized Additive Models (GAM) and builds on the MuMIn package.
#'
#' Advantages include the capacity to fit more predictors than there are replicates, automatic removal of models with correlated predictors, and model sets that include interactions between factors and smooth predictors, as all as smooth interactions with other smooths (via te).
#'
#' The function(s) takes a range of arguments that allow control over the model set being constructed, including specifying cyclic and linear continuous predictors, specification of the smoothing algorithm used and the maximum complexity allowed for smooth terms.
#'
#' The full subsets analysis can be carried out via one of two alternative methods allowed in the package.
#'
#' The first is through a direct call to full.subsets.gam (this is the original function).
#' This function both constructs and fits the complete model set, based on the user supplied input. This function requires that all model fits are saved, and is therefore
#' not suitable for extremely large models sets, as these will cause issues with memory. This method may be superceded in future versions of FSSgam, so for any new project please use the second method.
#'
#' The second method is via a call to generate.model.set followed by as second call to fit.model set. This pair of functions splits the process of generating the model set from actually fitting and extracting the relevant model data.
#' This method is useful for large model sets, because it allows the model set to be interrrogated before fitting and also optionally allows model fit data to not be saved, thus alleviating memory issues.
#'
#' The use of the function(s) is demonstrated via case studies that highlight how appropriate model sets can be easily constructed, and the broader utility of the approach for exploratory ecology.
#' Please see the case study files on github for usage examples at \url{https://github.com/beckyfisher/FSSgam}
#'
#' @name FSSgam-package
#' @aliases FSSgam FSSgam-package
#' @docType package
#' @author
#' Rebecca Fisher (Australian Institue of Marine Science)
#'
#' Maintainer: Rebecca Fisher \email{r.fisher@@aims.gov.au}
#' @references Fisher R, Wilson SK, Sin TM, Lee AC, Langlois TJ (2018) A simple function for full-subsets multiple regression in ecology with R. Ecology and Evolution
#' \url{https://onlinelibrary.wiley.com/doi/abs/10.1002/ece3.4134}
#' @examples
#' install.packages("FSSgam",dependencies=TRUE)
#' library(FSSgam)
NULL