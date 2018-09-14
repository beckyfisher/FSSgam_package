#' \tabular{ll}{
#' Package: \tab FSSgam\cr
#' Type: \tab Package\cr
#' Title: \tab FUll subsets multiple regression in R with gam(m4)\cr
#' Version: \tab 1.11\cr
#' Date: \tab 2018-09-14\cr
#' Author: \tab Rebecca Fisher\cr
#' Maintainer: \tab Rebecca Fisher\cr
#' License: \tab GPL (>= 3)\cr
#' LazyLoad: \tab yes\cr
#' Depends: \tab doSNOW, MuMIn, gamm4, mgcv, nnet\cr
#' }
#'
#' @details
#' Full subsets information theoretic approaches are becoming an increasingly popular tool for exploring predictive power and variable importance where a wide range of candidate predictors are being considered.
#'
#' This package provides a simple function that can be used to construct, fit and compare a complete model set of possible ecological or environmental predictors, given a response variable of interest. The function is based on Generalized Additive Models (GAM) and builds on the MuMIn package.
#'
#' Advantages include the capacity to fit more predictors than there are replicates, automatic removal of models with correlated predictors, and model sets that include interactions between factors and smooth predictors, as all as smooth interactions with other smooths (via te).
#'
#' The function takes a range of arguments that allow control over the model set being constructed, including specifying cyclic and linear continuous predictors, specification of the smoothing algorithm used and the maximum complexity allowed for smooth terms.
#'
#' The use of the function is demonstrated via case studies that highlight how appropriate model sets can be easily constructed, and the broader utility of the approach for exploratory ecology.
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