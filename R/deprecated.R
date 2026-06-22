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

#' generate.model.set
#'
#' Deprecated alias for \code{\link{generate_model_set}}. Retained because
#' \code{generate.model.set} is the function name cited in Fisher et al.
#' (2018, Ecology and Evolution) and is used by existing downstream code.
#' New code should call \code{\link{generate_model_set}} directly.
#'
#' @param ... Arguments passed on to \code{\link{generate_model_set}}.
#' @export
#' @return See \code{\link{generate_model_set}}.
#' @examples
#' library(mgcv)
#' data(case_study1)
#' use.dat <- case_study1
#' use.dat$site <- as.factor(use.dat$site)
#' test.fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
#'                  family = tw(), data = use.dat)
#' model.set <- generate.model.set(
#'   use.dat = use.dat,
#'   test.fit = test.fit,
#'   pred.vars.cont = c("complexity", "depth"),
#'   pred.vars.fact = "ZONE",
#'   null.terms = "s(site,bs='re')",
#'   max.predictors = 2,
#'   k = 3
#' )
generate.model.set <- function(...) {
  .Deprecated("generate_model_set")
  generate_model_set(...)
}

#' fit.model.set
#'
#' Deprecated alias for \code{\link{fit_model_set}}. Retained because
#' \code{fit.model.set} is the function name cited in Fisher et al.
#' (2018, Ecology and Evolution) and is used by existing downstream code.
#' New code should call \code{\link{fit_model_set}} directly.
#'
#' @param ... Arguments passed on to \code{\link{fit_model_set}}.
#' @export
#' @return See \code{\link{fit_model_set}}.
#' @examples
#' library(mgcv)
#' data(case_study1)
#' use.dat <- case_study1
#' use.dat$site <- as.factor(use.dat$site)
#' test.fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
#'                  family = tw(), data = use.dat)
#' model.set <- generate_model_set(
#'   use.dat = use.dat,
#'   test.fit = test.fit,
#'   pred.vars.cont = c("complexity", "depth"),
#'   pred.vars.fact = "ZONE",
#'   null.terms = "s(site,bs='re')",
#'   max.predictors = 2,
#'   k = 3
#' )
#' fit.model.set(model.set, parallel = FALSE)
fit.model.set <- function(...) {
  .Deprecated("fit_model_set")
  fit_model_set(...)
}

#' full.subsets.gam
#'
#' Deprecated alias for \code{\link{full_subsets_gam}}. Retained because
#' \code{full.subsets.gam} is used directly in the companion docs repository's
#' published case studies (\url{https://github.com/beckyfisher/FSSgam}), e.g.
#' case-study-2.Rmd. New code should call \code{\link{full_subsets_gam}}
#' directly.
#'
#' @param ... Arguments passed on to \code{\link{full_subsets_gam}}.
#' @export
#' @return See \code{\link{full_subsets_gam}}.
#' @examples
#' library(mgcv)
#' data(case_study1)
#' use.dat <- case_study1
#' use.dat$site <- as.factor(use.dat$site)
#' test.fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
#'                  family = tw(), data = use.dat)
#' full.subsets.gam(
#'   use.dat = use.dat,
#'   test.fit = test.fit,
#'   pred.vars.cont = c("complexity", "depth"),
#'   pred.vars.fact = "ZONE",
#'   null.terms = "s(site,bs='re')",
#'   max.predictors = 2,
#'   k = 3
#' )
full.subsets.gam <- function(...) {
  .Deprecated("full_subsets_gam")
  full_subsets_gam(...)
}
