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

#' full_subsets_gam
#'
#' Conducts a full subsets analysis based on gam(m4). In the most recent version of FSSgam this function is now a wrapper
#' for the two input functions, generate_model_set and fit_model_set. Input arguments are the same as these two underlying functions.
#' calling these underlying functions explicitly is the recommended method for running a full subsets analysis with the
#' FSSgam package because this enables the user to interrogate the candidate model set and the predictor correlation
#' matrix before actually running the analysis.
#'
#' @inheritParams generate_model_set
#'
#' @param parallel  A logical value indicating if parallel processing should be used. The default is FALSE.
#'
#' @param max.models The total number of models allowed to be fit and still save the model fits. Defaults to 200. If the candidate set is bigger than this value, a warning will be returned indicating that model fits will not be saved.
#'
#' @param save.model.fits Are the model fits to be saved in the output list? If TRUE this will be overwritten if the model candidate set is bigger than max.models. If FALSE only model output data are saved.
#'
#' @param n.cores An integer indicating the number of cores to use if parallel is TRUE. Defaults to 4.
#'
#' @param r2.type The value to extract from the gam model fit to use as the R squared value. Defaults to r2.lm.est which returns and estimated R squared value based on a linear regression between the observed and predicted values. r2 will return the adjusted R.sq as reported by gam, gamm or gamm4.dev will return the deviance explained as reported by gam or gamm. Note gamm4 does not currently return a deviance.
#'
#' @param  report.unique.r2 Should the r2.vals.unique column of mod.data.out be populated. Defaults to FALSE, which leaves it NA. When TRUE, the null model R2 is subtracted from each model R2 to give the variance explained beyond the terms supplied in null.terms. See \code{\link[=fit_model_set]{fit_model_set()}} for what the column is and is not.
#'
#' @param progress Should a text progress bar be written to the console while models are fitted. Defaults to interactive(), so the bar appears at the console but not in scripts, reports or checks.
#'
#' @param  VI.mods The set of models used to calculate summed variable importance scores. Defaults to 'min.n', which uses only the best n models for each variable (n being the minimum number of models any one predictor is present in). Set to 'all' to use all models in the candidate set instead.
#'
#' @param factor.interactions Deprecated. Superseded by factor.factor.interactions; retained only so older code does not break, and will warn if used.
#'
#' @param smooth.interactions Deprecated. Superseded by factor.smooth.interactions; retained only so older code does not break, and will warn if used.
#'
#' @param size Deprecated. Superseded by max.predictors; retained only so older code does not break, and will warn if used.
#'
#' @details The function constructs and fits a complete model set based on the supplied arguments.
#' for more information see Fisher R, Wilson SK, Sin TM, Lee AC, Langlois TJ (2018) A simple function for full-subsets multiple regression in ecology with R. Ecology and Evolution
#' https://onlinelibrary.wiley.com/doi/abs/10.1002/ece3.4134
#' @export
#' @return A list of the following output files:
#'
#' mod.data.out - A data.frame that contains the statistics associated with each model fit. This includes AICc and BIC, delta values (e.g. AICc-(min(AICc)), corresponding weight values (Burnham and Anderson 2003), an estimate of the model R2, and a column for each of the included predictor variables containing either 0 (variable not included in the model) or 1 (variable is present in the model).
#' A column r2.vals.unique is also present, and is NA unless report.unique.r2 is TRUE. This data.frame is the one fit_model_set() produced, passed through unaltered, so see \code{\link[=fit_model_set]{fit_model_set()}} for what the column is and is not.
#' Use of BIC in information theoretic approaches has been heavily criticised because of the inherent assumption of BIC that there is a true model that is represented in the candidate set (Anderson & Burnham 2002). Rather than decide a-priori which model selection tool users should adopt, we supply both as part of the function outputs.
#' To simplify output, only AICc and AICc based model weights, rather than AIC, are included as these are asymptotically equivalent at large sample sizes, and for small sample sizes AICc should be used in any case.
#' Calculating R2 values is non-trivial for mixed models, especially non-gaussian cases (and some argue should not be done at all). We have supplied a range of methods for estimating R2 (r2.type), as in our experience a single method rarely performs adequately across all scenarios.
#'
#' used.data - A data.frame which is identical to the data.frame initially supplied by the user, but with any hard coded interaction terms appended via cbind.
#'
#' predictor.correlations - The matrix of estimated predictor correlations returned by the function check_correlations and used for model exclusion based on cov.cutoff
#'
#' failed.models - A list containing the try-error catch associated with models that failed to fit. Ideally the list of failed models should be empty, but when this is not the case interrogating failed.models provides a useful means of troubleshooting. Users can examine which models are not fitting and explore the reasons for this by fitting the failed models outside the full_subsets_gam call based on the listed formula. When a large number of models fail to fit properly it usually indicates poor specification of the initial test.fit or other arguments in the call to full_subsets_gam (such as the inclusion of factor interactions when there are few data within each level of the factor), or that inappropriate variables are being included in the model set.
#'
#' success.models - A complete list of all successfully fitted models. This can be used for multimodel inference and creating model averaged predictions.
#'
#' variable.importance - A list containing importance scores for each included predictor.
#' To determine the relative importance of each predictor across the whole model set we summed the wi values for all models containing each variable. The higher the combined weights for an explanatory parameter, the more important it is in the analysis (Burnham & Anderson, 2002). An assumption of the use of summed model weights to infer variable importance is that the number of models in which the different predictors are present is uniform. As our function removes models with correlated predictors, this is not always the case. To overcome this issue, the summed variable.importance scores are the summed weights for the best n models, where n is equal to the minimum number of models any one predictor is present in.
#' @examples
#' library(mgcv)
#' data(case_study1)
#' use.dat <- case_study1
#' use.dat$site <- as.factor(use.dat$site)
#' test.fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
#'                  family = tw(), data = use.dat)
#' full_subsets_gam(
#'   use.dat = use.dat,
#'   test.fit = test.fit,
#'   pred.vars.cont = c("complexity", "depth"),
#'   pred.vars.fact = "ZONE",
#'   null.terms = "s(site,bs='re')",
#'   max.predictors = 2,
#'   k = 3
#' )

full_subsets_gam=function(use.dat,
                          test.fit,
                          pred.vars.cont=NA,
                          pred.vars.fact=NA,
                          cyclic.vars=NA,
                          linear.vars=NA,
                          factor.smooth.interactions=pred.vars.fact,
                          factor.factor.interactions=FALSE,
                          smooth.smooth.interactions=FALSE,
                          cov.cutoff=0.28,
                          cor.matrix=NA,
                          non.linear.correlations=FALSE,
                          max.predictors=3,
                          k=5,
                          bs.arg="'cr'",
                          null.terms="",
                          max.models=200,
                          save.model.fits=TRUE,
                          parallel=FALSE,
                          n.cores=4,
                          r2.type="r2.lm.est",
                          report.unique.r2=FALSE,
                          VI.mods='min.n',
                          progress=interactive(),
                          factor.interactions,
                          smooth.interactions,
                          size){
  # manage previous version arguments
  #
  # missing() rather than a "previous.arg" sentinel default. The sentinel had
  # to be compared with != or is.na(), and both fail on input the deprecated
  # arguments are documented to accept: a character vector of length > 1 gave
  # "the condition has length > 1", and NA gave "missing value where
  # TRUE/FALSE needed". missing() is correct for every type, length and NA.
  # full.subsets.gam() forwards through ..., so missing() still reports
  # correctly through the deprecated alias.
  if(!missing(factor.interactions)){
     factor.factor.interactions=factor.interactions
     warning('Argument factor.interactions has been replaced with factor.factor.interactions.
              Please update your code as usage of factor.interactions will not be supported in
              future versions.')
     }
  # NA no longer needs its own branch: it is simply assigned, which is its
  # documented meaning for factor.smooth.interactions ("If specified as NA no
  # factor-continuous predictor interactions will be included"), and it now
  # warns like any other use of the deprecated argument.
  if(!missing(smooth.interactions)){
     factor.smooth.interactions=smooth.interactions
     warning('Argument smooth.interactions has been replaced with factor.smooth.interactions.
              Please update your code as usage of smooth.interactions will not be supported in
              future versions.')
     }
  if(!missing(size)){
     max.predictors=size
     warning('Argument size has been replaced with max.predictors.
              Please update your code as usage of size will not be supported in
              future versions.')
     }

  # Validated here as well as in fit_model_set(), so an unrecognised value is
  # reported before the candidate set is built rather than after. Building the
  # set is not free: with several factors it fits a multinom() per ordered pair,
  # and with non.linear.correlations it fits a gam() per pair.
  r2.type <- match.arg(r2.type, c("r2.lm.est", "r2", "dev"))
  VI.mods <- match.arg(VI.mods, c("min.n", "all"))
  validate_progress(progress)

  model.set=generate_model_set(use.dat=use.dat,
                          test.fit=test.fit,
                          pred.vars.cont=pred.vars.cont,
                          pred.vars.fact=pred.vars.fact,
                          cyclic.vars= cyclic.vars,
                          linear.vars= linear.vars,
                          factor.smooth.interactions=factor.smooth.interactions,
                          factor.factor.interactions=factor.factor.interactions,
                          smooth.smooth.interactions=smooth.smooth.interactions,
                          cov.cutoff=cov.cutoff,
                          cor.matrix=cor.matrix,
                          non.linear.correlations=non.linear.correlations,
                          max.predictors=max.predictors,
                          k=k,
                          bs.arg=bs.arg,
                          null.terms=null.terms)

  out.dat=fit_model_set(model.set.list=model.set,
                          max.models=max.models,
                          save.model.fits=save.model.fits,
                          parallel=parallel,
                          n.cores=n.cores,
                          r2.type=r2.type,
                          report.unique.r2=report.unique.r2,
                          VI.mods=VI.mods,
                          progress=progress)

  # now return the list of outputs
  return(list(mod.data.out=out.dat$mod.data.out,
              used.data=model.set$used.data,
              predictor.correlations=model.set$predictor.correlations,
              #mod.formula=mod.formula,
              failed.models=out.dat$failed.models,
              success.models=out.dat$success.models,
              variable.importance=out.dat$variable.importance))
} #------------------ end function --------------------------------------------#