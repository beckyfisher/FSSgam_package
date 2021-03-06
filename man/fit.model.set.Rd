% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_fit_model_set.R
\name{fit.model.set}
\alias{fit.model.set}
\title{fit.model.set}
\usage{
fit.model.set(
  model.set.list,
  max.models = 200,
  save.model.fits = T,
  parallel = F,
  n.cores = 4,
  r2.type = "r2.lm.est",
  report.unique.r2 = F,
  VI.mods = "min.n"
)
}
\arguments{
\item{model.set.list}{A list as returned by generate.model.set()}

\item{max.models}{The total number of models allowed to be fit and the model fits to be saved during fitting. If the candidate set is bigger than this value, a warning message will be returned.}

\item{save.model.fits}{Should all successfully fitted models be saved to list success.models. Defaults to TRUE. Value is overwritten if the number of models n the set is bigger than max.models.}

\item{parallel}{A logical value indicating if parallel processing should be used. The default is FALSE.}

\item{r2.type}{The value to extract from the gam model fit to use as the R squared value. Defaults to r2.lm.est which returns and estimated R squared value based on a linear regression between the observed and predicted values. r2 will return the adjusted R.sq as reported by gam, gamm or gamm4.dev will return the deviance explained as reported by gam or gamm. Note gamm4 does not currently return a deviance.}

\item{report.unique.r2}{The estimated null model R2 is subtracted from each model R2 to give an idea of the unique variance explained. This can be useful where null terms are included in the model set.}
}
\value{
A list of the following output files:

mod.data.out - A data.frame that contains the statistics associated with each model fit. This includes AICc and BIC, delta values (e.g. AICc-(min(AICc)), corresponding weight values (Burnham and Anderson 2003), an estimate of the model R2, and a column for each of the included predictor variables containing either 0 (variable not included in the model) or 1 (variable is present in the model).
Use of BIC in information theoretic approaches has been heavily criticised because of the inherent assumption of BIC that there is a true model that is represented in the candidate set (Anderson & Burnham 2002). Rather than decide a-priori which model selection tool users should adopt, we supply both as part of the function outputs.
To simplify output, only AICc and AICc based model weights, rather than AIC, are included as these are asymptotically equivalent at large sample sizes, and for small sample sizes AICc should be used in any case.
Calculating R2 values is non-trivial for mixed models, especially non-gaussian cases (and some argue should not be done at all). We have supplied a range of methods for estimating R2 (r2.type), as in our experience a single method rarely performs adequately across all scenarios.

used.data - A data.frame which is identical to the data.frame initially supplied by the user, but with any hard coded interaction terms appended via cbind.

failed.models - A list of model formula that failed to fit. Ideally the list of failed models should be empty, but when this is not the case interrogating failed.models provides a useful means of troubleshooting. Users can examine which models are not fitting and explore the reasons for this by fitting the failed models outside the full.subsets.gam call based on the listed formula. When a large number of models fail to fit properly it usually indicates poor specification of the initial test.fit or other arguments in the call to full.subsets.gam (such as the inclusion of factor interactions when there are few data within each level of the factor), or that inappropriate variables are being included in the model set.

success.models - A complete list of all successfully fitted model formula. If models were saved, this can be used for multimodel inference and creating model averaged predictions.
otherwise the formula can be used to refit the top model set via a call to update of the test fit: update(test.fit,formula=mod.formula[[l]])

variable.importance - A list containing importance scores for each included predictor.
To determine the relative importance of each predictor across the whole model set we summed the wi values for all models containing each variable. 
The higher the combined weights for an explanatory parameter, the more important it is in the analysis (Burnham & Anderson, 2002). 
An assumption of the use of summed model weights to infer variable importance is that the number of models in which the different predictors are present is uniform. 
As our function removes models with correlated predictors, this is not always the case. 
To overcome this issue, the summed variable.importance scores are the summed weights for the best n models, where n is equal to the minimum number of models any one predictor is present in.
If you would like the variable importance scores to be based on all models in the set instead, set VI.mods="all".
}
\description{
Conducts a full subsets analysis based on gam(m4) using the list generated by a call to generate.model.set
}
\details{
The function constructs and fits a complete model set based on the supplied arguments.
for more information see Fisher R, Wilson SK, Sin TM, Lee AC, Langlois TJ (2018) A simple function for full-subsets multiple regression in ecology with R. Ecology and Evolution
https://onlinelibrary.wiley.com/doi/abs/10.1002/ece3.4134
}
