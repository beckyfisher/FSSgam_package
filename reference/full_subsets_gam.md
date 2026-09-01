# full_subsets_gam

Conducts a full subsets analysis based on gam(m4). In the most recent
version of FSSgam this function is now a wrapper for the two input
functions, generate_model_set and fit_model_set. Input arguments are the
same as these two underlying functions. calling these underlying
functions explicitly is the recommended method for running a full
subsets analysis with the FSSgam package because this enables the user
to interrogate the candidate model set and the predictor correlation
matrix before actually running the analysis.

## Usage

``` r
full_subsets_gam(
  use.dat,
  test.fit,
  pred.vars.cont = NA,
  pred.vars.fact = NA,
  cyclic.vars = NA,
  linear.vars = NA,
  factor.smooth.interactions = pred.vars.fact,
  factor.factor.interactions = FALSE,
  smooth.smooth.interactions = FALSE,
  cov.cutoff = 0.28,
  cor.matrix = NA,
  non.linear.correlations = FALSE,
  max.predictors = 3,
  k = 5,
  bs.arg = "'cr'",
  null.terms = "",
  max.models = 200,
  save.model.fits = TRUE,
  parallel = FALSE,
  n.cores = 4,
  r2.type = "r2.lm.est",
  report.unique.r2 = FALSE,
  VI.mods = "min.n",
  progress = interactive(),
  factor.interactions,
  smooth.interactions,
  size
)
```

## Arguments

- use.dat:

  A data.frame, with columns matching those included in pred.vars.cont
  and pred.vars.fact, the response variable to be analysed and any other
  fields required for the analysis (such as random effects, see
  test.fit). Note that any variables in use.dat that are used in model
  fits must not contain missing values, as this invalidates comparison
  via AICc/ BIC. If missing values occur among the predictor variables
  the function will return an error warning indicating that these rows
  need to be removed or interpolated.

- test.fit:

  A gam model fitted via a call to gam (mgcv) or uGamm (MuMIn). This can
  use any of the (preferably continuous) predictors in the call and will
  be used as a model to update in the fitting of the model set. The test
  fit must contain the appropriate random effects and call to family (if
  not gaussian) and if gamm4 should be used, or gamm in the case of a
  uGamm call (see ?uGamm). Both gamm from mgcv and gamm4 have slightly
  different features, as well as advantages and disadvantages, thus it
  is important that the full subsets function is able to deal with
  test.fit models based on either package. For example gamm4 is based on
  the lme4 package \[Bates, D.M. (2010) lme4: Mixed-Effects Modeling
  with R. Springer, New York\] which allows crossed random effects and
  avoids issues with PQL for non-gaussian model fits. On the other hand
  gamm (mgcv) is based on nlme which allows correlation structures
  \[Box, G.E.P., Jenkins, G.M., and Reinsel G.C. (1994) "Time Series
  Analysis: Forecasting and Control", 3rd Edition, Holden-Day\],
  variance structures \[Pinheiro, J.C. and Bates., D.M. (1996)
  "Unconstrained Parametrizations for Variance-Covariance Matrices",
  Statistics and Computing, 6, 289-296\], and a broader range of
  families that are not yet available in lmer (see ?family.mgcv). Models
  that have no random effects and are based only on gam (mgcv) are best
  fit via a direct call to gam, rather than using the uGamm wrapper.

- pred.vars.cont:

  A character vector indicating the continuous predictors to use. By
  default all continuous predictors will be fitted using a smoother (but
  see argument linear.vars). These must match column names in use.dat
  exactly. If NA is used the function can be run without any smooth
  predictors.

- pred.vars.fact:

  A character vector indicating the factor predictors to use. These must
  match column names in use.dat exactly. If NA is used the function can
  be run without any factor predictors.

- cyclic.vars:

  NA if there are no cyclic predictors, or if there are cyclic
  predictors, a character vector containing the names of any of the
  continuous predictors that should be modelled as cyclic variables.
  Note that these must also be contained in the pred.vars.cont character
  vector. Please also note there are issues with bs='cc' and model
  selection as this uses by default shrinkage. With shrinkage, variables
  are retained in models but with zero edf, which makes interpretation
  of AICc and BIC confusing. To account for this always select only the
  most parsimonious model (that with the fewest parameters), not just
  that with the lowest AICc. Reported estimated degrees of freedom (edf)
  in the model output table represent the sum of the edf of the smooth
  terms plus the number of parametric coefficients. When cyclic
  variables are included and shrinkage is used, any estimated edf of the
  smooth terms that are less than 1 are reset to 1 before summing to
  ensure the the total number of predictors in the model is captured
  properly.

- linear.vars:

  NA if there are no continuous predictors to be treated as linear (not
  fitted as smooths). Only use this where variables are clearly
  continuous in nature, but you are confident a linear relationship is
  valid. It may also be useful for continuous predictors that are not
  well distributed along the x-axis (ie, sampling was conducted in
  clumped distances from a feature of interest). Where this is
  necessary, transformations should be considered where they can be used
  to theoretically linearize response relationships. Does not need to be
  contained in vector pred.vars.cont

- factor.smooth.interactions:

  Default is the character vector pred.vars.fact, meaning that all
  factor predictors will be included as by arguments with all the
  continuous predictors. If factor.factor.interactions is TRUE, factor
  interactions variables will also be included as by arguments, yielding
  higher order interactions up to the specified model max.predictors. If
  a character vector is supplied, this must specify which factor
  variables should be included as "by" argument interaction terms with
  the continuous smooth predictors. If a list is supplied, this must be
  a named list containing the elements fact.vars, linear.vars,
  cont.vars, each a character vector indicating what predictors should
  be used to construct the factor smooth interactions. Note that
  specified factors, linear predictors and continuous predictors must
  also be included in their respective character vectors
  (pred.vars.fact, linear.vars, pred.vars.cont). If specified as NA no
  factor-continuous predictor interactions will be included.

- factor.factor.interactions:

  A logical value indicating if interactions between factors should be
  included, or only their main effects. Defaults to FALSE. Note that
  this can substantially increase the number of models in the candidate
  set. Not recommended when there are factors with many levels.
  Alternatively character vector specifying which factor predictors to
  include as interactions. These must be contained within
  pred.vars.fact. If factor.factor.interactions is set to TRUE the
  function automatically generates hard coded interaction variables up
  to the maximum number of predictors (see max.predictors below) using
  combn. New factors are generated by pasting the resulting unique
  combinations together. This method of generating interaction terms is
  necessary because smooth-factor interactions are specified as by
  arguments in calls to gam(m4). Because the full subsets function
  automatically checks for collinearity, there is no issue with
  constructing model sets with multiple factor arguments that are higher
  order factors of each other, as these are invariably collinear and
  subsequently removed (see cov.cutoff below). If a user defined
  cor.matrix is passed to the function (see cor.matrix) this must
  include these hard coded interactions.

- smooth.smooth.interactions:

  A logical value indicating if the function should include te smooths
  of second order continuous predictor interactions. If set to TRUE, all
  continuous predictors will be combined as bivariate calls to te.
  Alternatively character vector specifying which continuous predictors
  to include. These must be contained within pred.vars.cont.

- cov.cutoff:

  A numeric value between 0 and 1 indicating the correlation cutoff
  value to use for excluding collinear models, based on the cor.matrix
  (see below). The default value is 0.28 (see Graham MH (2003). It is
  highly recommended to keep this value low, as correlation among
  predictors can yield spurious results. Note that predictors with a
  correlation greater than the specified value will still appear in the
  model set but will never appear in the same model. Including highly
  correlated predictors can make interpreting variable importance values
  difficult.

- cor.matrix:

  A user-supplied pairwise correlation matrix, or NA (the default) to
  compute one from use.dat. When supplied it governs every stage that
  screens on correlation: which factor-factor interaction columns are
  built, which te smooth-smooth interaction terms are built, and which
  assembled candidate models are excluded. It must therefore carry a row
  and a column for every predictor, including any hard coded factor
  interactions that setting factor.factor.interactions causes to be
  created; missing names are reported by name. When supplied it replaces
  the automatic estimate rather than overriding it, so
  check_correlations is not called at all and a predictor of a class it
  does not accept can be used. By default predictor correlations are
  evaluated via a call to check_correlations, a function taking a
  data.frame (containing all predictors) as argument and generating a
  correlation matrix comprised of: 1) correlation coefficients between
  all continuous predictors via a call to cor; 2) approximate
  correlation values between continuous predictors and factors, as the
  square root of the R2 value obtained via a call to lm, where the
  continuous predictor is modelled as a response and the factor variable
  as a single fixed factor; and 3) approximate correlations values
  between factor predictors, as the square root of the R2 value obtained
  via a call multinom (from package nnet, Venables & Ripley 2002). Note
  that any user constructed pairwise matrix can be passed to the
  function and used for pairwise exclusion of variables from individual
  models.

- non.linear.correlations:

  Set this argument to TRUE if you would like to exclude continuous
  predictor combinations that are potentially "correlated" through
  non-linear relationships. See ?check_non_linear_correlations for more
  details.

- max.predictors:

  An integer indicating the maximum number of predictors to include in
  any one model.

- k:

  An integer indicating the dimension of the basis used to represent the
  smooth term (see ?s). The default value is 5. Higher values are not
  recommended unless a complex trend between the response variable and
  the continuous predictor variables is expected, and the data are
  sufficient to support this. k can be reduced to as low as 3 where
  there is trouble obtaining convergence, or sample size is low. Note
  that this must be set to override the default value, regardless of
  what k is used in the test.fit

- bs.arg:

  Specification of the smoother to use, see ?s for more information on
  smoother provided in gam (mgcv). Note that all continuous predictors
  specified in pred.vars.cont will be fitted using the same smooth,
  unless they are also specified as linear.terms or cyclic.vars. Note
  that any specification of bs in test.fit is discarded.

- null.terms:

  A character vector indicating the form of any re smooths to be
  included in gam \[e.g. s(site,bs=re)\] or any other fixed terms or
  smooths that the user wants to include in the null model. Use of bs=re
  is an alternative way of fitting simple random structures that avoids
  use of PQL and allows a the greater range of families available in
  gam.mgcv to be used. see ?s and links therein. Note: make sure you use
  gam instead of uGamm to make sure PQL is not used. to fit a
  correlation structure only (but no random effects) this must be
  achieved through a call the gamm via , with no random effects, fitted

- max.models:

  The total number of models allowed to be fit and still save the model
  fits. Defaults to 200. If the candidate set is bigger than this value,
  a warning will be returned indicating that model fits will not be
  saved.

- save.model.fits:

  Are the model fits to be saved in the output list? If TRUE this will
  be overwritten if the model candidate set is bigger than max.models.
  If FALSE only model output data are saved.

- parallel:

  A logical value indicating if parallel processing should be used. The
  default is FALSE.

- n.cores:

  An integer indicating the number of cores to use if parallel is TRUE.
  Defaults to 4.

- r2.type:

  The value to extract from the gam model fit to use as the R squared
  value. Defaults to r2.lm.est which returns and estimated R squared
  value based on a linear regression between the observed and predicted
  values. r2 will return the adjusted R.sq as reported by gam, gamm or
  gamm4.dev will return the deviance explained as reported by gam or
  gamm. Note gamm4 does not currently return a deviance.

- report.unique.r2:

  The estimated null model R2 is subtracted from each model R2 to give
  an idea of the unique variance explained. This can be useful where
  null terms are included in the model set.

- VI.mods:

  The set of models used to calculate summed variable importance scores.
  Defaults to 'min.n', which uses only the best n models for each
  variable (n being the minimum number of models any one predictor is
  present in). Set to 'all' to use all models in the candidate set
  instead.

- progress:

  Should a text progress bar be written to the console while models are
  fitted. Defaults to interactive(), so the bar appears at the console
  but not in scripts, reports or checks.

- factor.interactions:

  Deprecated. Superseded by factor.factor.interactions; retained only so
  older code does not break, and will warn if used.

- smooth.interactions:

  Deprecated. Superseded by factor.smooth.interactions; retained only so
  older code does not break, and will warn if used.

- size:

  Deprecated. Superseded by max.predictors; retained only so older code
  does not break, and will warn if used.

## Value

A list of the following output files:

mod.data.out - A data.frame that contains the statistics associated with
each model fit. This includes AICc and BIC, delta values (e.g.
AICc-(min(AICc)), corresponding weight values (Burnham and Anderson
2003), an estimate of the model R2, and a column for each of the
included predictor variables containing either 0 (variable not included
in the model) or 1 (variable is present in the model). Use of BIC in
information theoretic approaches has been heavily criticised because of
the inherent assumption of BIC that there is a true model that is
represented in the candidate set (Anderson & Burnham 2002). Rather than
decide a-priori which model selection tool users should adopt, we supply
both as part of the function outputs. To simplify output, only AICc and
AICc based model weights, rather than AIC, are included as these are
asymptotically equivalent at large sample sizes, and for small sample
sizes AICc should be used in any case. Calculating R2 values is
non-trivial for mixed models, especially non-gaussian cases (and some
argue should not be done at all). We have supplied a range of methods
for estimating R2 (r2.type), as in our experience a single method rarely
performs adequately across all scenarios.

used.data - A data.frame which is identical to the data.frame initially
supplied by the user, but with any hard coded interaction terms appended
via cbind.

predictor.correlations - The matrix of estimated predictor correlations
returned by the function check_correlations and used for model exclusion
based on cov.cutoff

failed.models - A list containing the try-error catch associated with
models that failed to fit. Ideally the list of failed models should be
empty, but when this is not the case interrogating failed.models
provides a useful means of troubleshooting. Users can examine which
models are not fitting and explore the reasons for this by fitting the
failed models outside the full_subsets_gam call based on the listed
formula. When a large number of models fail to fit properly it usually
indicates poor specification of the initial test.fit or other arguments
in the call to full_subsets_gam (such as the inclusion of factor
interactions when there are few data within each level of the factor),
or that inappropriate variables are being included in the model set.

success.models - A complete list of all successfully fitted models. This
can be used for multimodel inference and creating model averaged
predictions.

variable.importance - A list containing importance scores for each
included predictor. To determine the relative importance of each
predictor across the whole model set we summed the wi values for all
models containing each variable. The higher the combined weights for an
explanatory parameter, the more important it is in the analysis (Burnham
& Anderson, 2002). An assumption of the use of summed model weights to
infer variable importance is that the number of models in which the
different predictors are present is uniform. As our function removes
models with correlated predictors, this is not always the case. To
overcome this issue, the summed variable.importance scores are the
summed weights for the best n models, where n is equal to the minimum
number of models any one predictor is present in.

## Details

The function constructs and fits a complete model set based on the
supplied arguments. for more information see Fisher R, Wilson SK, Sin
TM, Lee AC, Langlois TJ (2018) A simple function for full-subsets
multiple regression in ecology with R. Ecology and Evolution
https://onlinelibrary.wiley.com/doi/abs/10.1002/ece3.4134

## Examples

``` r
library(mgcv)
data(case_study1)
use.dat <- case_study1
use.dat$site <- as.factor(use.dat$site)
test.fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr") + s(site, bs = "re"),
                 family = tw(), data = use.dat)
full_subsets_gam(
  use.dat = use.dat,
  test.fit = test.fit,
  pred.vars.cont = c("complexity", "depth"),
  pred.vars.fact = "ZONE",
  null.terms = "s(site,bs='re')",
  max.predictors = 2,
  k = 3
)
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
#> $used.data
#>         ZONE         site depth        SA complexity rugosity   LC   HC macro
#> 1     FISHED         MESA   3.2 321.83444        3.5 1.643836  0.0  0.0   0.0
#> 2     FISHED         MESA   3.9 309.10746        4.5 1.714286 20.0 16.0   2.0
#> 3     FISHED         MESA   3.9 118.71986        2.5 1.395349 38.0  0.0   0.0
#> 4     FISHED         MESA   1.7 252.36261        1.0 1.045296  0.0  0.0  67.0
#> 5     FISHED         MESA   3.2 132.35994        0.5 1.052632  0.0  0.0  69.0
#> 6  SANCTUARY     MANGROVE   2.1 149.24581        3.0 1.967213  4.0  6.0   6.0
#> 7  SANCTUARY     MANGROVE   2.5 344.89786        3.5 1.538462 24.0  2.0  24.0
#> 8  SANCTUARY     MANGROVE   2.9  61.93655        2.0 1.363636  4.0  1.0  20.0
#> 9  SANCTUARY     MANGROVE   3.1 380.20034        1.0 1.200000  2.0  6.0  75.0
#> 10 SANCTUARY     MANGROVE   3.7 188.89879        2.5 1.626016 48.6  0.0   8.0
#> 11 SANCTUARY     MANGROVE   4.5 436.68381        3.5 2.181818  3.0  5.6  18.4
#> 12    FISHED NTH MANGROVE   3.5  83.51676        2.5 1.846154  0.0 80.0   0.0
#> 13 SANCTUARY     MANGROVE   3.1  54.72278        3.5 1.363636 34.0  0.0   9.6
#> 14    FISHED NTH MANGROVE   3.0  88.19923        2.5 1.237113  0.0  5.6   0.0
#> 15    FISHED NTH MANGROVE   3.3 110.42952        2.0 1.239669  3.6 14.0   4.0
#> 16    FISHED NTH MANGROVE   4.1  62.63569        1.0 1.149425  0.0 40.0   0.0
#> 17    FISHED         MESA   2.0  68.95717        2.0 1.357466  7.0  3.0  17.0
#> 18    FISHED         MESA   3.0 136.83400        3.5 1.428571 40.0 14.0   2.0
#> 19    FISHED         MESA   2.0 517.59368        0.0 1.034483  0.0  0.0  22.0
#> 20    FISHED         MESA   2.1 405.58800        0.0 1.094891  0.0  5.0  74.0
#> 21    FISHED         MESA   3.2 136.11094        4.0 1.769912 48.0  0.0   0.0
#> 22 SANCTUARY        MANDU   2.1 164.04628        1.0 1.158301  0.0  6.0  46.0
#> 23 SANCTUARY        MANDU   2.2 257.49206        3.0 1.500000 28.0  4.0  34.0
#> 24 SANCTUARY        MANDU   1.7 296.18133        3.0 1.463415 80.0  0.0   0.0
#> 25 SANCTUARY        MANDU   1.4 288.81810        0.0 1.071429  0.0  0.0  26.0
#> 26    FISHED    NTH MANDU   3.2 191.68524        3.5 1.518987 26.0  6.0   4.0
#> 27    FISHED         MESA   2.1 115.56878        0.5 1.081081  0.0  0.0  22.0
#> 28 SANCTUARY       JURABI   3.8 167.13997        3.0 1.415094 18.0  6.0   0.0
#> 29 SANCTUARY       JURABI   3.7 235.32776        3.5 1.284797  0.0 28.0   0.0
#> 30 SANCTUARY       JURABI   3.6 224.16101        2.0 1.086957  0.0  0.0  52.0
#> 31 SANCTUARY       JURABI   3.4 475.01497        1.0 1.140684  0.0  0.0  34.0
#> 32 SANCTUARY       JURABI   4.1 100.49147        2.0 1.310044  6.0  2.0   0.0
#> 33 SANCTUARY       JURABI   3.6 109.62320        1.5 1.181102  0.0  0.0  22.0
#> 34    FISHED         MESA   3.7  74.20490        3.0 1.276596  6.0  4.0   4.0
#> 35    FISHED         MESA   2.9 135.40590        4.0 1.304348  8.0 10.0   0.0
#> 36    FISHED         MESA   2.1  98.84527        0.5 1.034483  0.0  0.0  44.0
#> 37 SANCTUARY     MANGROVE   2.6  32.99946        3.0 1.176471  6.0  2.0  10.0
#> 38 SANCTUARY     MANGROVE   1.8  71.59080        0.0 1.008403  0.0  0.0  30.0
#> 39 SANCTUARY     MANGROVE   2.2  76.09548        2.5 1.333333 12.0  2.0  24.0
#> 40 SANCTUARY     MANGROVE   2.2 159.26292        3.5 1.518987  8.0  0.0   2.0
#> 41 SANCTUARY     MANGROVE   1.8 339.83030        4.0 1.428571  0.0 12.0   8.0
#> 42    FISHED NTH MANGROVE   3.6  51.51412        2.0 1.071429  0.0  4.0  70.0
#> 43    FISHED NTH MANGROVE   2.7  94.46148        3.0 1.428571  0.0  0.0  16.0
#> 44    FISHED NTH MANGROVE   2.5 248.59391        3.0 1.538462  0.0  0.0  58.0
#> 45    FISHED NTH MANGROVE   2.0  90.19713        0.0 1.034483  0.0  0.0  62.0
#> 46 SANCTUARY        MANDU   1.8  23.64948        1.5 1.363636  0.0  0.0  52.0
#> 47 SANCTUARY        MANDU   2.4 142.37115        2.5 1.428571  7.0  8.0   8.0
#> 48 SANCTUARY        MANDU   2.0 146.10886        3.5 2.068966  0.0 30.0   0.0
#> 49    FISHED    NTH MANDU   2.1 131.59072        2.5 1.463415  0.0  4.0  30.0
#> 50    FISHED    NTH MANDU   2.1 245.86162        3.5 1.621622  0.0 14.0   4.0
#> 51    FISHED    NTH MANDU   1.8 104.98107        0.5 1.025641  0.0  0.0  42.0
#> 52 SANCTUARY       JURABI   3.5  68.67859        3.5 1.578947  6.0  2.0   0.0
#> 53 SANCTUARY       JURABI   3.1  93.17857        2.5 1.304348 24.0  0.0   4.0
#> 54 SANCTUARY       JURABI   3.9  16.65045        2.0 1.250000  6.0  1.0   0.0
#> 55 SANCTUARY        MANDU   2.2  35.70543        2.5 1.666667  0.0 32.0   0.0
#> 56 SANCTUARY        MANDU   1.8 189.27931        3.5 1.538462 16.0 12.0   0.0
#> 57 SANCTUARY        MANDU   2.1 125.48754        0.5 1.090909  0.0  0.0  58.0
#> 58 SANCTUARY        MANDU   1.7  74.33842        0.5 1.111111  0.0  4.0  30.0
#> 59 SANCTUARY        MANDU   2.4 121.48344        3.0 1.578947  0.0 22.0   0.0
#> 60    FISHED NTH MANGROVE   2.7  22.45597        1.5 1.071429  0.0  0.0  30.0
#> 61    FISHED NTH MANGROVE   4.0 128.19227        3.5 1.538462 20.0 24.0   0.0
#> 62    FISHED NTH MANGROVE   3.7 323.68014        3.0 1.463415  0.0  2.0   0.0
#> 63    FISHED NTH MANGROVE   3.5 181.77032        3.0 1.714286  0.0  6.0   0.0
#> 64    FISHED    NTH MANDU   1.8 178.87941        1.5 1.200000  0.0  0.0  64.0
#> 65    FISHED    NTH MANDU   2.1 106.65742        3.0 1.333333  0.0  0.0  36.0
#> 66    FISHED    NTH MANDU   2.2 247.92257        3.0 1.463415  4.0 10.0   0.0
#> 67    FISHED    NTH MANDU   2.9  51.58278        3.0 1.875000  4.0  0.0   0.0
#> 68    FISHED    NTH MANDU   3.1 385.04238        4.0 1.363636 22.0  0.0  34.0
#>     SCORE1   SCORE2 Herbivore.abundance Invertivore.abundance
#> 1   1.2200 -0.00782                  18                    54
#> 2   2.4400 -0.00349                  58                   252
#> 3   1.1000 -1.59000                  22                    27
#> 4  -2.7200 -0.15700                   3                    13
#> 5  -2.9600 -0.11600                   0                     4
#> 6   1.4600  0.24200                  85                   210
#> 7   1.0600 -1.02000                  46                   146
#> 8  -0.2760 -0.10400                  19                    25
#> 9  -2.2400  0.13500                  34                    64
#> 10  1.5500 -2.03000                  14                    50
#> 11  1.5900  0.19100                  48                    91
#> 12  2.6300  4.87000                  13                    66
#> 13  1.2000 -1.55000                  30                    44
#> 14  0.0904  0.33600                  38                    88
#> 15  0.0275  0.71400                  18                    50
#> 16 -0.1380  2.49000                  11                    27
#> 17 -0.1250 -0.10400                  19                    30
#> 18  1.8900 -0.93500                  52                    66
#> 19 -2.2400  0.09120                   0                    11
#> 20 -3.0500  0.22400                   6                    14
#> 21  2.5800 -2.09000                  39                    65
#> 22 -1.7800  0.31200                  32                    31
#> 23  0.6730 -1.07000                  49                    99
#> 24  2.3000 -3.42000                  48                    46
#> 25 -2.1900  0.09030                   3                    18
#> 26  1.5700 -0.80000                  54                   163
#> 27 -1.8500  0.06200                   3                    27
#> 28  1.0900 -0.41700                  20                    39
#> 29  1.1000  1.59000                  43                    68
#> 30 -1.8100 -0.18100                  20                    27
#> 31 -1.6900 -0.00734                   2                    27
#> 32  0.0951 -0.07070                  11                     9
#> 33 -1.0900  0.00160                  13                    28
#> 34  0.4150 -0.06610                  27                    45
#> 35  1.1700  0.13400                  63                    54
#> 36 -2.4900 -0.03330                   6                    13
#> 37 -0.0195 -0.23300                  33                    38
#> 38 -2.5100  0.05260                   0                     6
#> 39 -0.0260 -0.45400                  30                    25
#> 40  1.1400 -0.37800                  38                    71
#> 41  1.1400  0.59400                  93                   277
#> 42 -2.1700 -0.01280                   3                    47
#> 43  0.2880 -0.05580                  27                    55
#> 44 -0.4180 -0.18600                  30                    68
#> 45 -3.1000 -0.05160                   0                     9
#> 46 -1.2900 -0.05990                  12                    21
#> 47  0.5370  0.19400                  29                    34
#> 48  2.3100  1.83000                  69                    43
#> 49 -0.0907  0.18500                  47                    80
#> 50  1.3800  0.80700                  98                   169
#> 51 -2.4800 -0.02950                   1                     8
#> 52  1.2800 -0.15600                  44                    25
#> 53  0.5390 -1.02000                  46                    48
#> 54 -0.0707 -0.14500                   8                    19
#> 55  1.4400  1.99000                  30                    29
#> 56  1.6100  0.00253                  35                    68
#> 57 -2.5900 -0.06310                   2                     7
#> 58 -1.8400  0.28100                   2                     8
#> 59  1.3400  1.34000                  13                    27
#> 60 -1.6100 -0.06180                   0                    10
#> 61  1.9300  0.54400                  24                    55
#> 62  0.7380  0.12700                  61                   153
#> 63  1.2200  0.40500                  54                    72
#> 64 -1.9400 -0.14300                  10                    32
#> 65 -0.3420 -0.14700                  40                    44
#> 66  0.9780  0.43100                  71                    55
#> 67  1.3800 -0.10300                  20                    24
#> 68  0.6530 -1.17000                  86                   134
#>    Piscivore.abundance Planktivore.abundance Herbivore.biomass
#> 1                   30                   288        1221.26602
#> 2                   66                   196       12321.34367
#> 3                    1                   122        9475.76226
#> 4                    4                     0         473.19268
#> 5                    0                     0           0.00000
#> 6                    8                    80        8672.12624
#> 7                   65                   286        4438.89517
#> 8                    1                    47         461.30056
#> 9                    2                    64         589.29566
#> 10                   1                    72        1116.00484
#> 11                  54                   209        5330.11754
#> 12                  19                   550         659.49001
#> 13                   2                    75       13663.61201
#> 14                   4                   135       19267.44368
#> 15                  14                   136        1277.21068
#> 16                   2                    14         142.24824
#> 17                   5                     8         345.39049
#> 18                  59                   624        2094.78041
#> 19                   0                     1           0.00000
#> 20                   5                    11        1658.89489
#> 21                   3                   241        7390.80169
#> 22                   5                    23        1772.46304
#> 23                  18                   129       12806.41833
#> 24                   2                    29        8343.71297
#> 25                   3                     1         292.04715
#> 26                  38                   280        5886.49529
#> 27                   4                     0         351.29636
#> 28                   5                     0        2941.87591
#> 29                   6                    41        8358.85011
#> 30                   5                    42        2130.57399
#> 31                   9                    76         434.44048
#> 32                   1                    22         693.31316
#> 33                   4                    12        1109.14863
#> 34                  12                   262        2371.21813
#> 35                   6                   194        9025.61898
#> 36                  30                     3          89.99833
#> 37                   5                    11        4690.14049
#> 38                   0                     0           0.00000
#> 39                  12                   132        4977.87090
#> 40                  12                   199        7806.16908
#> 41                  46                   205       10043.70399
#> 42                   4                    12         915.29957
#> 43                   3                    89        3901.58471
#> 44                   6                    86        7624.26587
#> 45                   0                     0           0.00000
#> 46                   1                     2         697.54100
#> 47                   1                    50        1032.75346
#> 48                   4                    33       12744.96968
#> 49                   9                    14        5507.43485
#> 50                  24                    90       16834.03870
#> 51                   0                     0         534.62224
#> 52                   3                     3        9472.66874
#> 53                  35                   602       13787.90479
#> 54                   2                     4        3129.13641
#> 55                   4                    45        4958.27086
#> 56                  10                    55        2962.86207
#> 57                   2                     6         545.63788
#> 58                   3                     1         288.63227
#> 59                   4                    27        1345.25981
#> 60                   1                    33           0.00000
#> 61                  10                   150        4581.57478
#> 62                 229                   106       15379.43420
#> 63                   7                   107       15925.84843
#> 64                   9                    22        1049.28925
#> 65                  12                    61        3333.30344
#> 66                   2                   218        8010.40146
#> 67                   8                   105        1099.82292
#> 68                  12                   202       19124.63243
#>    Invertivore.biomass Piscivore.biomass Planktivore.biomass
#> 1           4724.72394       3364.995174        10444.547350
#> 2          22065.53828      11102.265070         2565.598694
#> 3          16687.42007         65.032060         2025.822123
#> 4            633.41995        437.724781            0.000000
#> 5            445.45148          0.000000            0.000000
#> 6          15458.46719        232.574849         3324.792220
#> 7          21960.41299       7486.252267         6996.160422
#> 8           1097.15729         43.188831          352.198147
#> 9           3415.47996         11.033823          212.875592
#> 10          5970.92133        251.552148          438.826833
#> 11         14939.55962       7528.882252         3271.764668
#> 12          4364.26506       1542.962813         4834.142689
#> 13         15135.23860        303.180173         1800.770512
#> 14         14044.66439        918.498976         1588.884565
#> 15          6743.62665       4581.361991         4630.976688
#> 16           603.13808        425.995509           40.903718
#> 17           838.23684        741.730954           77.157871
#> 18          4639.48603       8737.279021         8570.154587
#> 19           584.22670          0.000000            2.447276
#> 20          5272.54524       1248.245089          120.638318
#> 21          3966.23542        444.035931         5894.024441
#> 22          1332.31205        577.823779          121.797628
#> 23         11206.83537        787.674264         1123.240980
#> 24          3698.10591       1145.983808          381.785065
#> 25           191.01323        688.734817            2.093062
#> 26          9105.21989      72484.779400         7022.747589
#> 27           529.28767        287.208804            0.000000
#> 28          4535.50384        451.993676            0.000000
#> 29          4981.38080       1131.676338         1953.953348
#> 30          6023.06224       1457.164138          558.600000
#> 31          3860.00801       2670.895687          255.659034
#> 32           371.18500        124.827796          501.382200
#> 33          3050.82738        813.883554           41.751585
#> 34          2835.56622       1950.378852         8509.003663
#> 35          3144.97597        816.538835         2049.739582
#> 36            66.96136        118.567374           11.701137
#> 37          2915.73753        317.044032          272.158655
#> 38            83.79845          0.000000            0.000000
#> 39          1073.80203       1507.945419         3190.792956
#> 40          4955.72915        227.332930         3940.780294
#> 41         13920.16875       6945.401839         5342.955143
#> 42          1047.94820        412.032981          157.707413
#> 43         15067.24497        389.987748         2751.088325
#> 44          5038.40157        580.602071         3185.037271
#> 45            16.52421          0.000000            0.000000
#> 46           618.24909         27.129642            4.186124
#> 47          2020.89418          6.435908          154.957567
#> 48          5192.60000       1008.030507          234.478300
#> 49          2388.71761       1052.781185          129.659282
#> 50          7608.52332       2359.767985         1162.489716
#> 51           387.18156          0.000000            0.000000
#> 52          6984.73577        741.409371          123.156957
#> 53          8292.37163       4674.571751         6934.966326
#> 54          1162.33926        756.544266          121.333230
#> 55           137.22880       1138.818170          212.806942
#> 56          7473.90244        793.142385          718.255014
#> 57           848.03594        239.717133           29.667511
#> 58           933.93018        632.933705            1.347971
#> 59          3453.08026         32.687797          177.364700
#> 60            74.97125          7.537407          299.910347
#> 61          9062.73501       2176.892127          552.800944
#> 62         20539.10754      50338.413700         4229.137273
#> 63          2805.71613        560.008909         1251.178426
#> 64          2754.78323        724.264850          489.503660
#> 65          1789.92492        251.462539          603.152296
#> 66          4448.31141        449.248587         1083.167944
#> 67          1600.51803        896.055652         2111.212858
#> 68        144429.20740       1233.688056         4373.699001
#>    sqrt.Herbivore.abundance sqrt.Invertivore.abundance sqrt.Piscivore.abundance
#> 1                  4.242641                   7.348469                 5.477226
#> 2                  7.615773                  15.874508                 8.124038
#> 3                  4.690416                   5.196152                 1.000000
#> 4                  1.732051                   3.605551                 2.000000
#> 5                  0.000000                   2.000000                 0.000000
#> 6                  9.219544                  14.491377                 2.828427
#> 7                  6.782330                  12.083046                 8.062258
#> 8                  4.358899                   5.000000                 1.000000
#> 9                  5.830952                   8.000000                 1.414214
#> 10                 3.741657                   7.071068                 1.000000
#> 11                 6.928203                   9.539392                 7.348469
#> 12                 3.605551                   8.124038                 4.358899
#> 13                 5.477226                   6.633250                 1.414214
#> 14                 6.164414                   9.380832                 2.000000
#> 15                 4.242641                   7.071068                 3.741657
#> 16                 3.316625                   5.196152                 1.414214
#> 17                 4.358899                   5.477226                 2.236068
#> 18                 7.211103                   8.124038                 7.681146
#> 19                 0.000000                   3.316625                 0.000000
#> 20                 2.449490                   3.741657                 2.236068
#> 21                 6.244998                   8.062258                 1.732051
#> 22                 5.656854                   5.567764                 2.236068
#> 23                 7.000000                   9.949874                 4.242641
#> 24                 6.928203                   6.782330                 1.414214
#> 25                 1.732051                   4.242641                 1.732051
#> 26                 7.348469                  12.767145                 6.164414
#> 27                 1.732051                   5.196152                 2.000000
#> 28                 4.472136                   6.244998                 2.236068
#> 29                 6.557439                   8.246211                 2.449490
#> 30                 4.472136                   5.196152                 2.236068
#> 31                 1.414214                   5.196152                 3.000000
#> 32                 3.316625                   3.000000                 1.000000
#> 33                 3.605551                   5.291503                 2.000000
#> 34                 5.196152                   6.708204                 3.464102
#> 35                 7.937254                   7.348469                 2.449490
#> 36                 2.449490                   3.605551                 5.477226
#> 37                 5.744563                   6.164414                 2.236068
#> 38                 0.000000                   2.449490                 0.000000
#> 39                 5.477226                   5.000000                 3.464102
#> 40                 6.164414                   8.426150                 3.464102
#> 41                 9.643651                  16.643317                 6.782330
#> 42                 1.732051                   6.855655                 2.000000
#> 43                 5.196152                   7.416198                 1.732051
#> 44                 5.477226                   8.246211                 2.449490
#> 45                 0.000000                   3.000000                 0.000000
#> 46                 3.464102                   4.582576                 1.000000
#> 47                 5.385165                   5.830952                 1.000000
#> 48                 8.306624                   6.557439                 2.000000
#> 49                 6.855655                   8.944272                 3.000000
#> 50                 9.899495                  13.000000                 4.898979
#> 51                 1.000000                   2.828427                 0.000000
#> 52                 6.633250                   5.000000                 1.732051
#> 53                 6.782330                   6.928203                 5.916080
#> 54                 2.828427                   4.358899                 1.414214
#> 55                 5.477226                   5.385165                 2.000000
#> 56                 5.916080                   8.246211                 3.162278
#> 57                 1.414214                   2.645751                 1.414214
#> 58                 1.414214                   2.828427                 1.732051
#> 59                 3.605551                   5.196152                 2.000000
#> 60                 0.000000                   3.162278                 1.000000
#> 61                 4.898979                   7.416198                 3.162278
#> 62                 7.810250                  12.369317                15.132746
#> 63                 7.348469                   8.485281                 2.645751
#> 64                 3.162278                   5.656854                 3.000000
#> 65                 6.324555                   6.633250                 3.464102
#> 66                 8.426150                   7.416198                 1.414214
#> 67                 4.472136                   4.898979                 2.828427
#> 68                 9.273618                  11.575837                 3.464102
#>    sqrt.Planktivore.abundance log.Herbivore.biomass log.Invertivore.biomass
#> 1                   16.970563              3.087166                3.674468
#> 2                   14.000000              4.090693                4.343734
#> 3                   11.045361              3.976660                4.222415
#> 4                    0.000000              2.675955                2.802377
#> 5                    0.000000              0.000000                2.649774
#> 6                    8.944272              3.938176                4.189195
#> 7                   16.911535              3.647373                4.341660
#> 8                    6.855655              2.664924                3.040665
#> 9                    8.000000              2.771070                3.533579
#> 10                   8.485281              3.048055                3.776114
#> 11                  14.456832              3.726818                4.174367
#> 12                  23.452079              2.819866                3.640011
#> 13                   8.660254              4.135597                4.180018
#> 14                  11.618950              4.284847                4.147542
#> 15                  11.661904              3.106602                3.828958
#> 16                   3.741657              2.156089                2.781136
#> 17                   2.828427              2.539566                2.923885
#> 18                  24.979992              3.321346                3.666563
#> 19                   1.000000              0.000000                2.767324
#> 20                   3.316625              3.220081                3.722103
#> 21                  15.524175              3.868750                3.598488
#> 22                   4.795832              3.248822                3.124932
#> 23                  11.357817              4.107462                4.049522
#> 24                   5.385165              3.921411                3.568097
#> 25                   1.000000              2.466938                2.283331
#> 26                  16.733201              3.769931                3.959338
#> 27                   0.000000              2.546908                2.724512
#> 28                   0.000000              3.468772                3.656721
#> 29                   6.403124              3.922198                3.697437
#> 30                   6.480741              3.328700                3.779889
#> 31                   8.717798              2.638929                3.586701
#> 32                   4.690416              2.841555                2.570759
#> 33                   3.464102              3.045381                3.484560
#> 34                  16.186414              3.375155                3.452793
#> 35                  13.928388              3.955525                3.497755
#> 36                   1.732051              1.959033                1.832262
#> 37                   3.316625              3.671278                3.464897
#> 38                   0.000000              0.000000                1.928388
#> 39                  11.489125              3.697131                3.031328
#> 40                  14.106736              3.892494                3.695195
#> 41                  14.317821              4.001937                4.143676
#> 42                   3.464102              2.962037                3.020754
#> 43                   9.433981              3.591352                4.178063
#> 44                   9.273618              3.882255                3.702379
#> 45                   0.000000              0.000000                1.243639
#> 46                   1.414214              2.844192                2.791865
#> 47                   7.071068              3.014417                3.305758
#> 48                   5.744563              4.105373                3.715468
#> 49                   3.741657              3.741028                3.378347
#> 50                   9.486833              4.226214                3.881357
#> 51                   0.000000              2.728859                2.589035
#> 52                   1.732051              3.976518                3.844212
#> 53                  24.535688              4.139530                3.918731
#> 54                   2.000000              3.495563                3.065706
#> 55                   6.708204              3.695418                2.140599
#> 56                   7.416198              3.471858                3.873606
#> 57                   2.449490              2.737700                2.928926
#> 58                   1.000000              2.461847                2.970779
#> 59                   5.196152              3.129129                3.538332
#> 60                   5.744563              0.000000                1.880649
#> 61                  12.247449              3.661110                3.957307
#> 62                  10.295630              4.186969                4.312603
#> 63                  10.344080              4.202130                3.448198
#> 64                   4.690416              3.021309                3.440245
#> 65                   7.810250              3.523005                3.253077
#> 66                  14.764823              3.903708                3.648293
#> 67                  10.246951              3.041717                3.204532
#> 68                  14.212670              4.281616                5.159658
#>    log.Piscivore.biomass log.Planktivore.biomass
#> 1               3.527113               4.0189312
#> 2               4.045451               3.4093580
#> 3               1.819755               3.3068156
#> 4               2.642192               0.0000000
#> 5               0.000000               0.0000000
#> 6               2.368426               3.5218951
#> 7               3.874322               3.8449218
#> 8               1.645313               2.5480184
#> 9               1.080404               2.3301612
#> 10              2.402351               2.6432817
#> 11              3.876788               3.5149148
#> 12              3.188637               3.6844093
#> 13              2.483131               3.2556995
#> 14              2.963551               3.2013656
#> 15              3.661089               3.6657664
#> 16              2.630423               1.6222526
#> 17              2.870832               1.8929727
#> 18              3.941426               3.9330393
#> 19              0.000000               0.5374761
#> 20              3.096648               2.0850704
#> 21              2.648395               3.7704856
#> 22              2.762546               2.0891900
#> 23              2.896898               3.0508594
#> 24              3.059557               2.5829550
#> 25              2.838682               0.4903886
#> 26              4.860253               3.8465689
#> 27              2.459707               0.0000000
#> 28              2.656092               0.0000000
#> 29              3.054106               3.2911364
#> 30              3.163806               2.7478777
#> 31              3.426819               2.4093566
#> 32              2.099777               2.7010342
#> 33              2.911096               1.6309522
#> 34              3.290342               3.9299297
#> 35              2.912508               3.3119105
#> 36              2.077613               1.1038426
#> 37              2.502487               2.4364150
#> 38              0.000000               0.0000000
#> 39              3.178674               3.5040347
#> 40              2.358569               3.5956924
#> 41              3.841760               3.7278628
#> 42              2.615985               2.2005972
#> 43              2.592163               3.4396624
#> 44              2.764626               3.5032509
#> 45              0.000000               0.0000000
#> 46              1.449164               0.7148429
#> 47              0.871334               2.1930064
#> 48              3.003904               2.3719509
#> 49              3.022750               2.1161403
#> 50              3.373053               3.0657625
#> 51              0.000000               0.0000000
#> 52              2.870643               2.0939711
#> 53              3.669835               3.8411070
#> 54              2.879408               2.0875444
#> 55              3.056836               2.3300218
#> 56              2.899898               2.8568829
#> 57              2.381507               1.4866785
#> 58              2.802044               0.3706927
#> 59              1.527473               2.2513089
#> 60              0.931326               2.4784371
#> 61              3.338036               2.7433537
#> 62              4.701908               3.6263545
#> 63              2.748970               3.0976662
#> 64              2.860497               2.6906423
#> 65              2.402197               2.7811464
#> 66              2.653452               3.0350966
#> 67              2.952819               3.3247377
#> 68              3.091557               3.6409482
#> 
#> $predictor.correlations
#>            complexity      depth       ZONE
#> complexity 1.00000000 0.33364348 0.01410832
#> depth      0.33364348 1.00000000 0.08234682
#> ZONE       0.01410832 0.08234682 1.00000000
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
