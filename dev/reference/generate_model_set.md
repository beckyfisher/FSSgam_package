# generate_model_set

Generate a complete full subsets model set for analysis based on gam(m4)
via a call to fit_model_set.

## Usage

``` r
generate_model_set(
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
  null.cov.cutoff = 0.8
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
  gamm (mgcv), reached through MuMIn::uGamm since a bare mgcv::gamm fit
  cannot be refitted, is based on nlme which allows correlation
  structures \[Box, G.E.P., Jenkins, G.M., and Reinsel G.C. (1994) "Time
  Series Analysis: Forecasting and Control", 3rd Edition, Holden-Day\],
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
  subsequently removed (see cov.cutoff below). A user defined cor.matrix
  (see cor.matrix) need not include these hard coded interactions: rows
  and columns for any it does not name are computed from use.dat.

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
  and a column for every predictor named in pred.vars.cont,
  pred.vars.fact and linear.vars; any that are missing are reported by
  name. The hard coded factor interaction columns that setting
  factor.factor.interactions causes to be created are the exception,
  because which of them exist depends on the supplied matrix itself and
  so cannot be known in advance. Rows and columns for any of those the
  matrix does not carry are computed from use.dat and appended, leaving
  every supplied value as supplied. Each dimension is treated
  separately, so a name supplied as a column and not as a row keeps the
  column given for it and has only its row computed; because
  collinearity is screened in both directions, a value supplied in one
  dimension alone can tighten a screen but never loosen it. Computing
  them reads the data of every predictor, so this is the one case in
  which a predictor of a class check_correlations cannot classify has to
  be named in the supplied matrix along with the interaction columns.
  When supplied it replaces the automatic estimate rather than
  overriding it: except in the interaction case above,
  check_correlations is not called at all, and a predictor of a class it
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
  models, subject to two rules. It must be two-dimensional: a matrix, a
  data.frame, or any other object with two dimensions, so the value
  Matrix::Matrix returns is accepted while a length-one value of any
  other class is not, and NA is reserved for the default. And it must
  not contain NA between two terms that are actually screened against
  cov.cutoff, which is reported by naming the pairs; an NA on a pair no
  screen compares is accepted, so which cells matter depends on
  max.predictors and on the interaction arguments.

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
  gam instead of uGamm to make sure PQL is not used. To fit a
  correlation structure, pass MuMIn::uGamm(..., correlation = ...) as
  the test.fit. A fit produced by mgcv::gamm() directly cannot be used:
  it records no call, so the candidate models cannot be refitted from
  it, and generate_model_set stops saying so.

- null.cov.cutoff:

  The correlation above which a predictor is dropped for being
  correlated with a variable named in null.terms. Defaults to 0.8. Terms
  supplied through null.terms are forced into every candidate model and
  are outside the cov.cutoff screen, which covers pred.vars.cont,
  pred.vars.fact and linear.vars only, so without this a candidate could
  be arbitrarily strongly correlated with a forced term and still appear
  in every model in the set. That inflates the variance of the forced
  term's estimate, which is frequently the term the analysis exists to
  estimate, and nothing in the output would indicate it. A separate
  cutoff, with a much looser default than cov.cutoff, because the two
  screens answer different questions: cov.cutoff decides which
  predictors may appear together, and this decides which predictor is so
  nearly a restatement of a forced term that fitting both is not
  informative. Set it to 1 to admit every predictor whatever its
  correlation with a forced term. A variable inside a bs='re' smooth is
  exempt from this screen, however the argument is written – quoted
  either way, or wrapped in c(). A bs given as a variable rather than a
  literal is exempt too, since it cannot be read here and screening a
  grouping factor in error drops a legitimate predictor, while not
  screening a fixed term in error leaves the behaviour of earlier
  versions. A random-effect grouping factor is correlated with the
  predictors measured within it by construction – with null.terms =
  "s(Location,Site,bs='re')" and Status nested in Location their
  correlation is 1 – and that is the design of the study rather than
  collinearity to screen out. The screen applies to fixed forced terms,
  which compete with a candidate for the same variation. Correlations
  among the null.terms variables themselves are neither computed nor
  screened: those terms are forced in by your decision, and dropping one
  is what must not happen. Where the correlation estimate is asymmetric,
  as check_non_linear_correlations returns it, both directions are read
  and the larger is used, which is what the cov.cutoff screen already
  does. The correlations this screens on are returned by
  generate_model_set as null.term.correlations, each cell being the
  value screened on rather than one direction of it, whether or not
  anything is dropped, so they can be inspected even where no warning is
  raised; full_subsets_gam does not return them, its output being that
  of fit_model_set. A supplied cor.matrix is used for any forced term it
  names in either dimension, and the rest of the block is computed from
  use.dat, which is the ordinary case since a supplied matrix is indexed
  by predictor and a forced term is not one. Where that computation
  fails, which is what a predictor of a class check_correlations cannot
  classify causes, the screen is skipped with a warning rather than the
  call stopping. A variable named in null.terms that is not a column of
  use.dat, such as a function or a term written over several columns, is
  skipped, a correlation not being defined for it.

## Value

A list of the following output files:

n.mods - The number of candidate models generated, equal to
length(mod.formula).

predictor.correlations - The matrix of estimated predictor correlations
returned by the function check_correlations and used for model exclusion
based on cov.cutoff

null.term.correlations - A matrix of the correlations between each
variable named in null.terms and each candidate predictor, used for
model exclusion based on null.cov.cutoff. Rows are the null.terms
variables and columns the predictors. NULL where null.terms is empty,
where it names only random-effect grouping factors or no column of
use.dat, or where every predictor is itself a null.terms variable. The
element is always present in the returned list; it is its value that is
NULL. Correlations among the null.terms variables are not included,
being neither computed nor screened.

mod.formula - A named list containing the model formula that were
generated (and will be fitted by fit_model_set). The names are the
candidate model names used in the modname column of fit_model_set's
output.

used.data - A data.frame which is identical to the data.frame initially
supplied by the user, but with any hard coded interaction terms appended
via cbind.

test.fit - The test.fit supplied by the user, passed through so that
fit_model_set can update it.

included.vars - A character vector of the predictors included in the
model set, used to build the variable inclusion columns of
fit_model_set's output.

## Details

The function constructs a complete model set based on the supplied
arguments. for more information see Fisher R, Wilson SK, Sin TM, Lee AC,
Langlois TJ (2018) A simple function for full-subsets multiple
regression in ecology with R. Ecology and Evolution
https://onlinelibrary.wiley.com/doi/abs/10.1002/ece3.4134

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
```
