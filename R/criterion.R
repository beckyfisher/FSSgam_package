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

# The model selection criterion. Not exported.
#
# extract_mod_dat() reads AICc from MuMIn::AICc() and BIC from stats::BIC(),
# both of which resolve to the fitted family's own aic slot. Two families of
# test.fit break that route, in opposite ways:
#
#  - a quasi-likelihood has no log-likelihood, so the slot returns NA and every
#    candidate is recorded with an NA criterion (FSSgam_package#44);
#  - mgcv's censored families cnorm() and clog() return a number that is not
#    built from a censored log-likelihood, so the set is ranked on a criterion
#    that is wrong rather than missing (FSSgam_package#42).
#
# resolve_criterion() decides, once per model set and before any model is
# fitted, which of the two applies and what to do about it.

# A single finite number, and not the try-error a failed computation returns.
usable_scalar=function(x){
  !inherits(x,"try-error")&&is.numeric(x)&&length(x)==1&&is.finite(x)
}

# The family name of a fitted model, or NA_character_ where it cannot be read.
#
# stats::family() fails on a MuMIn::uGamm() fit -- it re-evaluates the recorded
# call and the first element of that call is not a function symbol -- so the
# gam component of a gamm/gamm4 fit is read directly instead. A dsm fit
# inherits from gam and is covered by the first branch.
fit_family_name=function(fit){
  fam=NULL
  if(inherits(fit,"gam")){fam=fit$family}
  if(is.null(fam)&&!is.null(fit$gam)){fam=fit$gam$family}
  if(is.null(fam)){
    fam=try(stats::family(fit),silent=TRUE)
    if(inherits(fam,"try-error")){fam=NULL}}
  if(is.null(fam)||is.null(fam$family)){return(NA_character_)}
  as.character(fam$family)[1]
}

# TRUE for mgcv's two censored continuous families, cnorm() and clog().
#
# Matched on the prefix rather than by equality: an extended family reports its
# estimated parameter as part of its name once fitted, so the value read from a
# fitted cnorm() model is "cnorm(0.524)" and not "cnorm". The trailing
# alternation keeps the test from also matching a family whose name merely
# begins with those letters.
is_censored_family=function(fam.name){
  !is.na(fam.name)&&grepl("^(cnorm|clog)($|\\()",fam.name)
}

# The censored log-likelihood of a fit from mgcv's cnorm() or clog().
#
# Computed from the fitted mean, the scale and the response's censoring
# coding, not from the family's aic slot, so it is unaffected by the three
# defects that slot has in mgcv 1.9-4 (FSSgam_package#42) and stays correct if
# they are repaired.
#
# mgcv codes a censored response as two columns. fit$y holds the first, and
# attr(fit$y, "censor") the second: equal to the first where the observation is
# uncensored, -Inf where it is left censored at the first, +Inf where it is
# right censored at it, and the upper bound where it is interval censored. A
# fit with no censored observation at all carries no censor attribute.
censored_loglik=function(fit){
  y=fit$y
  hi=attr(y,"censor")
  if(is.null(hi)){hi=y}
  mu=as.numeric(fit$fitted.values)
  # getTheta(TRUE) returns the scale itself, not its logarithm, despite theta
  # being documented as a log scale parameter -- exponentiating it inflates the
  # spread. A prior weight rescales the family's log scale by -log(wt)/2, so it
  # divides the scale by its square root.
  s=fit$family$getTheta(TRUE)/sqrt(fit$prior.weights)
  logistic=grepl("^clog($|\\()",fit$family$family)
  # log.p on the two tails rather than log() of the probability: a far
  # left-censored observation underflows pnorm() to zero and would contribute
  # -Inf to a sum that is finite.
  d.log=if(logistic){stats::dlogis(y,mu,s,log=TRUE)}else{stats::dnorm(y,mu,s,log=TRUE)}
  p.lower=if(logistic){stats::plogis(y,mu,s,log.p=TRUE)}else{stats::pnorm(y,mu,s,log.p=TRUE)}
  p.upper=if(logistic){
    stats::plogis(y,mu,s,lower.tail=FALSE,log.p=TRUE)}else{
    stats::pnorm(y,mu,s,lower.tail=FALSE,log.p=TRUE)}
  p.at=function(q) if(logistic){stats::plogis(q,mu,s)}else{stats::pnorm(q,mu,s)}

  out=d.log
  left=is.infinite(hi)&hi<0
  out[left]=p.lower[left]
  right=is.infinite(hi)&hi>0
  out[right]=p.upper[right]
  interval=is.finite(hi)&hi!=y
  if(any(interval)){
    out[interval]=log((p.at(hi)-p.at(y))[interval])}
  sum(out)
}

# AICc and BIC built from a supplied log-likelihood, at the degrees of freedom
# and sample size the default route would have used.
#
# Only the log-likelihood is replaced. mgcv's degrees of freedom convention is
# sum(edf2) + scale.estimated + family$n.theta, which is neither obvious nor
# the same as sum(edf) + scale.estimated -- on a cnorm fit measured here, 5.997
# against 4.919 -- so it is read from the fit rather than recomputed. Both
# formulas reproduce MuMIn::AICc() and stats::BIC() exactly when handed the
# log-likelihood those functions use.
#
# attr(logLik(fit), "df") is finite even for a quasi-likelihood fit, whose
# logLik() value is NA: mgcv::logLik.gam computes the degrees of freedom
# before it reads the family's aic slot. That is what makes logLik.fn a usable
# route for a quasi test.fit.
#
# Anything that cannot be computed here gives the all-NA criterion rather than
# stopping the run, which is how a candidate that fails to fit is already
# recorded: one unfittable candidate must not cost the whole model set.
# resolve_criterion() has already called logLik.fn on the test.fit, so a
# function that is broken outright is reported there, once, and only a genuinely
# model-specific failure reaches this.
criterion_from_loglik=function(mod.fit,logLik.fn){
  ll=try(logLik.fn(mod.fit),silent=TRUE)
  if(!usable_scalar(ll)){return(list(AICc=NA,BIC=NA))}
  k=try(attr(stats::logLik(mod.fit),"df"),silent=TRUE)
  n=try(stats::nobs(mod.fit),silent=TRUE)
  if(!usable_scalar(k)||!usable_scalar(n)){return(list(AICc=NA,BIC=NA))}
  list(AICc=-2*ll+2*k+2*k*(k+1)/(n-k-1),
       BIC=-2*ll+k*log(n))
}

# Decides which log-likelihood the model set is ranked on, once, before the
# fitting loop. Returns the function extract_mod_dat() is to use, or NULL for
# the default MuMIn::AICc()/stats::BIC() route.
#
# Reached from fit_model_set() rather than from generate_model_set(), because
# it is the argument logLik.fn that decides whether a test.fit whose family has
# no log-likelihood can be used at all, and generate_model_set() never sees it.
resolve_criterion=function(test.fit,logLik.fn){
  fam.name=fit_family_name(test.fit)
  named=if(is.na(fam.name)){"could not be read"}else{paste0("is ",fam.name)}

  if(!is.null(logLik.fn)){
    if(!is.function(logLik.fn)){
      stop("logLik.fn must be a function of one argument, a fitted model, returning a single log-likelihood value, or NULL to use MuMIn::AICc() and stats::BIC().")}
    # Called on test.fit here so that a function that cannot be applied to
    # these fits is reported once, naming the argument, rather than once per
    # candidate as an NA criterion column.
    ll=try(logLik.fn(test.fit),silent=TRUE)
    if(inherits(ll,"try-error")){
      stop(paste0("logLik.fn failed on the test.fit. It must take one argument, ",
           "a fitted model, and return a single log-likelihood value. The error ",
           "was: ",as.character(ll)))}
    if(!usable_scalar(ll)){
      stop("logLik.fn must return a single finite numeric log-likelihood value; on the test.fit it did not.")}
    return(logLik.fn)
  }

  if(is_censored_family(fam.name)){
    message(paste0("The fitted family ",fam.name," is one of mgcv's censored ",
            "families, whose reported AIC is not built from a censored ",
            "log-likelihood. AICc and BIC are computed from a censored ",
            "log-likelihood instead, so the ranking, the weights and variable ",
            "importance differ from those of FSSgam 1.1.0 and earlier ",
            "(FSSgam_package#42). Supply logLik.fn to use a log-likelihood of ",
            "your own."))
    return(censored_loglik)
  }

  # The general check that the criterion is usable at all. Read from the
  # test.fit, which is already fitted, so nothing is fitted to discover it. A
  # quasi-likelihood is what reaches it in practice, but the test is on the
  # value rather than on a list of family names, so any other family whose aic
  # slot is undefined is covered too.
  crit=suppressWarnings(try(MuMIn::AICc(test.fit),silent=TRUE))
  if(!usable_scalar(crit)){
    stop(paste0("The model set cannot be ranked: MuMIn::AICc() returns no usable ",
         "value for the test.fit, whose fitted family ",named,". A ",
         "quasi-likelihood has no log-likelihood and so no AIC or BIC, and ",
         "every candidate would be recorded with an NA criterion, leaving ",
         "delta.AICc, the Akaike weights and variable importance undefined ",
         "(FSSgam_package#44). Either refit the test.fit with a family that ",
         "has a likelihood -- nb() or tw() for overdispersed counts, betar() ",
         "for proportions -- or supply logLik.fn to fit_model_set(), a ",
         "function of one fitted model returning the log-likelihood the set ",
         "is to be ranked on."))}

  NULL
}
