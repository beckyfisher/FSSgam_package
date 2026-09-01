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

# Shared internal utilities. Not exported.

# Validates that every column of dat is factor/character/integer/numeric,
# then splits its column names into factor and continuous predictors. Used
# identically by check_correlations() and check_non_linear_correlations().
classify_correlation_predictors=function(dat){
  classes.dat=sapply(dat,class)
  classes.dat=lapply(classes.dat,FUN=paste,collapse=" ")
  valid.cols=which(match(unlist(classes.dat),c("factor","character", "integer","numeric"))>0)
  if(length(valid.cols)<ncol(dat)){
     invalid.cols=colnames(dat[-valid.cols])
     invalid.classes=classes.dat[invalid.cols]
     stop(
        paste("
        The predictor",
        invalid.cols,"is of class",invalid.classes,
        "which is not supported. Please check your input data.frame.")
     )}
  fact.vars=names(which(classes.dat=="factor" | classes.dat=="character"))
  cont.vars=names(which(classes.dat=="integer" | classes.dat=="numeric"))
  return(list(fact.vars=fact.vars,cont.vars=cont.vars))
}

# Confirms progress is a single, non-missing TRUE/FALSE. Without this the
# value reaches `if(progress)` inside the fitting helpers, where NA raises
# "missing value where TRUE/FALSE needed" and a string raises "argument is not
# interpretable as logical" -- neither of which names the argument at fault.
validate_progress=function(progress){
  if(!is.logical(progress) || length(progress)!=1 || is.na(progress)){
    stop("progress must be a single TRUE or FALSE.")
  }
}

# Resolves a fresh, independent family object to refit a candidate model
# from test.fit, for use by fit_mod_l() (functions_supporting.R) and
# fit_model_set()'s internal loops (fit-model-set.R). Re-evaluates the
# unfitted family expression captured in test.fit's own call -- in the
# environment where test.fit was originally created, not wherever the
# refit itself happens to run -- then hands the result to
# clone_independent_family() so repeated candidate refits never share
# mutable extended-family state. This combined approach is required to
# fix GitHub issues beckyfisher/FSSgam#10 and #12 together: #10 because a
# doSNOW worker process never has the calling session's variables in scope (so the
# expression must be evaluated on the calling process, not re-evaluated
# lazily inside update() on whatever process refits the candidate); #12
# because mgcv's extended families (nb(), tw(), ...) store their estimated
# extra parameter in mutable state that, if shared, gets warm-started from
# whatever any other refit (or test.fit itself) last left it at. Falls
# back to NULL (update()'s own no-override behaviour) if family can't be
# resolved this way for some test.fit type we haven't anticipated -- never
# worse than the pre-existing behaviour, just without this fix.
# Returns NULL if test.fit's call has no explicit family= (i.e. the
# default gaussian()).
resolve_candidate_family=function(test.fit){
  fam.expr=test.fit$call$family
  if(is.null(fam.expr)) return(NULL)
  fam=try(eval(fam.expr,envir=environment(stats::formula(test.fit))),silent=TRUE)
  if(inherits(fam,"try-error")) return(NULL)
  clone_independent_family(fam)
}

# Returns a copy of fam (a glm/gam family object) in which every mutable
# environment private to this particular family instance (e.g. the
# environment behind mgcv's extended-family getTheta()/putTheta() pair,
# used by nb(), tw(), ...) is duplicated into its own independent copy.
# Shared, non-instance environments (base/stats namespaces etc, identified
# via environmentName() being non-empty, since per-instance closures
# created by family generators like nb()/tw() are always anonymous) are
# left untouched, so functions still resolve any package-internal helpers
# they call via the original (preserved) parent environment chain. Plain
# (non-extended) families have no such mutable state and are returned
# as-is.
clone_independent_family=function(fam){
  if(!inherits(fam,"extended.family")) return(fam)
  is.fn=vapply(fam,is.function,logical(1))
  envs=lapply(fam[is.fn],environment)
  local.idx=vapply(envs,function(e) identical(environmentName(e),""),logical(1))
  local.envs=envs[local.idx]
  unique.envs=unique(local.envs)
  cloned.envs=lapply(unique.envs,function(e){
    list2env(as.list(e,all.names=TRUE),parent=parent.env(e))
  })
  # match() does not reliably compare environments by identity -- do it by hand.
  match.idx=vapply(local.envs,function(e){
    which(vapply(unique.envs,identical,logical(1),e))[1]
  },integer(1))
  fam[is.fn][local.idx]=Map(function(f,i){
    environment(f)=cloned.envs[[i]]
    f
  },fam[is.fn][local.idx],match.idx)
  fam
}
