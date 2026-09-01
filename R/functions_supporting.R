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

#' wi
#'
#' Supporting function for functions full_subsets_gam and fit_model_set. Not called directly.
#'
#' @param  AIC.vals vector of AICc, AIC or BIC values
#'
#' @details Calculates Akaike weight values from a vector of AICc, AIC or BIC values
#'
#' @export
#' @return A vector of Akaike weights
#' @examples
#' wi(c(100, 102, 105, 110))
wi <- function(AIC.vals){# This function calculate the Aikaike weights:
 # wi=(exp(-1/2*AICc.vals.adj))/Sum.wi=1 to r (exp(-1/2*AICc.vals.adj))
 AICc.vals.adj=AIC.vals-min(stats::na.omit(AIC.vals))
 wi.den=rep(NA,length(AICc.vals.adj))
 for(i in 1:length(AICc.vals.adj)){
  wi.den[i]=exp(-1/2*AICc.vals.adj[i])}
 wi.den.sum=sum(stats::na.omit(wi.den))
 wi=wi.den/wi.den.sum
 return(wi)}

#' extract_mod_dat
#'
#' Supporting function for functions full_subsets_gam and fit_model_set. Not called directly.
#'
#' @param  mod.fit A dsm, gam or uGamm fitted model object
#'
#' @param  r2.type. The type of r2 to extract. Passed through arguments supplied to fit_model_set
#'
#' @details Extracts model fit parameters from a dsm, gam or uGamm fitted model object
#'
#' @export
#' @return A list of model fit parameters
#' @examples
#' library(mgcv)
#' library(MuMIn)
#' data(case_study1)
#' fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
#'            family = tw(), data = case_study1)
#' extract_mod_dat(fit, r2.type. = "r2")
extract_mod_dat <- function(mod.fit,r2.type.="r2.lm.est"){
#x=mod.fit
 mod.dat <- list(AICc=NA,BIC=NA,r2.vals=NA,r2.vals.unique=NA,edf=NA,edf.less.1=NA)
 if(!inherits(mod.fit,"try-error")){
  # AIC and BIC
  mod.dat$AICc <- MuMIn::AICc(mod.fit)
  mod.dat$BIC <- stats::BIC(mod.fit)
  #R.sq
        tempOut=NA
        if(inherits(mod.fit,"gam") & r2.type.=="dev"){tempOut=summary(mod.fit)$dev.expl}
        if(inherits(mod.fit,"gam") & r2.type.=="r2"){tempOut=summary(mod.fit)$r.sq}
        if(inherits(mod.fit,"gam") & r2.type.=="r2.lm.est"){
           tempOut=summary(stats::lm(mod.fit$y~stats::predict(mod.fit)))$r.sq}
        if(inherits(mod.fit,"gamm4") & r2.type.=="dev"){
           tempOut=summary(mod.fit$gam)$dev.expl
           if(length(tempOut)==0){tempOut=NA}}
        if(inherits(mod.fit,"gamm4") & r2.type.=="r2"){tempOut=summary(mod.fit$gam)$r.sq}
        if(inherits(mod.fit,"gamm") & r2.type.=="r2"){tempOut=summary(mod.fit$gam)$r.sq}
        if(inherits(mod.fit,"gamm4") & r2.type.=="r2.lm.est"){
          if(stats::family(mod.fit$mer)[1]=="binomial"){
            y_dat <- attributes(mod.fit$mer)$frame$y
            y <- y_dat[,1]/(y_dat[,1] + y_dat[,2])
            x <- stats::predict(mod.fit[[1]],re.form=NA,type="response")

            tempOut=summary(stats::lm(y~x))$r.sq
           }else{
            tempOut=summary(stats::lm(attributes(mod.fit$mer)$frame$y~
                        stats::predict(mod.fit[[1]],re.form=NA,type="response")))$r.sq
           }
          }

           if(is.null(tempOut)){tempOut=NA}
  mod.dat$r2.vals=round(tempOut,5)
  # Summed edf
         if(inherits(mod.fit,"gam")){
          edf.m=summary(mod.fit)$edf
          p.coeff.m=summary(mod.fit)$p.coeff}else{
           #edf.m=summary(mod.fit$gam)$edf
           #p.coeff.m=summary(mod.fit$gam)$p.coeff
           edf.m=mod.fit$gam$edf
           p.coeff.m=mod.fit$gam$p.coeff
           }
        edf.m[which(edf.m<1)]=1 # any edf<0 are reset to 1 to ensure proper
                                # parameter count when there is shrinkage (bs='cc')
  mod.dat$edf=round(sum(c(edf.m,length(p.coeff.m))),2)
  # count the edf values less than 0.25 to check for serious shrinkage
         if(inherits(mod.fit,"gam")){
           edf.m <- summary(mod.fit)$edf}else{
             edf.m <- mod.fit$gam$edf}
  mod.dat$edf.less.1 <- length(which(edf.m<0.25))}
return(mod.dat)}

#' build_inclusion_mat
#'
#' Supporting function for functions full_subsets_gam and fit_model_set. Not called directly.
#'
#' @param  included.vars A character vector of variables included in the model set
#'
#' @param  formula.list A list of model formula, as obtained through generate_model_set
#'
#' @details Builds var.inclusion matrix based on the included variables and set of model formula
#'
#' @export
#' @return A matrix of variables included in the model set
#' @examples
#' included.vars <- c("depth", "complexity")
#' formula.list <- list(depth = ~1, "depth+complexity" = ~1)
#' build_inclusion_mat(included.vars, formula.list)
build_inclusion_mat <- function(included.vars,formula.list){
var.inclusions <- matrix(0,ncol=length(included.vars),length(formula.list))
colnames(var.inclusions) <- c(included.vars)

for(m in 1:length(formula.list)){
      pred.vars.m=unique(
        unlist(strsplit(unlist(strsplit(unlist(strsplit(unlist(strsplit(unlist(strsplit(unlist(strsplit(names(formula.list)[m],
        split="+",fixed=TRUE)),
        split=".by.",fixed=TRUE)),
        split=".I.",fixed=TRUE)),
        split="*",fixed=TRUE)),
        split=".t.",fixed=TRUE)),
        split=".te.",fixed=TRUE)))
      if(pred.vars.m[1]!="null"){var.inclusions[m,pred.vars.m]=1}}
return(var.inclusions)
}


#' fit_mod_l
#'
#' Supporting function for functions full_subsets_gam and fit_model_set. Not called directly.
#'
#' @param  formula.l A model formula
#'
#' @param  test.fit. A dsm, gam or uGamm fitted model object
#'
#' @param  use.dat the data used to fit test.fit#
#'
#' @param  family. The family to refit formula.l with. Defaults to a fresh,
#' independent re-evaluation of the family test.fit. itself used (see
#' resolve_candidate_family in R/utils.R), so repeated calls never share
#' mutable extended-family state (e.g. mgcv's nb()/tw() estimated theta).
#' fit_model_set() resolves this once per candidate up front, on the
#' calling process, and passes it in explicitly rather than relying on this
#' default -- see the comment below for why.
#'
#' @details Generates an updated model fit based on the supplied formula.
#' This wrapper was required to allow full_subsets_gam and fit_model_set to be applied to dsm models
#'
#' @export
#' @return An updated dsm, gam or uGamm fitted model object
#' @examples
#' library(mgcv)
#' data(case_study1)
#' base.fit <- gam(Herbivore.abundance ~ s(depth, k = 3, bs = "cr"),
#'                  family = tw(), data = case_study1)
#' fit_mod_l(formula.l = ~ s(complexity, k = 3, bs = "cr"),
#'           test.fit. = base.fit, use.dat = case_study1)
fit_mod_l <- function(formula.l,test.fit.,use.dat,family.=resolve_candidate_family(test.fit.)){
if(length(grep("dsm",class(test.fit.)))>0){
 mod.l=try(stats::update(test.fit.,formula=formula.l),
           silent=TRUE)}
if(length(grep("dsm",class(test.fit.)))==0){
 # family. is resolved via resolve_candidate_family() (R/utils.R), not by
 # omitting family= here and letting update() re-evaluate test.fit.'s
 # original family call itself. That was tried twice and breaks one of two
 # ways depending on which fix is "on":
 #  - Passing the already-FITTED family (family=stats::family(test.fit.))
 #    shares one mutable family object (theta etc.) across every candidate
 #    refit, warm-starting each from test.fit.'s unrelated estimate and
 #    destabilising mgcv's IRLS loop for most formulas (GitHub issue
 #    beckyfisher/FSSgam#12).
 #  - Omitting family= entirely so update() re-evaluates test.fit.'s
 #    original *expression* fixes that when the expression is a literal
 #    constructor call (family = tw()), but when family was supplied as a
 #    vector/list element (family = family.vec[[2]], GitHub issue
 #    beckyfisher/FSSgam#10) two
 #    things break: (a) under parallel = TRUE, update() re-evaluates that
 #    expression on a doSNOW worker process that never had family.vec
 #    exported to it, so every refit fails with "object not found"; and
 #    (b) even sequentially, re-evaluating an indexing expression doesn't
 #    construct a new family object at all -- it returns the same shared
 #    reference every time, silently reintroducing the #12 sharing problem.
 # resolve_candidate_family() avoids both: it evaluates the original family
 # expression once, in the environment test.fit. was actually created in
 # (not wherever the refit itself happens to run), and clones any mutable
 # state so every refit gets its own independent copy regardless of how
 # family was originally specified.
 if(is.null(family.)){
  mod.l=try(stats::update(test.fit.,formula=formula.l,data=use.dat),
            silent=TRUE)
 }else{
  mod.l=try(stats::update(test.fit.,formula=formula.l,data=use.dat,family=family.),
            silent=TRUE)
 }}
return(mod.l)
}





