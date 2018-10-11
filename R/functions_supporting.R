


#' wi
#'
#' Supporting function for functions full.subsets.gam and fit.model.set. Not called directly.
#'
#' @param  AIC.vals vector of AICc, AIC or BIC values
#'
#' @details Calculates Akaike weight values from a vector of AICc, AIC or BIC values
#'
#' @export
#' @return A vector of Akaike weights
#'
wi <- function(AIC.vals){# This function calculate the Aikaike weights:
 # wi=(exp(-1/2*AICc.vals.adj))/Sum.wi=1 to r (exp(-1/2*AICc.vals.adj))
 AICc.vals.adj=AIC.vals-min(na.omit(AIC.vals))
 wi.den=rep(NA,length(AICc.vals.adj))
 for(i in 1:length(AICc.vals.adj)){
  wi.den[i]=exp(-1/2*AICc.vals.adj[i])}
 wi.den.sum=sum(na.omit(wi.den))
 wi=wi.den/wi.den.sum
 return(wi)}

#' extract.mod.dat
#'
#' Supporting function for functions full.subsets.gam and fit.model.set. Not called directly.
#'
#' @param  mod.fit A dsm, gam or uGamm fitted model object
#'
#' @param  r2.type The type of r2 to extract. Passed through arguments supplied to fit.model.set
#'
#' @details Extracts model fit parameters from a dsm, gam or uGamm fitted model object
#'
#' @export
#' @return A list of model fit parameters
#'
extract.mod.dat <- function(mod.fit,r2.type.=r2.type){
#x=mod.fit
 mod.dat=list(AICc=NA,BIC=NA,r2.vals=NA,r2.vals.unique=NA,edf=NA,edf.less.1=NA)
 if(class(mod.fit)[[1]]!="try-error"){
  # AIC and BIC
  mod.dat$AICc=AICc(mod.fit)
  mod.dat$BIC=BIC(mod.fit)
  #R.sq
        tempOut=NA
        if(class(mod.fit)[1]=="gam" & r2.type.=="dev"){tempOut=summary(mod.fit)$dev.expl}
        if(class(mod.fit)[1]=="gam" & r2.type.=="r2"){tempOut=summary(mod.fit)$r.sq}
        if(class(mod.fit)[1]=="gam" & r2.type.=="r2.lm.est"){
           tempOut=summary(lm(mod.fit$y~predict(mod.fit)))$r.sq}
        if(class(mod.fit)[[1]]=="gamm4" & r2.type.=="dev"){
           tempOut=summary(mod.fit$gam)$dev.expl
           if(length(tempOut)==0){tempOut=NA}}
        if(class(mod.fit)[[1]]=="gamm4" & r2.type.=="r2"){tempOut=summary(mod.fit$gam)$r.sq}
        if(class(mod.fit)[[1]]=="gamm4" & r2.type.=="r2.lm.est"){
           tempOut=summary(lm(attributes(mod.fit$mer)$frame$y~
                        predict(mod.fit[[1]],re.form=NA,type="response")))$r.sq}
           if(is.null(tempOut)){tempOut=NA}
  mod.dat$r2.vals=round(tempOut)
  # Summed edf
         if(class(mod.fit)[1]=="gam"){
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
         if(class(mod.fit)[1]=="gam"){edf.m=summary(mod.fit)$edf}else{edf.m=mod.fit$gam$edf}
  mod.dat$edf.less.1=length(which(edf.m<0.25))}
return(mod.dat)}

#' build.inclusion.mat
#'
#' Supporting function for functions full.subsets.gam and fit.model.set. Not called directly.
#'
#' @param  included.vars A character vector of variables included in the model set
#'
#' @param  formula.list A list of model formula, as obtained through generate.model.set
#'
#' @details Builds var.inclusion matrix based on the included variables and set of model formula
#'
#' @export
#' @return A matrix of variables included in the model set
#'
build.inclusion.mat <- function(included.vars,formula.list){
var.inclusions=matrix(0,ncol=length(included.vars),length(formula.list))
colnames(var.inclusions)=c(included.vars)

for(m in 1:length(formula.list)){
      pred.vars.m=unique(
        unlist(strsplit(unlist(strsplit(unlist(strsplit(unlist(strsplit(unlist(strsplit(unlist(strsplit(names(formula.list)[m],
        split="+",fixed=T)),
        split=".by.",fixed=T)),
        split=".I.",fixed=T)),
        split="*",fixed=T)),
        split=".t.",fixed=T)),
        split=".te.",fixed=T)))
      if(pred.vars.m[1]!="null"){var.inclusions[m,pred.vars.m]=1}}
return(var.inclusions)
}


#' fit.mod.l
#'
#' Supporting function for functions full.subsets.gam and fit.model.set. Not called directly.
#'
#' @param  formula.l A model formula
#'
#' @param  test.fit A dsm, gam or uGamm fitted model object
#'
#' @param  use.dat the data used to fit test.fit#
#'
#' @details Generates an updated model fit based on the supplied formula.
#' This wrapper was required to allow full.subsets.gam and fit.model.set to be applied to dsm models
#'
#' @export
#' @return An updated dsm, gam or uGamm fitted model object
#'
fit.mod.l <- function(formula.l,test.fit.=test.fit,use.dat.=use.dat){
if(length(grep("dsm",class(test.fit.)))>0){
 mod.l=try(update(test.fit.,formula=formula.l),
           silent=T)}
if(length(grep("dsm",class(test.fit.)))==0){
 mod.l=try(update(test.fit.,formula=formula.l,data=use.dat.),
           silent=T)}
return(mod.l)
}







