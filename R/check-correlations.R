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

#' check_correlations
#'
#' generates a correlation matrix among all columns of a data.frame
#' @param dat the data.frame containing the columns for which a correlation
#' matrix is sought.
#' @param parallel a logical indicating if calaculation of the correlation matrix
#' should be done in parallel. Defaults to FALSE.
#' @param n.cores a numeric value indicating the number of cores to utilise if
#' parallel is TRUE.
#'
#' @details The function uses cor to calculate the Pearson correlation coefficient
#' among continuous variables, lm to approximate the correlation coefficient
#' among a continuous variable and a factor variable through the call lm(continuous~factor),
#' and nnet to approximate the correlation among factor variables using a multinomial
#' model fit.
#'
#' Missing values are handled pairwise: each pair of predictors is evaluated on
#' the rows for which both are present. For a factor-factor pair this applies to
#' the intercept-only model as well as the fitted one, so the two deviances the
#' estimate is a ratio of are always computed on the same rows.
#' @export
#' @return a correlation matrix
#' @examples
#' data(case_study1)
#' check_correlations(case_study1[, c("depth", "complexity", "ZONE")])
check_correlations=function(dat,parallel=FALSE,n.cores=4){
  vars=classify_correlation_predictors(dat)
  fact.vars=vars$fact.vars
  cont.vars=vars$cont.vars

  cor.mat=build_continuous_correlation_matrix(dat=dat,cont.vars=cont.vars)

  if(length(fact.vars)>0){
    out.cor.mat=build_factor_continuous_skeleton(dat=dat,fact.vars=fact.vars,
                          cont.vars=cont.vars,cor.mat=cor.mat)
    out.cor.mat=fill_factor_factor_correlations(dat=dat,fact.vars=fact.vars,
                          out.cor.mat=out.cor.mat,parallel=parallel,n.cores=n.cores)
  }else{
    out.cor.mat=cor.mat
  }
  return(out.cor.mat)
}

# Internal helpers for check_correlations(). Not exported.

# Pearson correlation matrix among the continuous predictors only.
build_continuous_correlation_matrix=function(dat,cont.vars){
  if(length(cont.vars)>1){
   cor.mat=stats::cor(dat[,cont.vars],use="pairwise.complete.obs")}else{
   cor.mat=matrix(1,ncol=1,nrow=1)
   colnames(cor.mat)=cont.vars
   rownames(cor.mat)=cont.vars}
  return(cor.mat)
}

# Builds the correlation matrix skeleton: continuous-continuous values from
# cor.mat, factor-continuous values estimated via lm(continuous~factor), and
# NA placeholders for factor-factor (filled in later by
# fill_factor_factor_correlations()).
build_factor_continuous_skeleton=function(dat,fact.vars,cont.vars,cor.mat){
   if(length(cont.vars)>0){
    lm.grid=expand.grid(list(fact.var=fact.vars,cont.var=cont.vars))
    r.estimates=cbind(lm.grid,apply(lm.grid,MARGIN=1,FUN=function(x){
        sqrt(summary(stats::lm(dat[,x[2]]~factor(dat[,x[1]])))$r.sq)}))

    fact.cont.upper.right=matrix(NA,ncol=length(fact.vars),nrow=length(cont.vars))
    colnames(fact.cont.upper.right)=fact.vars;rownames(fact.cont.upper.right)=cont.vars

    fact.cont.lower.left=matrix(NA,ncol=length(cont.vars),nrow=length(fact.vars))
    colnames(fact.cont.lower.left)=cont.vars;rownames(fact.cont.lower.left)=fact.vars

    fact.fact.lower.right=matrix(NA,ncol=length(fact.vars),nrow=length(fact.vars))
    colnames(fact.fact.lower.right)=fact.vars;rownames(fact.fact.lower.right)=fact.vars

    out.cor.mat=rbind(cbind(cor.mat,fact.cont.upper.right),
                      cbind(fact.cont.lower.left,fact.fact.lower.right))

    # assign the estimated r values to the upper right and lower left corners
    for(r in seq_len(nrow(r.estimates))){
       # upper right
       col.index=which(colnames(out.cor.mat)==r.estimates$fact.var[r])
       row.index=which(rownames(out.cor.mat)==r.estimates$cont.var[r])
       out.cor.mat[row.index,col.index]=r.estimates[r,3]
       # lower left
       col.index=which(colnames(out.cor.mat)==r.estimates$cont.var[r])
       row.index=which(rownames(out.cor.mat)==r.estimates$fact.var[r])
       out.cor.mat[row.index,col.index]=r.estimates[r,3]
     }
    }else{
        fact.fact.lower.right=matrix(NA,ncol=length(fact.vars),nrow=length(fact.vars))
    colnames(fact.fact.lower.right)=fact.vars;rownames(fact.fact.lower.right)=fact.vars
    out.cor.mat=fact.fact.lower.right}
  return(out.cor.mat)
}

# Estimates factor-factor correlations via nnet::multinom (parallel or
# sequential) and fills them into out.cor.mat.
fill_factor_factor_correlations=function(dat,fact.vars,out.cor.mat,parallel,n.cores){
  # estimate r values for fact-fact combinations
  #
  # Self-pairs are dropped here and the diagonal set to 1 at the end. Fitting
  # var ~ var gave a pseudo-R2 of about 0.999999 rather than exactly 1
  # (FSSgam_package#12), at the cost of one multinom() fit per factor for a
  # value that is known in advance.
  lm.grid=expand.grid(list(fact.var1=fact.vars,fact.var2=fact.vars))
  lm.grid=lm.grid[as.character(lm.grid$fact.var1)!=as.character(lm.grid$fact.var2),,drop=FALSE]

  # r.est is sqrt(1 - fit/null.fit), so both deviances have to be computed on
  # the same rows. `fit` below is fitted on na.omit() of the pair, so the null
  # has to be too (FSSgam_package#16). Fitting the null per pair
  # unconditionally would undo the reduction from n(n-1) null fits to n, so the
  # per-factor value is kept and refitted only for a pair that drops rows the
  # factor's own column does not.
  #
  # multinom() applies na.action itself, so the per-factor null is already
  # fitted on that factor's own complete cases, and complete.counts records how
  # many rows that is. The pair's rows are a subset of them, so equal counts
  # mean an identical row set and the per-factor value is exact.
  null.deviances=lapply(fact.vars,function(v){
    try(nnet::multinom(dat[,v] ~ 1,trace=FALSE)$deviance,silent=TRUE)})
  names(null.deviances)=fact.vars
  complete.counts=lapply(fact.vars,function(v){sum(!is.na(dat[,v]))})
  names(complete.counts)=fact.vars
  if(parallel==TRUE){
   cl=parallel::makeCluster(n.cores)
   doSNOW::registerDoSNOW(cl)
   out.cor.dat<-foreach::foreach(r = seq_len(nrow(lm.grid)),.packages=c('nnet'),.errorhandling='pass')%dopar%{
    var.1=as.character(lm.grid[r,1])
    var.2=as.character(lm.grid[r,2])
    dat.r=stats::na.omit(dat[,c(var.1,var.2)])
    # $deviance is taken from the fit itself rather than from summary(): the
    # summary method computes standard errors that are discarded here, and
    # sqrt()s a negative variance on a perfectly separated fit, which is where
    # the spurious "NaNs produced" warnings came from (FSSgam_package#10).
    fit <- try(nnet::multinom(dat.r[,var.1] ~ dat.r[,var.2],trace=FALSE)$deviance,silent=TRUE)
    null.fit=null.deviances[[var.1]]
    if(nrow(dat.r)<complete.counts[[var.1]]){
      null.fit=try(nnet::multinom(dat.r[,var.1] ~ 1,trace=FALSE)$deviance,silent=TRUE)}
    if(!inherits(fit,"try-error")&!inherits(null.fit,"try-error")){
       if(round(fit,4)==round(null.fit,4)){r.est=0}else{
      r.est=sqrt(1-(fit/null.fit))}
      c(var.1,var.2,r.est)}}
   parallel::stopCluster(cl)
   foreach::registerDoSEQ()
   }else{
    out.cor.dat=list()
    for(r in seq_len(nrow(lm.grid))){
          var.1=as.character(lm.grid[r,1])
          var.2=as.character(lm.grid[r,2])
          dat.r=stats::na.omit(dat[,c(var.1,var.2)])
          # see the matching comment in the parallel branch above
          fit <- try(nnet::multinom(dat.r[,var.1] ~ dat.r[,var.2],trace=FALSE)$deviance,silent=TRUE)
          null.fit=null.deviances[[var.1]]
          # see the matching comment on complete.counts above
          if(nrow(dat.r)<complete.counts[[var.1]]){
            null.fit=try(nnet::multinom(dat.r[,var.1] ~ 1,trace=FALSE)$deviance,silent=TRUE)}
          out=NA
          if(!inherits(fit,"try-error")&!inherits(null.fit,"try-error")){
           if(round(fit,4)==round(null.fit,4)){r.est=0}else{
                   r.est=sqrt(1-(fit/null.fit))}
                   out=c(var.1,var.2,r.est)}
      out.cor.dat=c(out.cor.dat,list(out))}
      }

    for(r in seq_along(out.cor.dat)){
       out.cor.mat[which(colnames(out.cor.mat)==out.cor.dat[[r]][1]),
                   which(rownames(out.cor.mat)==out.cor.dat[[r]][2])]=
                   as.numeric(out.cor.dat[[r]][3])}

  # Only the factor block's diagonal: the continuous block comes from cor(),
  # which already returns exactly 1 there.
  out.cor.mat[cbind(fact.vars,fact.vars)]=1
  return(out.cor.mat)
}
