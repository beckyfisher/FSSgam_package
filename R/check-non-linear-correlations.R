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

#' check_non_linear_correlations
#'
#' generates a correlation matrix among all columns of a data.frame
#' @param dat the data.frame containing the columns for which a correlation
#' matrix is sought.
#'
#' @details The function uses gam to estimate a correlation coefficient
#' among continuous variables (continuous~s(continuous), lm to approximate the correlation coefficient
#' between a continuous variable (as response) and a factor variable (as a predictor)
#' through the call lm(continuous~factor),
#' and nnet to apporoximate the correlation for factor variables as responses using a multnomial
#' model fit through a call to multinom(factor~factor) or (factor~continuous).
#' @export
#' @return an approximate correlation matrix
#' @note The resulting "correlation" matrix is assymetric as the row variable
#' is used as the "response" and the column variable is used as the "predictor".
#' The use of gam may be slightly oversensitive for continuous-continuous correlations
#' and users may wish to increase cor.cutoff. Inspect individual replationships manually.
#' Values are only approximate "correlations" and are in fact the sqrt of the
#' R-square values reported for each of the fitted relationships. Note that the function
#' assumes a gaussian distribution for continuous response variables. Substantial
#' deviations from this assumption will yield spurious "correlation" estimates.
#' @examples
#' data(case_study1)
#' check_non_linear_correlations(case_study1[, c("depth", "complexity", "ZONE")])
check_non_linear_correlations=function(dat){
  vars=classify_correlation_predictors(dat)
  fact.vars=vars$fact.vars
  cont.vars=vars$cont.vars

  test.mat=build_correlation_pair_grid(dat=dat)

  test.mat$r.sq=apply(test.mat,MARGIN=1,FUN=function(x){
    estimate_non_linear_correlation(response.var1=as.character(unlist(x["rows"])),
                          predictor.var2=as.character(unlist(x["cols"])),
                          dat=dat,fact.vars=fact.vars,cont.vars=cont.vars)})

  out.cor.mat=assemble_non_linear_correlation_matrix(dat=dat,test.mat=test.mat)

  return(out.cor.mat)
}

# Internal helpers for check_non_linear_correlations(). Not exported.

# All ordered (row, col) column-name pairs of dat, excluding the diagonal.
build_correlation_pair_grid=function(dat){
  test.mat=expand.grid(colnames(dat),colnames(dat))
  colnames(test.mat)=c("rows","cols")
  test.mat=test.mat[which(test.mat$rows!=test.mat$cols),]
  rownames(test.mat)=1:nrow(test.mat)
  return(test.mat)
}

# Estimates the asymmetric "correlation" of response.var1 (as response) on
# predictor.var2 (as predictor): a multinomial model if the response is a
# factor, lm() if the response is continuous and the predictor is a factor,
# or gam() if both are continuous.
estimate_non_linear_correlation=function(response.var1,predictor.var2,dat,fact.vars,cont.vars){
    r.est=NA
    class.response.var1=character()
    class.predictor.var2=character()
    if(length(stats::na.omit(match(response.var1,fact.vars)))==1){class.response.var1="factor"}
    if(length(stats::na.omit(match(response.var1,cont.vars)))==1){class.response.var1="continuous"}
    if(length(stats::na.omit(match(predictor.var2,fact.vars)))==1){class.predictor.var2="factor"}
    if(length(stats::na.omit(match(predictor.var2,cont.vars)))==1){class.predictor.var2="continuous"}

    dat.r=stats::na.omit(dat[,c(response.var1,predictor.var2)])
    dat.r$response.var1=dat.r[,response.var1]
    dat.r$predictor.var2=dat.r[,predictor.var2]

    # if the response.var1 variable is a factor use a multinomial model
    if(class.response.var1=="factor"){
      fit <- try(summary(nnet::multinom(response.var1 ~ predictor.var2,trace=FALSE,
          data=dat.r))$deviance,silent=TRUE)
      null.fit=try(summary(nnet::multinom(response.var1 ~ 1,trace=FALSE,data=dat.r))$deviance,silent=TRUE)
      if(!inherits(fit,"try-error")){
         if(round(fit,4)==round(null.fit,4)){r.est=0}else{
        r.est=sqrt(1-(fit/null.fit))}
        }
     }
    # if the response.var1 variable is continuous and the predictor.var2 is a factor, do an
    # anova
    if(class.response.var1=="continuous" & class.predictor.var2 == "factor"){
       r.est=sqrt(summary(stats::lm(response.var1~predictor.var2,data=dat.r))$r.sq)
     }
    # if both the response.var1 variable and the predictor.var2 are continuous, do a gam
    if(class.response.var1=="continuous" & class.predictor.var2 == "continuous"){
       nn=length(unique(dat.r$predictor.var2))
       k.use=4
       if(nn<k.use){k.use=nn}
       # s() must stay unqualified here: mgcv's formula parser matches smooth
       # terms by literal symbol name, so `mgcv::s(...)` is not recognised as a
       # smooth term and silently breaks model.frame() construction.
       fit=try(mgcv::gam(response.var1~s(predictor.var2,k=k.use),data=dat.r),silent=TRUE)
       if(!inherits(fit,"try-error")){
        r.sq=summary(fit)$r.sq
       if(r.sq>0){r.est=sqrt(r.sq)}else{r.est=0}}
     }
    return(r.est)
}

# Assembles the (rows, cols, r.sq) long-format test.mat into a square
# matrix with a unit diagonal.
assemble_non_linear_correlation_matrix=function(dat,test.mat){
  out.cor.mat=matrix(NA,nrow=ncol(dat),ncol=ncol(dat))
  rownames(out.cor.mat)=colnames(dat)
  colnames(out.cor.mat)=colnames(dat)
  diag(out.cor.mat)=1

  for(r in 1:nrow(test.mat)){
     out.cor.mat[test.mat$rows[r],test.mat$cols[r]]=test.mat$r.sq[r]
  }

  return(out.cor.mat)
}
