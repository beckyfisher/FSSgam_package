#' check.non.linear.correlations
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

check.non.linear.correlations=function(dat){
  classes.dat=sapply(dat,class)
  fact.vars=names(which(classes.dat=="factor" | classes.dat=="character"))
  cont.vars=names(which(classes.dat=="integer" | classes.dat=="numeric"))
  test.mat=expand.grid(colnames(dat),colnames(dat))
  colnames(test.mat)=c("rows","cols")
  test.mat=test.mat[which(test.mat$rows!=test.mat$cols),]
  rownames(test.mat)=1:nrow(test.mat)

  require(mgcv)
  if(length(fact.vars)>0){require(nnet)}

  #r=r+1
  #x=test.mat[r,]
  find.r.est.2=function(x){
    r.est=NA
    response.var1=as.character(unlist(x["rows"]))
    predictor.var2=as.character(unlist(x["cols"]))
    class.response.var1=character()
    class.predictor.var2=character()
    if(length(na.omit(match(response.var1,fact.vars)))==1){class.response.var1="factor"}
    if(length(na.omit(match(response.var1,cont.vars)))==1){class.response.var1="continuous"}
    if(length(na.omit(match(predictor.var2,fact.vars)))==1){class.predictor.var2="factor"}
    if(length(na.omit(match(predictor.var2,cont.vars)))==1){class.predictor.var2="continuous"}

    #return(c(class.response.var1,class.predictor.var2))}

    dat.r=na.omit(dat[,c(response.var1,predictor.var2)])
    dat.r$response.var1=dat.r[,response.var1]
    dat.r$predictor.var2=dat.r[,predictor.var2]

    # if the response.var1 variable is a factor use a multinomial model
    if(class.response.var1=="factor"){
      fit <- try(summary(multinom(response.var1 ~ predictor.var2,trace=F,
          data=dat.r))$deviance,silent=T)
      null.fit=try(summary(multinom(response.var1 ~ 1,trace=F,data=dat.r))$deviance,silent=T)
      if(class(fit)!="try-error"){
         if(round(fit,4)==round(null.fit,4)){r.est=0}else{
        r.est=sqrt(1-(fit/null.fit))}
        }
     }
    # if the response.var1 variable is continuous and the predictor.var2 is a factor, do an
    # anova
    if(class.response.var1=="continuous" & class.predictor.var2 == "factor"){
       r.est=sqrt(summary(lm(response.var1~predictor.var2,data=dat.r))$r.sq)
     }
    # if both the response.var1 variable and the predictor.var2 are continuous, do a gam
    if(class.response.var1=="continuous" & class.predictor.var2 == "continuous"){
       nn=length(unique(dat.r$predictor.var2))
       k.use=4
       if(nn<k.use){k.use=nn}
       fit=try(gam(response.var1~s(predictor.var2,k=k.use),data=dat.r),silent=T)
       if(class(fit)[[1]]!="try-error"){
        r.sq=summary(fit)$r.sq
       if(r.sq>0){r.est=sqrt(r.sq)}else{r.est=0}}
     }
    return(r.est)}

  out.cor.dat=apply(test.mat,MARGIN=1,FUN=find.r.est.2)
  test.mat$r.sq=out.cor.dat

  # now put the results in the correlation matrix
  out.cor.mat=matrix(NA,nrow=ncol(dat),ncol=ncol(dat))
  rownames(out.cor.mat)=colnames(dat)
  colnames(out.cor.mat)=colnames(dat)
  diag(out.cor.mat)=1

  for(r in 1:nrow(test.mat)){
     out.cor.mat[test.mat$rows[r],test.mat$cols[r]]=test.mat$r.sq[r]
  }

  return(out.cor.mat)
}






















