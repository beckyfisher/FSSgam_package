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
# identically by check.correlations() and check.non.linear.correlations().
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
