

# This script includes the basic steps required to create and
#    build the FFSgam package

library(devtools)
library(roxygen2)
library(knitr)
library(R.rsp)
library(digest)

setwd("C:/Users/rfisher/OneDrive - Australian Institute of Marine Science/Documents/AIMS/EcologicalRiskModelling/Full_subsets_methods/FSSgam_package")
# copy function files over to FFSgam/R

devtools::document()

use_package("doSNOW")
use_package("MuMIn")
use_package("gamm4")
use_package("mgcv")
use_package("nnet")
build()

#devtools::install_github("beckyfisher/FSSgam_package")
