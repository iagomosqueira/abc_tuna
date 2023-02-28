# utilities.R - DESC
# /home/mosquia/Active/ABC_tuna+iotc/abc_tuna/v2/utilities.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


logit <- function(x){
  return(log(x/(1-x)))
}

ilogit <- function(x){
  return(1/(1+exp(-x)))
}

