# model.R - DESC
# abc_tuna/albabc/model.R

# Copyright (c) WUR & CSIRO, 2023.
# Authors: Richard HILLARY (CSIRO) <rich.hillary@csiro.au>
#          Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(Rcpp)
source("utilities.R")

sourceCpp("src/init_pdyn.cpp")
sourceCpp("src/msy_pdyn.cpp")
sourceCpp("src/pdyn.cpp")
sourceCpp("src/pdyn_lfcpue.cpp")

# LOAD data (stk, inp)

load("data/data.RData")

# TEST sim()

system.time(
pro <- sim(R0=1e6, dep=0.5, h=0.75)
)

# COMPARE data inputs (inp) and simulator outputs (pro)

names(inp)
names(pro)

# I
dim(pro$I)
dim(inp$I)

# C
dim(pro$C)
dim(inp$C)

# LF
dim(pro$LF)
dim(inp$LF)

