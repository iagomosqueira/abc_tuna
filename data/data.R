# data.R - DESC
# /home/mosquia/Active/Doing/ABC_tuna+iotc/abc_tuna/data/data.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(FLCore)
library(ss3om)
library(data.table)

# LOAD SS output

load('run.RData')

# LOAD SS inputs

dat <- SS_readdat('sa/abt.dat')


# --- OUTPUT data matrices / arrays

# - samples

# catch per fleet, year and season: data.
catches <- data.table(dat$catch)[year >= 1954]

# length composition data: fleets 1-19
lencomp <- data.table(dat$lencomp)

# CPUE
cpue <- data.table(dat$CPUE)

# - stock

# CaA
caa <- catch.n(run)

# WaA
waa <- stock.wt(run)

# maturity
mat <- mat(run)

# M
m <- m(run)


# SAVE

save(catches, lencomp, cpue, caa, waa, mat, m,
    file="data.RData", compress="xz")
