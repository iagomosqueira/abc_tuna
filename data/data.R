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

# LOOKUP table: group LL & Other by area, keep PS, drop DN.

lookup <- data.table(
  fleet=c(seq(1, 16), 19, seq(20, 23)),
  unit=c(rep(seq(1, 4), 4), 5, rep(6, 4)))
setkey(lookup, "fleet")

# - EXTRACT catch per fleet, year and season

catches <- data.table(dat$catch)[year >= 1954]
setkey(catches, "fleet")

ncatches <- catches[lookup, on="fleet"]
setnames(ncatches, c('seas', 'catch'), c('season', 'data'))

ncatches <- ncatches[, .(data=sum(data, na.rm=TRUE)), by=.(year, season, unit)]

setorder(ncatches, year, season, unit)

ggplot(ncatches[unit %in% seq(1:4)],
  aes(x=ISOdate(year, season * 3, 1), y=data, group=unit)) +
  geom_line() + facet_wrap(~unit)+ geom_point()

ggplot(ncatches,
  aes(x=ISOdate(year, season * 3, 1), y=data, group=unit)) +
  geom_line() + facet_wrap(~unit, scales='free')+ geom_point()


# - EXTRACT length composition data

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
