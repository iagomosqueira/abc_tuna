# data.R - DESC
# abc_tuna/data/data.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(FLCore)
library(ggplotFL)
library(ss3om)
library(data.table)

# LOAD SS output

load('run.RData')

# LOAD SS inputs

dat <- SS_readdat('sa/abt.dat')


# --- OUTPUT data matrices / arrays

# LOOKUP table: group LL & Other by area, keep PS, drop DN.

# One fleet per season
lookup <- data.table(
  fleet=c(seq(1, 16), 19, seq(20, 23)),
  unit=c(rep(seq(1, 4), 4), 5, rep(6, 4)))
setkey(lookup, "fleet")

# One fleet per area
lookup <- data.table(
  fleet=c(seq(1, 16), 19, seq(20, 23)),
  unit=c(rep(seq(1, 4), each=4), 5, rep(6, 4)))
setkey(lookup, "fleet")


# - EXTRACT catch per fleet, year and season

catches <- data.table(dat$catch)[year >= 1954]
setkey(catches, "fleet")

# RECODE fleets
ncatches <- catches[lookup, on="fleet"]
setnames(ncatches, c('seas', 'catch'), c('season', 'data'))

catches <- ncatches[, .(data=sum(data, na.rm=TRUE)), by=.(year, season, unit)]

setorder(catches, year, season, unit)

ggplot(catches[unit %in% seq(1:4)],
  aes(x=ISOdate(year, season * 3, 1), y=data, group=unit)) +
  geom_line() + facet_wrap(~unit)+ geom_point() +
  ggtitle("LL fleets by area")


# - EXTRACT length composition data

lencomp <- data.table(dat$lencomp)

# SELECT 'f' columns
lencomp <- lencomp[, c(1,2,3,6,7:61)]
setnames(lencomp, c("year", "season", "fleet", "Nsamp", seq(30, 138, by=2)))

# RESHAPE to long
lencomp <- melt(lencomp, id=c("year", "season", "fleet", "Nsamp"),
  measure=seq(5, 59), variable.name = "length", value.name = "n")

# FIX season
lencomp[, season:=season - 1.5]

# A new fleet
lencomp <- lencomp[lookup[fleet < 20,], on="fleet"]

#
ggplot(lencomp[year > 2010], aes(x=length, y=n, group=season)) +
  geom_col() + facet_grid(year~unit, scales="free")


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
save(catches, lencomp, cpue, caa, waa, mat, m, file="data.RData",
  compress="xz")
