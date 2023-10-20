# om.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/abc_tuna/om/om.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# TODO: window(stk), keep full in observations(oem)


library(mse)
library(ss3om)

source("utilities.R")


# -- LOAD SS3 FLStock

stk <- readFLSss3('../data/base')
range(stk, c("minfbar", "maxfbar")) <- c(1, 12)

# -- LOAD abc4

run4 <- mget(load("../v2/runs/alb_abc_run4.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

# EXTRACT output for all iters
out4 <- mc.output(run4$mcvars)

# BUG: DIFF caa vs. catch.n(stk)
catch.n(stk) <- run4$caa

# -- CREATE FLStock for om: 2000-2020, no seasons, 500 iters

ostk <- window(stk, start=2000)
units(harvest(ostk)) <- "hr"

ostk <- simplify(ostk, 'season')
ostk <- propagate(ostk, 500)

# ASSIGN out4 FLStock

stk4 <- ostk

# m
m(stk4)[, ac(2000:2020)] <- seasonSums(out4$m)

# stock.n
stock.n(stk4)[, ac(2000:2020)] <- out4$stock.n[,,,1]

# Q4 y-1 * exp(-m) NOTE: rec Q1 after simplify, from Q4 y-1 + Z?
stock.n(stk4)[1, ac(2001:2020)] <- out4$stock.n[1, ac(2000:2019),, 4] /
  exp(-out4$m[1, ac(2000:2019),, 4])
stock.n(stk4)[1, ac(2000)] <- stock.n(stk)[1, ac(1999),, 4] /
  exp(-m(stk)[1, ac(1999),, 4])

# harvest
harvest(stk4)[, ac(2000:2020)] <- computeHarvest(stk4)[, ac(2000:2020)]

harvest(stk4) <- recomputeHarvest(stk4)

# spwn
m.spwn(stk4) <- harvest.spwn(stk4) <- 0.83

# EXTEND
fut4 <- fwdWindow(stk4, end=2040)

# ASSIGN selectivities
harvest(fut4)[, ac(2021:2040)] <- areaMeans(seasonMeans(out4$catch.sel))

# -- BUILD FLSR

# CORRECT srpars$R0 for Q1
srp <- out4$srpars
srp$R0 <- srp$R0 / exp(-c(seasonSums(out4$m[1,'2000','F', 1:3])))

srr4 <- FLSR(model=bevholtss3, params=srp)

# -- BUILD FLom

om4 <- FLom(stock=stk4, sr=srr4, refpts=out4$refpts)

om4 <- fwdWindow(om4, end=2040)

# DEVIANCES

rdev <- residuals(unitSums(rec(stk4)),
  predict(srr4, ssb=ssb(stk4)[,,'F']), 'log')

# BUG: deviances<- append fails
residuals(sr(om4)) <- rlnormar1(500, meanlog=0, sdlog=sqrt(yearVars(rdev)),
  rho=rho(rdev), years=2000:2040)

# test.R:fwd(F=0)

# -- BUILD FLoem

# stk
estk <- propagate(simplify(stk, c('season', 'unit')), 500)
estk[, ac(2000:2020)] <- simplify(stk4, 'unit')

# BUG: JUMP in ssb

plot(simplify(stk), estk)

# idx

LLCPUE1 <- FLIndexBiomass(
  index=out4$index.hat[,,,1],
  index.var=out4$index.hat[,,,1] %=% 0.2,
  catch.wt=unitMeans(stock.wt(stk4)),
  sel.pattern=expand(unitMeans(out4$catch.sel[,,,1,1]),
  year=2000:2020, fill=TRUE),
  range=c(startf=0, endf=0.25))

LLCPUE1 <- qapply(LLCPUE1, function(x) {
  dimnames(x)$unit <- 'unique'
  dimnames(x)$season <- '1'
  dimnames(x)$area <- 'unique'
  return(x)
  }
)

idx <- FLIndices(LLCPUE1=LLCPUE1)

# BUG: Your exploitation rate is not defined as F, cannot be added to M
index.q(idx[[1]]) <- computeQ(idx,
  simplify(stk4, 'unit', harvest=FALSE),
  FLQuants(LLCPUE1_Q1=out4$index.hat[,,,1]))[[1]]

# deviances

devs <- list(stk=FLQuants(catch.n=rlnorm(500, unitSums(catch.n(om4) %=% 0),
  0.3)))

oem4 <- FLoem(observations=list(stk=fwdWindow(estk, end=2040),
  idx=fwdWindow(idx, end=2040)),
  deviances=devs, method=sampling.oem)

# TODO: verify(oem, om)

# -- SAVE
save(om4, oem4, file='om4.rda', compress='xz')
