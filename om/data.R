# om.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/abc_tuna/om/om.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(mse)
library(ss3om)

source("utilities.R")


# -- LOAD SS3 FLStock & FLIndices

stk <- readFLSss3('../data/base')
# BUG: PICK UP from starter.ss
range(stk, c("minfbar", "maxfbar")) <- c(1, 12)

# FLIndices
ids <- readFLIBss3('../data/base')


# -- LOAD abc4

run <- mget(load("../v2/runs/alb_abc_run4b.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

# EXTRACT output for all iters
out <- mc.output(run$mcvars, run$C)


# -- CREATE FLStock for om: 2000-2020, no seasons, 500 iters

stko <- propagate(simplify(window(stk, start=2000), 'season'), 500)

# ADD m, sum accross seasons
m(stko)[, ac(2000:2020)] <- seasonSums(out$m)

# ADD stock.n Q1
stock.n(stko)[, ac(2000:2020)] <- out$stock.n[,,,1]

# age 0 Q1 = Q4 * exp(m Q1-3)
stock.n(stko)[1, ac(2000:2020)] <- out$stock.n[1, ac(2000:2020),, 4] *
  exp(seasonSums(out$m[1, ac(2000:2020),, 1:3]))

# ADD harvest rate, mean over fleets (area) & seasons

# hy = 1-(1-h1)*(1-h2)*(1-h3)*(1-h4)
hrya <- 1 - (1 - out$hra[,,,1]) * (1 - out$hra[,,,2]) *
  (1 - out$hra[,,,3]) * (1 - out$hra[,,,4])

harvest(stko)[, ac(2000:2020)] <- areaMeans(hrya)
units(harvest(stko)) <- 'hr'

# spwn, mid Q4
m.spwn(stko) <- harvest.spwn(stko) <- 0.83

# EXTEND
futo <- fwdWindow(stko, end=2040)

# ASSIGN selectivities
# TODO: catch.sel(FLStock/FLom) <- FLQuant +flr
harvest(futo)[, ac(2021:2040)] <- areaMeans(seasonMeans(out$catch.sel))


# -- BUILD FLSR

srp <- out$srpars

# CORRECT srpars$R0 for Q1
srp$R0 <- srp$R0 * exp(c(seasonSums(out$m[1, 1, 'F', 1:3])))

# spr0y(stk) / rps$B0
corr <- yearMeans((FLSRTMB::spr0y(stko) / 1000)) / (srp$v/srp$R0)
srp$v <- srp$v * exp(c(seasonSums(out$m[1, 1,'F', 1:3]))) * c(corr)

srr <- FLSR(model=bevholtss3, params=srp)

# -- BUILD refpts

rps <- out$refpts
rps$R0 <- srp$R0
rps$B0 <- srp$v

# TODO: fwd(F=0), rec = R0 / srr

tes0 <- fwd(fwdWindow(stko, end=2100),
  sr=FLQuant(rep(srp$R0, each=80), dimnames=list(age='0',
    year=2021:2100, iter=1:500)),
  fbar=FLQuant(0, dimnames=list(year=2021:2100)))

tesR <- fwd(fwdWindow(stko, end=2100), sr=srr,
  fbar=FLQuant(0, dimnames=list(year=2021:2100)))

plot(tes0, tesR)

plot(ssb(tes0)[,,'F'] / rps$B0, ssb(tesR)[,,'F'] / rps$B0) +
  geom_hline(yintercept=1)

# PLOT on recalculated B0
plot(ssb(tes0)[,,'F'] %/% (yearMeans((FLSRTMB::spr0y(stko) / 1000)) * srp$R0)) +
  geom_hline(yintercept=1)

# PREDICT R0
predict(srr, ssb=FLQuant(c(srp$v), dimnames=list(year=1, iter=1:500)))
srp$R0

# -- BUILD FLom

om <- fwdWindow(FLom(stock=stko, sr=srr, refpts=rps), end=2040)

# DEVIANCES

rdev <- residuals(unitSums(rec(stko)),
  predict(srr, ssb=ssb(stko)[,,'F']), 'log')

# BUG: deviances<- append fails
residuals(sr(om)) <- rlnormar1(500, meanlog=0, sdlog=sqrt(yearVars(rdev)),
  rho=rho(rdev), years=2000:2040)


# -- BUILD FLoem

# stk
estk <- propagate(simplify(stk), 500)

# ADD OM stk in 2000:2020, single sex
estk[, ac(2000:2020)] <- simplify(stko, 'unit')

# idx - LLCPUE1_Q1-Q4

# SET index.q
index.q <- Reduce(sbind, lapply(1:4, function(i)
  log(window(index(ids[[i]]), start=2000) / out$index.hat[,,,i])))

# PICK index
ind <- propagate(ids[['LLCPUE1_Q1']], 500)
range(ind, c('startf', 'endf')) <- c(0, 0.25)

# BUG:
effort(ind) <- unitSums(effort(ind))
sel.pattern(ind) <- unitMeans(sel.pattern(ind))
catch.n(ind) <- unitMeans(catch.n(ind))

# ADD bias-corrected index.q
index.q(ind)[, ac(2000:2020)] <- 
  exp(index.q[,,,1] - 0.5 * sqrt(iterVars(index.q[,,,1])) ^ 2)

# ADD catch.wt Q1
catch.wt(ind) <- catch.wt(stk)[,,,1][, ac(1975:2020)]

# ADD sel.pattern: Q1, fleet (area) 1, Q1
sel.pattern(ind)[, ac(2000:2020)] <- expand(out$catch.sel[,,,1,1],
  year=2000:2020, fill=TRUE)

idx <- FLIndices(LLCPUE1=ind)

survey(stock(om), idx)

# deviances
devs <- list(stk=FLQuants(catch.n=rlnorm(500, unitSums(catch.n(om) %=% 0),
  0.3)))

oem <- FLoem(observations=list(stk=fwdWindow(estk, end=2040),
  idx=fwdWindow(idx, end=2040)),
  deviances=devs, method=sampling.oem)

# TODO: verify(oem, om)

# -- SAVE

save(om, oem, file='data/om4b.rda', compress='xz')
