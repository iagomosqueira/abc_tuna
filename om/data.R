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
range(stk, c("minfbar", "maxfbar")) <- c(1, 12)

# FLIndices
ids <- readFLIBss3('../data/base')


# -- LOAD abc4

run4 <- mget(load("../v2/runs/alb_abc_run4.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

# EXTRACT output for all iters
out4 <- mc.output(run4$mcvars)


# -- CREATE FLStock for om: 2000-2020, no seasons, 500 iters

stk4 <- propagate(simplify(window(stk, start=2000), 'season'), 500)

# ADD m, sum accross seasons
m(stk4)[, ac(2000:2020)] <- seasonSums(out4$m)

# ADD stock.n Q1
stock.n(stk4)[, ac(2000:2020)] <- out4$stock.n[,,,1]

# age 0 Q1 = Q4 * exp(m Q1-3)
stock.n(stk4)[1, ac(2000:2020)] <- out4$stock.n[1, ac(2000:2020),, 4] *
  exp(seasonSums(out4$m[1, ac(2000:2020),, 1:3]))

# ADD harvest rate, mean over fleets (area) & seasons
harvest(stk4)[, ac(2000:2020)] <- areaMeans(seasonMeans(out4$hra))

# spwn, mid Q4
m.spwn(stk4) <- harvest.spwn(stk4) <- 0.83

# TODO: CHECK

# EXTEND
fut4 <- fwdWindow(stk4, end=2040)

# ASSIGN selectivities
harvest(fut4)[, ac(2021:2040)] <- areaMeans(seasonMeans(out4$catch.sel))


# -- BUILD FLSR

srp <- out4$srpars

# CORRECT srpars$R0 for Q1
srp$R0 <- srp$R0 * exp(c(seasonSums(out4$m[1,'2000','F', 1:3])))
srp$v <- srp$v * exp(c(seasonSums(out4$m[1,'2000','F', 1:3])))

srr4 <- FLSR(model=bevholtss3, params=srp)

# -- BUILD refpts

rps <- out4$refpts
rps$R0 <- srp$R0
rps$B0 <- srp$v


# TODO: fwd(F=0), rec=R0

tes <- fwd(fut4,
  sr=FLQuant(c(srp$R0), dimnames=list(age='0', year=2021:2040,
    iter=seq(500))),
  fbar=FLQuant(0, dimnames=list(year=2021:2040)))

plot(tes)
plot(ssb(tes)[,,'F'] / rps$B0)

# TODO: fwd(F=0), rec=SRR

tes <- fwd(fut4, sr=srr4, fbar=FLQuant(0, dimnames=list(year=2021:2040)))

plot(tes)
plot(ssb(tes)[,,'F'] / rps$B0)

# TODO: 

harvest(fut4)[, ac(2000:2020)] <- recomputeHarvest(fut4)[, ac(2000:2020)]

tes <- fwd(fut4, sr=srr4, fbar=FLQuant(0, dimnames=list(year=2021:2040)))

plot(tes)
plot(ssb(tes)[,,'F'] / rps$B0)


# -- BUILD FLom

om4 <- FLom(stock=stk4, sr=srr4, refpts=out4$refpts)

om4 <- fwdWindow(om4, end=2040)

# DEVIANCES

rdev <- residuals(unitSums(rec(stk4)),
  predict(srr4, ssb=ssb(stk4)[,,'F']), 'log')

# BUG: deviances<- append fails
residuals(sr(om4)) <- rlnormar1(500, meanlog=0, sdlog=sqrt(yearVars(rdev)),
  rho=rho(rdev), years=2000:2040)


# -- BUILD FLoem

# PROPAGATE no unit/season stk
estk <- propagate(simplify(stk), 500)

# ADD OM stk in 2000:2020
# NOTE: SHOULD oem be single sex?
estk[, ac(2000:2020)] <- simplify(stk4, 'unit')

# PLOT
plot(FLStocks(SS3=simplify(stk), OEM=estk)) +
  ggtitle("SS3 SA vs. OEM")

# idx - LLCPUE1

# SET index.q
index.q <- Reduce(sbind, lapply(1:4, function(i)
  log(window(index(ids[[i]]), start=2000) / out4$index.hat[,,,i])))

ll1 <- propagate(ids[[1]], 500)
range(ll1, c('startf', 'endf')) <- c(0, 0.25)

index(ll1)[, ac(2000:2020)] <- exp(index.q[,,,1]) * out4$index.hat[,,,1]

plot(FLIndices(SS=ll1, ABC=survey(stk4, ll1[, ac(2000:2020)])))

idx <- FLIndices(LLCPUE1=idx)


# LLCPUE1 <- FLIndexBiomass(
#   index=out4$index.hat[,,,1],
#   index.q=index.q[,,,1],
#   index.var=out4$index.hat[,,,1] %=% 0.2,
#   catch.wt=unitMeans(stock.wt(stk4)),
#   sel.pattern=expand(unitMeans(out4$catch.sel[,,,1,1]),
#   year=2000:2020, fill=TRUE),
#   range=c(startf=0, endf=0.25))
# 
# LLCPUE1 <- qapply(LLCPUE1, function(x) {
#   dimnames(x)$unit <- 'unique'
#   dimnames(x)$season <- '1'
#   dimnames(x)$area <- 'unique'
#   return(x)
#   }
# )


# deviances
devs <- list(stk=FLQuants(catch.n=rlnorm(500, unitSums(catch.n(om4) %=% 0),
  0.3)))

oem4 <- FLoem(observations=list(stk=fwdWindow(estk, end=2040),
  idx=fwdWindow(idx, end=2040)),
  deviances=devs, method=sampling.oem)

# TODO: verify(oem, om)

# -- SAVE

save(om4, oem4, file='om4.rda', compress='xz')
