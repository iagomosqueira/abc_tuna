# om.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/abc_tuna/om/om.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(mse)
library(ss3om)

source("utilities.R")


# -- LOAD SS3

# FLStock
stk <- readFLSss3('boot/data/base')
range(stk, c("minfbar", "maxfbar")) <- c(1, 12)

# FLIndices
ids <- readFLIBss3('boot/data/base')

# iol + FLFisheries
out <- readOutputss3('boot/data/base')
bfs <- buildFLBFss330(out)
bio <- bfs$biol
fis <- bfs$fisheries

# -- LOAD ABC (5b)

run <- mget(load("model/alb_abc_run5b.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

# EXTRACT output for all iters
out <- mc.output(run$mcvars, run$C)


# -- UPDATE FLBiol + FLFisheries

# - bio
bio <- propagate(window(bio, start=2000), 500)

# ADD m, sum accross seasons
m(bio) <- out$m

# ADD stock.n
n(bio) <- out$stock.n

# spwn, mid Q4
spwn(bio)[,,,4] <- 0.5

# FIX mat BUG: CHECK readFLSss3 +ss3om 
mat(bio)[,,'F',1:4] <- mat(bio)[,,'F',1]

# - fisheries

# landings,n, landings.wt, catch.sel, catch.q

# SPLIT catch.n and selex by fishery
can <- divide(out$catch.n, 5)
sel <- divide(out$catch.sel %*% (out$catch.n %=% 1), 5)

# NOTE: effort from partial Fs
pfs <- parallel::mclapply(can, function(x) harvest(n(bio), x, m(bio)),
  mc.cores=2)
efs <- Map(function(x, y) quantMeans(unitMeans(x %/% y)), x=pfs, y=sel)

efs <- list(FLQuant(0, dimnames=list(quant="effort", year=2000:2020)))

# BUILD fisheries
fis <- FLFisheries(Map(function(x, y, ef) FLFishery(effort=ef,
  ALB=FLCatch(landings.n=x, landings.wt=wt(bio),
    discards.n=x %=% 0 , discards.wt=wt(bio), catch.sel=y)),
  x=can, y=sel, ef=efs))
names(fis) <- c(paste0("LL", 1:4), "PS", "Other")

# BUILD OM
om <- FLombf(biols=FLBiols(ALB=bio),
  fisheries=FLFisheries(lapply(window(fis, start=2000), propagate, 500)))

#

save(om, file='model/om.rda')


# ---------

rr <- predictModel(model=rec~a, params=FLPar(c(-1,0,0,1.343e4 * 2),
  dimnames=list(params='a', season=1:4, iter=1)))

tes <- fwd(futo, sr=rr, catch=unitSums(catch(futo)[, ac(2021:2035)]) %=% 0.0001)

rec(tes)[, ac(2019:2023), 'F']
ssb(tes)[, ac(2021:2023)]

plot(window(stko, end=2020), window(tes, end=2035))
plot(window(tes, end=2035))

plot(ssb(tes)[,,,4] / srp$v)

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
