# data.R - LOAD SS3 model and add ABC output
# abc_tuna/om/data.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(mse)

source("utilities.R")

# LOAD SS3 base run

load("output/base.rda")

# --- LOAD abc run

run4 <- mget(load("model/alb_abc_run4.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

run5b <- mget(load("model/alb_abc_run5b.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

run6 <- mget(load("model/alb_abc_run6.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

names(run4$mcvars[[1]])
names(run5b$mcvars[[1]])
names(run6$mcvars[[1]])

run <- mget(load("model/alb_abc_run5b.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

run <- mget(load("model/alb_abc_run4a.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

run <- mget(load("model/alb_abc_run6.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

run <- mget(load("model/alb_abc_run6a.rda", verbose=FALSE,
  envir=(.NE <- new.env())), envir=.NE)

# EXTRACT output for all iters
out <- mc.output(mcvars, C)


# --- CREATE om

its <- dims(out$m)$iter

# - FLBiol (stock.n, m)

bio <- propagate(window(sbio, start=2000), its)

n(bio) <- out$stock.n
m(bio) <- out$m

sr(bio) <- predictModel(model=bevholtss3()$model, params=out$srpars)

# spwn, mid Q4
spwn(bio)[,,,4] <- 0.5

# FIX mat BUG: CHECK readFLSss3 +ss3om 
mat(bio)[,,'F',1:4] <- mat(bio)[,,'F',1]

# TODO: DEVIANCES

deviances(bio)

# BUG:
predict(sr(bio), ssb=out$ssb)
seasonSums(unitSums(rec(bio)))

# - FLFisheries (catch.n, catch.sel)

# FLCatch(es)
cas <- Map(function(x, y) FLCatch(landings.n=x, landings.wt=wt(bio),
  catch.sel=y, discards.n=x %=% 0, discards.wt=wt(bio)),
  x=divide(out$catch.n, 5), y=divide(out$catch.sel %*% (out$catch.n %=% 1), 5))

# FLFisheries
fis <- FLFisheries(lapply(cas, function(x)
  FLFishery(effort=unitSums(catch(cas[[1]])) %=% 0, ALB=x)))

names(fis) <- c(paste0("LL", 1:4), "PS", "Other")

# PLOT
plot(fis)

om <- FLombf(biols=FLBiols(ALB=bio), fisheries=fis,
  refpts=FLPars(ALB=out$refpts))

# ADD hr
# TODO: SET harvest as slot
attr(om, 'harvest') <- FLQuants(ALB=expand(n(biol(om)), area=1:6) %=% as.numeric(NA))

# TEST om


# - BUILD oem

# idx: FLIndexBiomass by season, with sel.pattern by sex

NW <- FLIndexBiomass(index=out$index.hat %*% out$index.q,
  index.q=expand(out$index.q, year=2000:2020),
  sel.pattern=expand(out$sel, year=2000:2020),
  catch.wt=wt(biol(om)),
  range=c(startf=0.5, endf=0.5))

# stk: no units
oem <- FLoem(
    observations=list(ALB=list(idx=FLIndices(NW=NW),
      stk=simplify(stock(om)[[1]], 'unit'))),
    method=sampling.oem)

survey(observations(oem)$ALB$stk, observations(oem)$ALB$idx)



# FLStock


# -- BUILD FLoem

# stk

# idx - LLCPUE1_Q1-Q4

# TODO: verify(oem, om)


# -- SAVE

save(om, oem, file='output/om5b.rda', compress='xz')
