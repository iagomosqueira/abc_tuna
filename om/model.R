# test.R - DESC
# abc_tuna/om/test.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# BUG: PLOT to use FLStock plot quants

library(mse)

library(doFuture)
plan(multicore, workers=3)

load('data/om4b.rda')

mets <- list(Rec=function(x) unitSums(rec(x)),
  SB=function(x) unitSums(ssb(x)),
  Catch=function(x) unitSums(catch(x)),
  F=function(x) unitMeans(fbar(x)),
  HR=function(x) unitMeans(quantMeans(harvest(x)[1:12,])))

# - fwd(F=0)

tf0 <- fwd(om, control=fwdControl(year=2021:2040, quant='fbar', value=0))

plot(stock(tf0), metrics=mets) +
  geom_vline(xintercept=2020, linetype=3) +
  ggtitle("F=0")

# NOTE: SHOULDN'T ssb recover to B0?
plot(ssb(stock(tf0))[,,'F'] / refpts(om)$B0) +
  geom_hline(yintercept=1) +
  ggtitle("F=0, relative")


# - fwd(F=runif(0.1, 0.3))

tfu <- fwd(om4, control=fwdControl(year=2021:2040, quant='fbar',
  value=runif(20, 0.1, 0.3)))

plot(stock(tfu), metrics=mets) +
  geom_vline(xintercept=2020, linetype=3) +
  ggtitle("F=0")

plot(FLStocks(OM=observations(oem4)$stk, FU=stock(tfu)), metrics=mets) +
  ggtitle("F=runif()")

# - mp(perfect, fixedF)

control <- mpCtrl(list(
  est=mseCtrl(method=perfect.sa),
  hcr=mseCtrl(method=fixedF.hcr, args=list(ftrg=0.10))
))

mrun <- mp(om, ctrl=control, args=list(iy=2020, frq=3))

plot(FLStocks(OM=window(stock(om), end=2020), RUN=stock(mrun))) +
  ggtitle("perfect.sa + fixedF.hcr(0.1)")

# - mp(perfect, hockeystick)

control <- mpCtrl(list(
  est = mseCtrl(method=perfect.sa),
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0.10, trigger=0.40, target=mean(refpts(om)$MSY) * 0.90,
      metric=depletion, output="catch", dlow=0.85, dupp=1.15))))

prun <- mp(om, ctrl=control, args=list(iy=2020, frq=3))

plot(FLStocks(OM=window(stock(om), end=2020), RUN=stock(prun)),
  metrics=mets) +
  ggtitle("perfect.sa + hockeystick.hcr(0.90*MSY)")

# - mp(cpue)

control <- mpCtrl(list(
  est = mseCtrl(method=cpue.ind, args=list(index=1)),
  hcr = mseCtrl(method=cpue.hcr, args=list(k1=0.2, k2=0.2, k3=0.2, k4=0.2,
    dlow = NA, dupp = NA, target=6))))
      
crun <- mp(om, oem, ctrl=control, args=list(iy=2020, fy=2030, frq=3))

plot(FLStocks(OM=window(stock(om), end=2020), RUN=stock(crun))) +
  ggtitle("cpue.ind + cpue.hcr")

plan(multisession, workers=2)

cruns <- mps(om, oem, ctrl=control, args=list(iy=2020, fy=2030, frq=3),
  hcr=list(target=seq(0.1, 0.5, length=8)))

data(statistics)
tperiod <- list(2020 + seq(11, 15))

# P(kobe=green) = 0.5
jabba_hcst_05_targ <- tunebisect(om, oem=oem, control=control, 
  args=list(iy=2020, frq=3), metrics=mets[c("SB", "F")],
  statistic=statistics["green"], years=tperiod,
  tune=list(target=c(0.1, 5)), prob=0.5, tol=0.02,
  maxit=12, parallel=TRUE)

# BUG: refpts for tuning
refpts(om)$SBMSY<-refpts(om)$B0*0.40
refpts(om)$FMSY<-refpts(om)$SBMSY
refpts(om)$FMSY[]<-0.15

# SAVE
save(tf0, tfu, mrun, prun, crun, file="model/runs4.rda", compress="xz")

# TODO: CHECK

# sel.pattern for 1 sex
sel.pattern(observations(oem4)$idx[[1]])  <- unitMeans(sel.pattern(observations(oem4)$idx[[1]]))

# catch.wt
catch.wt(observations(oem4)$idx[[1]]) <- 
  unitMeans(catch.wt(observations(oem4)$idx[[1]]))
catch.wt(observations(oem4)$idx[[1]]) <- 
  c(stock.wt(observations(oem4)$stk)[, ac(1975:2040)])

plot(FLStocks(OM=window(stock(om4), end=2020), RUN=stock(crun)))
