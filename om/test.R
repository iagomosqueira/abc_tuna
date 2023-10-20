# test.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/abc_tuna/om/test.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

library(mse)

load('om4.rda')


# - fwd(F=0)

tf0 <- fwd(om4, control=fwdControl(year=2021:2040, quant='fbar', value=0))

# BUG: PLOT to use FLStock plot quants
plot(stock(tf0))

plot(FLStocks(OM=observations(oem4)$stk, F0=stock(tf0)))

# NOTE: SHOULDN'T ssb recover to B0?
plot(ssb(stock(tf0))[,,'F'] / refpts(om4)$B0) +
  geom_hline(yintercept=1)

# - fwd(F=runif(0.1, 0.3))

fu <- fwd(om4, control=fwdControl(year=2021:2040, quant='fbar',
  value=runif(20, 0.1, 0.3)))

plot(stock(tfu))

# - mp(perfect, fixedF)

control <- mpCtrl(list(
  est=mseCtrl(method=perfect.sa),
  hcr=mseCtrl(method=fixedF.hcr, args=list(ftrg=0.1))
))

mrun <- mp(om4, oem4, ctrl=control, args=list(iy=2020, frq=3))

plot(FLStocks(OM=window(stock(om4), end=2020), RUN=stock(mrun)))

# - mp(perfect, hockeystick)

control <- mpCtrl(list(
  est = mseCtrl(method=perfect.sa),
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0.10, trigger=0.40, target=mean(refpts(om4)$MSY) * 0.90,
      metric=depletion, output="catch", dlow=0.85, dupp=1.15))))

prun <- mp(om4, oem4, ctrl=control, args=list(iy=2020, frq=3, fy=2025))

plot(FLStocks(OM=window(stock(om4), end=2020), RUN=stock(prun)))

# - mp(cpue)

control <- mpCtrl(list(
  est = mseCtrl(method=cpue.ind, args=list(index=1)),
  hcr = mseCtrl(method=cpue.hcr, args=list(k1=0.2, k2=0.2, k3=0.2, k4=0.2,
    dlow = NA, dupp = NA, target=6))))
      
crun <- mp(om4, oem4, ctrl=control, args=list(iy=2020, frq=3, fy=2032))


# sel.pattern for 1 sex
sel.pattern(observations(oem4)$idx[[1]])  <- unitMeans(sel.pattern(observations(oem4)$idx[[1]]))

# catch.wt
catch.wt(observations(oem4)$idx[[1]]) <- 
  unitMeans(catch.wt(observations(oem4)$idx[[1]]))
catch.wt(observations(oem4)$idx[[1]]) <- 
  c(stock.wt(observations(oem4)$stk)[, ac(1975:2040)])

plot(FLStocks(OM=window(stock(om4), end=2020), RUN=stock(crun)))
