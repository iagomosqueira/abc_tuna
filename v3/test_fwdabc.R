# fwdabc.R - DESC
# abc_tuna/v3/fwdabc.R

# Copyright (c) WUR, 2024.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

library(Rcpp)
library(FLCore)
library(ggplotFL)
library(parallel)
library(mvtnorm)
library(mse)
source("utilities.R")

sourceCpp("utilities/init_pdyn.cpp")
sourceCpp("utilities/msy_pdyn.cpp")
sourceCpp("utilities/pdyn_lfcpue.cpp")

# LOAD pcbar (catch proportions by fleet & season)
pcbar <- as.matrix(fread('data/pcbar.dat'))

# LOAD ALK [len, age, season, sex]
load('data/pla.rda')

# LOAD om
load('data/om5b.rda')

# ASSIGN args to om@projection
args(projection(om)) <- list(pla=pla, pcbar=pcbar)

# ADD HR attribute.  TODO: SET as slot
attr(om, 'harvest') <- FLQuants(ALB=setunits(expand(n(biol(om)), area=1:6) %=% as.numeric(NA), 'hr'))

attr(om, 'hrbar') <- FLQuants(ALB=setunits(expand(quantSums(n(biol(om))), area=1:6) %=% as.numeric(NA), 'hr'))

# --- TESTS

its <- seq(100)
iy <- 2017

# - TEST fwd 2011:2020 with historical catch
# BUG: differences maybe due to recruitment deviances

ctrl <- fwdControl(year=iy:2020, quant="catch",
  value=c(unitSums(seasonSums(Reduce('+', catch(om))))[, ac(iy:2020),,,,1]))

tes2 <- fwdabc.om(iter(om, its), ctrl, pcbar=pcbar, pla=pla)

quantSums(unitSums(seasonSums(Reduce('+', catch(fisheries(tes2$om))))))

(plot(biol(iter(tes2$om, its))) + geom_vline(xintercept=ISOdate(iy,1,1))) +
(plot(biol(iter(om, its))) + geom_vline(xintercept=ISOdate(iy,1,1)))

Reduce('+', lapply(fisheries(om), function(x)
  unitSums(seasonSums(landings(x[[1]])[,2017:2020,,,,]))))

ref <- Reduce('+', lapply(fisheries(om), function(x) unitSums(seasonSums(catch(x)[[1]][, ac(2017:2020)]))))

tes <- Reduce('+', lapply(fisheries(tes2$om), function(x) unitSums(seasonSums(catch(x)[[1]][, ac(2017:2020)]))))

com <- ref/tes
all.equal(c(com[, '2017']), rep(1, 500), tolerance=0.001)

# - TEST fwd iy:2020 with no catch

ctrl <- fwdControl(year=iy:2020, quant="catch", value=0)

tes1 <- fwdabc.om(iter(om, its), ctrl, pcbar=pcbar, pla=pla)

(plot(biol(tes1$om)) + geom_vline(xintercept=ISOdate(iy,1,1))) +
(plot(biol(iter(om, its))) + geom_vline(xintercept=ISOdate(iy,1,1)))



# TODO: RUN with actual rec

# sr(biol(om)) <- predictModel(model=rec~a, params=FLPar(c(rec(om)[[1]][,,'F',]),
#   dimnames=list(params='a', year=2000:2020, season=1:4, iter=1:10)))

# - TEST fwd 2011:2020 for different catch levels

ctrl <- fwdControl(lapply(2011:2020, function(x)
  list(year=x, quant="catch", value=seq(50, 50000, length=10))))

tes3 <- fwdabc.om(iter(om, its), ctrl, pcbar=pcbar, pla=pla)

ggplot(unitSums(ssb(biol(tes1[[1]]))), aes(x=date, y=data, group=iter)) +
  geom_line(aes(colour=iter)) + facet_wrap(~as.factor(iter))

# - TEST mp in past

method(oem) <- sampling.oem

projection(om) <- mseCtrl(method=fwdabc.om, args=list(pcbar=pcbar, pla=pla))

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=perfect.sa),
  # HCR
  hcr = mseCtrl(method=fixedC.hcr, args=list(ctrg=rep(35000, 10)))))

tes <- mp(om, oem=oem, ctrl=ctrl, args=list(iy=2015, fy=2020))

ctrl <- mpCtrl(list(
  # EST
  est = mseCtrl(method=cpue.ind),
  # HCR
  hcr = mseCtrl(method=cpue.hcr, args=list(target=2))))

tes <- mp(om, oem=oem, ctrl=ctrl, args=list(iy=2015, fy=2020))

plot(om, tes)


# TODO: verify(FLombf), verify(FLBiol), ...

# CHECK harvest
attr(om(tes), 'harvest')
harvest(om(tes))

tcatch <- function(x) unitSums(seasonSums(catch(x)[[1]]))
plot(window(tcatch(om), end=2015), tcatch(om(tes))) + ylim(0, NA)

# - TEST fwd in future

omf <- fwdWindow(om, end=2040)
attr(omf, 'harvest') <- window(attr(om, 'harvest'), end=2040)

ctrl <- fwdControl(year=2021:2040, quant="catch", value=rep(0, 200))
ctrl <- fwdControl(year=2021:2040, quant="catch", value=rep(0))

tes3 <- fwdabc.om(omf, ctrl, pcbar=pcbar, pla=pla)

plot(om, tes3)

unitSums(seasonSums(quantSums(n(biol(tes3$om)) * wt(biol(tes3$om)))))
unitSums(seasonSums(quantSums(n(biol(tes3$om))[1,])))

# HR metric
# hmsy

# H / hmsy


