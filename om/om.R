# om.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/abc_tuna/om/om.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# TODO: window(stk), keep full in observations(oem)

# REFERENCE
#
# mcpars: c(log(R0),logit(dep),epsr,log(as.vector(selpars)))
#
# mcvars: "N"  "Rtot" "SSB" "dep" "dbmsy" "Cmsy" "hmsyrat" "H" "Ihat" "LFhat"


library(mse)
library(ss3om)

library(patchwork)

# LOAD SS3 FLStock

stk <- readFLSss3('../data/base')

# ABC dimnames
dmns <- dimnames(stk)
dmns$year <- dmns$year[seq(71 - 20, 71)]


# --- LOAD abc5

run5 <- mget(load("../v2/alb_abc_run5.rda", verbose=FALSE,
  envir=(NE. <- new.env())), envir=NE.)

stk5 <- propagate(stk, run5$nits)

# SET mat
mat(stk5)[,,,4] <- mat(stk5)[,,,1]
mat(stk5)[,,,1] <- 0

# LOAD stock.n

nabc <- Reduce(combine, lapply(run5$mcvars, function(x)
  FLQuant(aperm(x$N, c(2,1,4,3)), dimnames=dmns) / 1000
))

stock.n(stk5)[, dmns$year] <- nabc

stock(stk5) <- computeStock(stk5)

# TODO: CHECK recruitment timing

# LOAD catch.n

# LOAD refpts: FMSY, BMSY, B0

# LOAD SRR: B0, R0, h, spr0

# double B0 = R0*spr0;
# double alp = 4.*hh/(spr0*(1.-hh));
# double bet = (5.*hh-1.)/(B0*(1.-hh));
# Rtot = (alp*S[ysp][spwn]/(1.+bet*S[ysp][spwn]))*exp(epsr(y-1)); 

r0 <- exp(run5$mcpars[,1])

# LOAD other results

dep <- run5$mcpars[,2]
sigmar <- exp(run5$mcpars[,3]) / (1 + exp(run5$mcpars[,3]))


# TODO: COMPARE ssbs

nssb <- Reduce(combine, lapply(run5$mcvars, function(x)
  FLQuant(x$SSB, dimnames=list(age='all', year=dmns$year, unit='F', season=4))
))


# --- PLOT

pdf(file="abc5_ss.pdf")

plot(append(window(stock.n(stk), end=1999), nabc)) +
  geom_vline(xintercept=ISOdate(2000, 1, 1)) +
  ggtitle("SS3 + ABC Ns")

plot(stock.n(stk)) + ggtitle("SS3 Ns")

plot(stock.n(stk)['0',, 'F', 4], stock.n(stk5)['0',, 'F', 4]) +
  ggtitle("SS3 + ABC age 0")
plot(stock.n(stk)['1',, 'F', 1], stock.n(stk5)['1',, 'F', 1]) +
  ggtitle("SS3 + ABC age 1")

plot(ssb(stk5)[, dmns$year, 'F', 4], nssb, ssb(stk)[, dmns$year, 'F', 1]) +
  ggtitle("SS3 + ABC SSB") + ylim(c(0,NA))

dev.off()
