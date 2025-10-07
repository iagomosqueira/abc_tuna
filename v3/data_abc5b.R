# abc5.R - DESC
# /home/mosquia/Active/Doing/ABC_tuna+iotc/abc_tuna/v2/abc.R
# resample (h,M) from joint distribution
# resample sigmaR using informative conjugate prior

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# data

#' @param laa Length at ageby sex
#' @param wal Weight at length by sex
#' @param mal Maturity at length for females
#' @param cpues CPUE indices by fishery, area [fishery, year]
#' @param lenfreq Length-frequency in catch by fishery [len, year, fishery]
#' @param iALK, VB

# arguments

#' @param priors List for h, K
#' @param pvars Proposal variances
#' @param Fstatus F/FMSY prior

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

# NC by fleet [y, s, f]
# 

load("data/alb_abcdata.rda")
load("data/hmuprior.rda")
load("data/base.rda")

# fnscale 
fnscale <- 1

# years and values for stock status priors

yof <- 1:ny # apply over-fishing prior penalty
sdof <- 0.5 # idea P(hmsyrat > 2) <= 0.05 essentially
ybmsy <- c(ny-1,ny)
mubmsy <- c(2.25,2)
sdbmsy <- c(0.35,0.35)
ydep <- 1
mudep <- 0.5
sddep <- 0.1

# rescale weight to tonnes

wta <- wta[]*1e-3
 
pobs <- apply(LFfits,2,function(x){x <- x/sum(x)}) # what we will fit to

# --- simulator

# - arguments
#   - parameters: R0, dep, h, epsr, selpars
#   - biology
#   - fishery

R0 <- 14e6
dep <- 0.5
h <- 0.8

# number of fisheries

nf <- dim(C)[3]

# number of distinct selectivity groups

nselg <- 5

# selectivity for each fishery

selidx <- c(1,2,3,4,5,5)

# fisheries with LF data

flf <- c(1,2,3,4,6)
nflf <- length(flf)
pobs <- apply(LFfits,2,function(x){x <- x/sum(x)}) # what we will fit to
pobs <- pobs[,flf]

# catch distro targets

pctarg <- C[1,,] / sum(C[1,,])

# set up selectivity parameters (all double normal)

smax <- c(120,125,85,85,115)
sL <- c(20,20,7,7,10)
sR <- c(35,30,25,15,30)
selpars <- cbind(smax,sL,sR)

# recruitment variations (ny-1)

epsr <- rep(0,ny-1)

# informative priors (roughly between 0.2-0.5 mean of 0.3)

musigma2R <- 0.3^2
cvx <- 0.4
mux <- 1/musigma2R
vx <- (mux*cvx)^2
alpR <- mux^2/vx
betR <- mux/vx
psigmaR <- sqrt(1/rgamma(1000,alpR,betR))
round(quantile(psigmaR,c(0.025,0.5,0.975)),2)

# recruitment season

srec <- 4

# sex ratio at birth (fiddy:fiddy)

psi <- 0.5

# dimensions

dms <- c(na,ns,nf,nselg)

## MCMC algorithm 

# set up MCMC controls for unconditional sampling of (h,M)

acphmu <- 0.25 # force acceptance rate at "optimal" MCMC value

# Gibbs sampling parameter groupings
# 1. B0 and dep
# 2. recruitment deviates
# 3. selectivity 

npar <- 2+ny-1+3*nselg
ngibbs <- 3
paridx <- list()
paridx[[1]] <- 1:2
paridx[[2]] <- 3:(ny+1)
paridx[[3]] <- (ny+2):npar
lidx <- unlist(lapply(paridx,length))

# SD in CPUE index

fcpue <- 1
scpue <- 1:4
sd.cpue <- rep(NA,length(scpue))

sdcpue <- mean(sd.cpue)

fcpue <- 1
scpue <- 1:4
sd.cpue <- rep(NA,length(scpue))
for(s in scpue) {
  idf <- data.frame(t=yrs,y=log(I[,s,fcpue]))
  ires <- loess(y~t,idf)
  sd.cpue[s] <- sd(residuals(ires))
}

sdcpue <- mean(sd.cpue)

# KLmax

KLmax <- 0.8 # consistent with minimum Neff = 20 multinomial

# seasonal q for CPUE (T or F)

seasonq <- TRUE

# catchability trend in q

qtrend <- FALSE

# burn-in and thinning factor
 
burn <- 10
thin <- 1

###################
# run the sampler #
###################

# set up initial h and M

hold <- hmu
Mold <- Mmu
sigmarold <- sqrt(musigma2R)

# set up initial guess parameter vector

parvecold <- c(log(R0),logit(dep),epsr,log(as.vector(selpars)))

# RW variance by Gibbs grouping

rwsd <- rep(0,npar)
rwsd[paridx[[1]]] <- c(0.1,0.05)
rwsd[paridx[[2]]] <- 0.08
rwsd[paridx[[3]]] <- 0.025

# TEST
nits1 <- 10 # total number of retained samples
system.time(zzz <- mcmc5.abc(nits1))
zzz$acp/nits1

# parallelised efficient version

parvecold <- zzz$pars[nits1,1:npar]
hold <- zzz$pars[nits1,npar+1]
Mold <- zzz$pars[nits1,npar+2]
sigmarold <- zzz$pars[nits1,npar+3]
nits <- 500
ncore <- 5
thin <- 100
mcnits <- floor(nits/ncore)

# SAVE image
save.image(file="data/image/abc5b.rda", compress="xz")

# RUN
system.time(mczzz <- mclapply(rep(mcnits,ncore),mcmc5.abc,mc.cores=ncore))

# EXTRACT
mcpars <- do.call(rbind, lapply(mczzz, '[[', 'pars'))
mcvars <- get.mcmc2.vars(mcpars)

# SAVE
save(mczzz, file="data/mcmc/abc5b.rda", compress="xz")
save(mcvars, C, file="data/mcvars/abc5b.rda", compress="xz")

# --- CREATE om {{{

load('data/mcmc/abc5b.rda')
load('data/mcvars/abc5b.rda')

# EXTRACT output for all iters
out <- mc.output(mcvars, C)

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

# DEVIANCES
deviances(bio)[,,,4] <- residuals(
  rec(bio)[,,,4],
  expand(predict(sr(bio), ssb=out$ssb) / 2, unit=c('F', 'M'))
)

# - FLFisheries (catch.n, catch.sel)

# FLCatch(es)
cas <- Map(function(x, y) FLCatch(landings.n=x, landings.wt=wt(bio),
  catch.sel=y, discards.n=x %=% 0, discards.wt=wt(bio)),
  x=divide(out$catch.n, 5), y=divide(out$catch.sel %*% (out$catch.n %=% 1), 5))

# FLFisheries
fis <- FLFisheries(lapply(cas, function(x)
  FLFishery(effort=unitSums(catch(cas[[1]])) %=% 0, ALB=x)))

names(fis) <- c(paste0("LL", 1:4), "PS", "Other")

om <- FLombf(biols=FLBiols(ALB=bio), fisheries=fis,
  refpts=FLPars(ALB=out$refpts))


# TEST: total catch
unitSums(seasonSums(areaSums(Reduce('+', catch(fisheries(om))))))


# TEST: fwdabc.om hindcast


# - BUILD oem

# idx: FLIndexBiomass by season, with sel.pattern by sex

# BUG: mc.output to output FLQuants by fiosheru, not 'area'
sp <- expand(divide(out$catch.sel, 5)[[1]], year=2000:2020)
dimnames(sp)$area <- 'unique'

NW <- FLIndexBiomass(index=out$index.hat %*% out$index.q,
  index.q=expand(out$index.q, year=2000:2020),
  sel.pattern=sp,
  catch.wt=wt(biol(om)),
  range=c(startf=0.5, endf=0.5))

# stk: no units
oem <- FLoem(observations=list(ALB=list(idx=FLIndices(NW=NW),
  stk=simplify(stock(om)[[1]], 'unit'))), method=sampling.oem)

survey(observations(oem)$ALB$stk, observations(oem)$ALB$idx)

# TODO: verify(oem, om)

# SAVE
save(om, oem, file='data/om5b.rda', compress='xz')

# }}}








