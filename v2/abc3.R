# abc3.R - DESC
# /home/mosquia/Active/Doing/ABC_tuna+iotc/abc_tuna/v2/abc.R
# INITIAL DEPLETION NOW A DISTRIBUTION  
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
source("utilities.R")

sourceCpp("init_pdyn.cpp")
sourceCpp("msy_pdyn.cpp")
#sourceCpp("pdyn.cpp")
sourceCpp("pdyn_lfcpue.cpp")

# NC by fleet [y, s, f]
# 

#load('../data/data.RData')
load("alb_abcdata.rda")

# fnscale 
fnscale <- 1

# years and values for stock status priors

#yfmsy <- c(1,2,ny-1,ny) # first two, last two 
#mufmsy <- c(0.6,0.55,0.57,0.6)
#sdfmsy <- c(0.125,0.125,0.175,0.18)
ybmsy <- c(ny-1,ny)
mubmsy <- c(2.25,2)
sdbmsy <- c(0.35,0.35)
ydep <- 1
mudep <- 0.5
sddep <- 0.05
# solve for alpdep and betdep of beta distribution

fn <- function(x) {

  vdep <- sddep^2
  expr <- (x+x*(1-mudep)/mudep)^2*(x+x*(1-mudep)/mudep+1)*vdep-x^2*(1-mudep)/mudep

  return(expr)

}

alpdep <- uniroot(fn,interval=c(40,70))$root
betdep <- alpdep*(1-mudep)/mudep

# recale weight to tonnes

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

# recruitment season

srec <- 4

# sex ratio at birth (fiddy:fiddy)

psi <- 0.5

# dimensions

dms <- c(na,ns,nf,nselg)

## MCMC algorithm 

# Gibbs sampling parameter groupings
# 1. R0 
# 2. recruitment deviates
# 3. selectivity 

npar <- 1+ny-1+3*nselg
ngibbs <- 3
paridx <- list()
paridx[[1]] <- 1
paridx[[2]] <- 2:(ny)
paridx[[3]] <- (ny+1):npar
lidx <- unlist(lapply(paridx,length))

# SD in CPUE index

fcpue <- 1
scpue <- 1:4
sd.cpue <- rep(NA,length(scpue))
par(mfrow=c(2,2))
for(s in scpue) {

  idf <- data.frame(t=yrs,y=log(I[,s,fcpue]))
  ires <- loess(y~t,idf)
  plot(yrs,ires$fitted,type='l')
  points(yrs,log(I[,s,fcpue]))
  sd.cpue[s] <- sd(residuals(ires))

}

sdcpue <- mean(sd.cpue)

# sigmaR

sigmar <- 0.3 # from assessment

# KLmax

KLmax <- 0.8 # consistent with minimum Neff = 20 multinomial

# seasonal q for CPUE (T or F)

seasonq <- FALSE

# burn-in and thinning factor
 
burn <- 10
thin <- 1

###################
# run the sampler #
###################

# initial depletion

depx <- 0.5

# set up initial guess parameter vector

parvecold <- c(log(R0),epsr,log(as.vector(selpars)))

# RW variance by Gibbs grouping

rwsd <- rep(0,npar)
rwsd[paridx[[1]]] <- 0.01
rwsd[paridx[[2]]] <- 0.02
rwsd[paridx[[3]]] <- 0.005

nits1 <- 10 # total number of retained samples
system.time(zzz <- mcmcb.abc(nits1))
zzz$acp/nits1
boxplot(zzz$pars,outline=F,col='magenta')

# parallelised efficient version

depx <- zzz$dep[nits1]
parvecold <- zzz$pars[nits1,]
nits <- 500
ncore <- 10
thin <- 50
mcnits <- floor(nits/ncore)
system.time(mczzz <- mclapply(rep(mcnits,ncore),mcmcb.abc,mc.cores=ncore))

mcacp <- apply(matrix(unlist(lapply(mczzz,function(x){x <- x$acp})),ncol=ngibbs,byrow=T),2,sum)/(nits*thin)
mcacp
mcpars <- mczzz[[1]]$pars
mcdep <- mczzz$dep
for(i in 2:ncore) mcpars <- rbind(mcpars,mczzz[[i]]$pars)
boxplot(mcpars,outline=F,col='magenta') 

mcvars <- get.mcmc.vars(mcpars)

plot.mcmc.vars(mcvars,'dep')
plot.mcmc.vars(mcvars,'bmsy')
plot.mcmc.vars(mcvars,'rec')
plot.mcmc.vars(mcvars,'cpue')
plot.mcmc.vars(mcvars,'lf')
plot.mcmc.sel(mcpars)

save.image("alb_abc_run3.rda")


