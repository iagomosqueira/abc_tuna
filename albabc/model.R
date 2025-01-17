# model.R - DESC
# abc_tuna/albabc/model.R

# Copyright (c) WUR & CSIRO, 2023.
# Authors: Richard HILLARY (CSIRO) <rich.hillary@csiro.au>
#          Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(Rcpp)
source("utilities.R")
library(parallel)
library(mvtnorm)

# SOURCE Cpp

sourceCpp("utilities/init_pdyn.cpp")
sourceCpp("utilities/msy_pdyn.cpp")
sourceCpp("utilities/pdyn_lfcpue.cpp")

# LOAD data
load("data/alb_abcdata.rda")
load("data/hmuprior.rda")


# --- SETUP

# 
fnscale <- 1
#
ybmsy <- c(ny-1, ny)
#
mubmsy <- c(2.25, 2)
#
sdbmsy <- c(0.35, 0.35)
#
ydep <- 1
#
mudep <- 0.5
#
sddep <- 0.1

# rescale weight to tonnes
wta <- wta[]*1e-3
 
#
pobs <- apply(LFfits,2,function(x){x <- x/sum(x)}) # what we will fit to

# - parameters: R0, dep, h, epsr, selpars
R0 <- 14e6
dep <- 0.5
h <- 0.8

# dims
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

# burn-in and thinning factor
burn <- 10
thin <- 1

# run the sampler #

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

nits1 <- 10 # total number of retained samples

system.time(zzz <- mcmc3.abc(nits1))


# parallelised efficient version

parvecold <- zzz$pars[nits1,1:npar]
hold <- zzz$pars[nits1,npar+1]
Mold <- zzz$pars[nits1,npar+2]
sigmarold <- zzz$pars[nits1,npar+3]
nits <- 200
ncore <- 5
thin <- 10
mcnits <- floor(nits/ncore)
system.time(mczzz <- mclapply(rep(mcnits,ncore),mcmc3.abc,mc.cores=ncore))

mcacp <- apply(matrix(unlist(lapply(mczzz,function(x){x <- x$acp})),ncol=ngibbs,byrow=T),2,sum)/(nits*thin)
mcacp
mcpars <- mczzz[[1]]$pars
for(i in 2:ncore) mcpars <- rbind(mcpars,mczzz[[i]]$pars)
boxplot(mcpars,outline=F,col='magenta') 

mcvars <- get.mcmc2.vars(mcpars)

plot.mcmc.vars(mcvars,'dep')
plot.mcmc.vars(mcvars,'bmsy')
plot.mcmc.vars(mcvars,'rec')
plot.mcmc.cpue(mcvars)
plot.mcmc.lf(mcvars)
plot.mcmc.sel(mcpars)

# prior vs posterior sigmaR

sigrpost <- mcpars[,npar+3]
dpost <- density(sigrpost)
dprior <- density(psigmaR)
p1 <- dpost$y
p1 <- p1/sum(p1)
x1 <- dpost$x
p2 <- dprior$y
p2 <- p2/sum(p2)
x2 <- dprior$x
pmax <- max(c(max(p1),max(p2)))
xmin <- min(c(min(x1),min(x2)))
xmax <- max(c(max(x1),max(x2)))
plot(x1,p1,type='l',xlim=c(xmin,xmax),ylim=c(0,pmax),xlab=expression(sigma[R]),ylab='density')
lines(x2,p2,lty=2,col='purple')
legend("topright",lty=c(1,2),col=c("black","purple"),legend=c("Posterior","Prior"),bty='n')

save.image("alb_abc_run5.rda")

