###########################################################
# v2 projection code with MCMC iters embedded in C++ ######
###########################################################
# R. Hillary & I. Mosqueira (2025) ########################
###########################################################

library(Rcpp)
library(FLCore)
library(ggplotFL)
library(parallel)
library(mvtnorm)
source("utilities.R")

load("alb_abc_run5a.rda")

##############################
# all-in-one projection code #
##############################

sourceCpp("projv2.cpp")

# projection control file

prj.ctrl <- list(TAC=30000,nyproj=10)
nyprj <- prj.ctrl$nyproj

nitsx <- 500
dms2 <- c(nyprj,ns,na,nf,nitsx)

yref <- 2017:2020
cref <- C[as.character(yref),,]
cbar <- apply(cref,c(2,3),mean)
pcbar <- cbar/sum(cbar)
Cproj <- array(dim=c(nyprj,ns,nf,nitsx))
for(y in 1:nyprj) Cproj[y,,,] <- prj.ctrl$TAC*pcbar

# create arrays needed

sela.mc <- array(dim=c(na,ns,2,nf,nitsx))
Ninit <- array(dim=c(na,ns,2,nitsx))
R0.mc <- B0.mc <- hh.mc <- M.mc <- spr0.mc <- sigmar.mc <- rep(NA,nitsx)
q.mc <- array(dim=c(ns,nitsx))
dep.mc <- array(dim=c(ny,nitsx))

# extract MCMC variables - only have to do this once!

for(n in 1:nitsx) {

  tmp <- mcvars[[n]]
  parv <- mcpars[n,] 
  selvx <- exp(parv[(ny+2):npar])
  selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
  sela.mc[,,,,n] <- get.sel.age(6,5,selidx,selparsx)
  hh.mc[n] <- parv[npar+1]
  M.mc[n] <- parv[npar+2]
  sigmar.mc[n] <- parv[npar+3]
  R0.mc[n] <- exp(parv[1])
  B0.mc[n] <-  tmp$B0
  spr0.mc[n] <- B0.mc[n]/R0.mc[n]
  Ninit[,,,n] <- mcvars[[n]]$N[ny,,,]
  dep.mc[,n] <- tmp$dep
  resq <- log(I[,,fcpue]/mcvars[[n]]$Ihat)
  if(seasonq) q.mc[,n] <- exp(apply(resq,2,mean))
  if(!seasonq) q.mc[n] <- mean(resq)

  if(n %% 50 == 0) cat("Iteration ",n,"of",nitsx,"\n")

}

system.time(resv2 <- projv2(dms2,srec,R0.mc,hh.mc,psi,sigmar.mc,spr0.mc,M.mc,as.vector(mata),as.vector(wta),as.vector(sela.mc),as.vector(Ninit),as.vector(Cproj),as.vector(q.mc),fcpue))

# reconstruction of projection variables

Nhat <- array(resv2$N,dim=c(nyprj,na,ns,2,nitsx))
Shat <- array(resv2$S,dim=c(nyprj,ns,nitsx))
Ihat <- array(resv2$I,dim=c(nyprj,ns,nitsx))
Hhat <- array(resv2$H,dim=c(nyprj,ns,nf,nitsx))

# join together historical and projected

depprj.mc <- t(apply(Shat[,srec-1,],1,function(x,B0.mc){x <- x/B0.mc},B0.mc))
yrs.prj <- (yrs[1]):(yrs[length(yrs)]+nyprj-1)
depall.mc <- rbind(dep.mc,depprj.mc[-1,])
dq <- apply(depall.mc,1,quantile,c(0.025,0.5,0.975))
plot(yrs.prj,dq[2,],ylim=c(0,max(dq[3,])),type='l',col='blue',lwd=1.5)
lines(yrs.prj,dq[1,],lty=2,lwd=2,col='blue')
lines(yrs.prj,dq[3,],lty=2,lwd=2,col='blue')

