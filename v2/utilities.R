# utilities.R - DESC
# /home/mosquia/Active/ABC_tuna+iotc/abc_tuna/v2/utilities.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


logit <- function(x){
  return(log(x/(1-x)))
}

ilogit <- function(x){
  return(1/(1+exp(-x)))
}

get.sel.age <- function(nf=6,nselg=5,selidx=c(1,2,3,4,5,5),selpars)
{

  seltmp <- rep(NA,20)
  sela <- array(dim=c(na,ns,2,nf))
  for(g in 1:2) {
    for(s in 1:4) {
      for(a in 1:na) {

        lmin <- max(0,mula[a,s,g]-sdla[a,s,g]*1.96)
        lmax <- mula[a,s,g]+sdla[a,s,g]*1.96
        lref <- seq(lmin,lmax,length=20)
        dl <- dnorm(lref,mula[a,s,g],sdla[a,s,g])
        dl <- dl/sum(dl)
        for(f in 1:nf) {
        
          fref <- selidx[f]
          for(l in 1:20) seltmp[l] <- ifelse(lref[l] < selpars[fref,1],2^{-((lref[l]-selpars[fref,1])/selpars[fref,2])^2},2^{-((lref[l]-selpars[fref,1])/selpars[fref,3])^2})

          sela[a,s,g,f] <- sum(seltmp*dl)
        }
      }
    }
  }

  return(sela)
}

objfn.init <- function(theta,targv,sela) {

  hxinit <- 1 / (1 + exp(-theta))
  resx <- initpdyn(c(ns, na, nf), srec, psi, M, as.vector(mata),
  as.vector(wta), as.vector(sela), hxinit) 

  cx <- resx$C
  px <- cx/sum(cx)
  tmpv <- c(resx$rho, as.vector(px))
  objv <- logit(tmpv)

  return(fnscale * (sum((objv - targv) ^ 2)))

}

robjfn.init <- function(theta,targv,sela) {

  hxinit <- 1 / (1 + exp(-theta))
  resx <- rinitpdyn(c(ns, na, nf),srec,psi,M,mata,wta,sela,hxinit) 

  cx <- resx$C
  px <- cx/sum(cx)
  tmpv <- c(resx$rho, as.vector(px))
  objv <- logit(tmpv)

  return(fnscale * (sum((objv - targv) ^ 2)))
 
}

msyfn <- function(H,ph,sela) {

  hx <- H * ph
  resx <- msypdyn(c(ns,na,nf),srec,R0,h,psi,M,as.vector(mata),as.vector(wta),as.vector(sela),hx)
    
  return(sum(resx$C))
}

rmsyfn <- function(H,ph,sela) {

  hx <- H * ph
  resx <- rmsypdyn(c(ns,na,nf),srec,R0,h,psi,M,mata,wta,sela,hx)
    
  return(sum(resx$C))
}

# vanilla R initial popn dynamics

# rinitpdyn {{{

rinitpdyn <- function(dms,srec,psi,M,mata,wta,sela,hstart) {

  spwn <- ifelse(srec == 1,ns,srec-1)

  ns <- dms[1]
  na <- dms[2]
  nf <- dms[3]
  H <- array(hstart,dim=c(ns,nf))

  # unfished level

  neqm <- array(dim=c(na,ns,2))
  ceqm <- array(dim=c(ns,nf))
  heqm <- array(dim=c(na,ns,2,nf))

  for(s in 1:ns) {

    if(s < srec) neqm[1,s,] <- 0
    if(s == srec) neqm[1,s,] <- c(psi,1-psi)
    if(s > srec) neqm[1,s,] <- neqm[1,s-1,]*exp(-M)
    
  }

  for(a in 2:na) {
    for(s in 1:ns) {

      if(s == 1) neqm[a,1,] <- neqm[a-1,ns,]*exp(-M)
      if(s > 1) neqm[a,s,] <- neqm[a,s-1,]*exp(-M)

    }
  }

  spr0 <- sum(neqm[,spwn,1]*mata[,spwn,1]*wta[,spwn,1])

  # exploited eqm

  for(s in 1:ns) {
    for(g in 1:2) {
      for(f in 1:nf) {
        
        heqm[,s,g,f] <- H[s,f]*sela[,s,g,f]

      }
    }
  }

  hsum <- apply(heqm,1:3,sum)

  for(s in 1:ns) {

    if(s < srec) neqm[1,s,] <- 0
    if(s == srec) neqm[1,s,] <- c(psi,1-psi)
    if(s > srec) {

      for(g in 1:2) neqm[1,s,g] <- neqm[1,s-1,g]*exp(-M)*(1-hsum[1,s-1,g])

    }
  }

  for(a in 2:na) {
    for(s in 1:ns) {
      for(g in 1:2) {
      
        htmp <- min(0.9,hsum[a,s,g])
        if(s == 1) neqm[a,1,g] <- neqm[a-1,ns,g]*exp(-M)*(1-htmp)
        if(s > 1) neqm[a,s,g] <- neqm[a,s-1,g]*exp(-M)*(1-htmp)

      }
    }
  }
 
  sprf <- sum(neqm[,spwn,1]*mata[,spwn,1]*wta[,spwn,1]) 
    
  # catch biomass by fleet

  for(s in 1:ns) {
    for(f in 1:nf) {

      ceqm[s,f] <- 0
      for(g in 1:2) ceqm[s,f] <- ceqm[s,f]+sum(neqm[,s,g]*wta[,s,g]*sela[,s,g,f]*H[s,f])

    }
  }

  # SPR ratio

  rho <- sprf/spr0
  return(list(rho=rho,C=ceqm,N=neqm))
}
# }}}

# {{{ rmsypdyn
rmsypdyn <- function(dms,srec,R0,h,psi,M,mata,wta,sela,hx) {

  spwn <- ifelse(srec == 1,ns,srec-1)

  ns <- dms[1]
  na <- dms[2]
  nf <- dms[3]
  H <- array(hx,dim=c(ns,nf))

  # unfished level

  neqm <- array(dim=c(na,ns,2))
  ceqm <- array(dim=c(ns,nf))
  heqm <- array(dim=c(na,ns,2,nf))

  for(s in 1:ns) {

    if(s < srec) neqm[1,s,] <- 0
    if(s == srec) neqm[1,s,] <- c(psi,1-psi)
    if(s > srec) neqm[1,s,] <- neqm[1,s-1,]*exp(-M)
    
  }

  for(a in 2:na) {
    for(s in 1:ns) {

      if(s == 1) neqm[a,1,] <- neqm[a-1,ns,]*exp(-M)
      if(s > 1) neqm[a,s,] <- neqm[a,s-1,]*exp(-M)

    }
  }

  spr0 <- sum(neqm[,spwn,1]*mata[,spwn,1]*wta[,spwn,1])

  # exploited eqm

  for(s in 1:ns) {
    for(g in 1:2) {
      for(f in 1:nf) {
        
        heqm[,s,g,f] <- H[s,f]*sela[,s,g,f]

      }
    }
  }

  hsum <- apply(heqm,1:3,sum)

  for(s in 1:ns) {

    if(s < srec) neqm[1,s,] <- 0
    if(s == srec) neqm[1,s,] <- c(psi,1-psi)
    if(s > srec) {

      for(g in 1:2) neqm[1,s,g] <- neqm[1,s-1,g]*exp(-M)*(1-hsum[1,s-1,g])

    }
  }

  for(a in 2:na) {
    for(s in 1:ns) {
      for(g in 1:2) {
      
        htmp <- min(0.9,hsum[a,s,g])
        if(s == 1) neqm[a,1,g] <- neqm[a-1,ns,g]*exp(-M)*(1-htmp)
        if(s > 1) neqm[a,s,g] <- neqm[a,s-1,g]*exp(-M)*(1-htmp)

      }
    }
  }
 
  sprf <- sum(neqm[,spwn,1]*mata[,spwn,1]*wta[,spwn,1]) 

  # S-R bit

  rho <- sprf/spr0
  B0 <- R0*spr0
  alp <- 4*h/(spr0*(1-h))
  bet <- (5*h-1)/(B0*(1-h))
  Bbar <- ifelse((alp*sprf-1.)/bet > 0,(alp*sprf-1)/bet,0)
  Rbar <- Bbar/sprf

  # scale population by Rbar

  neqm[] <- neqm[]*Rbar

  # catch biomass by fleet

  for(s in 1:ns) {
    for(f in 1:nf) {

      ceqm[s,f] <- 0
      for(g in 1:2) ceqm[s,f] <- ceqm[s,f]+sum(neqm[,s,g]*wta[,s,g]*sela[,s,g,f]*H[s,f])

    }
  } 
 
  return(list(rho=rho,C=ceqm,N=neqm,spr0=spr0,Bmsy=Bbar,Rmsy=Rbar))

}
# }}}

# rpdynlfcpue {{{
rpdynlfcpue <- function(dms,srec,R0,h,psi,epsr,spr0,M,mata,wta,sela,Ninit,C,pla,fcpue) {

  ny <- dms[1]
  ns <- dms[2]
  na <- dms[3]
  nl <- dms[4]
  nf <- dms[5]

  spwn <- ifelse(srec == 1,ns,srec-1)

  B0 <- R0*spr0
  alp <- 4*h/(spr0*(1-h))
  bet <- (5*h-1)/(B0*(1-h)) 

  N <- array(dim=c(ny,na,ns,2))
  H <- array(dim=c(ny,ns,nf))
  S <- I <- array(dim=c(ny,ns))
  LF <- array(dim=c(ny,nl,ns,nf))

  # initial population

  N[1,,,] <- Ninit
  for(s in 1:ns) S[1,s] <- sum(N[1,,s,1]*mata[,s,1]*wta[,s,1])
  for(f in 1:nf) {
    for(s in 1:ns) {

      xsum <- 0
      for(g in 1:2) xsum <- xsum+sum(N[1,,s,g]*sela[,s,g,f]*wta[,s,g])
      H[1,s,f] <- min(C[1,s,f]/xsum,0.9)
      if(f == fcpue) I[1,s] <- xsum

      # LF data
 
      for(l in 1:nl) {

        LF[1,l,s,f] <- 0 
        for(g in 1:2) LF[1,l,s,f] <- LF[1,l,s,f]+sum(N[1,,s,g]*sela[,s,g,f]*pla[l,,s,g])
      }
    }
  }

  # annual loop

  for(y in 2:ny) {

    ysp <- ifelse(srec == 1,y-1,y)

    for(s in 1:ns) {

      if(s < srec) N[y,1,s,] <- 0
      if(s == srec) {

        Rtot = (alp*S[ysp,spwn]/(1+bet*S[ysp,spwn]))*exp(epsr[y-1])
        N[y,1,s,] <- Rtot*c(psi,1-psi)

      }

      if(s > srec) {

        for(g in 1:2) {

          hsum <- min(sum(H[y,s-1,]*sela[1,s-1,g,]),0.9)
          N[y,1,s,g] <- N[y ,1,s-1,g]*exp(-M)*(1-hsum)

        }
      }

      # loop through ages

      for(a in 2:na) {
        for(g in 1:2) {

          if(s == 1) {

            hsum <- min(sum(H[y-1,s-1,]*sela[a-1,ns,g,]),0.9)
            N[y,a,s,g] <- N[y-1,a-1,ns,g]*exp(-M)*(1-hsum) 
 
          } else {
  
            hsum <- min(sum(H[y,s-1 ,]*sela[a,s-1,g,]),0.9) 
            N[y,a,s,g] <- N[y,a,s-1,g]*exp(-M)*(1-hsum) 

          }
        }
      }

      # female SSB
  
      S[y,s] <- sum(N[y,,s,1]*mata[,s,1]*wta[,s,1]) 

      # harvest rates, CPUE and LF 
  
      for(f in 1:nf) {

        xsum <- 0
        for(g in 1:2) xsum <- xsum+sum(N[y,,s,g]*sela[,s,g,f]*wta[,s,g])
        H[y,s,f] <- min(C[y,s,f]/xsum,0.9)
        if(f == fcpue) I[y,s] <- xsum

        # LF data
 
        for(l in 1:nl) {

          LF[y,l,s,f] <- 0 
          for(g in 1:2) LF[y,l,s,f] <- LF[y,l,s,f]+sum(N[y,,s,g]*sela[,s,g,f]*pla[l,,s,g])
        } 
      }
    }
  }

  return(list(S=S,N=N,H=H,I=I,LF=LF))
}
# }}}

# rmcmc.abc {{{
rmcmc.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar)
  acp <- rep(0,ngibbs)

  # get initial guess discrepancy

  xx <- rsim(R0,dep,h,M,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- apply(resq,2,mean)
    resq <- t(apply(resq,1,function(x,lnq){x <- x-lnq},lnq))

  } else {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))
  }  

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmar,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    for(gg in 1:ngibbs) {

      epsrw <- rnorm(lidx[gg],0,rwsd[paridx[[gg]]])
      parvecnew <- parvecold
      parvecnew[paridx[[gg]]] <- parvecnew[paridx[[gg]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- rsim(R0x,depx,h,M,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- apply(resq,2,mean)
        resq <- resq-lnq

      } else {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 

      # parameter priors

      pprior <- sum(dnorm(epsrx,0,sigmar,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[gg] <- acp[gg]+1

      }
    }

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- parvecold

  }

  return(list(pars=theta.mcmc,acp=acp))
} 
# }}}

# mcmc.abc {{{
mcmc.abc <- function(nits) {

  theta.mcmc <- matrix(nrow=nits,ncol=npar)
  acp <- rep(0,ngibbs)

  # get initial guess discrepancy

  xx <- sim(R0,dep,h,M,selpars,epsr,dms,pctarg,selidx) 

  # LF discrepancy

  phat <- xx$LF
  kllf <- pobs*log(pobs/phat)
  dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

  # CPUE discrepancy

  if(seasonq) {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- apply(resq,2,mean)
    resq <- resq-lnq

  } else {

    resq <- log(I[,,fcpue]/xx$I)
    lnq <- mean(resq)
    resq <- resq-lnq

  }

  dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

  ## priors (parameters + stock status)

  # status priors

  Bmsyrat <- xx$S[,3]/xx$Bmsy
  Bmsyrat <- Bmsyrat[ybmsy]
  dSSB <- xx$S[,srec-1]/xx$B0
  dSSB <- dSSB[ydep]

  if(length(ybmsy) == 1) {

    sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

  } else {

    sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

  }

  if(length(ydep) == 1) {

    sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

  } else {

    sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

  }  

  # parameter priors

  pprior <- sum(dnorm(epsr,0,sigmar,TRUE))

  # starting discrepancy

  dtotold <- dcpue+sprior+pprior-dlf

  for(n in 1:(burn+thin*nits)) {

    for(g in 1:ngibbs) {

      epsrw <- rnorm(lidx[g],0,rwsd[paridx[[g]]])
      parvecnew <- parvecold
      parvecnew[paridx[[g]]] <- parvecnew[paridx[[g]]]+epsrw
      R0x <- exp(parvecnew[1])
      depx <- ilogit(parvecnew[2])
      epsrx <- parvecnew[3:(ny+1)]
      selvx <- exp(parvecnew[(ny+2):npar])
      selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
      xx <- sim(R0x,depx,h,M,selparsx,epsrx,dms,pctarg,selidx)

      # LF discrepancy

      phat <- xx$LF
      kllf <- pobs*log(pobs/phat)
      dlf <- sum(apply(kllf,2,function(x){sum(x[!is.nan(x)])}))

      # CPUE discrepancy

      if(seasonq) {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- apply(resq,2,mean)
        resq <- resq-lnq

      } else {

        resq <- log(I[,,fcpue]/xx$I)
        lnq <- mean(resq)
        resq <- resq-lnq

      }

      dcpue <- sum(dnorm(resq,0,sdcpue,TRUE))

      ## priors (parameters + stock status)

      # status priors

      Bmsyrat <- xx$S[,3]/xx$Bmsy
      Bmsyrat <- Bmsyrat[ybmsy]
      dSSB <- xx$S[,srec-1]/xx$B0
      dSSB <- dSSB[ydep]

      if(length(ybmsy) == 1) {

        sprior <- dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE)

      } else {

        sprior <- sum(dnorm(Bmsyrat,mubmsy,sdbmsy,TRUE))

      }

      if(length(ydep) == 1) {

        sprior <- sprior+dnorm(dSSB,mudep,sddep,TRUE)

      } else {

        sprior <- sprior+sum(dnorm(dSSB,mudep,sddep,TRUE))

      } 

      # parameter priors

      pprior <- sum(dnorm(epsr,0,sigmar,TRUE))

      ## ABC accept/reject:
      # 1. KL(LF data) < KL_max or reject immediately
      # 2. If 1 is true accept/reject given remaining discrepancy

      if(dlf < KLmax) {

        dtotnew <- dcpue+sprior+pprior-dlf
        pirat <- min(dtotnew-dtotold,0)
        uvar <- log(runif(1,0,1))
        accpt <- ifelse(pirat>uvar,TRUE,FALSE)

      } else {

        accpt <- FALSE

      }

      if(accpt) {

        parvecold <- parvecnew
        dtotold <- dtotnew
        if(n > burn) acp[g] <- acp[g]+1

      }
    }

    # outputs
  
    if(n > burn & (n-burn) %% thin == 0) theta.mcmc[(n-burn)/thin,] <- parvecold

    if(n %% 100 == 0) cat("Iteration",n,"of",burn+nits*thin,"\n")
 
  }

  return(list(pars=theta.mcmc,acp=acp))
} 
# }}}

# sim {{{
sim <- function(R0=1e6, dep=0.5, h=0.75, M=0.075, selpars, epsr, dms, pctarg,selidx) {

  na <- dms[1]
  ns <- dms[2]
  nf <- dms[3]
  nselg <- dms[4]

  # SPR ratio at exploited eqm given steepness and depletion
  # (based on derivation: dep = (4*h*rho+h-1)/(5*h-1))
  rhotarg <- (dep*(5*h-1)+1-h)/(4*h)

  # create selectivity-at-age

  sela <- get.sel.age(nf,nselg,selidx,selpars) 
    
  # target vector (rho+pc)
  targv <- logit(c(rhotarg,as.vector(pctarg)))

  # wrapper objective function to solve (minimise)

  # Estimate initial Fs
  
  hinit <- array(0.15*pctarg, dim=c(ns, nf)) 

  # scalar to give objective function some bite

  theta <- logit(as.vector(hinit))

  zz <- optim(theta,objfn.init,targv=targv,sela=sela,method=("L-BFGS-B"),control=list(trace=0))

  hinit <- array(ilogit(zz$par), dim=c(ns, nf))
  resinit <- initpdyn(c(ns, na, nf), srec, psi, M, as.vector(mata),
    as.vector(wta), as.vector(sela), as.vector(hinit)) 

  # relative H-split for MSY calcs
  ph <- as.vector(hinit[] / sum(hinit))

  msy <- optimise(msyfn,interval=c(0,0.7),ph=ph,sela=sela,maximum=TRUE)
  Hmsy <- msy$maximum
  Cmsy <- msy$objective
  resmsy <- msypdyn(c(ns,na,nf), srec, R0, h, psi, M, as.vector(mata),
    as.vector(wta), as.vector(sela), Hmsy * ph)
  
  Bmsy <- resmsy$Bmsy
  spr0 <- resinit$spr0
  B0 <- R0*spr0
  alp <- 4*h/(spr0*(1-h))
  bet <- (5*h-1)/(B0*(1-h))

  Bratio <- Bmsy/B0
  Rratio <- resmsy$Rmsy/R0
  hmsyv <- apply(array(Hmsy*ph,dim=c(ns,nf)),1,sum)

  if(!all.equal(Rratio, (4*h*Bratio)/(h*(5*Bratio-1)+1-Bratio)))
    warning("B-H invariant check - should be same as Rratio")

  # set up initial numbers-at-age for input to population dynamics

  Rinit <- R0*(4*h*dep)/(h*(5*dep-1)+1-dep)
  Ninit <- array(resinit$N,dim=c(na,ns,2))
  Ninit[] <- Ninit[]*Rinit
  nvec <- as.vector(Ninit)

  # expected catch @ hinit
  zinit <- msypdyn(c(ns,na,nf),srec,R0,h,psi,M,as.vector(mata),
    as.vector(wta),as.vector(sela),as.vector(hinit))
  Cinit <- array(zinit$C,dim=c(ns,nf))

  # main population stuff

  cvec <- as.vector(C)

  ## generating predicted LF and CPUE 
    
  # fishery for CPUE generation

  fcpue <- 1
  resp2 <- pdynlfcpue(c(ny,ns,na,nbins,nf),srec,R0,h,psi,epsr,spr0,M,
    as.vector(mata),as.vector(wta),as.vector(sela),nvec,cvec,as.vector(pla),fcpue)

  N <- array(resp2$N,dim=c(ny,na,ns,2))
  S <- array(resp2$S,dim=c(ny,ns))
  H <- array(resp2$H,dim=c(ny,ns,nf))
  LFhat <- array(resp2$LF,dim=c(ny,nbins,ns,nf))
  Ihat <- array(resp2$I,dim=c(ny,ns))

  # predicted LF distro for relevant fisheries

  LFhat <- apply(LFhat,c(2,4),sum)
  phat <- LFhat[,flf]
  phat <- apply(phat,2,function(x){x <- x/sum(x)})

  return(list(N=N,S=S,H=H,LF=phat,I=Ihat,Bmsy=Bmsy,Cmsy=Cmsy,Hmsy=hmsyv,B0=B0))

}
# }}}

# rsim {{{
rsim <- function(R0=1e6, dep=0.5, h=0.75, M=0.075, selpars, epsr, dms, pctarg,selidx) {

  na <- dms[1]
  ns <- dms[2]
  nf <- dms[3]
  nselg <- dms[4]

  # SPR ratio at exploited eqm given steepness and depletion
  # (based on derivation: dep = (4*h*rho+h-1)/(5*h-1))
  rhotarg <- (dep*(5*h-1)+1-h)/(4*h)

  # create selectivity-at-age

  sela <- get.sel.age(nf,nselg,selidx,selpars) 
    
  # target vector (rho+pc)
  targv <- logit(c(rhotarg,as.vector(pctarg)))

  # wrapper objective function to solve (minimise)

  # Estimate initial Fs
  
  hinit <- 0.15*pctarg[]

  # scalar to give objective function some bite

  theta <- logit(as.vector(hinit))

   zz <- optim(theta,robjfn.init,targv=targv,sela=sela,method=("L-BFGS-B"),control=list(trace=0))

  hinit <- array(ilogit(zz$par), dim=c(ns, nf))
  resinit <- rinitpdyn(c(ns,na,nf),srec,psi,M,mata,wta,sela,hinit) 

  # relative H-split for MSY calcs
  ph <- as.vector(hinit[] / sum(hinit))

  msy <- optimise(rmsyfn,interval=c(0,0.7),ph=ph,sela=sela,maximum=TRUE)
  Hmsy <- msy$maximum
  Cmsy <- msy$objective
  resmsy <- rmsypdyn(c(ns,na,nf),srec,R0,h,psi,M,mata,wta,sela,Hmsy*ph)
  
  Bmsy <- resmsy$Bmsy
  spr0 <- resmsy$spr0
  B0 <- R0*spr0
  alp <- 4*h/(spr0*(1-h))
  bet <- (5*h-1)/(B0*(1-h))
  Bratio <- Bmsy/B0
  Rratio <- resmsy$Rmsy/R0
  hmsyv <- apply(array(Hmsy*ph,dim=c(ns,nf)),1,sum)

  #if(!all.equal(Rratio, (4*h*Bratio)/(h*(5*Bratio-1)+1-Bratio)))
    #warning("B-H invariant check - should be same as Rratio")

  # set up initial numbers-at-age for input to population dynamics

  Rinit <- R0*(4*h*dep)/(h*(5*dep-1)+1-dep)
  Ninit <- array(resinit$N,dim=c(na,ns,2))
  Ninit[] <- Ninit[]*Rinit

  # expected catch @ hinit
  zinit <- rmsypdyn(c(ns,na,nf),srec,R0,h,psi,M,mata,wta,sela,hinit)
  Cinit <- array(zinit$C,dim=c(ns,nf))

  # main population stuff

  ## generating predicted LF and CPUE 
    
  # fishery for CPUE generation

  fcpue <- 1
  resp2 <- rpdynlfcpue(c(ny,ns,na,nbins,nf),srec,R0,h,psi,epsr,spr0,M,mata,wta,sela,Ninit,C,pla,fcpue)

  N <- resp2$N
  S <- resp2$S
  H <- resp2$H
  LFhat <- resp2$LF
  Ihat <- resp2$I

  # predicted LF distro for relevant fisheries

  LFhat <- apply(LFhat,c(2,4),sum)
  phat <- LFhat[,flf]
  phat <- apply(phat,2,function(x){x <- x/sum(x)})

  return(list(N=N,S=S,H=H,LF=phat,I=Ihat,Bmsy=Bmsy,Cmsy=Cmsy,Hmsy=hmsyv,B0=B0))

}
# }}}

# get.mcmc.vars {{{

get.mcmc.vars <- function(parsmat) {

  varlist <- list()                    
  nnits <- dim(parsmat)[1]
  for(nn in 1:nnits) {

    R0x <- exp(mcpars[nn,1])
    depx <- ilogit(mcpars[nn,2])
    epsrx <- mcpars[nn,3:(ny+1)]
    selvx <- exp(mcpars[nn,(ny+2):npar])
    selparsx <- cbind(selvx[1:nselg],selvx[(nselg+1):(2*nselg)],selvx[(2*nselg+1):(3*nselg)])
    xx <- rsim(R0x,depx,h,M,selparsx,epsrx,dms,pctarg,selidx)

    varlist[[nn]] <- list()
    varlist[[nn]][['Rtot']] <- apply(xx$N[,1,srec,],1,sum)
    varlist[[nn]][['SSB']] <- xx$S[,srec-1]
    varlist[[nn]][['dep']] <- xx$S[,srec-1]/xx$B0
    varlist[[nn]][['dbmsy']] <- xx$S[,srec-1]/xx$Bmsy
    varlist[[nn]][['Cmsy']] <- xx$Cmsy
    varlist[[nn]][['Ihat']] <- xx$I
    varlist[[nn]][['LFhat']] <- xx$LF 

    if(nn %% 100 == 0) cat("Iteration",nn,"of",nnits,"\n")
  }

  return(varlist)

}
# }}}

# plot.mcmc.vars {{{
plot.mcmc.vars <- function(varlist,type='dep') {

  nnits <- length(varlist)

  if(type == 'dep') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$dep
    vmin <- 0
    vmax <- 1
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975))
    plot(yrs,vq[2,],ylim=c(vmin,vmax),xlab='year',ylab='SSB depletion',col='blue',type='l')
    lines(yrs,vq[1,],lty=2,col='blue')
    lines(yrs,vq[3,],lty=2,col='blue') 

  }

  if(type == 'bmsy') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$dbmsy
    vmin <- 0
    vmax <- max(vv)
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975))
    plot(yrs,vq[2,],ylim=c(vmin,vmax),xlab='year',ylab='Bmsy ratio',col='blue',type='l')
    lines(yrs,vq[1,],lty=2,col='blue')
    lines(yrs,vq[3,],lty=2,col='blue') 

  } 

  if(type == 'rec') {

    vv <- matrix(nrow=nnits,ncol=ny)
    for(nn in 1:nnits) vv[nn,] <- varlist[[nn]]$Rtot
    vmin <- 0
    vmax <- max(vv)
    vq <- apply(vv,2,quantile,c(0.025,0.5,0.975))
    plot(yrs,vq[2,],ylim=c(vmin,vmax),xlab='year',ylab='Recruitment',col='blue',type='l')
    lines(yrs,vq[1,],lty=2,col='blue')
    lines(yrs,vq[3,],lty=2,col='blue') 

  }

  if(type == 'cpue') {

    vv <- array(dim=c(nnits,ny,ns))
    for(nn in 1:nnits) {
      
      tmpv <- varlist[[nn]]$Ihat
      iobs <- I[,,fcpue]
      if(seasonq) {

        resq <- log(iobs/tmpv)
        lnq <- apply(resq,2,mean)
        vv[nn,,] <- t(apply(tmpv,1,function(x,lnq){x <- x*exp(lnq)},lnq))*rlnorm(ny*ns,0,sdcpue)

      } else {

        resq <- log(iobs/tmpv)
        lnq <- mean(resq)
        vv[nn,,] <- tmpv*exp(lnq)*rlnorm(ny*ns,0,sdcpue)  

      }
    }   
    
    vq <- apply(vv,c(2,3),quantile,c(0.025,0.5,0.975))

    # ggplot the sumbitch

    vdf <- expand.grid(year=yrs,season=1:ns,obs=NA,hat=NA,lq=NA,uq=NA)
    vdf$obs <- as.vector(iobs)
    vdf$hat <- as.vector(vq[2,,])
    vdf$lq <- as.vector(vq[1,,]) 
    vdf$uq <- as.vector(vq[3,,]) 
    ggplot(vdf)+geom_line(aes(x=year,y=hat),colour='blue')+geom_line(aes(x=year,y=lq),colour='blue',linetype='dashed')+geom_line(aes(x=year,y=uq),colour='blue',linetype='dashed')+geom_point(aes(x=year,y=obs),colour='magenta')+facet_wrap(~season)+ylab("CPUE")

  }

  if(type == 'lf') {

    vv <- array(dim=c(nnits,nbins,nselg)) 
    for(nn in 1:nits) vv[nn,,] <- varlist[[nn]]$LF
    vq <- apply(vv,c(2,3),quantile,c(0.025,0.5,0.975))
    vdf <- expand.grid(length=mulbins,fishery=1:nselg,obs=NA,hat=NA,lq=NA,uq=NA) 
    vdf$obs <- as.vector(pobs)
    vdf$hat <- as.vector(vq[2,,])
    vdf$lq <- as.vector(vq[1,,]) 
    vdf$uq <- as.vector(vq[3,,]) 
    ggplot(vdf)+geom_line(aes(x=length,y=hat),colour='blue')+geom_line(aes(x=length,y=lq),colour='blue',linetype='dashed')+geom_line(aes(x=length,y=uq),colour='blue',linetype='dashed')+geom_point(aes(x=length,y=obs),colour='magenta')+facet_wrap(~fishery)+ylab("Length frequency")

  }

}
# }}}
