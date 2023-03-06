# utilities.R - DESC
# abc_tuna/albabc/utilities.R

# Copyright (c) WUR & CSIRO, 2023.
# Authors: Richard HILLARY (CSIRO) <rich.hillary@csiro.au>
#          Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# logit & ilogit {{{

logit <- function(x){
  return(log(x/(1-x)))
}

ilogit <- function(x){
  return(1/(1+exp(-x)))
}
# }}}

# readss3abc {{{
readss3abc <- function(dat) {

  # - EXTRACT catch per fleet, year and season
  catches <- data.table(dat$catch)[year != -999]
  setkey(catches, "fleet")
  
  # LOOKUP table: group LL & Other by area, keep PS, drop DN.
  lookup <- data.table(
    fleet=c(seq(1, 16), 19, seq(20, 23)),
    unit=c(rep(seq(1, 4), each=4), 5, rep(6, 4)))
  setkey(lookup, "fleet")

  # RECODE fleets
  ncatches <- catches[lookup, on="fleet"]
  setnames(ncatches, c('seas', 'catch'), c('season', 'data'))

  catches <- ncatches[, .(data=sum(data, na.rm=TRUE)), by=.(year, season, unit)]
  setorder(catches, year, season, unit)

  # - EXTRACT length composition data
  lencomp <- data.table(dat$lencomp)

  # SELECT 'f' columns
  lencomp <- lencomp[, c(1,2,3,6,7:61)]
  lnames <- ac(seq(30, 138, by=2))
  setnames(lencomp, c("year", "season", "fleet", "Nsamp", lnames))

  # COERCE ALL length columns to double
  lencomp[, (lnames) := lapply(.SD, as.double), .SDcols=lnames]

  # RESHAPE to long
  lencomp <- melt(lencomp, id=c("year", "season", "fleet", "Nsamp"),
    measure=seq(5, 59), variable.name = "length", value.name = "n")

  # FIX season
  lencomp[, season:=(season + 0.5) / 3]

  # A new fleet
  lencomp <- lencomp[lookup[fleet < 20,], on="fleet"]

  # - EXTRACT CPUE

  cpue <- data.table(dat$CPUE)
  setnames(cpue, c("seas"), c("season"))

  lookup <- data.table(
    index=seq(24, 38),
    unit=rep(seq(1, 4), each=4)[-16])
  setkey(lookup, "index")

  cpue <- cpue[lookup, on="index"]

  cpue[,season:=(season + 0.5) / 3]
  cpue <- cpue[, .(obs=mean(obs)), by=.(year, season, unit)]

  return(list(cpue=cpue, lencomp=lencomp, catches=catches))
}
# }}}

# setabc {{{
setabc <- function(file, stk, ymin=2000, ymax=max(inp$catches$year, na.rm=TRUE)) {

  # LOAD dat file
  dat <- SS_readdat(file)

  # READ input data from SS3 dat
  inp <- readss3abc(dat)

  # DIMENSIONS
  yrs <- seq(ymin, ymax)
  ny <- length(yrs)

  nf <- length(unique(inp$catches$unit))
  ns <- length(unique(inp$catches$season))

  # REFORMAT catches (C)

  cdf <- inp$catches[year >= ymin,]
  C <- array(0, dim=c(ny, ns, nf))
  dimnames(C)[[1]] <- as.character(yrs)
  for(f in 1:nf) {
    for(s in 1:ns) {
      zz <- subset(cdf,season == s & unit == f)
      yz <- as.character(zz$year)
      C[yz,s,f] <- zz$data
    }
  }

  # REFORMAT CPUE (I)

  idf <- inp$cpue[year >= ymin,]
  cpuef <- 1:4 
  # fleets for which we consider abundance indices for
  nfcpue <- length(cpuef)
  I <- array(dim=c(ny,ns,nfcpue))
  dimnames(I)[[1]] <- as.character(yrs)
  sx <- unique(idf$season)
  for(f in cpuef) {
    for(s in 1:ns) {

      zz <- subset(idf,season == sx[s] & unit == f)
      yz <- as.character(zz$year)
      I[yz,s,f] <- zz$obs
    }
  }

  # REFORMAT LF (LF)
  lencomp <- inp$lencomp
  ldf <- lencomp[year >= ymin, ]
  fref <- 19
  ldf$l <- as.numeric(as.character(ldf$length))
  lmin <- min(as.numeric(lencomp$length))
  lmax <- max(as.numeric(lencomp$length))
  ldel <- 4 # 4cm length bins
  lorig <- as.numeric(unique(ldf$l))
  lbins <- seq(lmin,lmax,by=ldel)
  lagg <- as.numeric(lbins)
  nbins <- length(lbins)-1
  mulbins <- 0.5*(lbins[-1]+lbins[-(nbins+1)])
  LF <- array(dim=c(ny,nbins,ns,nf))

  for(y in 1:ny) {
    for(f in 1:nf) {

      if(f == 1) fref <- 1:4
      if(f == 2) fref <- 5:8 
      if(f == 3) fref <- 9:12 
      if(f == 4) fref <- 13:16 
      
      if(f <= 4) {
  
        zz <- subset(ldf,year == yrs[y] & fleet %in% fref)
        for(ff in 1:4) {
          
          zzz <- subset(zz,fleet  == fref[ff])
          if(dim(zzz)[1] > 0) {
  
            nn <- zzz$n
            names(nn) <- lorig
            nx <- rep(0,nbins)
            ndf <- data.frame(l=lorig,n=nn)
            for(ll in 1:nbins) nx[ll] <- sum(subset(ndf,l >= lagg[ll] & l < lagg[ll+1])$n)
            LF[y,,ff,f] <- nx
  
          }
        }
      }
  
      if(f == 6) {
  
        fref <- 19
        zz <- subset(ldf,year == yrs[y] & fleet == fref) 
        for(s in 1:ns) {
  
          zzz <- subset(zz,season == sx[s])
          if(dim(zzz)[1] > 0) {
  
            nn <- zzz$n
            names(nn) <- lorig
            nx <- rep(0,nbins)
            ndf <- data.frame(l=lorig,n=nn)
            for(ll in 1:nbins) nx[ll] <- sum(subset(ndf,l >= lagg[ll] & l < lagg[ll+1])$n) 
            LF[y,,s,f] <- nx 
  
          }
        }
      }
    }
  }
  
  # year sums (probably best use of data)
  
  LFagg <- apply(LF,c(2,3,4),sum,na.rm=T)
  LFfits <- apply(LF,c(2,4),sum,na.rm=T)

  # SETUP biology (mata, wta)
  ages <- seq(0, dat$Nages)
  na <- length(ages)
  M <- 0.075 # per quarter i.e. seasonal

  mata <- wta <- array(dim=c(na,ns,2))
  for(s in 1:ns) {
    for(g in 1:2) {
      mata[,s,g] <- unname(mat(stk)@.Data[,1,1,1,1,1])
      wta[,s,g] <- unname(stock.wt(stk)@.Data[,1,g,s,1,1])
    }
  }

  # p(l|a,s,g) (pla)

  pla <- array(dim=c(nbins,na,ns,2))

  # mean length-at-age (Schnute LvB) (cva)

  a1 <- c(1,1)
  a2 <- c(10,10)
  l1 <- c(52.6038,52.036)
  l2 <- c(103.8,110.6)
  k <- c(0.38,0.34)
  cvjuv <- c(0.06,0.06)
  cvadu <- c(0.025,0.025) 
  cva <- matrix(nrow=na,ncol=2)
  amax <- max(ages)
  for(s in 1:2) {
    for(a in 2:na) {
    
      atrue <- ages[a]
      cva[a,s] <- cvjuv[s]+(atrue-1)*(cvjuv[s]-cvadu[s])/(1-amax)
    }
  cva[1,s] <- cva[2,s]
  }

  # (mula, sdla)

  mula <- sdla <- array(dim=c(na,ns,2))
  for(g in 1:2) {
    for(s in 1:ns) {
      aa <- ages+(s-1)*0.25
      mula[,s,g] <- l1[g]+(l2[g]-l1[g])*(1-exp(-k[g]*(aa-a1[g])))/(1-exp(-k[g]*(a2[g]-a1[g])))
      sdla[,s,g] <- cva[,g]*mula[,s,g]
    }
  }

  for(g in 1:2) {
    for(s in 1:ns) {
      for(a in 1:na) {
        dx <- dnorm(mulbins,mula[a,s,g],sdla[a,s,g],FALSE)
        dx <- dx/sum(dx)
        pla[,a,s,g] <- dx 
      }
    }
  }

  return(list(C=C, I=I, LF=LF, pla=pla, cva=cva, mula=mula, sdla=sdla))
}
# }}}

# sim {{{
sim <- function(R0=1e6, dep=0.5, h=0.75) {

  # DIMENSIONS
  # seasons
  ns <- 4
  # ages
  na <- 20
  # TODO: First age is 0 or 1? 
  ages <- seq(1, na)
  # fisheries
  nf <- 2
  # sexes
  ng <- 2
  # length bins
  nl <- 20

  # VARIABLES
  # sex ratio at birth
  psi <- 0.5
  # recruitment season
  srec <- 2

  # growth-by-sex
  k <- c(0.2, 0.2)
  linf <- c(100, 100)
  l0 <- c(10, 10)
  sdla <- c(0.1, 0.1)

  # maturity-at-length
  ml50 <- 60
  ml95 <- 80
  
  # selectivity-at-length
  sl50 <- c(45, 45)
  sl95 <- c(65, 65)
  
  # natural mortality (by season)
  # TODO: Set as vector at age?
  M <- 0.1
  
  # weight-at-length
  alw <- 2e-6
  blw <- 3

  # INITIALIZE output objects
  # C
  C <- array(10000, dim=c(ns, nf))

  # BIOLOGY

  # mean length-at-age by season/sex
  
  mula <- array(dim=c(na, ns, ng))

  for(g in seq(ng)) {
    for(s in seq(ns)) {
    a <- ages + (s - 1) * 0.25
    mula[, s, g] <- l0[g] + (linf[g] - l0[g]) * (1 - exp(-k[g] * a))
    }
  }

  # maturity-at-age, weight-at-age and selectivity-at-age
  mata <- wta <- array(dim=c(na, ns, 2))
  sela <- array(dim=c(na, ns, 2, nf))

  for(g in seq(ng)) {
    for(s in seq(ns)) {
      for(a in ages) {

      lmin <- max(0, mula[a, s, g] * (1 - sdla[g] * 1.96))
      lmax <- mula[a, s, g] * (1 + sdla[g] * 1.96)
      lref <- seq(lmin, lmax, length=nl)
      dl <- dlnorm(lref, log(mula[a,s,g]), sdla[g])
      dl <- dl / sum(dl)
      mlref <- 1 / (1 + 19 ^ (-(lref - ml50) / (ml95-ml50)))
      wlref <- alw * lref ^ blw
      mata[a, s, g] <- sum(mlref*dl)
      wta[a, s, g] <- sum(wlref*dl)

        for(f in seq(nf)) {
          slref <- 1 / (1 + 19 ^ (-(lref - sl50[f]) / (sl95[f] - sl50[f])))
          sela[a, s, g, f] <- sum(slref * dl)
        }
      }
    }
  }
  
  # SPR ratio at exploited eqm given steepness and depletion
  # (based on derivation: dep = (4*h*rho+h-1)/(5*h-1))
  rhotarg <- (dep * (5 * h - 1) + 1 - h)/(4 * h)

  # catch fraction by season 
  pctarg <- C / sum(C)

  # Estimate initial Fs
  hinit <- array(0.025, dim=c(ns, nf))
  res <- initpdyn(c(ns, na, nf), srec, psi, M, as.vector(mata),
    as.vector(wta), as.vector(sela), as.vector(hinit))

  # target vector (rho+pc)
  targv <- logit(c(rhotarg,pctarg))

  # wrapper objective function to solve (minimise)

  fnscale <- 1 # scalar to give objective function some bite

  objfn.init <- function(theta) {

    hxinit <- 1 / (1 + exp(-theta))
    resx <- initpdyn(c(ns, na, nf), srec, psi, M, as.vector(mata),
      as.vector(wta), as.vector(sela), hxinit) 

    px <- resx$C / sum(resx$C)
    tmpv <- c(resx$rho, as.vector(px))
    objv <- logit(tmpv)

    return(fnscale * (sum((objv - targv) ^ 2)))

  }

  theta <- logit(as.vector(hinit))

  system.time(zz <- optim(theta, objfn.init, method=("L-BFGS-B"),
    control=list(trace=0)))

  hinit <- array(ilogit(zz$par), dim=c(ns, nf))
  resinit <- initpdyn(c(ns, na, nf), srec, psi, M, as.vector(mata),
    as.vector(wta), as.vector(sela), as.vector(hinit)) 

  # MSY estimation
  res <- msypdyn(c(ns, na, nf), srec, R0, h, psi, M,
    as.vector(mata), as.vector(wta), as.vector(sela), hinit)

  # relative H-split for MSY calcs
  ph <- as.vector(hinit[] / sum(hinit))

  # MSY wrapper
  msyfn <- function(H) {

    hx <- H * ph
    resx <- msypdyn(c(ns,na,nf),srec,R0,h,psi,M,as.vector(mata),
      as.vector(wta),as.vector(sela),hx)
    
    return(sum(resx$C))
  }

  msy <- optimise(msyfn, interval=c(0, 0.5), maximum=TRUE)

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

  ny <- 10
  epsr <- rep(0,ny-1)
  Cb <- array(dim=c(ny,ns,nf))
  for(y in 1:ny)
    Cb[y,,] <- Cinit
  cvec <- as.vector(Cb)

  resp <- pdyn(c(ny,ns,na,nf),srec,R0,h,psi,epsr,spr0,M,as.vector(mata),
    as.vector(wta),as.vector(sela),nvec,cvec)

  N <- array(resp$N,dim=c(ny,na,ns,2))
  S <- array(resp$S,dim=c(ny,ns))
  H <- array(resp$H,dim=c(ny,ns,nf))

  # generating predicted LF and CPUE #
  # p(l | a, s, g)

  lbins <- seq(0,120,by=6)
  nbins <- length(lbins)-1
  mulbins <- 0.5*(lbins[-1]+lbins[-length(lbins)])

  pla <- array(dim=c(nbins,na,ns,2))

  for(g in 1:2) {
    for(s in 1:ns) {
      for(a in 1:na) {
        dx <- dnorm(log(mulbins),log(mula[a,s,g]),sdla[g],FALSE)
        dx <- dx/sum(dx)
        pla[,a,s,g] <- dx 
      }
    }
  }

  # fishery for CPUE generation

  fref <- 1

  resp2 <- pdynlfcpue(c(ny,ns,na,nl,nf),srec,R0,h,psi,epsr,spr0,M,
    as.vector(mata),as.vector(wta),as.vector(sela),nvec,cvec,as.vector(pla),
    fref)

  N <- array(resp2$N,dim=c(ny,na,ns,ng))
  S <- array(resp2$S,dim=c(ny,ns))
  H <- array(resp2$H,dim=c(ny,ns,nf))
  LF <- array(resp2$LF,dim=c(ny,nbins,ns,nf))
  I <- array(resp2$I,dim=c(ny,ns))

  return(list(N=N, S=S, H=H, LF=LF, I=I))

}
# }}}
