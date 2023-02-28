# abc.R - DESC
# /home/mosquia/Active/Doing/ABC_tuna+iotc/abc_tuna/v2/abc.R

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

# --- simulator

# - arguments
#   - parameters: R0, dep, h
#   - biology
#   - fishery

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
      mata[a, s, g] <- sum(mlref  *dl)
      wta[a, s, g] <- sum(wlref, dl)

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
    as.vector(mata),as.vector(wta),as.vector(sela),nvec,cvec,as.vector(pla),fref)

  N <- array(resp2$N,dim=c(ny,na,ns,ng))
  S <- array(resp2$S,dim=c(ny,ns))
  H <- array(resp2$H,dim=c(ny,ns,nf))
  LF <- array(resp2$LF,dim=c(ny,nbins,ns,nf))
  I <- array(resp2$I,dim=c(ny,ns))

  return(list(N=N, S=S, H=H, LF=LF, I=I))

}
# }}}

#' @examples
#' system.time(
#'   d <- sim(R0=1e6, dep=0.5, h=0.75)
#' )
#' naa <- FLQuant(aperm(d$N, c(2,1,4,3)), dimnames=list(age=seq(0,19),
#'  year=2001:2010, unit=c("F", "M"), season=1:4))
