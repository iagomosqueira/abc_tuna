# utilities.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/abc_tuna/om/utilities.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


# N - stock.n
# Rtot - rec
# SSB - ssb
# dep - dep
# dbmsy - dep@BMSY
# Cmsy - refpts$MSY
# hmsyrat - ratio hr/hrmsy year
# H - hr
# Ihat - (index.hat)
# LFhat - (lf.hat)
# B0 - refpts$B0
# R0 - refpts$R0
# M - m
# h - srpars$s
# sela - catch.sel

# mc.output {{{

mc.output <- function(x, C) {

  nits <- length(x)
  dmns <- list(age=0:14, year=2000:2020, season=1:4, unit=c('F', 'M'))

  # - FLQuants at age

  # N - stock.n (y, a, s, u)
  stock.n <- Reduce(combine, lapply(x, function(i)
   FLQuant(aperm(i$N, c(2,1,4,3)), dimnames=dmns, units='1000') / 1000
  ))

  # M - m
  m <- expand(FLQuant(unlist(lapply(x, '[[', 'M')),
    quant='age',dim=c(1,1,1,1,1,nits)),
    age=0:14, year=2000:2020, season=1:4, unit=c('F', 'M'))

  # Ihat - index.hat
  index.hat <- Reduce(combine, lapply(x, function(i)
   FLQuant(c(i$Ihat), dimnames=list(age='all', year=2000:2020, season=1:4))
  ))

  # H [y, s, f] - hr
  hr <- FLQuant(unlist(lapply(x, '[[', 'H')),
    dimnames=list(age='all', year=2000:2020, season=1:4, area=1:6,
    iter=seq(nits)), units='hr')

  # sela - catch.sel (a, s, u, f)
  catch.sel <- Reduce(combine, lapply(x, function(i) {
    res <- FLQuant(dimnames=list(age=0:14, year=2000, unit=c('F', 'M'),
      season=1:4, area=1:6), units='')
    res[] <- aperm(i$sela, c(1,3,2,4))
    return(res %/% apply(res, 2:6, max))
    }
  ))
 
  # HR @age[a,y,s,f,u]
  hra <- expand(hr, age=0:14, unit=c('F', 'M'),
    fill=TRUE) * expand(catch.sel, year=2000:2020, fill=TRUE)

  # catches (y, s, f)
  caf <- FLQuant(dimnames=list(year=2018:2020, season=1:4, area=1:6))
  caf[] <- C[19:21,,]
  cap <- caf %/% areaSums(caf)
  sel <- yearMeans(areaMeans(catch.sel %*% cap))
  sel <- sel %/% apply(sel, 2:6, max)

  # - FLPar

  # B0
  B0 <- unlist(lapply(x, '[[', 'B0'))

  # R0, value in thousands
  R0 <- unlist(lapply(x, '[[', 'R0')) / 1000

  # h
  h <- unlist(lapply(x, '[[', 'h'))

  # Cmsy
  cmsy <- unlist(lapply(x, '[[', 'Cmsy'))

  # srpars
  srpars <-  FLPar(v=B0, R0=R0, s=h)

  # refpts
  refpts <- FLPar(R0=R0, B0=B0, MSY=cmsy)

  # - FLQuant

  # SSB
  ssb <- Reduce(combine, lapply(x, function(i)
   FLQuant(i$SSB, dimnames=list(age='all', year=2000:2020))
  ))

  # Rtot
  rec <- Reduce(combine, lapply(x, function(i)
   FLQuant(i$Rtot, dimnames=list(age='0', year=2000:2020))
  ))

  # dep
  dep <- Reduce(combine, lapply(x, function(i)
   FLQuant(i$dep, dimnames=list(age='all', year=2000:2020))
  ))

  # TODO: Rtot does not match unitSums(N[1,,,4])

  return(list(stock.n=stock.n, m=m, catch.sel=catch.sel, ssb=ssb, dep=dep,
    srpars=srpars, refpts=refpts, hr=hr, rec=rec, index.hat=index.hat,
    hra=hra))
}
# }}}

# buildOM (stk,vars) {{{

buildOM <- function(stk, vars) {

  # WINDOW
  stk <- window(stk, start=2000)

  # ASSIGN stock.n and m
  stock.n(stk) <- vars$stock.n
  m(stk) <- vars$m

  # SIMPLIFY

}
# }}}
