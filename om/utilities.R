# utilities.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/abc_tuna/om/utilities.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# XX {{{
# }}}

  # "Rtot", "dbmsy", "Cmsy", "hmsyrat" "H", "Ihat", "LFhat"

mc.output <- function(x) {

  nits <- length(x)
  dmns <- list(age=0:14, year=2000:2020, season=1:4, unit=c('F', 'M'))

  # - FLQuants at age

  # N - stock.n (y, a, s, u)
  stock.n <- Reduce(combine, lapply(x, function(i)
   FLQuant(aperm(i$N, c(2,1,4,3)), dimnames=dmns) / 1000
  ))

  # M - m
  m <- expand(FLQuant(unlist(lapply(x, '[[', 'M')),
    quant='age',dim=c(1,1,1,1,1,nits)),
    age=0:14, year=2000:2020, season=1:4, unit=c('F', 'M'))

  # sela - catch.sel (a, s, u, f)
  catch.sel <- Reduce(combine, lapply(x, function(i) {
    res <- FLQuant(dimnames=list(age=0:14, year=2000, unit=c('F', 'M'),
      season=1:4, area=1:6))
    res[] <- aperm(i$sela, c(1,3,2,4))
    return(res %/% apply(res, 2:6, max))
    }
  ))

  # catches (y, s, f)
  caf <- FLQuant(dimnames=list(year=2018:2020, season=1:4, area=1:6))
  caf[] <- run4$C[19:21,,]
  cap <- caf %/% areaSums(caf)
  sel <- yearMeans(areaMeans(catch.sel %*% cap))
  sel <- sel %/% apply(sel, 2:6, max)

  ggplot(iter(sel, 1), aes(x=age,y=data, group=unit)) +
    geom_line() +
    facet_wrap(~season)


  # - FLPar

  # B0
  B0 <- unlist(lapply(x, '[[', 'B0'))

  # R0
  R0 <- unlist(lapply(x, '[[', 'R0'))

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

  # dep
  dep <- Reduce(combine, lapply(x, function(i)
   FLQuant(i$dep, dimnames=list(age='all', year=2000:2020))
  ))

  # TODO: Rtot does not match unitSums(N[1,,,4])

  return(list(stock.n=stock.n, m=m, catch.sel=catch.sel, ssb=ssb, dep=dep,
    srpars=srpars, refpts=refpts))
}


buildOM <- function(stk, vars) {

  # WINDOW
  stk <- window(stk, start=2000)

  # ASSIGN stock.n and m
  stock.n(stk) <- vars$stock.n
  m(stk) <- vars$m

  # SIMPLIFY

}
