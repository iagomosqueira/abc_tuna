# test_data.R - TEST hr, catch and catch match
# abc_tuna/v3/test.R

# Copyright (c) WMR, 2025.
# Author: Iago MOSQUEIRA <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(mse)

source("utilities.R")

load('data/om5b.rda')

# --- TEST hr and catch.n {{{

# hr {{{
hrf <- function(om) {

  # C_f / sum(S_f * N * W)
  lapply(fisheries(om), function(x) {
    # ADD catch across sexes
    res <- unitSums(catch(x[[1]])) /
      # ADD N * WT * S_f across ages and units
      unitSums(quantSums(n(biol(om)) * wt(biol(om)) * catch.sel(x[[1]])))
    # SET to max = 0.9
    res <- ifelse(res > 0.9, 0.9, res)

    units(res) <- 'hr'
    
    return(res)
  })
}

# }}} 

# COMPUTE hr_fys
hom <- hrf(om)

dim(hom[[1]]) # [1]   1  21   1   4   1 500

# COMPUTE catch.n by age, year, unit, season, iter
can <- lapply(hom, function(x)
  expand(x, age=0:14, unit=c('F', 'M')) * n(biol(om)))

# MULTIPLY catch.n by catch.wt by fishery
caw <- Map('*', e1=can, e2=lapply(fisheries(om), function(x) landings.wt(x[[1]])))

# COMPARE computed from HR, ...
# INFO: ADD over ages(sexes(seasons(list elements in caw
quantSums(unitSums(seasonSums(Reduce('+', caw))))

# and calculated from catch.n * catch.wt in om
# INFO: ADD over ages(sexes(seasons(elements of catch by fishery list
quantSums(unitSums(seasonSums(Reduce('+', catch(fisheries(om))))))

# DIRECT calculation on F1

# TODO: fishery(om, 1)
print(seasonSums(unitSums(catch(fisheries(om)[[1]][[1]]))))

# H_f = C_f / sum(S_f * N * W)
hr1 <- unitSums(catch(fisheries(om)[[1]][[1]])) /
  quantSums(unitSums(catch.sel(fisheries(om)[[1]][[1]]) * n(biol(om)) * wt(biol(om))))

# COMPUTE catch.n by age, year, unit, season, iter
ca1 <- expand(hr1, age=0:14, unit=c('F', 'M')) * n(biol(om))

# COMPARE computed from HR, ...
quantSums(unitSums(seasonSums(ca1 * wt(biol(om)))))

# and calculated from catch.n * catch.wt in om
unitSums(seasonSums(catch(fisheries(om)[[1]][[1]])))

# }}}

# --- TEST fwdabc.om {{{

library(Rcpp)
library(parallel)
library(mvtnorm)

sourceCpp("utilities/init_pdyn.cpp")
sourceCpp("utilities/msy_pdyn.cpp")
sourceCpp("utilities/pdyn_lfcpue.cpp")

# LOAD pcbar (catch proportions by fleet & season)
pcbar <- as.matrix(fread('data/pcbar.dat'))

# LOAD ALK [len, age, season, sex]
load('data/pla.rda')

# LOAD om
load('data/om5b.rda')

its <- seq(50)
iy <- 2017

# - TEST fwd iy:2020 with historical catch

ctrl <- fwdControl(year=iy:2020, quant="catch",
  value=c(unitSums(seasonSums(Reduce('+', catch(om))))[, ac(iy:2020),,,,1]))

tes1 <- fwdabc.om(iter(om, its), ctrl, pcbar=pcbar, pla=pla)

# GET catch
quantSums(unitSums(seasonSums(Reduce('+', catch(fisheries(tes1$om))))))
quantSums(unitSums(seasonSums(Reduce('+', catch(fisheries(om))))))

# GET hr
plot(Reduce('+', unitMeans(hrf(tes1$om))),
  Reduce('+', unitMeans(hrf(om))))

# }}}
