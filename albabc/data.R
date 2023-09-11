# data.R - DESC
# /home/mosquia/Active/ABC_tuna+iotc/abc_tuna/albabc/data.R

# Copyright (c) WUR & CSIRO, 2023.
# Authors: Richard HILLARY (CSIRO) <rich.hillary@csiro.au>
#          Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(ss3om)
library(data.table)

source("utilities.R")

# LOAD SS3 inputs & outputs

stk <- readFLSss3("bootstrap/data/base")
inp <- setabc("bootstrap/data/base/abt.dat", stk, ymin=2000)
out <- readOutputss3("boot/data/base")

# DIMENSIONS

yrs <- seq(2000, 2020)
ny <- length(yrs)
nf <- 6
ns <- 4


# --- EXTRACT data: C, I, LF, pla, cva, mula, sdla

# C, [ny,ns,nf]

C <- data.table(out$catch)

# I, [ny, ns, ni]

# "LF"   "pla"  "cva"  "mula" "sdla"


# --- SETUP

  fnscale <- 1
  ybmsy <- c(ny-1,ny)
  mubmsy <- c(2.25,2)
  sdbmsy <- c(0.35,0.35)
  ydep <- 1
  mudep <- 0.5
  sddep <- 0.1




# SAVE

save(stk, inp, file="data/data.RData", compress="xz")
