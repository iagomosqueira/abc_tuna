# data.R - DESC
# /home/mosquia/Active/ABC_tuna+iotc/abc_tuna/albabc/data.R

# Copyright (c) WUR & CSIRO, 2023.
# Authors: Richard HILLARY (CSIRO) <rich.hillary@csiro.au>
#          Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(ss3om)

source("utilities.R")

# LOAD FLStock and abc input set

stk <- readFLSss3("bootstrap/data/base")

inp <- setabc("bootstrap/data/base/abt.dat", stk, ymin=2000)


# SAVE

save(stk, inp, file="data/data.RData", compress="xz")
