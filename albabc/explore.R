# explore.R - DESC
# /home/mosquia/Active/ABC_tuna+iotc/abc_tuna/albabc/explore.R

# Copyright (c) WUR, 2023.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2


library(ss3om)

# LOAD SS_output and FLStock

out <- readOutputss3("bootstrap/data/base")

stk <- readFLSss3("bootstrap/data/base")

# 
