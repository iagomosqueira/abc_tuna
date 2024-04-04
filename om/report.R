# report.R - DESC
# /home/mosqu003/Active/ABC_tuna+iotc/abc_tuna/om/report.R

# Copyright (c) WUR, 2024.
# Author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
#
# Distributed under the terms of the EUPL-1.2

# XX {{{
# }}}

pdf('report/figures.pdf')

# data.R

load('data/om.rda')

# HR ~ F

plot(unitMeans(harvest(stko)), unitSums(hrya)[, ac(2000:2020)]) +
  facet_wrap(~age, scales='free') +
  ggtitle("F ABC (red), HR OM (blue)")

plot(unitMeans(harvest(stko)), unitMeans(harvest(stky))) +
  facet_grid(~age, scales='free') +
  ggtitle("F ABC (red), SS3 (blue)")

# ABC vs. SS3

plot(FLStocks(ABC=stko, SS3=stky), metrics=total)

# FWD(F=0) 
plot(tes) + ggtitle("fwd(F=0), ABC SRR")

# FWD compored
plot(FLStocks(CORRECT=ates, ORIG=tes), metrics=total)

# REFPTS
(plot(remap(refpts(rps))) + ggtitle("Aggregated")) |
(plot(remap(refpts(rpsf))) + ggtitle("Female"))

# 
plot(unitSums(ssb(ates)) / (rep$SB0)) + ggtitle("SB/SB0") +
  geom_hline(yintercept=1)
plot(unitSums(ssb(ates)) / (rep$SB0 / 2)) + ggtitle("SB/(SB0 * 0.5)") +
  geom_hline(yintercept=1)

dev.off()
