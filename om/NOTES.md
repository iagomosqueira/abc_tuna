---
title: 
author: Iago MOSQUEIRA (WMR) <iago.mosqueira@wur.nl>
tags: 
created: 22-10-2023T23:00:56+02:00
updated: 22-10-2023T23:00:56+02:00
---

# OM4

- LOAD SS3 base case run
- LOAD mcvars
 [1] "stock.n"   "m"         "catch.sel" "ssb"       "dep"       "srpars"    "refpts"    "hr"        "rec"       "index.hat"
- BUG: DIFF caa vs. catch.n(stk)
- SIMPLIFY season and start in 2000
- ADD from out mcvars
  - M
  - stock.n @ season 1
- SET stock.n[1,] as Q4 / exp(-M)
- CALCULATE F
- SET spwn at 0.83


