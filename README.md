# abc_tuna

# FILES

- data/: initial datasets
  - sa/: input SS3 files 2022 WPTmT SA
  - run.RData: output from base case runs
    - out: SS_output list
    - run: FLStock
    - retro: SS_retrosummaries
  - data.RData: list containing

# SS3 model structure

- 1954 - 2020
- 4 seasons
- 2 sexes
- 4 fishing areas
- 40 fleets, separated by gear, area and quarter
  - 23 fishing fleets (code_gear+area_quarter)
    - 1   F1_LL1_Q1
    - 2   F2_LL1_Q2
    - 3   F3_LL1_Q3
    - 4   F4_LL1_Q4
    - 5   F5_LL2_Q1
    - 6   F6_LL2_Q2
    - 7   F7_LL2_Q3
    - 8   F8_LL2_Q4
    - 9   F9_LL3_Q1
    - 10 F10_LL3_Q2
    - 11 F11_LL3_Q3
    - 12 F12_LL3_Q4
    - 13 F13_LL4_Q1
    - 14 F14_LL4_Q2
    - 15 F15_LL4_Q3
    - 16 F16_LL4_Q4
    - 17    F17_DN3
    - 18    F18_DN4
    - 19    F19_PS1
    - 20 F20_Other1
    - 21 F21_Other2
    - 22 F22_Other3
    - 23 F23_Other4
  - 17 survey fleets (gear+CPUE+area_quarter)
    - 24 LLCPUE1_Q1
    - 25 LLCPUE1_Q2
    - 26 LLCPUE1_Q3
    - 27 LLCPUE1_Q4
    - 28 LLCPUE2_Q1
    - 29 LLCPUE2_Q2
    - 30 LLCPUE2_Q3
    - 31 LLCPUE2_Q4
    - 32 LLCPUE3_Q1
    - 33 LLCPUE3_Q2
    - 34 LLCPUE3_Q3
    - 35 LLCPUE3_Q4
    - 36 LLCPUE4_Q1
    - 37 LLCPUE4_Q2
    - 38 LLCPUE4_Q3
    - 39 LLCPUE4_Q4
    - 40    DNCPUE4


# src
  - init_pdyn.cpp: generate initial population-per-recruit for ABC simulator.
  - msy_pdyn.cpp: generate eqm popn for MSY estimation IOTC ABC simulator.
  - pdyn.cpp: dynamic non-eqm population dynamics for  ABC simulator.
  - pdyn_lfcpue.cpp: generate predicted LF and CPUE

# test
  - test_mdarr.cpp
  - test_mdclass.cpp

# old
  - mdclass.cpp
  - abc_pdyn.cpp

# include
  - mdarrays.h
  - mdarrclass.h

# R
  - abc.R
  - init_pdyn.R
  - test.R
  - test_season.R
