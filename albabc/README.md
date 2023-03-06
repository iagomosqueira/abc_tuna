
- LOAD and format data.

```
setabc() {
  SS_readdat()
  readss3abc()
}

- C[ny, ns, nf], catch by year, season & fleet.
- I[ny, ns, nf(LL)], CPUE index by year, season & fleet (only LL).
- LF[ny, na, ns, nf], length-frequency data by year, age, season and fleet.
- pla, cva, mula, sdla

```

- GENERATE simulated population from priors and biology.

```
sim() {

  initpdyn(dm, srec, psi, M, mata, wta, sela, hinit)
    - out: rho, C, N, spr0
  msypdyn(dm, srec, R0, hh, psi, M, mata, wta, sela, hinit)
    - out: rho, C, spr0, Bmsy, Rmsy
  pdyn(dm, srec, R0, hh, psi, epsr, spr0, M, mata, wta, sela, Ninit, Cb)
    - out: S, N, H
  pdynlfcpue(dm, srec, R0, hh, psi, epsr, spr0, M, mata, wta, sela, Ninit, Cb,
    pla, fref)
    - out: S, N, H, I, LF
}
```

- MAIN


```
main(data, priors) {

  sim(new_params) 
    
  get_observations(popn)

  calc_discrepancy(data, popn)

  accept_reject()

  iterator++
}
```
