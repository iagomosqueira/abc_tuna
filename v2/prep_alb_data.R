######################################################
# prep IOTC ALB stock assessment data ################
######################################################
# R. Hillary & I. Mosqueira 2023 #####################
######################################################

library(FLCore)
library(ggplot2)
library(ggplotFL)

load("../data/data.RData")

# years of analyses

ymin <- 2000
ymax <- 2020
yrs <- ymin:ymax
ny <- length(yrs)
nf <- 6
ns <- 4

## catches

cdf <- subset(catches,year >= ymin)
ggplot(cdf)+geom_line(aes(x=year,y=data))+facet_grid(unit~season)
C <- array(0,dim=c(ny,ns,nf))
dimnames(C)[[1]] <- as.character(yrs)
for(f in 1:nf) {
  for(s in 1:ns) {

    zz <- subset(cdf,season == s & unit == f)
    yz <- as.character(zz$year)
    C[yz,s,f] <- zz$data

  }
}

## CPUE 

idf <- subset(cpue,year >= ymin)
ggplot(idf)+geom_line(aes(x=year,y=obs))+facet_grid(unit~season)

cpuef <- 1:4 # fleets for which we consider abundance indices for
nfcpue <- length(cpuef)
I <- array(dim=c(ny,ns,nfcpue))
dimnames(I)[[1]] <- as.character(yrs)
sx <- unique(idf$season)
for(f in cpuef) {
  for(s in 1:ns) {

    zz <- subset(idf,season == sx[s] & unit == f)
    yz <- as.character(zz$year)
    I[yz,s,f] <- zz$obs

  }
}

## LF

ldf <- subset(lencomp,year >= ymin)
fref <- 19
ggplot(subset(ldf,fleet==fref))+geom_point(aes(x=l,y=n))+facet_grid(season~year)
ldf$l <- as.numeric(as.character(ldf$length))

lmin <- 30
lmax <- 138
ldel <- 4 # 4cm length bins
lbins <- seq(lmin,lmax,by=ldel)
nbins <- length(lbins)-1
mulbins <- 0.5*(lbins[-1]+lbins[-(nbins+1)])
LF <- array(dim=c(ny,nbins,ns,nf))

## biology

ages <- as.numeric(dimnames(waa)$age)
na <- length(ages)
M <- 0.075 # per quarter i.e. seasonal

mata <- wta <- array(dim=c(na,ns,2))
for(s in 1:ns) {
  for(g in 1:2) {

    mata[,s,g] <- unname(mat@.Data[,1,1,1,1,1])
    wta[,s,g] <- unname(waa@.Data[,1,g,s,1,1])

  }
}

# p(l|a,s,g)

pla <- array(dim=c(nbins,na,ns,2))

# mean length-at-age (Schnute LvB)

a1 <- c(1,1)
a2 <- c(10,10)
l1 <- c(52.6038,52.036)
l2 <- c(103.8,110.6)
k <- c(0.38,0.34)
cvjuv <- c(0.06,0.06)
cvadu <- c(0.025,0.025) 
cva <- matrix(nrow=na,ncol=2)
amax <- max(ages)
for(s in 1:2) {
  for(a in 2:na) {
    
    atrue <- ages[a]
    cva[a,s] <- cvjuv[s]+(atrue-1)*(cvjuv[s]-cvadu[s])/(1-amax)

  }

  cva[1,s] <- cva[2,s]

}

mula <- sdla <- array(dim=c(na,ns,2))
for(g in 1:2) {
  for(s in 1:ns) {

    aa <- ages+(s-1)*0.25
    mula[,s,g] <- l1[g]+(l2[g]-l1[g])*(1-exp(-k[g]*(aa-a1[g])))/(1-exp(-k[g]*(a2[g]-a1[g])))
    sdla[,s,g] <- cva[,g]*mula[,s,g]

  }
}

for(g in 1:2) {
  for(s in 1:ns) {
    for(a in 1:na) {

      dx <- dnorm(mulbins,mula[a,s,g],sdla[a,s,g],FALSE)
      dx <- dx/sum(dx)
      pla[,a,s,g] <- dx 

    }
  }
}

# save it

save.image("alb_abcdata.rda")

