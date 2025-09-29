# summarise the OM outputs

library(ggplot2)

# run

load("alb_abc_run5a.rda")

source("utilities.R")

# summaries

get.mcmc.vars(mcvars,'dep')
get.mcmc.vars(mcvars,'bmsy')
get.mcmc.vars(mcvars,'hmsy')

# plots

plot.mcmc.vars(mcvars,'dep')
plot.mcmc.vars(mcvars,'bmsy')
plot.mcmc.vars(mcvars,'rec')
plot.mcmc.vars(mcvars,'hmsy')
plot.mcmc.cpue(mcvars)
plot.mcmc.lf(mcvars)
plot.mcmc.sel(mcpars)

# sigmaR stuff

sigrpost <- mcpars[,npar+3]
dpost <- density(sigrpost)
dprior <- density(psigmaR)
p1 <- dpost$y
p1 <- p1/sum(p1)
x1 <- dpost$x
p2 <- dprior$y
p2 <- p2/sum(p2)
x2 <- dprior$x
pmax <- max(c(max(p1),max(p2)))
xmin <- min(c(min(x1),min(x2)))
xmax <- max(c(max(x1),max(x2)))
plot(x1,p1,type='l',xlim=c(xmin,xmax),ylim=c(0,pmax),xlab=expression(sigma[R]),ylab='density',main='OM scenario R2a')
lines(x2,p2,lty=2,col='purple')
legend("topright",lty=c(1,2),col=c("black","purple"),legend=c("Posterior","Prior"),bty='n')

round(quantile(psigmaR,c(0.025,0.5,0.975)),2)
round(quantile(sigrpost,c(0.025,0.5,0.975)),2)



