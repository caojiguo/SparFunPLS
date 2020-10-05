rm(list=ls())

# Please change the directory of this file in your computer
setwd("~/Desktop/SFPLS_Supplementary_Materials_Oct2020/SFPLS")

install.packages("psych")
install.packages("pls")
install.packages("plsdof")
install.packages("fda")
install.packages("MASS")
install.packages("glmnet")

library(psych)
library(pls)
library(plsdof)
library(fda)
library(MASS)
library(glmnet)
source("functions.R")



# GENERATE DATA USING B-SPLINES BASIS with 50 evenly-spaced knots
# number of simulations = 100
# sample size = 500
nknots  = 50
norder  = 4
snr     = 5
nsim    = 100
n       = 500
ntest   = 5000
domain  = rng = c(0,1)


# ************** scenario I: beta1(t) ************** #
betaind = 1
Y = array(NA,c(n,nsim))
X = list()

for(itersim in 1:nsim)
{
  dat = data.generator(n=n, nknots=nknots, norder=norder, domain = c(0, 1), snr=snr, betaind=betaind)
  Y[,itersim]  = dat$Y
  X[[itersim]] = dat$X
}


Ytest = array(NA,c(ntest,nsim))
Xtest = list()

for(itersim in 1:nsim)
{
  dat = data.generator(n=ntest, nknots=nknots, norder=norder, domain = c(0, 1), snr=snr, betaind=betaind)
  Ytest[,itersim]  = dat$Y
  Xtest[[itersim]] = dat$X
}

# FPCR method
rng     = c(0,1)
Mobs    = 500
tobs    = seq(rng[1],rng[2],length.out = Mobs+1)
d       = 3
K       = 30
lambda  = 10^(-10:6)

lam100opt = lam500opt = K100opt = K500opt = pmse100 = pmse500 = rep(NA,nsim)
beta100 = beta500 = array(NA,c(Mobs+1,nsim))

for(itersim in 1:nsim)
{
  # sample size n = 100
  nsample = 100
  x100simcoef = X[[itersim]]$coefs[,1:nsample]
  x100simbasis = X[[itersim]]$basis
  x100sim = fd(coef=x100simcoef,basisobj=x100simbasis)
  FPCRfit = fit.FPCR.BIC(y=Y[1:nsample,itersim],x=x100sim,K=K,Mobs=Mobs,lambda=lambda,rng=rng)
  K100opt[itersim]   = FPCRfit$Kopt
  lam100opt[itersim] = FPCRfit$lamopt
  beta100[,itersim]  = FPCRfit$betaobs
  
  # predict on test set
  x100testsimcoef  = Xtest[[itersim]]$coefs[,1:nsample]
  x100testsimbasis = Xtest[[itersim]]$basis
  x100testsim      = fd(coef=x100testsimcoef,basisobj=x100testsimbasis)
  predict100       = predict.FPCR(y=Ytest[1:nsample,itersim], x=x100testsim, rng=rng, Mobs=Mobs,fpcr_Rfit=FPCRfit)
  pmse100[itersim] = predict100$PMSE
  
  
  
  # sample size n = 500
  nsample = 500
  FPCRfit = fit.FPCR.BIC(y=Y[,itersim],x=X[[itersim]],K=K,Mobs=Mobs,lambda=lambda,rng=rng)
  K500opt[itersim]   = FPCRfit$Kopt
  lam500opt[itersim] = FPCRfit$lamopt
  beta500[,itersim]  = FPCRfit$betaobs
  
  # predict on test set
  predict500       = predict.FPCR(y=Ytest[,itersim], x=Xtest[[itersim]], rng=rng, Mobs=Mobs,fpcr_Rfit=FPCRfit)
  pmse500[itersim] = predict500$PMSE
  
  print(itersim)
}

beta0 = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
plot(tobs, beta100[,1],type="l",ylim=c(-1,1.5))
matlines(tobs,beta100)
lines(tobs,beta0,lwd=2)

plot(tobs, beta500[,1],type="l",ylim=c(-1,1.5))
matlines(tobs,beta500)
lines(tobs,beta0,lwd=2)


# ise
nullrg = 251:501
nnrg   = 1:250

ise100nn   = apply(beta100[nnrg,],2,ise2,b0=beta0[nnrg],rng=c(0,0.5))
isenull100 = apply(beta100[nullrg,],2,ise,b0=beta0[nullrg],rng=c(0.5,1))  
ise500nn   = apply(beta500[nnrg,],2,ise2,b0=beta0[nnrg],rng=c(0,0.5))
isenull500 = apply(beta500[nullrg,],2,ise,b0=beta0[nullrg],rng=c(0.5,1))

boxplot(isenull100)
boxplot(isenull500)
boxplot(ise100nn)
boxplot(ise500nn)
boxplot(pmse100)
boxplot(pmse500)







# ************** scenario II: beta2(t) ************** #
betaind = 2
Y = array(NA,c(n,nsim))
X = list()

for(itersim in 1:nsim)
{
  dat = data.generator(n=n, nknots=nknots, norder=norder, domain = c(0, 1), snr=snr, betaind=betaind)
  Y[,itersim]  = dat$Y
  X[[itersim]] = dat$X
}


Ytest = array(NA,c(ntest,nsim))
Xtest = list()

for(itersim in 1:nsim)
{
  dat = data.generator(n=ntest, nknots=nknots, norder=norder, domain = c(0, 1), snr=snr, betaind=betaind)
  Ytest[,itersim]  = dat$Y
  Xtest[[itersim]] = dat$X
}

# FPCR method
rng     = c(0,1)
Mobs    = 500
tobs    = seq(rng[1],rng[2],length.out = Mobs+1)
d       = 3
K       = 30
lambda  = 10^(-10:6)

lam100opt = lam500opt = K100opt = K500opt = pmse100 = pmse500 = rep(NA,nsim)
beta100 = beta500 = array(NA,c(Mobs+1,nsim))

for(itersim in 1:nsim)
{
  # sample size n = 100
  nsample = 100
  x100simcoef = X[[itersim]]$coefs[,1:nsample]
  x100simbasis = X[[itersim]]$basis
  x100sim = fd(coef=x100simcoef,basisobj=x100simbasis)
  FPCRfit = fit.FPCR.BIC(y=Y[1:nsample,itersim],x=x100sim,K=K,Mobs=Mobs,lambda=lambda,rng=rng)
  K100opt[itersim]   = FPCRfit$Kopt
  lam100opt[itersim] = FPCRfit$lamopt
  beta100[,itersim]  = FPCRfit$betaobs
  
  # predict on test set
  x100testsimcoef  = Xtest[[itersim]]$coefs[,1:nsample]
  x100testsimbasis = Xtest[[itersim]]$basis
  x100testsim      = fd(coef=x100testsimcoef,basisobj=x100testsimbasis)
  predict100       = predict.FPCR(y=Ytest[1:nsample,itersim], x=x100testsim, rng=rng, Mobs=Mobs,fpcr_Rfit=FPCRfit)
  pmse100[itersim] = predict100$PMSE
  
  
  
  # sample size n = 500
  nsample = 500
  FPCRfit = fit.FPCR.BIC(y=Y[,itersim],x=X[[itersim]],K=K,Mobs=Mobs,lambda=lambda,rng=rng)
  K500opt[itersim]   = FPCRfit$Kopt
  lam500opt[itersim] = FPCRfit$lamopt
  beta500[,itersim]  = FPCRfit$betaobs
  
  # predict on test set
  predict500       = predict.FPCR(y=Ytest[,itersim], x=Xtest[[itersim]], rng=rng, Mobs=Mobs,fpcr_Rfit=FPCRfit)
  pmse500[itersim] = predict500$PMSE
  
  print(itersim)
}

beta0 = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
plot(tobs, beta100[,1],type="l",ylim=c(-1,1.5))
matlines(tobs,beta100)
lines(tobs,beta0,lwd=2)

plot(tobs, beta500[,1],type="l",ylim=c(-1,1.5))
matlines(tobs,beta500)
lines(tobs,beta0,lwd=2)


# ise
nullrg = 251:501
nnrg   = 1:250

ise100nn   = apply(beta100[nnrg,],2,ise2,b0=beta0[nnrg],rng=c(0,0.5))
isenull100 = apply(beta100[nullrg,],2,ise,b0=beta0[nullrg],rng=c(0.5,1))  
ise500nn   = apply(beta500[nnrg,],2,ise2,b0=beta0[nnrg],rng=c(0,0.5))
isenull500 = apply(beta500[nullrg,],2,ise,b0=beta0[nullrg],rng=c(0.5,1))

boxplot(isenull100)
boxplot(isenull500)
boxplot(ise100nn)
boxplot(ise500nn)
boxplot(pmse100)
boxplot(pmse500)









# ************** scenario III: beta3(t) ************** #
betaind = 3
Y = array(NA,c(n,nsim))
X = list()

for(itersim in 1:nsim)
{
  dat = data.generator(n=n, nknots=nknots, norder=norder, domain = c(0, 1), snr=snr, betaind=betaind)
  Y[,itersim]  = dat$Y
  X[[itersim]] = dat$X
}


Ytest = array(NA,c(ntest,nsim))
Xtest = list()

for(itersim in 1:nsim)
{
  dat = data.generator(n=ntest, nknots=nknots, norder=norder, domain = c(0, 1), snr=snr, betaind=betaind)
  Ytest[,itersim]  = dat$Y
  Xtest[[itersim]] = dat$X
}

# FPCR method
rng     = c(0,1)
Mobs    = 500
tobs    = seq(rng[1],rng[2],length.out = Mobs+1)
d       = 3
K       = 30
lambda  = 10^(-10:6)

lam100opt = lam500opt = K100opt = K500opt = pmse100 = pmse500 = rep(NA,nsim)
beta100 = beta500 = array(NA,c(Mobs+1,nsim))

for(itersim in 1:nsim)
{
  # sample size n = 100
  nsample = 100
  x100simcoef = X[[itersim]]$coefs[,1:nsample]
  x100simbasis = X[[itersim]]$basis
  x100sim = fd(coef=x100simcoef,basisobj=x100simbasis)
  FPCRfit = fit.FPCR.BIC(y=Y[1:nsample,itersim],x=x100sim,K=K,Mobs=Mobs,lambda=lambda,rng=rng)
  K100opt[itersim]   = FPCRfit$Kopt
  lam100opt[itersim] = FPCRfit$lamopt
  beta100[,itersim]  = FPCRfit$betaobs
  
  # predict on test set
  x100testsimcoef  = Xtest[[itersim]]$coefs[,1:nsample]
  x100testsimbasis = Xtest[[itersim]]$basis
  x100testsim      = fd(coef=x100testsimcoef,basisobj=x100testsimbasis)
  predict100       = predict.FPCR(y=Ytest[1:nsample,itersim], x=x100testsim, rng=rng, Mobs=Mobs,fpcr_Rfit=FPCRfit)
  pmse100[itersim] = predict100$PMSE
  
  
  
  # sample size n = 500
  nsample = 500
  FPCRfit = fit.FPCR.BIC(y=Y[,itersim],x=X[[itersim]],K=K,Mobs=Mobs,lambda=lambda,rng=rng)
  K500opt[itersim]   = FPCRfit$Kopt
  lam500opt[itersim] = FPCRfit$lamopt
  beta500[,itersim]  = FPCRfit$betaobs
  
  # predict on test set
  predict500       = predict.FPCR(y=Ytest[,itersim], x=Xtest[[itersim]], rng=rng, Mobs=Mobs,fpcr_Rfit=FPCRfit)
  pmse500[itersim] = predict500$PMSE
  
  print(itersim)
}

beta0 = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
plot(tobs, beta100[,1],type="l",ylim=c(-1,5.5))
matlines(tobs,beta100)
lines(tobs,beta0,lwd=2)

plot(tobs, beta500[,1],type="l",ylim=c(-1,5.5))
matlines(tobs,beta500)
lines(tobs,beta0,lwd=2)


# ise
ise100nn   = apply(beta100,2,ise2,b0=beta0,rng=c(0,1))
ise500nn   = apply(beta500,2,ise2,b0=beta0,rng=c(0,1))


boxplot(ise100nn)
boxplot(ise500nn)
boxplot(pmse100)
boxplot(pmse500)