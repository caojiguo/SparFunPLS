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
domain  = c(0,1)


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

# FPLS method
M       = 2:15
Mobs    = 500
tobs    = seq(domain[1], domain[2], length.out = Mobs+1)
d       = 3
K       = 15 

M100opt = M500opt = rep(NA,nsim)

# SFPLS method
M.spa    = 100
d.spa    = 3
norder.spa   = d+1
knots.spa    = seq(0,1, length.out=M.spa+1)
nknots.spa   = length(knots.spa)
nbasis.spa   = length(knots.spa) + norder.spa - 2 # i+2
basis.spa    = create.bspline.basis(knots.spa,nbasis.spa,norder.spa)
beta.basis.spa = basis.spa
W        = slos.compute.weights(beta.basis.spa)


kappa   = 0.1
lambda  = exp(seq(-11,-6,length.out = 15))
lambda500 = exp(seq(-10,-5,length.out = 10))
delta   = 10^(-1)
gamma   = 10^(seq(-5,-4,length.out = 3))
Maxiter = 100
absTol  = 10^(-15)
Cutoff  = 10^(-8)
k       = 5

lam100opt.spa = lam500opt.spa = gam100opt.spa = gam500opt.spa = K100opt.spa = K500opt.spa = pmse100.spa = pmse500.spa = rep(NA,nsim)
beta100.spa = beta500.spa = array(NA,c(Mobs+1,nsim))


for(itersim in 1:nsim)
{
  # sample size n = 100
  nsample = 100
  x100simcoef = X[[itersim]]$coefs[,1:nsample]
  x100simbasis = X[[itersim]]$basis
  x100sim = fd(coef=x100simcoef,basisobj=x100simbasis)
  
  # FPLS fit
  FunPLSfit = fit.FunPLS(y=Y[1:nsample,itersim],xind="fd", x=x100sim, domain=domain, M=M, d=norder-1, K=K, Mobs=Mobs,mod.select="BIC")
  M100opt[itersim] = FunPLSfit$Mopt
  
  # SFPLS fit
  SpaFunPLSfit = SpaFunPLS.tune.BIC(Y=Y[1:nsample,itersim],xind="fd", X=x100sim, k=k, basis=beta.basis.spa, W=W, M.FunPLS=M100opt[itersim],kappa=kappa, d=norder-1,lambda=lambda, delta=delta, gamma=gamma, domain=domain)
  lam100opt.spa[itersim]  = SpaFunPLSfit$lamopt
  gam100opt.spa[itersim]  = SpaFunPLSfit$gamopt
  K100opt.spa[itersim]  = SpaFunPLSfit$Kopt
  beta100.spa[,itersim] = SpaFunPLSfit$beta
  
  # SFPLS predict on test set
  x100testsimcoef  = Xtest[[itersim]]$coefs[,1:nsample]
  x100testsimbasis = Xtest[[itersim]]$basis
  x100testsim      = fd(coef=x100testsimcoef,basisobj=x100testsimbasis)
  
  predict100       = predict.SpaFunPLS(y=Ytest[1:nsample,itersim], x=x100testsim, domain=domain, Mobs=Mobs,sfplsfit.beta = beta100.spa[,itersim])
  pmse100.spa[itersim] = predict100$PMSE
  
 
  
  # sample size n = 500
  nsample = 500
 
  # FPLS fit
  FunPLSfit = fit.FunPLS(y=Y[,itersim],xind="fd", x=X[[itersim]], domain=domain, M=M, d=norder-1, K=K, Mobs=Mobs,mod.select="BIC")
  M500opt[itersim] = FunPLSfit$Mopt
  
  # SFPLS fit
  SpaFunPLSfit = SpaFunPLS.tune.BIC(Y=Y[,itersim],xind="fd", X=X[[itersim]], k=k, basis=beta.basis.spa, W=W, M.FunPLS=M500opt[itersim],kappa=kappa, d=norder-1,lambda=lambda500, delta=delta, gamma=gamma, domain=domain)
  lam500opt.spa[itersim]  = SpaFunPLSfit$lamopt
  gam500opt.spa[itersim]  = SpaFunPLSfit$gamopt
  K500opt.spa[itersim]  = SpaFunPLSfit$Kopt
  beta500.spa[,itersim] = SpaFunPLSfit$beta
  
  # SFPLS predict on test set
  predict500       = predict.SpaFunPLS(y=Ytest[,itersim], x=Xtest[[itersim]], domain=domain, Mobs=Mobs,sfplsfit.beta = beta500.spa[,itersim])
  pmse500.spa[itersim] = predict500$PMSE
  
  print(itersim)
}

beta0 = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
plot(tobs, beta100.spa[,1],type="l",ylim=c(-1,1.5))
matlines(tobs,beta100.spa)
lines(tobs,beta0,lwd=2)

plot(tobs, beta500.spa[,1],type="l",ylim=c(-1,1.5))
matlines(tobs,beta500.spa)
lines(tobs,beta0,lwd=2)


# ise
nullrg = 251:501
nnrg   = 1:250

ise100nn   = apply(beta100.spa[nnrg,],2,ise2,b0=beta0[nnrg],rng=c(0,1))
isenull100 = apply(beta100.spa[nullrg,],2,ise,b0=beta0[nullrg],rng=c(0.5,1))  
ise500nn   = apply(beta500.spa[nnrg,],2,ise2,b0=beta0[nnrg],rng=c(0,1))
isenull500 = apply(beta500.spa[nullrg,],2,ise,b0=beta0[nullrg],rng=c(0.5,1))


boxplot(isenull100)
boxplot(isenull500)
boxplot(ise100nn)
boxplot(ise500nn)
boxplot(pmse100.spa)
boxplot(pmse500.spa)











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

# FPLS method
M       = 3:10
Mobs    = 500
tobs    = seq(domain[1], domain[2], length.out = Mobs+1)
d       = 3
K       = 10 

M100opt = M500opt = rep(NA,nsim)

# SFPLS method
M.spa    = 100
d.spa    = 3
norder.spa   = d+1
knots.spa    = seq(0,1, length.out=M.spa+1)
nknots.spa   = length(knots.spa)
nbasis.spa   = length(knots.spa) + norder.spa - 2 # i+2
basis.spa    = create.bspline.basis(knots.spa,nbasis.spa,norder.spa)
beta.basis.spa = basis.spa
W        = slos.compute.weights(beta.basis.spa)


kappa   = 0.1
delta   = 10^(-1)
lambda  = exp(seq(-12.5,-7,length.out = 10))
lambda500 = exp(seq(-10,-5,length.out = 10))
gamma   = 10^(seq(-6,-4,length.out = 3))
Maxiter = 100
absTol  = 10^(-15)
Cutoff  = 10^(-8)
k       = 5

lam100opt.spa = lam500opt.spa = gam100opt.spa = gam500opt.spa = K100opt.spa = K500opt.spa = pmse100.spa = pmse500.spa = rep(NA,nsim)
beta100.spa = beta500.spa = array(NA,c(Mobs+1,nsim))


for(itersim in 1:nsim)
{
  # sample size n = 100
  nsample = 100
  x100simcoef = X[[itersim]]$coefs[,1:nsample]
  x100simbasis = X[[itersim]]$basis
  x100sim = fd(coef=x100simcoef,basisobj=x100simbasis)
  
  # FPLS fit
  FunPLSfit = fit.FunPLS(y=Y[1:nsample,itersim],xind="fd", x=x100sim, domain=domain, M=M, d=norder-1, K=K, Mobs=Mobs,mod.select="BIC")
  M100opt[itersim] = FunPLSfit$Mopt
  
  # SFPLS fit
  SpaFunPLSfit = SpaFunPLS.tune.BIC(Y=Y[1:nsample,itersim],xind="fd", X=x100sim, k=k, basis=beta.basis.spa, W=W, M.FunPLS=M100opt[itersim],kappa=kappa, d=norder-1,lambda=lambda, delta=delta, gamma=gamma, domain=domain)
  lam100opt.spa[itersim]  = SpaFunPLSfit$lamopt
  gam100opt.spa[itersim]  = SpaFunPLSfit$gamopt
  K100opt.spa[itersim]  = SpaFunPLSfit$Kopt
  beta100.spa[,itersim] = SpaFunPLSfit$beta
  
  # SFPLS predict on test set
  x100testsimcoef  = Xtest[[itersim]]$coefs[,1:nsample]
  x100testsimbasis = Xtest[[itersim]]$basis
  x100testsim      = fd(coef=x100testsimcoef,basisobj=x100testsimbasis)
  
  predict100       = predict.SpaFunPLS(y=Ytest[1:nsample,itersim], x=x100testsim, domain=domain, Mobs=Mobs,sfplsfit.beta = beta100.spa[,itersim])
  pmse100.spa[itersim] = predict100$PMSE
  
  
  
  # sample size n = 500
  nsample = 500
  
  # FPLS fit
  FunPLSfit = fit.FunPLS(y=Y[,itersim],xind="fd", x=X[[itersim]], domain=domain, M=M, d=norder-1, K=K, Mobs=Mobs,mod.select="BIC")
  M500opt[itersim] = FunPLSfit$Mopt
  
  # SFPLS fit
  SpaFunPLSfit = SpaFunPLS.tune.BIC(Y=Y[,itersim],xind="fd", X=X[[itersim]], k=k, basis=beta.basis.spa, W=W, M.FunPLS=M500opt[itersim],kappa=kappa, d=norder-1,lambda=lambda500, delta=delta, gamma=gamma, domain=domain)
  lam500opt.spa[itersim]  = SpaFunPLSfit$lamopt
  gam500opt.spa[itersim]  = SpaFunPLSfit$gamopt
  K500opt.spa[itersim]  = SpaFunPLSfit$Kopt
  beta500.spa[,itersim] = SpaFunPLSfit$beta
  
  # SFPLS predict on test set
  predict500       = predict.SpaFunPLS(y=Ytest[,itersim], x=Xtest[[itersim]], domain=domain, Mobs=Mobs,sfplsfit.beta = beta500.spa[,itersim])
  pmse500.spa[itersim] = predict500$PMSE
  
  print(itersim)
}

beta0 = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
plot(tobs, beta100.spa[,1],type="l",ylim=c(-1,1.5))
matlines(tobs,beta100.spa)
lines(tobs,beta0,lwd=2)

plot(tobs, beta500.spa[,1],type="l",ylim=c(-1,1.5))
matlines(tobs,beta500.spa)
lines(tobs,beta0,lwd=2)


# ise
nullrg = 251:501
nnrg   = 1:250

ise100nn   = apply(beta100.spa[nnrg,],2,ise2,b0=beta0[nnrg],rng=c(0,1))
isenull100 = apply(beta100.spa[nullrg,],2,ise,b0=beta0[nullrg],rng=c(0.5,1))  
ise500nn   = apply(beta500.spa[nnrg,],2,ise2,b0=beta0[nnrg],rng=c(0,1))
isenull500 = apply(beta500.spa[nullrg,],2,ise,b0=beta0[nullrg],rng=c(0.5,1))

boxplot(isenull100)
boxplot(isenull500)
boxplot(ise100nn)
boxplot(ise500nn)
boxplot(pmse100.spa)
boxplot(pmse500.spa)











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

# FPLS method
M       = seq(3,15,by=3)
Mobs    = 500
tobs    = seq(domain[1], domain[2], length.out = Mobs+1)
d       = 3
K       = 10 

M100opt = M500opt = rep(NA,nsim)

# SFPLS method
M.spa    = 100
d.spa    = 3
norder.spa   = d+1
knots.spa    = seq(0,1, length.out=M.spa+1)
nknots.spa   = length(knots.spa)
nbasis.spa   = length(knots.spa) + norder.spa - 2 # i+2
basis.spa    = create.bspline.basis(knots.spa,nbasis.spa,norder.spa)
beta.basis.spa = basis.spa
W        = slos.compute.weights(beta.basis.spa)


kappa   = 0.1
delta   = 10^(-1)
lambda  = exp(seq(-100,-10,length.out = 3))
lambda500 = exp(seq(-100,-10,length.out = 3))
gamma   = 10^(seq(-5,-4,length.out = 3))
Maxiter = 100
absTol  = 10^(-15)
Cutoff  = 10^(-8)
k       = 15

lam100opt.spa = lam500opt.spa = gam100opt.spa = gam500opt.spa = K100opt.spa = K500opt.spa = pmse100.spa = pmse500.spa = rep(NA,nsim)
beta100.spa = beta500.spa = array(NA,c(Mobs+1,nsim))


for(itersim in 1:nsim)
{
  # sample size n = 100
  nsample = 100
  x100simcoef = X[[itersim]]$coefs[,1:nsample]
  x100simbasis = X[[itersim]]$basis
  x100sim = fd(coef=x100simcoef,basisobj=x100simbasis)
  
  # FPLS fit
  FunPLSfit = fit.FunPLS(y=Y[1:nsample,itersim],xind="fd", x=x100sim, domain=domain, M=M, d=norder-1, K=K, Mobs=Mobs,mod.select="BIC")
  M100opt[itersim] = FunPLSfit$Mopt
  
  # SFPLS fit
  SpaFunPLSfit = SpaFunPLS.tune.BIC(Y=Y[1:nsample,itersim],xind="fd", X=x100sim, k=k, basis=beta.basis.spa, W=W, M.FunPLS=M100opt[itersim],kappa=kappa, d=norder-1,lambda=lambda, delta=delta, gamma=gamma, domain=domain)
  lam100opt.spa[itersim]  = SpaFunPLSfit$lamopt
  gam100opt.spa[itersim]  = SpaFunPLSfit$gamopt
  K100opt.spa[itersim]  = SpaFunPLSfit$Kopt
  beta100.spa[,itersim] = SpaFunPLSfit$beta
  
  # SFPLS predict on test set
  x100testsimcoef  = Xtest[[itersim]]$coefs[,1:nsample]
  x100testsimbasis = Xtest[[itersim]]$basis
  x100testsim      = fd(coef=x100testsimcoef,basisobj=x100testsimbasis)
  
  predict100       = predict.SpaFunPLS(y=Ytest[1:nsample,itersim], x=x100testsim, domain=domain, Mobs=Mobs,sfplsfit.beta = beta100.spa[,itersim])
  pmse100.spa[itersim] = predict100$PMSE
  
  
  
  # sample size n = 500
  nsample = 500
  
  # FPLS fit
  FunPLSfit = fit.FunPLS(y=Y[,itersim],xind="fd", x=X[[itersim]], domain=domain, M=M, d=norder-1, K=K, Mobs=Mobs,mod.select="BIC")
  M500opt[itersim] = FunPLSfit$Mopt
  
  # SFPLS fit
  SpaFunPLSfit = SpaFunPLS.tune.BIC(Y=Y[,itersim],xind="fd", X=X[[itersim]], k=k, basis=beta.basis.spa, W=W, M.FunPLS=M500opt[itersim],kappa=kappa, d=norder-1,lambda=lambda500, delta=delta, gamma=gamma, domain=domain)
  lam500opt.spa[itersim]  = SpaFunPLSfit$lamopt
  gam500opt.spa[itersim]  = SpaFunPLSfit$gamopt
  K500opt.spa[itersim]  = SpaFunPLSfit$Kopt
  beta500.spa[,itersim] = SpaFunPLSfit$beta
  
  # SFPLS predict on test set
  predict500       = predict.SpaFunPLS(y=Ytest[,itersim], x=Xtest[[itersim]], domain=domain, Mobs=Mobs,sfplsfit.beta = beta500.spa[,itersim])
  pmse500.spa[itersim] = predict500$PMSE
  
  print(itersim)
}

beta0 = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
plot(tobs, beta100.spa[,1],type="l",ylim=c(-1,5.5))
matlines(tobs,beta100.spa)
lines(tobs,beta0,lwd=2)

plot(tobs, beta500.spa[,1],type="l",ylim=c(-1,5.5))
matlines(tobs,beta500.spa)
lines(tobs,beta0,lwd=2)


# ise
ise100  = apply(beta100.spa,2,ise,b0=beta0,rng=c(0,1))
ise500  = apply(beta500.spa,2,ise,b0=beta0,rng=c(0,1))


boxplot(ise100nn)
boxplot(ise500nn)
boxplot(pmse100.spa)
boxplot(pmse500.spa)