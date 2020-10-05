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


# slos estimator
# B-splines basis used by slos method
n       = 500
ntest   = 5000
nsim    = 100
Mobs    = 500
domain  = c(0,1)
tobs    = seq(domain[1],domain[2],length.out = Mobs+1)
h = (domain[2]-domain[1])/Mobs
cef = c(1, rep(c(4,2), (Mobs-2)/2), 4, 1)

# B-splines basis used by slos method
M = 100
d = 3
norder   = d+1
nknots   = M+1
knots    = seq(domain[1],domain[2],length.out = nknots)
nbasis   = nknots + norder - 2
basis    = create.bspline.basis(knots,nbasis,norder)
basismat = eval.basis(tobs, basis) # 101 103
V = eval.penalty(basis,int2Lfd(2))


Xs     = array(NA, c(n,Mobs+1,nsim))
Xstest = array(NA, c(ntest,Mobs+1,nsim))
ustest = array(NA,c(ntest,nbasis,nsim))
for(itersim in 1:nsim)
{
  ptm = proc.time()
  Xs[,,itersim]     = t(eval.fd(tobs,X[[itersim]]))
  Xstest[,,itersim] = t(eval.fd(tobs,Xtest[[itersim]]))
  ustest[,,itersim] = h/3*Xstest[,,itersim]%*%diag(cef)%*%basismat
  timeint = as.numeric((proc.time() - ptm)[3])
  print(paste("n =", itersim, "completed", timeint))
}




lambda  = exp(seq(-20,-11, length.out = 10))
gamma   = 10^(-8:-6)
Maxiter = 100
absTol  = 10^(-10)
Cutoff  = 10^(-6)
beta100 = beta500 = array(NA,c(Mobs+1,nsim))
b100    = b500    = array(NA,c(M+d,nsim))
pmse100 = pmse500 = rep(NA, nsim)
Optgamma100 = Optgamma500 = Optlambda100 = Optlambda500 = rep(NA,nsim)

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  slosfit = slos.fit(Y=Y[1:n,itersim],t(Xs[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff,tobs=tobs)
  beta100[,itersim] = slosfit$beta
  b100[,itersim]    = slosfit$b
  Optgamma100[itersim]  = slosfit$Optgamma
  Optlambda100[itersim] = slosfit$Optlambda
  y100testhat  = ustest[,,itersim]%*%b100[,itersim]    
  pmse100[itersim] = t(y100testhat-Ytest[,itersim])%*%(y100testhat-Ytest[,itersim])/ntest
  
  # sample size n = 500
  n=500
  slosfit = slos.fit(Y=Y[1:n,itersim],t(Xs[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff,tobs=tobs)
  beta500[,itersim] = slosfit$beta
  b500[,itersim]    = slosfit$b
  Optgamma500[itersim]  = slosfit$Optgamma
  Optlambda500[itersim] = slosfit$Optlambda
  y500testhat  = ustest[,,itersim]%*%b500[,itersim]
  pmse500[itersim] = t(y500testhat-Ytest[,itersim])%*%(y500testhat-Ytest[,itersim])/ntest
  
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


# slos estimator
# B-splines basis used by slos method
n       = 500
ntest   = 5000
nsim    = 100
Mobs    = 500
domain  = c(0,1)
tobs    = seq(domain[1],domain[2],length.out = Mobs+1)
h = (domain[2]-domain[1])/Mobs
cef = c(1, rep(c(4,2), (Mobs-2)/2), 4, 1)

# B-splines basis used by slos method
M = 100
d = 3
norder   = d+1
nknots   = M+1
knots    = seq(domain[1],domain[2],length.out = nknots)
nbasis   = nknots + norder - 2
basis    = create.bspline.basis(knots,nbasis,norder)
basismat = eval.basis(tobs, basis) # 101 103
V = eval.penalty(basis,int2Lfd(2))


Xs     = array(NA, c(n,Mobs+1,nsim))
Xstest = array(NA, c(ntest,Mobs+1,nsim))
ustest = array(NA,c(ntest,nbasis,nsim))
for(itersim in 1:nsim)
{
  ptm = proc.time()
  Xs[,,itersim]     = t(eval.fd(tobs,X[[itersim]]))
  Xstest[,,itersim] = t(eval.fd(tobs,Xtest[[itersim]]))
  ustest[,,itersim] = h/3*Xstest[,,itersim]%*%diag(cef)%*%basismat
  timeint = as.numeric((proc.time() - ptm)[3])
  print(paste("n =", itersim, "completed", timeint))
}




lambda  = exp(seq(-20,-12, length.out = 10))
gamma   = 10^(-8:-6)
Maxiter = 100
absTol  = 10^(-10)
Cutoff  = 10^(-6)
beta100 = beta500 = array(NA,c(Mobs+1,nsim))
b100    = b500    = array(NA,c(M+d,nsim))
pmse100 = pmse500 = rep(NA, nsim)
Optgamma100 = Optgamma500 = Optlambda100 = Optlambda500 = rep(NA,nsim)

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  slosfit = slos.fit(Y=Y[1:n,itersim],t(Xs[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff,tobs=tobs)
  beta100[,itersim] = slosfit$beta
  b100[,itersim]    = slosfit$b
  Optgamma100[itersim]  = slosfit$Optgamma
  Optlambda100[itersim] = slosfit$Optlambda
  y100testhat  = ustest[,,itersim]%*%b100[,itersim]    
  pmse100[itersim] = t(y100testhat-Ytest[,itersim])%*%(y100testhat-Ytest[,itersim])/ntest
  
  # sample size n = 500
  n=500
  slosfit = slos.fit(Y=Y[1:n,itersim],t(Xs[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff,tobs=tobs)
  beta500[,itersim] = slosfit$beta
  b500[,itersim]    = slosfit$b
  Optgamma500[itersim]  = slosfit$Optgamma
  Optlambda500[itersim] = slosfit$Optlambda
  y500testhat  = ustest[,,itersim]%*%b500[,itersim]
  pmse500[itersim] = t(y500testhat-Ytest[,itersim])%*%(y500testhat-Ytest[,itersim])/ntest
  
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


# slos estimator
# B-splines basis used by slos method
n       = 500
ntest   = 5000
nsim    = 100
Mobs    = 500
domain  = c(0,1)
tobs    = seq(domain[1],domain[2],length.out = Mobs+1)
h = (domain[2]-domain[1])/Mobs
cef = c(1, rep(c(4,2), (Mobs-2)/2), 4, 1)

# B-splines basis used by slos method
M = 100
d = 3
norder   = d+1
nknots   = M+1
knots    = seq(domain[1],domain[2],length.out = nknots)
nbasis   = nknots + norder - 2
basis    = create.bspline.basis(knots,nbasis,norder)
basismat = eval.basis(tobs, basis) # 101 103
V = eval.penalty(basis,int2Lfd(2))


Xs     = array(NA, c(n,Mobs+1,nsim))
Xstest = array(NA, c(ntest,Mobs+1,nsim))
ustest = array(NA,c(ntest,nbasis,nsim))
for(itersim in 1:nsim)
{
  ptm = proc.time()
  Xs[,,itersim]     = t(eval.fd(tobs,X[[itersim]]))
  Xstest[,,itersim] = t(eval.fd(tobs,Xtest[[itersim]]))
  ustest[,,itersim] = h/3*Xstest[,,itersim]%*%diag(cef)%*%basismat
  timeint = as.numeric((proc.time() - ptm)[3])
  print(paste("n =", itersim, "completed", timeint))
}




lambda  = exp(seq(-60,-50, length.out = 2))
gamma   = 10^(-8:-6)
Maxiter = 100
absTol  = 10^(-10)
Cutoff  = 10^(-6)
beta100 = beta500 = array(NA,c(Mobs+1,nsim))
b100    = b500    = array(NA,c(M+d,nsim))
pmse100 = pmse500 = rep(NA, nsim)
Optgamma100 = Optgamma500 = Optlambda100 = Optlambda500 = rep(NA,nsim)

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  slosfit = slos.fit(Y=Y[1:n,itersim],t(Xs[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff,tobs=tobs)
  beta100[,itersim] = slosfit$beta
  b100[,itersim]    = slosfit$b
  Optgamma100[itersim]  = slosfit$Optgamma
  Optlambda100[itersim] = slosfit$Optlambda
  y100testhat  = ustest[,,itersim]%*%b100[,itersim]    
  pmse100[itersim] = t(y100testhat-Ytest[,itersim])%*%(y100testhat-Ytest[,itersim])/ntest
  
  # sample size n = 500
  n=500
  slosfit = slos.fit(Y=Y[1:n,itersim],t(Xs[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff,tobs=tobs)
  beta500[,itersim] = slosfit$beta
  b500[,itersim]    = slosfit$b
  Optgamma500[itersim]  = slosfit$Optgamma
  Optlambda500[itersim] = slosfit$Optlambda
  y500testhat  = ustest[,,itersim]%*%b500[,itersim]
  pmse500[itersim] = t(y500testhat-Ytest[,itersim])%*%(y500testhat-Ytest[,itersim])/ntest
  
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
ise100nn   = apply(beta100,2,ise,b0=beta0,rng=c(0,1))
ise500nn   = apply(beta500,2,ise,b0=beta0,rng=c(0,1))


boxplot(ise100nn)
boxplot(ise500nn)
boxplot(pmse100)
boxplot(pmse500)