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




# Data
data = as.matrix(read.table('OSB.csv',sep=',',header=TRUE))[,-1]

Y = data[,1]
Y = Y - mean(Y)

spectra = data[,-1]

n      = length(Y)  # sample size
wave   = 350:2500
domain = c(wave[1],wave[length(wave)])

# random process X matrix to fd object
bbasis    = create.bspline.basis(c(domain[1],domain[2]), nbasis = 50, norder = 4)
fdParzt   = fdPar(bbasis, 2, 100)
zfd       = smooth.basis(domain[1]:domain[2], y = t(spectra), fdParzt)$fd
zcd       = center.fd(zfd)

waveobs = seq(domain[1],domain[2],length.out = 501)
X0 = t(eval.fd(waveobs, zfd))
X = t(eval.fd(waveobs, zcd))


par(mar=c(4.5,4.5,2,1))
# different proportions of sound wood
plot(wave, spectra[87,],type = "l", ylim = c(0,1.4), yaxt='n',lwd=1,cex.lab=1.6,cex.axis=1.2,cex.main=1.5,
     xlab='Wavelength (nm)',ylab='Spectra',main="(b)") # 100% sound
axis(2, at=seq(0,1.4,by=0.2), labels=seq(0,1.4,by=0.2),las=2,cex.axis=1.2)
lines(wave,spectra[77,],lty=2,col=2,lwd=1) # 80% sound
lines(wave,spectra[59,],lty=3,col=4,lwd=1) # 40% sound 
lines(wave,spectra[3,],lty=4,col="brown",lwd=1) # 0% sound 


# Plot smoothed centered spectra curves
set.seed(333)
seto = seq(1,length(Y),by=3)
set1 = sample(seto, 10, replace = F)
plot(waveobs, t(X[set1[1],]),type = "l", 
     ylim = c(-0.3,0.6), yaxt='n',cex.lab=1.2,cex.axis=1.2,cex.main=1.2,
     xlab='Wavelength',ylab='Spectra', main = "(a)")
axis(2, at=-3:6/10, labels=-3:6/10,las=2,cex.axis=1.2)
matlines(waveobs,t(X[set1,]), lty = 1, col =1)





#--------------- FPLS method ---------------#
M       = 3:20
d       = 3
K       = 10
Mobs    = dim(X)[2]-1

FunPLSfit = fit.FunPLS(y=Y, xind = "matrix", x=X, domain=domain, M=M, d=d, K=K, Mobs=Mobs,mod.select = "CV",ncv=5)
Mopt = FunPLSfit$Mopt
Kopt = FunPLSfit$Kopt
beta.FunPLS = FunPLSfit$betahatobs


#--------------- SFPLS method ---------------#
M.spa    = 100
d.spa    = 3
norder.spa   = d+1
knots.spa    = seq(domain[1],domain[2], length.out=M.spa+1)
nknots.spa   = length(knots.spa)
nbasis.spa   = length(knots.spa) + norder.spa - 2 # i+2
basis.spa    = create.bspline.basis(knots.spa,nbasis.spa,norder.spa)
beta.basis.spa = basis.spa
W        = slos.compute.weights(beta.basis.spa)

kappa   = 0.1
lambda  = exp(seq(-4,8,length.out = 10))
delta   = 10^(-1)
gamma   = 10^(10)
Maxiter = 1000
absTol  = 10^(-15)
Cutoff  = 10^(-8)
k       = 3 # the number of components chosen by the FunPLS method is 3, so we set k is no larger than 3

SpaFunPLSfit = SpaFunPLS.tune.BIC(Y=Y,xind="matrix", X=X, k=k, basis=beta.basis.spa, W=W, M.FunPLS=Mopt,kappa=kappa, d=d,lambda=lambda, delta=delta, gamma=gamma, domain=domain)
lam.spa   = SpaFunPLSfit$lamopt
Kopt.spa  = SpaFunPLSfit$Kopt
beta.spa  = SpaFunPLSfit$beta

par(mar=c(4,4.5,2,1))
plot(waveobs,beta.spa, yaxt='n', type="l",lwd=2,xlab='Second', ylim=c(-0.02,0.07),ylab=expression(beta(t)), main = "(b)")
axis(2, at=seq(-0.04,0.08,by=0.02), labels=seq(-0.04,0.08,by=0.02),las=2)
lines(waveobs,beta.FunPLS,lty=2,col="grey",lwd=2)
abline(h=0)