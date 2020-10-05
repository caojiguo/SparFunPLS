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
data = as.matrix(read.table('Truck.csv',sep=',',header=FALSE))

# Scale Y up by 1000
Y = data[,1]*1e3
Y = Y - mean(Y) # centered Y
# Covariates are everything but the first column
Z = t(data[,-1])

# Basis expansion
bbasis = create.bspline.basis(c(0,90),norder=4,breaks = seq(0,90,by=3))

# Smooth for covariates, centered
zfd = smooth.basis((0:89)+0.5,Z,fdPar(bbasis,2))
zct = center.fd(zfd$fd)

tt = seq(0,90,by=0.25)
X  = -eval.fd(tt,zct,1)
nsec  = 241

X = X[1:nsec,]



# Plot accelrations
set.seed(22)
set1 = sample(1:108, 10, replace = F)
tt = seq(0, 60, length.out = nsec)
par(mar=c(4,4.5,2,1))
matplot(tt,X[,set1],type='l',lty=1,col=1,xlab='Second',ylab='Acceleration', main = "(a)")
plot(tt, X[,set1[1]],type = "l", ylim = c(-5, 4), yaxt='n',
     xlab='Second',ylab='Acceleration', main = "(a)")
axis(2, at=seq(-6,4,by=2), labels=seq(-6,4,by=2),las=2)
matlines(tt,X[,set1], lty = 1, col =1)

# plot all acceleration curves
matplot(tt,X,type='l',lty=1,col=1,xlab='Second',yaxt='n',
        ylim=c(-8,4),ylab='Acceleration', main = "(a)")
axis(2, at=seq(-8,4,by=2), labels=seq(-8,4,by=2),las=2)





#--------------- FPLS method ---------------#
n       = length(Y) # 108
domain  = c(0,60)
Mobs    = length(tt)-1
M       = 3:20
d       = 3
K       = 15     

FunPLSfit = fit.FunPLS(y=Y, xind = "matrix", x=t(X), domain=domain, M=M, d=d, K=K, Mobs=Mobs,mod.select = "BIC")
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
lambda  = exp(seq(5,9,length.out = 10))
delta   = 10^(1)
gamma   = 10^(9)
Maxiter = 1000
absTol  = 10^(-15)
Cutoff  = 10^(-8)
k       = 8

SpaFunPLSfit = SpaFunPLS.tune.BIC(Y=Y,xind="matrix", X=t(X), k=k, basis=beta.basis.spa, W=W, M.FunPLS=Mopt,kappa=kappa, d=d,lambda=lambda, delta=delta, gamma=gamma, domain=domain)
lam.spa   = SpaFunPLSfit$lamopt
Kopt.spa  = SpaFunPLSfit$Kopt
beta.spa  = SpaFunPLSfit$beta


par(mar=c(4,4.5,2,1))
plot(tt,beta.spa, yaxt='n', type="l",lwd=2,xlab='Second', ylim=c(-0.4,1.2),ylab=expression(beta(t)), main = "(b)")
axis(2, at=seq(-0.4,1.2,by=0.4), labels=seq(-0.4,1.2,by=0.4),las=2)
lines(tt,beta.FunPLS,lty=2,col="grey",lwd=2)
abline(h=0)