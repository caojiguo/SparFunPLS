rm(list=ls())

# Please change the directory of this file in your computer
setwd("~/Desktop/SFPLS_Supplementary_Materials_Oct2020/SFPLS")

install.packages("psych")
install.packages("pls")
install.packages("plsdof")
install.packages("fda")
install.packages("MASS")
install.packages("glmnet")
install.packages("nlme")


library(psych)
library(pls)
library(plsdof)
library(fda)
library(MASS)
library(glmnet)
library(nlme)
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




#--------------- FPLS method ---------------#
Mweight = seq(3,30,by=3)
Mobs    = 500
dsim    = 3
Acomp   = 5
Doft    = array(NA,c(Acomp,length(Mweight),nsim))

tobs     = seq(domain[1],domain[2],length.out = Mobs+1)
h        = (domain[2]-domain[1])/Mobs
cef      = c(1, rep(c(4,2), (Mobs-2)/2), 4, 1)

BICpls100 = BICpls500 = array(NA,c(Acomp,length(Mweight),nsim))
Mfpls100  = Afpls100 =  pmse100  = Mfpls500  = Afpls500 =  pmse500  = array(NA,c(Acomp,nsim))
beta_fpls_hat100 = beta_fpls_hat500 = array(NA,c(length(tobs),Acomp,nsim))
etahat = list()

for(itersim in 1:nsim)
{
  # n = 100
  nsample = 100
  y100sim      = Y[1:nsample,itersim] 
  x100simcoef  = X[[itersim]]$coefs[,1:nsample]
  x100simbasis = X[[itersim]]$basis
  x100sim      = fd(coef=x100simcoef,basisobj=x100simbasis)
  x100obs      = t(eval.fd(tobs,x100sim))

  # n = 500
  nsample = 500
  y500sim      = Y[1:nsample,itersim] 
  x500simcoef  = X[[itersim]]$coefs[,1:nsample]
  x500simbasis = X[[itersim]]$basis
  x500sim      = fd(coef=x500simcoef,basisobj=x500simbasis)
  x500obs      = t(eval.fd(tobs,x500sim))
  
  ytest     = Ytest[,itersim] 
  xtest     = Xtest[[itersim]]
  
  for(iterM in 1:length(Mweight))
  {
    norder   = dsim+1
    nknots   = Mweight[iterM]+1
    knots    = seq(domain[1],domain[2], length.out = nknots)
    nbasis   = nknots + norder - 2 
    basis.weight = create.bspline.basis(knots,nbasis,norder)
    basis    = basis.weight
    
    basismat = eval.basis(tobs,basis)
    u100     = h/3*x100obs%*%diag(cef)%*%basismat
    u500     = h/3*x500obs%*%diag(cef)%*%basismat
    
    Aind = which(1:Acomp<=Mweight[iterM]+dsim)
    
    fitpls100  = linear.pls.fit(X=u100,y=y100sim,m=Acomp,compute.jacobian = T)
    fitpls500  = linear.pls.fit(X=u500,y=y500sim,m=Acomp,compute.jacobian = T)
    
    for(iterk in Aind)
    {
      # n = 100
      RSS     = fitpls100$RSS[iterk+1]
      Dof     = fitpls100$DoF[iterk+1]
      sigma   = fitpls100$sigmahat[iterk+1]
      BICpls100[iterk,iterM,itersim]  = RSS+log(n)*sigma^2*Dof
      
      # n = 500
      RSS     = fitpls500$RSS[iterk+1]
      Dof     = fitpls500$DoF[iterk+1]
      sigma   = fitpls500$sigmahat[iterk+1]
      BICpls500[iterk,iterM,itersim]  = RSS+log(n)*sigma^2*Dof
    }
  }
  
  for(iterk in 1:Acomp)
  {
    # n = 100
    BICind100  = which(BICpls100[iterk,,itersim]==min(BICpls100[iterk,,itersim],na.rm = TRUE),arr.ind = TRUE)
    Afpls100[iterk,itersim]   = iterk
    Mfpls100[iterk,itersim]   = Mweight[BICind100]
    
    norder   = dsim+1
    nknots   = Mfpls100[iterk,itersim]+1
    knots    = seq(domain[1],domain[2], length.out = nknots)
    nbasis   = nknots + norder - 2 
    basis.weight = create.bspline.basis(knots,nbasis,norder)
    basis    = basis.weight
    basismat = eval.basis(tobs,basis)
    u100sim  = h/3*x100obs%*%diag(cef)%*%basismat
    
    fplsfit = convFPLS.mat(y=y100sim,u=u100sim,basis.weight=basis.weight,A=Afpls100[iterk,itersim])
    etahat  = fplsfit$etahat
    RA      = fplsfit$RA
    W0      = fplsfit$W0
    basis.weight = fplsfit$basis.weight
    betacoef     = ginv(W0)%*%RA%*%etahat
    betahat      = fd(coef=betacoef,basisobj=basis.weight)
    beta_fpls_hat100[,iterk,itersim] = eval.fd(tobs,betahat)
    
    u100test = get.U(xtest,basis = basis.weight,rng=domain,Mdiv=Mobs)
    y100testhat = u100test%*%ginv(W0)%*%RA%*%etahat
    pmse100[iterk,itersim] = sum((y100testhat - ytest)^2)/ntest
    
    
    # n = 500
    BICind500  = which(BICpls500[iterk,,itersim]==min(BICpls500[iterk,,itersim],na.rm = TRUE),arr.ind = TRUE)
    Afpls500[iterk,itersim]   = iterk
    Mfpls500[iterk,itersim]   = Mweight[BICind500]
    
    nknots   = Mfpls500[iterk,itersim]+1
    knots    = seq(domain[1],domain[2], length.out = nknots)
    nbasis   = nknots + norder - 2 
    basis.weight = create.bspline.basis(knots,nbasis,norder)
    basis    = basis.weight
    basismat = eval.basis(tobs,basis)
    u500sim  = h/3*x500obs%*%diag(cef)%*%basismat
    
    fplsfit = convFPLS.mat(y=y500sim,u=u500sim,basis.weight=basis.weight,A=Afpls500[iterk,itersim])
    etahat  = fplsfit$etahat
    RA      = fplsfit$RA
    W0      = fplsfit$W0
    basis.weight = fplsfit$basis.weight
    betacoef     = ginv(W0)%*%RA%*%etahat
    betahat      = fd(coef=betacoef,basisobj=basis.weight)
    beta_fpls_hat500[,iterk,itersim] = eval.fd(tobs,betahat)
    
    u500test    = get.U(xtest,basis = basis.weight,rng=domain,Mdiv=Mobs)
    y500testhat = u500test%*%ginv(W0)%*%RA%*%etahat
    pmse500[iterk,itersim] = sum((y500testhat - ytest)^2)/ntest
  }
  print(itersim)
}


nullrg = 251:501
nnrg = 1:250

beta0 = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
isenn1.n100.FPLS = isenull1.n100.FPLS = pmse1.n100.FPLS = array(NA, c(Acomp,nsim))
isenn1.n500.FPLS = isenull1.n500.FPLS = pmse1.n500.FPLS = array(NA, c(Acomp,nsim))
for(iter in 1:Acomp)
{
  isenn1.n100.FPLS[iter,]    = apply(beta_fpls_hat100[nnrg,iter,],2,ise2,b0=beta0[nnrg],rng=c(0,0.5)) 
  isenull1.n100.FPLS[iter,]  = apply(beta_fpls_hat100[nullrg,iter,],2,ise,b0=beta0[nullrg],rng=c(0.5,1)) 

  isenn1.n500.FPLS[iter,]    = apply(beta_fpls_hat500[nnrg,iter,],2,ise2,b0=beta0[nnrg],rng=c(0,0.5)) 
  isenull1.n500.FPLS[iter,]  = apply(beta_fpls_hat500[nullrg,iter,],2,ise,b0=beta0[nullrg],rng=c(0.5,1)) 
}

pmse1.n100.FPLS  = pmse100
pmse1.n500.FPLS  = pmse500

boxplot(t(isenn1.n100.FPLS))
boxplot(t(isenull1.n100.FPLS))
boxplot(t(pmse1.n100.FPLS))
boxplot(t(isenn1.n500.FPLS))
boxplot(t(isenull1.n500.FPLS))
boxplot(t(pmse1.n500.FPLS))






#--------------- SFPLS method ---------------#
M        = 100
d        = 3
norder   = d+1
knots    = seq(0,1, length.out=M+1)
nknots   = length(knots)
nbasis   = length(knots) + norder - 2 # i+2
basis    = create.bspline.basis(knots,nbasis,norder)
beta.basis = basis
W        = slos.compute.weights(beta.basis)

kappa   = 0.1
lambda  = exp(seq(-12,-5,length.out = 15))
delta   = 10^(-1)
gamma   = 10^(-6:-4)
Maxiter = 100
absTol  = 10^(-15)
Cutoff  = 10^(-8)
k       = 5
BICs100 = BICs500 = array(NA,c(k,length(lambda),length(gamma),nsim))
lamselect100 = Kselect100 = gamselect100 = pmse100 = array(NA,c(k,nsim))
lamselect500 = Kselect500 = gamselect500 = pmse500 = array(NA,c(k,nsim))
beta.all.eval100 = beta.all.eval500 = array(NA,c(Mobs+1,k,nsim))

for(itersim in 1:nsim)
{
  # n = 100
  nsample = 100
  y100sim      = Y[1:nsample,itersim] 
  x100simcoef  = X[[itersim]]$coefs[,1:nsample]
  x100simbasis = X[[itersim]]$basis
  x100sim      = fd(coef=x100simcoef,basisobj=x100simbasis)
  x100obs      = t(eval.fd(tobs,x100sim))
  
  # n = 500
  nsample = 500
  y500sim      = Y[1:nsample,itersim] 
  x500simcoef  = X[[itersim]]$coefs[,1:nsample]
  x500simbasis = X[[itersim]]$basis
  x500sim      = fd(coef=x500simcoef,basisobj=x500simbasis)
  x500obs      = t(eval.fd(tobs,x500sim))
  
  ytest     = Ytest[,itersim] 
  xtest     = Xtest[[itersim]]
  
  for(iterk in 1:k)
  {
    # n = 100
    M2        = Mfpls100[iterk,itersim]
    knots2    = seq(domain[1],domain[2], length.out=M2+1)
    nknots2   = length(knots2)
    nbasis2   = length(knots2) + norder - 2 # i+2
    beta.basis.FPLS = create.bspline.basis(knots2,nbasis2,norder)
    
    for(iterlam in 1:length(lambda))
    {
      for(itergam in 1:length(gamma))
      {
        sfplsfit = try(fit.SpaFunPLS(y=y100sim,xind="fd",x=x100sim,beta.basis=basis,W=W,beta.basis.FPLS = beta.basis.FPLS,Mobs=Mobs,k=iterk,kappa=kappa,lambda=lambda[iterlam],delta=delta,gamma=gamma[itergam],Maxiter=Maxiter,absTol=absTol,Cutoff=Cutoff))
        if(isTRUE(class(sfplsfit)=="try-error")) {BICs100[iterk,iterlam,itergam,itersim] =  NA} else {BICs100[iterk,iterlam,itergam,itersim] = sfplsfit$BIC}
      }
    }
    
    
    # n = 500
    M2        = Mfpls500[iterk,itersim]
    knots2    = seq(domain[1],domain[2], length.out=M2+1)
    nknots2   = length(knots2)
    nbasis2   = length(knots2) + norder - 2 # i+2
    beta.basis.FPLS = create.bspline.basis(knots2,nbasis2,norder)
    
    for(iterlam in 1:length(lambda))
    {
      for(itergam in 1:length(gamma))
      {
        sfplsfit = try(fit.SpaFunPLS(y=y500sim,xind="fd",x=x500sim,beta.basis=basis,W=W,beta.basis.FPLS = beta.basis.FPLS,Mobs=Mobs,k=iterk,kappa=kappa,lambda=lambda[iterlam],delta=delta,gamma=gamma[itergam],Maxiter=Maxiter,absTol=absTol,Cutoff=Cutoff))
        if(isTRUE(class(sfplsfit)=="try-error")) {BICs500[iterk,iterlam,itergam,itersim] =  NA} else {BICs500[iterk,iterlam,itergam,itersim] = sfplsfit$BIC}
      }
    }
  }

  
  for(iterk in 1:k)
  {
    # n = 100
    BICind100 = which(BICs100[iterk,,,itersim]==min(BICs100[iterk,,,itersim],na.rm = TRUE),arr.ind = TRUE)
  
    if(dim(BICind100)[1]==0){pmse100[iterk,itersim]=NA}
    if(dim(BICind100)[1]!=0){
      
      if(dim(BICind100)[1]>1){
        BICind100=BICind100[1,]}
      
      Kselect100[iterk,itersim]   = iterk
      lamselect100[iterk,itersim] = lambda[BICind100[1]]
      gamselect100[iterk,itersim] = gamma[BICind100[2]]
      
      M2        = Mfpls100[iterk,itersim]
      knots2    = seq(domain[1],domain[2], length.out=M2+1)
      nknots2   = length(knots2)
      nbasis2   = length(knots2) + norder - 2 # i+2
      beta.basis.FPLS = create.bspline.basis(knots2,nbasis2,norder)
      
      sfplsfit = fit.SpaFunPLS(y=y100sim,xind="fd",x=x100sim,beta.basis=basis,W=W,beta.basis.FPLS = beta.basis.FPLS,Mobs=Mobs,k=Kselect100[iterk,itersim],kappa=kappa,lambda=lamselect100[iterk,itersim],delta=delta,gamma=gamselect100[iterk,itersim],Maxiter=Maxiter,absTol=absTol,Cutoff=Cutoff)
      beta.all.eval100[,iterk,itersim] = sfplsfit$beta[,Kselect100[iterk,itersim]]
      
      x100testobs = t(eval.fd(tobs,xtest))
      y100testhat = h/3*x100testobs%*%diag(cef)%*% beta.all.eval100[,iterk,itersim]
      pmse100[iterk,itersim] = sum((ytest-y100testhat)^2)/ntest}
    
    
    # n = 500
    BICind500 = which(BICs500[iterk,,,itersim]==min(BICs500[iterk,,,itersim],na.rm = TRUE),arr.ind = TRUE)
    
    if(dim(BICind500)[1]==0){pmse500[iterk,itersim]=NA}
    if(dim(BICind500)[1]!=0){
      
      if(dim(BICind500)[1]>1){
        BICind500=BICind500[1,]}
      
      Kselect500[iterk,itersim]   = iterk
      lamselect500[iterk,itersim] = lambda[BICind500[1]]
      gamselect500[iterk,itersim] = gamma[BICind500[2]]
      
      M2        = Mfpls500[iterk,itersim]
      knots2    = seq(domain[1],domain[2], length.out=M2+1)
      nknots2   = length(knots2)
      nbasis2   = length(knots2) + norder - 2 # i+2
      beta.basis.FPLS = create.bspline.basis(knots2,nbasis2,norder)
      
      sfplsfit = fit.SpaFunPLS(y=y500sim,xind="fd",x=x500sim,beta.basis=basis,W=W,beta.basis.FPLS = beta.basis.FPLS,Mobs=Mobs,k=Kselect500[iterk,itersim],kappa=kappa,lambda=lamselect500[iterk,itersim],delta=delta,gamma=gamselect500[iterk,itersim],Maxiter=Maxiter,absTol=absTol,Cutoff=Cutoff)
      beta.all.eval500[,iterk,itersim] = sfplsfit$beta[,Kselect500[iterk,itersim]]
      
      x500testobs = t(eval.fd(tobs,xtest))
      y500testhat = h/3*x500testobs%*%diag(cef)%*% beta.all.eval500[,iterk,itersim]
      pmse500[iterk,itersim] = sum((ytest-y500testhat)^2)/ntest}
  }
  print(itersim)
}


isenn1.n100.SFPLS = isenull1.n100.SFPLS = pmse1.n100.SFPLS = array(NA, c(k,nsim))
isenn1.n500.SFPLS = isenull1.n500.SFPLS = pmse1.n500.SFPLS = array(NA, c(k,nsim))
for(iter in 1:k)
{
  isenn1.n100.SFPLS[iter,]    = apply(beta.all.eval100[nnrg,iter,],2,ise2,b0=beta0[nnrg],rng=c(0,0.5)) 
  isenull1.n100.SFPLS[iter,]  = apply(beta.all.eval100[nullrg,iter,],2,ise,b0=beta0[nullrg],rng=c(0.5,1)) 
  
  isenn1.n500.SFPLS[iter,]    = apply(beta.all.eval500[nnrg,iter,],2,ise2,b0=beta0[nnrg],rng=c(0,0.5)) 
  isenull1.n500.SFPLS[iter,]  = apply(beta.all.eval500[nullrg,iter,],2,ise,b0=beta0[nullrg],rng=c(0.5,1)) 
}

pmse1.n100.SFPLS  = pmse100
pmse1.n500.SFPLS  = pmse500

boxplot(t(isenn1.n100.SFPLS))
boxplot(t(isenull1.n100.SFPLS))
boxplot(t(pmse1.n100.SFPLS))
boxplot(t(isenn1.n500.SFPLS))
boxplot(t(isenull1.n500.SFPLS))
boxplot(t(pmse1.n500.SFPLS))






#--------------- FPCR method ---------------#
Acomp  = 5
lambda = 10^(-10:6)
BICfpca100 = BICfpca500 = array(NA,c(length(lambda),Acomp,nsim))
beta_BIC = list()
lam100_BIC  = A100_BIC = array(NA, c(Acomp,nsim))
lam500_BIC  = A500_BIC = array(NA, c(Acomp,nsim))
pmse100_BIC = pmse500_BIC = array(NA,c(Acomp,nsim))
betahat100_BIC  = betahat500_BIC  = list()
betaeval100_BIC = betaeval500_BIC = array(NA,c(length(tobs),Acomp,nsim))

for(itersim in 1:nsim)
{
  # n = 100
  nsample = 100
  y100sim      = Y[1:nsample,itersim] 
  x100simcoef  = X[[itersim]]$coefs[,1:nsample]
  x100simbasis = X[[itersim]]$basis
  x100sim      = fd(coef=x100simcoef,basisobj=x100simbasis)
  x100obs      = t(eval.fd(tobs,x100sim))
  
  # n = 500
  nsample = 500
  y500sim      = Y[1:nsample,itersim] 
  x500simcoef  = X[[itersim]]$coefs[,1:nsample]
  x500simbasis = X[[itersim]]$basis
  x500sim      = fd(coef=x500simcoef,basisobj=x500simbasis)
  x500obs      = t(eval.fd(tobs,x500sim))
  
  ytest     = Ytest[,itersim] 
  xtest     = Xtest[[itersim]]
  
  for(iterlambda in 1:length(lambda))
  {
    fplsfit100 = FPCA(y=y100sim,x=x100sim,A=Acomp,lambda = lambda[iterlambda],rng=domain)
    fplsfit500 = FPCA(y=y500sim,x=x500sim,A=Acomp,lambda = lambda[iterlambda],rng=domain)
    
    for(iterA in 1:Acomp)
    {
      BICfpca100[iterlambda,iterA,itersim] = fplsfit100$BIC[iterA]
      BICfpca500[iterlambda,iterA,itersim] = fplsfit500$BIC[iterA]
    }
  }

  
  for(iterk in 1:Acomp)
  {
    # n = 100
    idx_BIC = which(BICfpca100[,iterk,itersim]==min(BICfpca100[,iterk,itersim],na.rm=TRUE),arr.ind = TRUE)
    A100_BIC[iterk,itersim]   = iterk
    lam100_BIC[iterk,itersim] = lambda[idx_BIC]
    
    fplsfit.BIC = FPCA(y=y100sim,x=x100sim,A=A100_BIC[iterk,itersim],lambda = lam100_BIC[iterk,itersim],rng=domain)
    betaeval100_BIC[,iterk,itersim] = eval.fd(tobs,fplsfit.BIC$betaA)
    phihat = fplsfit.BIC$fpcbasis
    phijhat=list()
    bcoefA = fplsfit.BIC$bcoef[[A100_BIC[iterk,itersim]]]
    u100test = array(NA,c(ntest,A100_BIC[iterk,itersim]))
    for(j in 1:A100_BIC[iterk,itersim])
    {
      phijhat[[j]]  = fd(coef=phihat[[j]]$coefs,basisobj = phihat[[j]]$basis)
      u100test[,j]  = inprod(xtest, phijhat[[j]], rng=domain)
    }
    y100testhat_BIC = u100test%*%bcoefA
    pmse100_BIC[iterk,itersim] = sum((y100testhat_BIC - ytest)^2)/ntest
    
    
    # n = 500
    idx_BIC = which(BICfpca500[,iterk,itersim]==min(BICfpca500[,iterk,itersim],na.rm=TRUE),arr.ind = TRUE)
    A500_BIC[iterk,itersim]   = iterk
    lam500_BIC[iterk,itersim] = lambda[idx_BIC]
    
    fplsfit.BIC = FPCA(y=y500sim,x=x500sim,A=A500_BIC[iterk,itersim],lambda = lam500_BIC[iterk,itersim],rng=domain)
    betaeval100_BIC[,iterk,itersim] = eval.fd(tobs,fplsfit.BIC$betaA)
    phihat = fplsfit.BIC$fpcbasis
    phijhat=list()
    bcoefA = fplsfit.BIC$bcoef[[A500_BIC[iterk,itersim]]]
    u500test = array(NA,c(ntest,A500_BIC[iterk,itersim]))
    for(j in 1:A500_BIC[iterk,itersim])
    {
      phijhat[[j]]  = fd(coef=phihat[[j]]$coefs,basisobj = phihat[[j]]$basis)
      u500test[,j]  = inprod(xtest, phijhat[[j]], rng=domain)
    }
    y500testhat_BIC = u500test%*%bcoefA
    pmse500_BIC[iterk,itersim] = sum((y500testhat_BIC - ytest)^2)/ntest
  }
  print(itersim)
}

isenn1.n100.FPCR = isenull1.n100.FPCR = pmse1.n100.FPCR = array(NA, c(k,nsim))
isenn1.n500.FPCR = isenull1.n500.FPCR = pmse1.n500.FPCR = array(NA, c(k,nsim))
for(iter in 1:k)
{
  isenn1.n100.FPCR[iter,]    = apply(betaeval100_BIC[nnrg,iter,],2,ise2,b0=beta0[nnrg],rng=c(0,0.5)) 
  isenull1.n100.FPCR[iter,]  = apply(betaeval100_BIC[nullrg,iter,],2,ise,b0=beta0[nullrg],rng=c(0.5,1)) 
  
  isenn1.n500.FPCR[iter,]    = apply(betaeval100_BIC[nnrg,iter,],2,ise2,b0=beta0[nnrg],rng=c(0,0.5)) 
  isenull1.n500.FPCR[iter,]  = apply(betaeval100_BIC[nullrg,iter,],2,ise,b0=beta0[nullrg],rng=c(0.5,1)) 
}

pmse1.n100.FPCR  = pmse100_BIC
pmse1.n500.FPCR  = pmse500_BIC

boxplot(t(isenn1.n100.FPCR))
boxplot(t(isenull1.n100.FPCR))
boxplot(t(pmse1.n100.FPCR))
boxplot(t(isenn1.n500.FPCR))
boxplot(t(isenull1.n500.FPCR))
boxplot(t(pmse1.n500.FPCR))

