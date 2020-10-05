# This file includes the following 
# *** basic functions
# *** functions for the SFPLS method
# *** functions for the FPLS method
# *** functions for the FPCR method
# *** functions for the FPLS_R method (Reiss and Ogden 2007)
# *** functions for the SLoS method (Lin et al. 2017)




#--------------- basis function ---------------#
inner.prod = function(f,basis,j)
{
  
  rng    = getbasisrange(basis)
  knots  = c(rng[1],basis$params,rng[2])
  nbasis = basis$nbasis
  norder = basis$nbasis - length(knots) + 2
  
  a = rng[1]
  if(j-norder > 0) a = knots[j-norder+1]
  
  b = rng[2]
  if (j <= nbasis-norder) b = knots[j+1]
  
  bfun = function(t)
  {
    mat = eval.basis(t,basis)
    z = t(mat[,j])
  }
  
  y = integrate(function(t) {f(t)*bfun(t)},a,b)
  y$value
  
}



# INput:
# x: fd object
# basis: basis
get.U = function(x,basis,rng,Mdiv)
{
  tobs     = seq(rng[1],rng[2],length.out = Mdiv+1)
  h        = (rng[2]-rng[1])/Mdiv
  cef      = c(1, rep(c(4,2), (Mdiv-2)/2), 4, 1)
  xobs     = eval.fd(tobs,x)
  basismat = eval.basis(tobs,basis)
  
  U = h/3*t(xobs)%*%diag(cef)%*%basismat
  return(U)
}



get.W0 = function(basis)
{
  V0 = eval.penalty(basis,int2Lfd(0))
  eig = eigen(V0) 
  eig$values[eig$values < 0] = 0 
  W0 = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)
  return(W0)
}


ise = function(b,b0,rng) # input b and b0 is a vector, domain is c(start, end) and is c(0, 1) here
{
  lgb = length(b) 
  l = rng[2]-rng[1]
  h = (rng[2]-rng[1])/(lgb-1)
  cef = c(1, rep(c(4,2), (lgb-3)/2), 4, 1)
  difb2=(b-b0)^2
  isetp = 1/l*h/3*cef%*%difb2 
  return(isetp)
}



ise2 = function(b,b0,rng)
{
  # rng  = rngnull1
  #  b  = betahat.FPLS[nullrg1,1]
  #  b0 = beta0[nullrg1]
  
  lgb   = length(b)
  l     = rng[2]-rng[1]
  h     = (rng[2]-rng[1])/(lgb)
  difb2 = (b-b0)^2*h
  isetp = 1/l*sum(difb2) 
  #return(list(ise=isetp,l=l))
  return(isetp)
}




#--------------- slope functions in the simulation studies ---------------#
beta_fun = function(t, ii)
{
  if(ii == 1){bf = ifelse(t<0.5 && t>=0,1,0)}
  else if(ii == 2){bf = sin(2*pi*t)*ifelse(t<0.5 && t>=0,1,0)}
  else if(ii == 3){bf = 3*t+exp(t^2)*cos(3*pi*t)+1}
  else{print("model does not exit")}
  return(bf)
}




#--------------- Data generation ---------------#
# n: a numeric value, sample size
# nknots and norder: numeric values, number of knots and order of bspline basis used to generate x
# snr: a numeric value, signal to noise ratio
# betaind: a numeric value of 1, 2 or 3, indicating the scenarios
data.generator = function(n, nknots, norder, domain = c(0, 1), snr, betaind) 
{
  knots  = seq(domain[1], domain[2], length.out = nknots)
  nbasis = nknots + norder - 2
  basis  = create.bspline.basis(knots, nbasis, norder)
  
  cMat1 = matrix(rnorm(n*nbasis),n,nbasis)
  x = fd(coef=t(cMat1),basisobj=basis)
  
  G1 = matrix(0,nbasis,1)
  beta.func = function(t){beta_fun(t,ii=betaind)}
  for(j in 1:nbasis) G1[j] = inner.prod(beta.func,basis,j)
  y0 = cMat1 %*% G1
  
  eps0 = sd(y0)
  
  y = y0 + rnorm(n, mean = 0, sd = eps0/sqrt(snr))
  return(list(X = x, Y = y))
}
# OUTPUT: X: an fd object, Y: a vector










#--------------- FPLS method ---------------#
# INPUTS 
# y: vector of n
# u: matrix
# basis.weight: bspline basis
# A: a numeric value, indicating the max number of components
convFPLS.mat = function(y,u,basis.weight,A)
{
  n = length(y)
  
  W0      = get.W0(basis.weight)
  ustar   = u%*%ginv(W0)
  
  # fit PLS regression
  fitpls  = plsr(y~ustar, ncomp = A,scale=TRUE)
  
  yhat = array(NA,c(n,A))
  betacoef=list()
  
  WA      = fitpls$loading.weights
  PA      = fitpls$loadings
  RA      = WA%*%ginv((t(PA)%*%WA))
  
  DA      = ustar%*%RA
  Hat     = DA%*%ginv(t(DA)%*%DA)%*%t(DA)
  yhat    = Hat%*%y
  etahat  = ginv(t(DA)%*%DA)%*%t(DA)%*%y
  
  return(list(etahat = etahat,betacoef = betacoef,basis.weight=basis.weight,W0=W0,RA=RA,DA=DA,yhat=yhat))  
}



#--- Fit FPLS regression model ---#
# INPUTS 
# y: vector of n
# x: fd object
# domain: vector of 2, domain of x
# M: vector, bspline basis functions have M+1 equally spaced knots on domain
# d: degree of the bspline basis functions
# K: a numeric value, indicating the max number of components
# Mobs: a numeric value, evaluate bspline basis functions over Mobs+1 points on domain
# mod.select: string of eighter "BIC" or "CV", model selection criteria
fit.FunPLS = function(y, xind, x, domain, M, d, K, Mobs,mod.select, ncv=NULL)
{
  n    = length(y)
  tobs = seq(domain[1], domain[2], length.out = Mobs+1)
  h    = (domain[2]-domain[1])/Mobs
  cef  = c(1, rep(c(4,2), (Mobs-2)/2), 4, 1)
  
  if(xind=="fd"){xobs = t(eval.fd(tobs,x))}
  if(xind=="matrix"){xobs = x}
  
  if(mod.select=="BIC"){BICpls = array(NA, c(K,length(M)))
  
  for(iterM in 1:length(M))
  {
    norder   = d + 1
    nknots   = M[iterM] + 1
    knots    = seq(domain[1], domain[2], length.out = nknots)
    nbasis   = nknots + norder - 2 
    basis    = create.bspline.basis(knots, nbasis, norder)
    basismat = eval.basis(tobs,basis)
    u        = h/3*xobs%*%diag(cef)%*%basismat
    
    Kind = which(1:K<=M[iterM]+d)
    
    fitpls  = linear.pls.fit(X=u,y=y,m=K,compute.jacobian = T)
    
    for(iterk in Kind)
    {
      RSS     = fitpls$RSS[iterk+1]
      Dof     = fitpls$DoF[iterk+1]
      sigma   = fitpls$sigmahat[iterk+1]
      BICpls[iterk,iterM] = RSS+log(n)*sigma^2*Dof
    }
  }
  
  BICnsim = min(BICpls,na.rm = TRUE)
  BICind  = which(BICpls==min(BICpls,na.rm = TRUE),arr.ind = TRUE)
  
  Kfpls   = BICind[1]
  Mfpls   = M[BICind[2]]
}
  
  if(mod.select=="CV"){

    neachcv = floor(n/ncv)
    msepcv  = rep(NA, ncv)
    msepKM  = array(NA,c(K,length(M)))
    
    for(iterk in 1:K)
    {
      Mind = which(M+d >= iterk)
      for(iterM in Mind)
      {
        norder   = d + 1
        nknots   = M[iterM] + 1
        knots    = seq(domain[1], domain[2], length.out = nknots)
        nbasis   = nknots + norder - 2 
        basis    = create.bspline.basis(knots, nbasis, norder)
        basismat = eval.basis(tobs,basis)
        u        = h/3*xobs%*%diag(cef)%*%basismat
        
        for(itercv in 1:ncv)
        {
          if(itercv < ncv){test   = ((itercv-1)*neachcv+1):(itercv*neachcv)}else{test   = ((itercv-1)*neachcv+1):n}
          
          ntest  = length(test)
          ntrain = n-ntest
          ytrain  =  y[-test]
          ytest   =  y[test]
          
          utrain = u[-test,]
          utest  = u[test,]
          
          fplsfit = convFPLS.mat(y=ytrain,u=utrain,basis.weight=basis,A=iterk)
          etahat = fplsfit$etahat
          basis.weight = fplsfit$basis.weight
          W0 = fplsfit$W0
          RA = fplsfit$RA
          
          utestsstar = utest%*%ginv(W0)%*%RA
          ytesthat = utestsstar%*%etahat
          msepcv[itercv] = sum((ytest-ytesthat)^2)/ntest
        }
        msepKM[iterk,iterM] = mean(msepcv)
    }
  }
    idx     = which(msepKM==min(msepKM,na.rm=TRUE),arr.ind = TRUE)
    Kfpls   = idx[1]
    Mfpls   = M[idx[2]]
  }    
 
    norder   = d+1
    nknots   = Mfpls+1
    knots    = seq(domain[1],domain[2], length.out = nknots)
    nbasis   = nknots + norder - 2 
    basis    = create.bspline.basis(knots,nbasis,norder)
    basismat = eval.basis(tobs,basis)
    u        = h/3*xobs%*%diag(cef)%*%basismat
    
    fplsfit = convFPLS.mat(y=y,u=u,basis.weight=basis,A=Kfpls)
    etahat  = fplsfit$etahat
    RA      = fplsfit$RA
    W0      = fplsfit$W0
    basis.weight  = fplsfit$basis.weight
    betacoef      = ginv(W0)%*%RA%*%etahat
    betahat       = fd(coef=betacoef,basisobj=basis.weight)
    beta_fpls_hat = eval.fd(tobs,betahat)   
  
  return(list(Kopt=Kfpls, Mopt=Mfpls, betahat=betahat, betahatobs = beta_fpls_hat, etahat = etahat, basis.weight=basis.weight, W0=W0, RA=RA))
}


predict.FunPLS = function(y, x, domain, Mobs,fplsfit)
{
  n = length(y)
  etahat  = fplsfit$etahat
  RA      = fplsfit$RA
  W0      = fplsfit$W0
  basis.weight = fplsfit$basis.weight

  u = get.U(x,basis = basis.weight,rng=domain,Mdiv=Mobs)
  yhat = u%*%ginv(W0)%*%RA%*%etahat
  pmse = sum((y - yhat)^2)/n
  
  return(list(Yhat = yhat, PMSE = pmse))
}









#--------------- SFPLS method ---------------#
# y: vector, n, centered response
# x: matrix, n by p, centered 
# Maxiter: a numeric number, the maximum number of iterations for convergence of beta
# lambda: a positive numeric number, tuning parameter for fSCAD penalty
# gamma: a positive numeric number, tuning parameter for the roughness penalty
# beta.basis: basis for beta(t), the coefficient function 
# absTol: a numeric number, of max(norm(bHat)) is smaller than absTol, we stop another iteration
# Cutoff: a numeric number, if bHat is smaller than Cutoff, set it to zero to avoide numerical unstable
slosLQA_fpls= function(U,y,V,V0,bHat,W,gamma,lambda,delta,Maxiter,M,L,L2NNer,absTol,a)
{
  
  betaNormj = c(0,M)
  bZeroMat  = rep(FALSE,L)
  betaNorm  = Inf
  n         = length(y)
  
  it = 1
  while(it <= Maxiter)
  {
    betaNormOld = betaNorm
    betaNorm    = sqrt(sum(bHat^2))
    
    change = (betaNormOld-betaNorm)^2
    if(change < absTol || betaNorm==0) break
    
    lqaW = NULL
    lqaWk = matrix(0,L,L)
    for(j in 1:M)
    {
      index = j:(j+d)
      betaNormj[j] = t(bHat[j:(j+d)])%*%W[,,j]%*%bHat[j:(j+d)]
      cjk = Dpfunc(betaNormj[j]*L2NNer,lambda,a)  
      
      if(cjk != 0)
      {
        if(betaNormj[j] < absTol) {bZeroMat[index] = TRUE}else{lqaWk[index,index] = lqaWk[index,index] + cjk*(L2NNer/betaNormj[j])*W[,,j]}
      }
    }
    
    lqaW = lqaWk
    lqaW = lqaW / 2
    
    bZeroVec    = bZeroMat
    bNonZeroVec = !bZeroVec
    
    if(sum(bNonZeroVec)!=0){
      UtU = t(U[,bNonZeroVec,drop=F])%*%U[,bNonZeroVec,drop=F]
      Ut  = t(U[,bNonZeroVec,drop=F])
      Vp  = gamma*V[bNonZeroVec,bNonZeroVec]
      Vq  = delta*V0[bNonZeroVec,bNonZeroVec]
      
      theta = solve(UtU+Vp+Vq+lqaW[bNonZeroVec,bNonZeroVec,drop=F],Ut %*% y)
      bHat = matrix(0,length(bNonZeroVec),1)
      bHat[bNonZeroVec] = theta}else{bHat = rep(0,L)}
    
    it = it + 1
  }
  tempa = as.numeric(t(bHat)%*%V0%*%bHat)
  if(tempa == 0){bHat = rep(0,M+d)}else{bHat = sqrt(1/as.numeric(t(bHat)%*%V0%*%bHat))*bHat}
  return(list(bHat = bHat)) 
}


slos.compute.weights = function(basis)
{
  L       = basis$nbasis
  rng     = getbasisrange(basis)
  breaks  = c(rng[1],basis$params,rng[2])
  M       = length(breaks) - 1
  norder  = L-M+1
  W       = array(0,dim=c(norder,norder,M))
  
  for (j in 1:M)
  {
    temp = inprod(basis,basis,rng=c(breaks[j],breaks[j+1]))
    W[,,j] = temp[j:(j+norder-1),j:(j+norder-1)]
  }
  W
}

#- SCAD function
Dpfunc = function(u,lambda,a)
{
  if(u<=lambda) Dpval = lambda
  else if(u<a*lambda) Dpval = -(u-a*lambda)/(a-1)
  else Dpval = 0
  Dpval
}


#- calculate U
GetDesignMatrix = function(xind, x,beta.basis,M0)
{
  rng      = getbasisrange(beta.basis)
  tobs     = seq(rng[1],rng[2],length.out = M0+1)
  basismat = eval.basis(tobs,beta.basis)
  if(xind =="fd"){xobs  = eval.fd(tobs,x)}
  if(xind =="matrix"){xobs  = t(x)}
  h        = (rng[2]-rng[1])/M0
  cef      = c(1, rep(c(4,2), (M0-2)/2), 4, 1)
  U        = h/3*t(xobs)%*%diag(cef)%*%basismat
  return(U)
}


# fit SFPLS regression
# INPUT: y: a vector of n x 1 in a matrix form, a centered univariate response
#        x: fd object
#        beta.basis: B-spline basis for expanding beta and direction functions
#        k: a number, number of PLS basis
#        lambda1: a number, representing tuning parameter of the fSCAD penalty
#        lambda2: a number, representing tuning parameter of the L2 penalty
#        gamma: a number, representing tuning parameter of the roughness penalty
#        beta.basis.FPLS is from conventioanl FPLS fit selected M
fit.SpaFunPLS = function(y,xind,x,beta.basis,W,beta.basis.FPLS,Mobs,k,kappa,lambda,delta,gamma,Maxiter,absTol,Cutoff)
{
  n = length(y)
  
  # beta b-spline basis
  rng     = getbasisrange(beta.basis)
  breaks  = c(rng[1],beta.basis$params,rng[2])
  tobs    = seq(rng[1],rng[2],length.out = Mobs+1)
  L       = beta.basis$nbasis
  M       = length(breaks) - 1
  norder  = L-M+1
  d       = L-M
  
  h       = (rng[2]-rng[1])/(Mobs-1)
  cef     = c(1, rep(c(4,2), (Mobs-2)/2), 4, 1)
  L2NNer  = sqrt(M/(rng[2]-rng[1]))  
  
  # calculate design matrix U and roughness penalty matrix V
  U   = GetDesignMatrix(xind=xind, x=x,beta.basis=beta.basis,M0=Mobs)
  V0  = eval.penalty(beta.basis,int2Lfd(0))
  V   = eval.penalty(beta.basis,int2Lfd(2))
  VV0 = delta*V0
  VV  = gamma*V
  
  eig = eigen(V0) 
  eig$values[eig$values < 0] = 0 
  W0  = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)
  
  y1   = y
  
  # step 1
  beta_pls0 = rep(0,Mobs+1)
  if(xind=="fd"){ xobs = t(eval.fd(tobs,x))}
  if(xind=="matrix"){xobs = x}
  betafd    = list()
  betaeval  = array(NA,c(Mobs+1,k))
  
  iterfunk = 1
  while(iterfunk <= k)
  {
    P      = t(U)%*%y1
    
    # initial estimate of b
    Ztilde = ginv(W0)%*%P
    Ztilde = Ztilde/sqrt(sum(Ztilde^2))
    bHat   = ginv(W0)%*%Ztilde
    # cHat  = (1-kappa)*solve((1-kappa)*P%*%t(P)+VV0+VV)%*%P%*%t(P)%*%bHat
    
    # (a) find w_hat
    Up        = sqrt(1-kappa)*t(P)
    yp        = sqrt(1-kappa)*t(P)%*%bHat
    slosfit   = slosLQA_fpls(U=as.matrix(Up,ncol=1),y=yp,V=V,V0=V0,bHat=bHat,W=W,gamma=gamma,lambda=lambda,delta=delta,Maxiter=Maxiter,M=M,L=L,L2NNer=L2NNer,absTol=absTol,a=3.7)
    chatemp   = slosfit$bHat
    chat      = chatemp/sqrt(sum(chatemp^2))
    
    # if(sum(abs(chat))==0) break
    
    # (b) update set A
    basismat  = eval.basis(tobs,beta.basis)
    weightfun = fd(coef=chat,basisobj = beta.basis)
    
    wt = basismat%*%chat
    Aset = union(which(wt!=0),which(beta_pls0!=0))
    
    # (c) fit PLS with x_A using kk number of latent components
    if(length(Aset) == 0){xA = array(0,c(n,Mobs+1))}else if(length(Aset) == Mobs+1){xA = xobs}else{
      xA = xobs
      xA[,-Aset] = 0} 
    
    
    
    basismat.fpls  = eval.basis(tobs,beta.basis.FPLS)
    UA.SFPLS       = h/3*xA%*%diag(cef)%*%basismat.fpls

    fplsfit        = convFPLS.mat(y=y,u=UA.SFPLS,basis.weight=beta.basis.FPLS,A=iterfunk)
    
    etahat.SFPLS    = fplsfit$etahat
    RA.SFPLS        = fplsfit$RA
    DA.SFPLS        = fplsfit$DA
    W0.SFPLS        = fplsfit$W0
    
    
    betacoef        = ginv(W0.SFPLS)%*%RA.SFPLS%*%etahat.SFPLS
    #betafd[[iterk]] = fd(coef=betacoef,basisobj=beta.basis.FPLS)
    
    if(length(Aset) == 0){basismat.fpls = array(0,c(Mobs+1,beta.basis.FPLS$nbasis))}else{
      basismat.fpls[-Aset,] = 0}
    betaeval[,iterfunk] = beta_pls0 = basismat.fpls%*%betacoef
    #   yhat = UA.SFPLS%*%betacoef
    yhat = h/3*xA%*%diag(cef)%*%betaeval[,iterfunk]
    
    # update y1
    y1 = y - yhat
    
    iterfunk = iterfunk + 1
  }
  
  fitpls  = linear.pls.fit(X=UA.SFPLS,y=y,m=k,compute.jacobian = T)
  RSS     = fitpls$RSS[k+1]
  Dof     = fitpls$DoF[k+1]
  sigma   = fitpls$sigmahat[k+1]
  BICpls  = RSS+log(n)*sigma^2*Dof
  
  return(list(beta = betaeval,etahat = etahat.SFPLS,yhat = yhat,DoF = Dof,BIC=BICpls))
}  


#- Tune the SFPLS fit using BIC
SpaFunPLS.tune.BIC = function(Y, xind,X, k, basis, W, M.FunPLS,kappa, d,lambda, delta, gamma, domain)
{
  knots    = seq(domain[1],domain[2], length.out= M.FunPLS+1)
  nknots   = length(knots)
  norder   = d+1
  nbasis   = M.FunPLS + d # i+2
  beta.basis.FPLS = create.bspline.basis(knots,nbasis,norder)
  
  Kind = which(M.FunPLS+d >= 1:k)
  BICs = array(NA,c(k,length(lambda),length(gamma)))
  for(iterk in Kind)
  {
    for(iterlam in 1:length(lambda))
    {
      for(itergam in 1:length(gamma))
      {
        sfplsfit = try(fit.SpaFunPLS(y=Y,xind=xind,x=X,beta.basis=basis,W=W,beta.basis.FPLS = beta.basis.FPLS,Mobs=Mobs,k=iterk,kappa=kappa,lambda=lambda[iterlam],delta=delta,gamma=gamma[itergam],Maxiter=Maxiter,absTol=absTol,Cutoff=Cutoff))
        if(isTRUE(class(sfplsfit)=="try-error")) { BICs[iterk,iterlam,itergam] =  NA} else {BICs[iterk,iterlam,itergam] = sfplsfit$BIC}
      }
    }
  }
  
  BICnsim = min(BICs,na.rm = TRUE)
  BICind = which(BICs==min(BICs,na.rm = TRUE),arr.ind = TRUE)
  
  if(dim(BICind)[1]>1){
    BICind=BICind[1,]}
  
  Kopt   = BICind[1]
  lamopt = lambda[BICind[2]]
  if(length(gamma)==1){gamopt = gamma[1]}else{gamopt = gamma[BICind[3]]}
  
  
  sfplsfit = fit.SpaFunPLS(y=Y,xind=xind,x=X,beta.basis=basis,W=W,beta.basis.FPLS = beta.basis.FPLS,Mobs=Mobs,k=Kopt,kappa=kappa,lambda=lamopt,delta=delta,gamma=gamopt,Maxiter=Maxiter,absTol=absTol,Cutoff=Cutoff)
  beta.spa = sfplsfit$beta[,Kopt]
  
  return(list(Kopt = Kopt, lamopt = lamopt, gamopt = gamopt, beta = beta.spa,BIC=BICs))
}


# SFPLS method predict 
predict.SpaFunPLS = function(y, x, domain, Mobs,sfplsfit.beta)
{
  n = length(y)
 
  tobs   = seq(domain[1],domain[2],length.out = Mobs+1)
  hobs   = (domain[2]-domain[1])/Mobs
  cefobs = c(1, rep(c(4,2), (Mobs-2)/2), 4, 1)
  xobs   = t(eval.fd(tobs,x))
  
  yhat = hobs/3*xobs%*%diag(cefobs)%*%sfplsfit.beta
  pmse = sum((y-yhat)^2)/n
 
  return(list(Yhat = yhat, PMSE = pmse))
}










#--------------- FPCR method ---------------#
# y: n x 1 response vector
# x: an fd object
# M: M+1 is the number of knots of Bspline basis functions to expand a weight function
# d: degree of Bspline basis functions to expand a weight function
# Mobs: number of observed t to evaluate the weight function
# A: number of components
# rng: support of x(t)
FPCA = function(y,x,A,lambda,rng)
{
  # y = ytrain[,1,3]
  # x = xtrain[[1]]
  # A = 2
  # lambda = 0.01
  # rng = c(0,1)
  # 
  
  n = length(y)
  # estimate g
  xcoef    = x$coef
  xbasis   = x$basis
  
  # obtain thetahat and estimate gjhat
  tempfdPar = fdPar(xbasis,2,lambda)
  fpcfit    = pca.fd(x,nharm=A,2, harmfdPar=tempfdPar)
  thetahat  = fpcfit$values
  phihat    = fpcfit$harmonics
  phijhat   = list()
  
  Umtrx = array(NA,c(n,A))
  bcoef = Hat = list()
  AICpca = BICpca = df = rep(NA,A)
  for(j in 1:A)
  {
    phijhat[[j]]  = fd(coef=phihat$coefs[,j],basisobj = phihat$basis)
    Umtrx[,j]     = inprod(x, phijhat[[j]], rng=rng)
    Uuptoj        = Umtrx[,1:j]
    bcoef[[j]]    = ginv(t(Uuptoj)%*%Uuptoj)%*%t(Uuptoj)%*%y
    
    # beta[[j]] = sum bcoef[[j]]*phijhat[[j]]
    
    Hat[[j]]      = Uuptoj%*%ginv(t(Uuptoj)%*%Uuptoj)%*%t(Uuptoj)
    df[j]         = tr(Hat[[j]])
    yhat          = Hat[[j]]%*%y
    
    resids        = t(y-yhat)%*%(y-yhat)
    AICpca[j]     = n*log(resids/n)+2*df[j]
    BICpca[j]     = n*log(resids/n)+log(n)*df[j]
  }
  
  
  betaAcoef = t(t(bcoef[[A]])%*%t(phihat$coefs))
  betaAbasis = phihat$basis 
  return(list(fpcbasis = phijhat,betaA=fd(coef=betaAcoef,basisobj = betaAbasis), AIC = AICpca,BIC=BICpca,bcoef = bcoef))  
}
# OUTPUT
# fpcbasis: list with j fd objects
# AIC: a vector with length A
# BIC: a vector with length A
# bcoef: list with j objects, the jth object is a vector of length j
fit.FPCR = function(y,x,K,Mobs,lambda,rng,ncv)
{
  n        = length(y)
  tobs     = seq(rng[1],rng[2],length.out = Mobs+1)
  neachcv  = floor(n/ncv)
  msepcv   = array(NA,c(K,length(lambda),ncv))
  msepKlam = array(NA,c(K,length(lambda)))

  for(itercv in 1:ncv)
  {
    if(itercv < ncv){test   = ((itercv-1)*neachcv+1):(itercv*neachcv)}else{test   = ((itercv-1)*neachcv+1):n}
    
    ntest    = neachcv
    ntrain   = n-ntest
    ytrain   =  y[-test]
    ytest    =  y[test]
    
    xbasis     = x$basis
    xtraincoef = x$coef[,-test]
    xtestcoef  = x$coef[,test]
    xtrain  = fd(coef=xtraincoef,basisobj=xbasis)
    xtest   = fd(coef=xtestcoef,basisobj=xbasis)
    
    for(iterlambda in 1:length(lambda))
    { 
      fplsfit = FPCA(y=ytrain,x=xtrain,A=K,lambda = lambda[iterlambda],rng=rng)
      phihat  = fplsfit$fpcbasis
      utest0 = array(NA,c(ntest,K))
      for(j in 1:K)
      {
        utest0[,j]  = inprod(xtest, phihat[[j]], rng=rng)
      }
  
      for(iterk in 1:K)
      {
        bcoefA      = fplsfit$bcoef[[iterk]]
        utest       = utest0[,1:iterk]
        ytesthat    = utest%*%bcoefA
        msepcv[iterk,iterlambda,itercv] = sum((ytest-ytesthat)^2)/ntest
      }
     }
    }
  
    for(iterk in 1:K)
    {
      for(iterlambda in 1:length(lambda))
      {
        msepKlam[iterk,iterlambda] = mean(msepcv[iterk,iterlambda,])
      }
    }
  
  idx    = which(msepKlam==min(msepKlam,na.rm=TRUE),arr.ind = TRUE)
  Kopt   = idx[1]
  lamopt = lambda[idx[2]]
  
  fplsfit  = FPCA(y=y,x=x,A=Kopt,lambda = lamopt,rng=rng)
  betahat  = fplsfit$betaA
  betaeval = eval.fd(tobs,betahat)
  phihat   = fplsfit$fpcbasis
  bcoefA   = fplsfit$bcoef
  
  return(list(beta = betahat, Kopt=Kopt, lamopt = lamopt,betaobs = betaeval, phi = phihat,bcoef = bcoefA))
}



fit.FPCR.BIC = function(y,x,K,Mobs,lambda,rng)
{
  n        = length(y)
  tobs     = seq(rng[1],rng[2],length.out = Mobs+1)
  BICfpca  = array(NA,c(length(lambda),K))
  
  for(iterlambda in 1:length(lambda))
  { 
    fplsfit = FPCA(y=y,x=x,A=K,lambda = lambda[iterlambda],rng=rng)
    for(iterA in 1:K)
    {
      BICfpca[iterlambda,iterA] = fplsfit$BIC[iterA]
    }
  }
  
  idx_BIC  = which(BICfpca==min(BICfpca,na.rm=TRUE),arr.ind = TRUE)
  Kopt     = idx_BIC[2]
  lamopt   = lambda[idx_BIC[1]]
  
  fplsfit  = FPCA(y=y,x=x,A=Kopt,lambda = lamopt,rng=rng)
  betahat  = fplsfit$betaA
  betaeval = eval.fd(tobs,betahat)
  phihat   = fplsfit$fpcbasis
  bcoefA   = fplsfit$bcoef
  
  return(list(beta = betahat, Kopt=Kopt, lamopt = lamopt,betaobs = betaeval, phi = phihat,bcoef = bcoefA))
}




predict.FPCR = function(y, x, rng, Mobs,fpcr_Rfit)
{
  n         = length(y)
  phihat    = fpcr_Rfit$phi
  bcoefA    = fpcr_Rfit$bcoef
  K         = fpcr_Rfit$Kopt
  
  u = array(NA,c(n,K))
  for(j in 1:K)
  {
    u[,j] = inprod(x, phihat[[j]], rng=rng)
  }
  
  yhat = u%*%bcoefA[[K]]
  pmse = sum((yhat - y)^2)/n
  
  return(list(Yhat = yhat, PMSE = pmse))
}










#--------------- FPLS_R method by Reiss and Ogden 2007 ---------------#
# y: n x 1 response vector
# x: fd object
# M: M+1 is the number of knots of B-splines
# d: the degree of B-spline 
# Mobs: number of points to evaluate U, the integral between x and basis functions
# A: number of components
# rng: domain of x
ReissFPLS.REML = function(y,x,M,d,Mobs,A,rng)
{
  n = length(y)
  norder   = d+1
  nknots   = M+1
  knots    = seq(rng[1],rng[2], length.out = nknots)
  nbasis   = nknots + norder - 2 
  basis    = create.bspline.basis(knots,nbasis,norder)
  u        = get.U(x,basis,rng,Mdiv=Mobs)
  
  
  fitpls  = plsr(y~u, ncomp = A)
  WA      = fitpls$loading.weights
  PA      = fitpls$loadings
  RA      = WA%*%ginv((t(PA)%*%WA))
  
  V = eval.penalty(basis,int2Lfd(2))
  ustar = u%*%RA
  vstar = t(RA)%*%V%*%RA
  
  
  # estimate lambda using REML
  svdvstar = svd(vstar)
  dstar = svdvstar$d
  svdUstar = svdvstar$u
  dstar[dstar < 0] = 0 
  s = length(dstar[dstar>0]) 
  
  if(s==A){D1 = diag(sqrt(dstar),nrow=A, ncol=A)}else{D1 = diag(c(sqrt(dstar[1:s]),rep(1,A-s)),nrow=A,ncol=A)}
  ginvD1 = ginv(D1)
  XZ = (ustar%*%svdUstar%*%ginvD1)[,1:s]
  
  if(s==A){
    id          = 1:n
    #data       = data.frame(Y=y,XZ=XZ,id=id)
    dummyID     = factor(rep(1,n))
    fit         = lme(y ~ -1,random = list(dummyID = pdIdent(~-1+XZ)))
    sigsqepsHat = fit$sigma^2
    sigsquHat   = as.numeric(VarCorr(fit)[1,1])
    lambda      = sigsqepsHat/sigsquHat}else{
      
      XF          = (ustar%*%svdustar%*%ginvD1)[,(s+1):A]
      id          = 1:n
      # data        = data.frame(Y=y,XF = XF, XZ=XZ,id=id)
      dummyID     = factor(rep(1,n))
      fit         = lme(y ~ -1+XF,random = list(dummyID = pdIdent(~-1+XZ)))
      sigsqepsHat = fit$sigma^2
      sigsquHat   = as.numeric(VarCorr(fit)[1,1])
      lambda      = sigsqepsHat/sigsquHat}
  
  
  etahat       = ginv(t(ustar)%*%ustar+n*lambda*vstar)%*%t(ustar)%*%y
  betacoef     = RA%*%etahat
  
  return(list(beta=fd(coef = betacoef,basisobj = basis),basis=basis,lambda=lambda,RA=RA,eta = etahat))
}

# M: a numeric value
# K: a numeric value, number of components
fit.FPLS_R = function(y,x,M,d,Mobs,K,rng,ncv)
{
  n       = length(y)
  tobs    = seq(rng[1],rng[2],length.out = Mobs+1)
  neachcv = floor(n/ncv)
  msepcv  = rep(NA, ncv)
  msepK   = rep(NA, K)
  
  for(iterk in 1:K)
  {
    for(itercv in 1:ncv)
    {
      if(itercv < ncv){test   = ((itercv-1)*neachcv+1):(itercv*neachcv)}else{test   = ((itercv-1)*neachcv+1):n}
      
      ntest    = neachcv
      ntrain   = n-ntest
      ytrain   =  y[-test]
      ytest    =  y[test]
      
      xbasis     = x$basis
      xtraincoef = x$coef[,-test]
      xtestcoef  = x$coef[,test]
      xtrain     = fd(coef=xtraincoef,basisobj=xbasis)
      xtest      = fd(coef=xtestcoef,basisobj=xbasis)
      
      fit_REML = ReissFPLS.REML(y=ytrain,x=xtrain,M=M,d=d,Mobs=Mobs,A=iterk,rng=rng)
      
      etahat_R = fit_REML$eta
      RA_R     = fit_REML$RA
      basis_R  = fit_REML$basis

      # calculate test error
      utest_R = get.U(xtest,basis = basis_R,rng=rng,Mdiv=Mobs)
      uteststar_R = utest_R%*%RA_R
      
      ytesthat_R =  uteststar_R%*%etahat_R
      msepcv[itercv] = sum((ytest-ytesthat_R)^2)/ntest
    }
    msepK[iterk] = mean(msepcv)
  }
  
  Kopt  = which.min(msepK)
  
  fit_REML   = ReissFPLS.REML(y=y,x=x,M=M,d=d,Mobs=Mobs,A=Kopt,rng=rng)
  betahat_R  = fit_REML$beta
  beta_R_hat = eval.fd(tobs,betahat_R)
  
  lambdaR   = fit_REML$lambda
  etahat_R  = fit_REML$eta
  RA_R      = fit_REML$RA
  basis_R   = fit_REML$basis

  return(list(betahat = betahat_R, betaobs = beta_R_hat,lamopt = lambdaR,Kopt=Kopt,eta = etahat_R, R=RA_R,basis=basis_R))
}


predict.FPLS_R = function(y, x, domain, Mobs,fpls_Rfit)
{
  n         = length(y)
  utest     = get.U(x,basis = fpls_Rfit$basis,rng=domain,Mdiv=Mobs)
  uteststar = utest%*%fpls_Rfit$R
  yhat      =  uteststar%*%fpls_Rfit$eta
  pmse      = sum((y-yhat)^2)/n
  
  return(list(Yhat = yhat, PMSE = pmse))
}





#--------------- SLoS method ---------------#
# y: vector, n, centered response
# x: matrix, n by p, centered 
# Maxiter: a numeric number, the maximum number of iterations for convergence of beta
# lambda: a positive numeric number, tuning parameter for fSCAD penalty
# gamma: a positive numeric number, tuning parameter for the roughness penalty
# beta.basis: basis for beta(t), the coefficient function 
# absTol: a numeric number, of max(norm(bHat)) is smaller than absTol, we stop another iteration
# Cutoff: a numeric number, if bHat is smaller than Cutoff, set it to zero to avoide numerical unstable
slos_temp = function(Y, X, Maxiter,lambda,gamma,beta.basis,absTol,Cutoff)
{
  n = length(Y)
  m = dim(X)[2]
  
  # beta b-spline basis
  rng = getbasisrange(beta.basis)
  breaks = c(rng[1],beta.basis$params,rng[2])
  L = beta.basis$nbasis
  M = length(breaks) - 1
  norder = L-M+1
  d = L-M
  
  L2NNer = sqrt(M/(rng[2]-rng[1]))  
  
  # calculate design matrix U and roughness penalty matrix V
  U = GetDesignMatrix2(X=X,beta.basis=beta.basis)
  V = eval.penalty(beta.basis,int2Lfd(2))
  VV = n*gamma*V
  
  # calculate W
  W = slos.compute.weights(beta.basis)
  
  # initial estimate of b
  bHat = solve(t(U)%*%U+VV)%*%t(U)%*%Y
  
  bTilde = bHat
  
  if(lambda > 0)
  {
    changeThres = absTol
    bTilde = slosLQA(U,Y,V,bHat,W,gamma,lambda,Maxiter,M,L,L2NNer,absTol,a=3.7)
    bZero = (abs(bTilde) < Cutoff)
    bTilde[bZero] = 0
    
    bNonZero = !bZero
    
    U1 = U[,bNonZero]
    V1 = VV[bNonZero,bNonZero]
    bb = solve(t(U1)%*%U1+V1,t(U1)%*%Y)
    bTilde = matrix(0,dim(U)[2],1)
    bTilde[bNonZero,1] = matrix(bb,length(bb),1)
  }
  
  bNonZero = as.vector((bTilde != 0))
  
  projMat = U1 %*% solve(t(U1)%*%U1+V1,t(U1))
  result = list(beta=NULL,projMat=projMat,intercept=0,fitted.values=projMat%*%Y)
  
  betaobj = list()
  bfd = fd(coef=bTilde,basisobj=beta.basis)
  betaobj = bfd
  
  result$b = bTilde
  result$beta = betaobj
  result$U = U
  class(result) = 'slos'
  result
}


GetDesignMatrix2 = function(X,beta.basis)
{
  rng     = getbasisrange(beta.basis)
  breaks  = c(rng[1],beta.basis$params,rng[2])
  Mobs    = dim(X)[2]-1
  tobs    = seq(rng[1],rng[2],length.out = Mobs+1)
  
  beta.basismat = eval.basis(tobs,beta.basis)
  hDM           = (rng[2]-rng[1])/Mobs
  cefDM         = c(1, rep(c(4,2), (Mobs-2)/2), 4, 1)
  U             = hDM/3*X%*%diag(cefDM)%*%beta.basismat
  return(U=U)
}




slosLQA = function(U,Y,V,bHat,W,gamma,lambda,Maxiter,M,L,L2NNer,absTol,a)
{
  betaNormj = c(0,M)
  bZeroMat = rep(FALSE,L)
  betaNorm = Inf
  n = length(Y)
  
  it = 1
  while(it <= Maxiter)
  {
    betaNormOld = betaNorm
    betaNorm = sqrt(sum(bHat^2))
    
    change = (betaNormOld-betaNorm)^2
    if(change < absTol) break
    
    lqaW = NULL
    lqaWk = matrix(0,L,L)
    for(j in 1:M)
    {
      index = j:(j+d)
      betaNormj[j] = t(bHat[j:(j+d)])%*%W[,,j]%*%bHat[j:(j+d)]
      cjk = Dpfunc(betaNormj[j]*L2NNer,lambda,a)  
      
      if(cjk != 0)
      {
        if(betaNormj[j] < absTol) {bZeroMat[index] = TRUE}else{lqaWk[index,index] = lqaWk[index,index] + cjk*(L2NNer/betaNormj[j])*W[,,j]}
      }
    }
    
    lqaW = lqaWk
    lqaW = lqaW / 2
    
    bZeroVec = bZeroMat
    bNonZeroVec = !bZeroVec
    
    UtU = t(U[,bNonZeroVec])%*%U[,bNonZeroVec]
    Ut = t(U[,bNonZeroVec])
    Vp = n*gamma*V[bNonZeroVec,bNonZeroVec]
    
    theta = solve(UtU+Vp+n*lqaW[bNonZeroVec,bNonZeroVec,drop=F],Ut %*% Y)
    bHat = matrix(0,length(bNonZeroVec),1)
    bHat[bNonZeroVec] = theta
    
    it = it + 1
  }
  bHat 
}




BIC.fun = function(Y, b, U, V, n, alpha)
{
  sparse.idx   = which(b == 0)
  if(length(sparse.idx) == 0)
  {
    ula = U
    vla = V
  }
  else{
    ula  = U[, -sparse.idx]
    vla  = V[-sparse.idx, -sparse.idx]
  }
  
  hat2 = ula%*%solve(t(ula)%*%ula + n*alpha*vla)%*%t(ula)
  df2  = tr(hat2)
  Yhat = U%*%b
  RSS.gbric  = t(Y - Yhat)%*%(Y - Yhat)
  BIC2temp.gbr = n*log(RSS.gbric/n) + log(n)*df2
  return(list(BIC2 = BIC2temp.gbr))
}




slos.fit = function(Y, X, V, Maxiter,lambda,gamma,M,d,domain,absTol,Cutoff,tobs)
{
  
  # Y=Y[1:n,itersim]
  # X = t(X[1:n,,itersim])
  # 
  Yn = Y
  Xn = t(X)
  Ymean=mean(Yn)
  Xmean = apply(Xn,2,mean) # column mean of the covariate curves
  n = length(Y)
  m = dim(X)[2]
  
  norder   = d+1
  knots    = seq(domain[1],domain[2], length.out=M+1)
  nknots   = length(knots)
  nbasis   = length(knots) + norder - 2 # i+2
  beta.basis  = create.bspline.basis(knots,nbasis,norder)
  
  beta.slos.all = array(NA,c(length(tobs),length(lambda),length(gamma)))
  b.slos.all = array(NA,c(M+d,length(lambda),length(gamma)))
  BIC_slos  = array(NA,c(length(lambda),length(gamma)))
  
  for(ii in 1:length(lambda))
  {
    for(jj in 1:length(gamma))
    {
      cYn = Yn - Ymean
      cXn = Xn - matrix(rep(Xmean,dim(Xn)[1]),byrow=TRUE,nrow=dim(Xn)[1])
      
      slosresult = slos_temp(Y=cYn, X=cXn, Maxiter=Maxiter,lambda=lambda[ii],gamma=gamma[jj],beta.basis=beta.basis,absTol=absTol,Cutoff=Cutoff)
      beta.temp =  slosresult$beta
      beta.slos.all[,ii,jj] =  eval.fd(tobs,beta.temp)
      b.slos.all[,ii,jj] = slosresult$b
      
      BIC_temp = BIC.fun(Y=Yn,b=slosresult$b,U=slosresult$U,V=V, n=length(Yn), alpha=gamma[jj])
      BIC_slos[ii,jj]  = BIC_temp$BIC2
    }
  }
  
  idx = which(BIC_slos == min(BIC_slos), arr.ind = TRUE)
  idx = idx[1,]
  
  Optlambda = lambda[idx[1]]
  Optgamma  = gamma[idx[2]]
  
  beta_slos  = beta.slos.all[,idx[1],idx[2]]
  b_slos     = b.slos.all[,idx[1],idx[2]]
  
  result = list(beta=beta_slos,b = b_slos,Optgamma=Optgamma,Optlambda=Optlambda)
  class(result) = 'slos'
  result 
}

