"
    Directly download from https://www.ida.liu.se/~olesy12/include/soft.en.shtml
"

library(fANCOVA)
library(limSolve)

K1<-function(X1,X2,p) {
  return(1/(abs(X1-X2)^p))  
}

K2<-function(X,i,j,p) {
  return(K1(X[i],X[j],p))
}


ARsmooth=function(X,Y,W,lambda,Block, sB, Ys,Ws,correction=T) {
  n=length(Y)
  if (sB==1) return(Yf=sum(Y*W)/sum(W))
  A.da=rep(0, length=sB-1)
  A.db=rep(0, length=sB-1)
  A.d=rep(0, length=sB)
  Yb=rep(0, length=sB)
  P=rep(0, length=sB)
  Q=rep(0, length=sB)
  NB=rep(0, length=sB)
  Yb=rep(0,length=sB)
  r=1:sB
  
  ctr=0; 
  p=1
  
  while(p!=n+1) {
    ctr=ctr+1;
    p=P[ctr]=Block[p,1]
    q=Q[ctr]=Block[p,2]
    NB[ctr]=Ws[p]
    Yb[ctr]=Ys[p]/Ws[p]
    p=q+1
  }
  ind1=1:(sB-1)
  ind2=2:sB
  ind3=2:(sB-1)
  A.da=-lambda*Kern2(X, Q[ind1],Q[ind1]+1)/NB[ind1]
  A.db=-lambda*Kern2(X,P[ind2], P[ind2]-1)/NB[ind2]
  A.d[1]=1+lambda*Kern2(X,Q[1],Q[1]+1)/NB[1]
  ss=lambda*Kern2(X,P[sB],P[sB]-1)/NB[sB]
  if (length(ss)==0) browser()
  A.d[sB]=1+lambda*Kern2(X,P[sB],P[sB]-1)/NB[sB]
  if (sB>2) A.d[ind3]=1+lambda*(Kern2(X,Q[ind3],Q[ind3]+1)+Kern2(X,P[ind3],P[ind3]-1))/NB[ind3]

  Yf=Solve.tridiag(A.db,A.d,A.da,Yb)
  if(correction) {
    e=rep(0, length=sB)
    e[1]=1/NB[1]
    e[sB]=-1/NB[sB]
    ehat=Solve.tridiag(A.db,A.d,A.da,e)
    beta=mean((Yb-Yf)*ehat)/mean(ehat*ehat)
    Yf=Yf+beta*ehat
  }
  
  return(Yf)
}


ARInt<-function(X,Y,W,lambda, correction){
  r=length(X)
  Yf=ARsmooth(X,Y, W,lambda, cbind(1:r,1:r),r, Y, W, correction)

  return(Yf)
}

restoreY=function(Block,res,n){
  Yr=c(0,length=n)
  ctr1=1
  ctr=0
  while(ctr1<=n){
    ctr=ctr+1
    ctr1=Block[ctr1,1]
    ctr2=Block[ctr1,2]
    Yr[ctr1:ctr2]=res[ctr]
    ctr1=ctr2+1
  }
  return(Yr)
}

SPAVInt=function(X,Y,W, lambda, correction) {
  n=length(X)
  Knots=1:n
  Block=cbind(1:n,1:n)
  sB=n
  Yb=rep(0,n)
  Yr=rep(0,n)
  Ys=Y*W
  Ws=W
  foundMonot=T
  nIter=0
  
  res=ARsmooth(X,Y,W,lambda, Block, sB, Ys, Ws, correction)


  while(foundMonot) {
    p=1    
    ctr1=1
    ctr=0
    while(ctr1<=n){
      ctr=ctr+1
      ctr1=Block[ctr1,1]
      ctr2=Block[ctr1,2]
      Yb[ctr1]=res[ctr]
      Yb[ctr2]=res[ctr]
      ctr1=ctr2+1
    }
    p=Block[p,1]
    q=Block[p,2]
    foundMonot=F
    p1=q+1
    while((p1<=n)) {
      #search for first nonmonot
      if (Yb[p1-1]<Yb[p1]){
        p=q+1;
        q=Block[p,2]
        p1=q+1
      }  else{
        sB=sB-1
        Block[p,2]=Block[p1,2]
        Ys[p]=Ys[p]+Ys[p1]
        Ws[p]=Ws[p]+Ws[p1]
        Block[p1,1]=p
        foundMonot=T
        p1=Block[p1,2]+1
        q=Block[p,2]
      }
    }
        
    res=ARsmooth(X,Y,W,lambda, Block, sB, Ys, Ws,correction) 
    nIter=nIter+1
  }  
  
  ctr1=1
  while(ctr1<=n){
    ctr1=Block[ctr1,1]
    ctr2=Block[ctr1,2]
    Yr[ctr1:ctr2]=Yb[ctr1]
    ctr1=ctr2+1
  }
  
  return(Yr)
}


predict.SMR=function(smrobj, newdata=NULL) {
  Yhat=smrobj$fitted
  X=smrobj$X  
  n=length(X)
  
  if (is.null(newdata)) Xnew=X else Xnew=newdata
  n1=length(Xnew)
  Ynew=vector(length=n1)
  
  for(i in 1:n1) {
    Z=Xnew[i]<X
    indMi=which.max(Z)
    if (Z[indMi]==F) indMi=n else indMi=indMi-1
    Z=!Z
    indMa=which.min(Z)
    if (Z[indMa]==T) indMa=0
    
    if(indMi==0) {
      Ynew[i]=Yhat[1]
    } else if (indMa==0) {
      Ynew[i]=Yhat[n]
    } else if (X[indMi]==Xnew[i]){
      Ynew[i]=Yhat[i]
    }else {
      Ynew[i]=(Kern1(X[indMi],Xnew[i])*Yhat[indMi]+Kern1(Xnew[i],X[indMa])*Yhat[indMa])/(Kern1(X[indMi],Xnew[i])+Kern1(Xnew[i],X[indMa]))
    }
  }
  return(Ynew)
  
}
interpolate=function(X,Yhat,fold) {

  n=length(X)
  sz=blocks[1,fold]
  bls=blocks[2:sz, fold]
  Xnew=X[bls]
  Yh=numeric(n)
  Yh[-bls]=Yhat
  XMi=X[bL[bls]]
  XMa=X[bU[bls]]
  
  YMi=Yh[bL[bls]]
  YMa=Yh[bU[bls]]
  
  Ynew=(Kern1(XMi,Xnew)*YMi+Kern1(Xnew,XMa)*YMa)/(Kern1(XMi,Xnew)+Kern1(Xnew,XMa))

  return(Ynew)
  
}

prepareFolds=function(n, nfolds, correction){
  S=floor(n/nfolds)
  block=(1:n) %% nfolds
  perm=sample(block,n)
  blocks=matrix(0,ncol=nfolds, nrow=2*S)
  blocks[1,]=1
  bU=numeric(n)
  bL=numeric(n)
  bl=0
  blS=0
  for(i in 1:n){
    if (bl==0){
      bl=perm[i]+1
    } else if (bl==perm[i]+1){
    } else{
      if(blS==0){
        bL[1:(i-1)]=i
        bU[1:(i-1)]=i
        blS=i
        blE=i
        bl=perm[i]+1
      } else{
        bL[blS:(i-1)]=blS-1
        bU[blS:(i-1)]=i
        blS=i
        bl=perm[i]+1
      }
    }
    if(i==n){
      bL[blS:i]=blS-1
      bU[blS:i]=blS-1
    }

  }
  
  for (i in 1:n){
    bl=perm[i]+1
    blocks[blocks[1,bl]+1,bl]=i
    blocks[1,bl]=blocks[1,bl]+1
    
  }
  
  
  assign("blocks", blocks, envir=.GlobalEnv)
  assign("bU", bU, envir=.GlobalEnv)
  assign("bL", bL, envir=.GlobalEnv)

}

cvSSE=function(funct,lambda, X,Y,W, nfolds, correction){
  if(lambda<0) return(Inf)
    n=length(X)
  
  funct1=function(fold){
    len=blocks[1,fold]
    range=blocks[2:len,fold]
    Xs=X[-range]
    Ws=W[-range]
    Wp=W[range]
    Ys=Y[-range]
    Yp=Y[range]
    r=length(Xs)
    Yf=funct(Xs,Ys, Ws, lambda, correction)

    Yp0=interpolate(X,Yf, fold)
    SSE= sum(Wp*(Yp0-Yp)^2)
    return(SSE)
  }

  SSEs =sapply(1:nfolds, FUN = funct1)

  return(sum(SSEs))
}

cvLambda=function(X,Y,W,lambdarange, nfolds,type, correction) {
  ln=length(lambdarange)
  n=length(X)
  
  if(type==1) {
    func=ARInt
  }  else func=SPAVInt
  
  SSE=rep(0, ln)
  for(i in 1:ln) {
    SSE[i]=cvSSE(func,lambdarange[i],X,Y,W,nfolds, correction)
  }
  
  lambda=lambdarange[which.min(SSE)]
  return(list(lambda=lambda, lambdarange=lambdarange,SSE=SSE))
  
}

findLambda=function(X,Y, W, nfolds,type,correction,power) {
  if(type==1) {
    func=ARInt
  }  else func=SPAVInt
  
  r=optimize(f=cvSSE, interval=c(0,100/length(X)^(1/8*power^2)),funct=func,X=X, Y=Y, W=W,nfolds=nfolds, correction=correction)
  return(list(lambda=r$minimum))
}

Kern1=function(){return()}
Kern2=function(){return()}


originalSPAV<-function(X0,Y0,W0=numeric(length(X0))+1,method="SPAV", cv="GCV", lambda=NULL, nfolds=10, kernel="linear", correction=T){
  N=length(Y0)
  X=X0
  Y=Y0
  W=W0
  res<-list(X=X0,Y=Y0,fitted=rep(0,N),lambda=Inf)
  class(res)<-append(class(res), "SMR")
  

  if(class(kernel)=="character" && kernel=="linear") {
    Kern1<-function(X1,X2) {
      return(K1(X1,X2,1))
    }
    Kern2<-function(X, i,j) {
      return(K2(X,i,j,1))
    }
  } else if(class(kernel)=="character" && kernel=="quadratic") {
    Kern1<-function(X1,X2) {
      return(K1(X1,X2,2))
    }
    Kern2<-function(X,i,j) {
      return(K2(X,i,j,2))
    }
  } else if(class(kernel)=="function")
  {
    Kern1<-kernel
    Kern2<-function(X,i,j) {
      return(kernel(X[i],X[j]))
    }
  }
    else{
      stop("Kernel function is misspecified")
 }
 assign("Kern1", Kern1, envir=globalenv())
 assign("Kern2", Kern2, envir=globalenv())
 
  type=Inf
 
 if (is.null(cv)){} 
 else {
   if(cv=="GCV") {
    type=1
  } else if(cv=="CV"){
    type=2    
  } else stop("Wrong type of cross validation")
  
  if (cv=="GCV"||cv=="CV") {
    prepareFolds(N, nfolds, correction)
    if (is.null(lambda)) {
      result=findLambda(X,Y,W,nfolds,type,correction, ifelse(kernel=="quadratic",2,1))
      lambda=result$lambda
    } else {
      result=cvLambda(X,Y,W,lambda,nfolds,type,correction)
      lambda=result$lambda
      res$lambdarange=result$lambdarange
      res$SSE=result$SSE
    }
  } 
 }
  
  if(method=="SPAV"){
      Yr=SPAVInt(X,Y,W,lambda, correction)
      res$fitted=Yr
      res$lambda=lambda                         
    
  } else if (method=="AR"){
    Yr=ARInt(X,Y,W,lambda,correction)
    res$fitted=Yr
    res$lambda=lambda
    
   } else stop("Wrongly specified method")

  return(res)
}


plot.SMR=function(SMRobject,panel=1,...) {
  if(panel==1) {
    ind=order(SMRobject$X)
    plot(SMRobject$X,SMRobject$Y, xlab="Predictor", ylab="response")
    points(SMRobject$X[ind], SMRobject$fitted[ind], type="l",...)
  } else if(panel==2 && !is.null(SMRobject$SSE)){
    plot(SMRobject$lambdarange, SMRobject$SSE, type="b",...)
  }
}
print.SMR=function(SMRobject){
  print("Fitted values:")
  print(SMRobject$fitted)
  print("")
  print("Penalty factor is:")
  print(SMRobject$lambda)
}
  


