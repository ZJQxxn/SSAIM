"
Description:  
    Implementation of Isotonic Lasso, SLIM and Cyclic SPAV. 

Programmer: 
    Jiaqi Zhang

Reference:
    ``Isotonic regression meets LASSO''
    -- https://projecteuclid.org/euclid.ejs/1550632213

    ``Sparse Linear Isotonic Models''
    -- http://proceedings.mlr.press/v84/chen18d/chen18d.pdf

    ``A smoothed monotonic regression via L2 regularization''
    -- https://link.springer.com/article/10.1007/s10115-018-1201-2

Functions:
    uniqueProject: Project vector for uniqueness.
    sampleCorr: Calculate the Kendall's tau correlation matrix and vector.

    **IsoLasso**: Main function for Isotonic Lasso.
    **SLIM**: Main function for solving SLIM.
    **CyclicSPAV**: Main function for Cyclic SPAV.

Update history: 
    -- 2010/04/07: Finish 4 functions, 'uniqueProject', 'sampleCorr', 'IsoLasso', 'SLIM', 'CyclicSPAV'.
"

library('VGAM')
library('MTS')
library('glmnet')
library('fastclime')
library('foreach')
source('SPAV.r')
source('originalSPAV.r')

uniqueProject<-function(vec){
    "
    Project estimated F.
    "
    vec=vec-sum(vec)/length(vec)
    vec=min(c(sqrt(length(vec))/norm(vec,type="2"),1))*vec
    return (vec)
}

sampleCorr<-function(X,Y){
    "
    Compute the Kendall's tau correlation matrix and vector.
    "
    varNum=length(X[1,])
    sampleNum=length(Y)
    corr_matrix=matrix(0,varNum,varNum)
    # Compute sample Kendall's tau correlation matrix (p times p)
    for (row in  1:varNum){
        for (col in 1:varNum){
            corr_matrix[row,col]=sin((pi/2)*kendall.tau(X[,row],X[,col]))
        }
    # Compute sample Kendall's tau correlation vector (p)
    corr_vector=matrix(0,1,varNum)
    }
    for (index in 1:varNum){
        corr_vector[index] = sin((pi/2)*kendall.tau(X[,index],Y))
    }

    return (list(mat=corr_matrix,vec=corr_vector))
}

IsoLasso<-function(X,Y,sparsePar=0.1){
    "
    Main function for Isotonic Lasso.
    =================================
    Reference:
        ``Isotonic regression meets LASSO''
        -- https://projecteuclid.org/euclid.ejs/1550632213
    "
    varNum=length(X[1,])
    sampleNum=length(Y)
    #start=proc.time()
    sampCov=cov(X,X)
    # Catch the rank-deficient problem in high-dim settings
    out=tryCatch(
        {
            rootCov=msqrt(sampCov)$mtxsqrt
        },
        error=function(cond) {
            message("Error:Can't resolve the root of sample covariance!")
            return(0)
        },
        warning=function(cond) {
            message("Warning:Can't resolve the root of sample covariance!")
            return(0)
        }
    )
    if(length(out)==1){
        return (0)
    }
    else{
        rootCov=out
    }
    # Lasso to obtain theta
    lassoModel=glmnet(X,Y,intercept=FALSE)
    lassoTheta=as.matrix(coef(lassoModel,s=sparsePar)[1:varNum+1])
    #print(proc.time()-start)
    estimateTheta=matrix(0,varNum,1)
    for (i in 1:varNum){
        estimateTheta[i]=lassoTheta[i]
    }
    estimateTheta=estimateTheta/norm(rootCov %*% estimateTheta,type="F")
    # PAV step
    latentX=X %*% estimateTheta
    reorder=sort(latentX,index.return=TRUE)$ix
    isoModel=isoreg(Y[reorder])
    #plot(isoModel)
    estimateY=isoModel$yf
    reorderY=matrix(0,1,sampleNum)
    for (i in 1:sampleNum){
        reorderY[i]=estimateY[reorder[i]]
    }
    return (list(theta=estimateTheta,estimatedY=reorderY,order=reorder))
}

SLIM<-function(X,Y,sparsePar=0.1,iter=100){
    "
    Main function for solving SLIM.
    ================================
    Reference:
        ``Sparse Linear Isotonic Models''
        -- http://proceedings.mlr.press/v84/chen18d/chen18d.pdf
    "
    varNum=length(X[1,])
    sampleNum=length(Y)
    # Kendall's tau correlation
    #start=proc.time()
    kendall=sampleCorr(X,Y)
    corrMat=kendall$mat
    corrVec=kendall$vec
    # Dantzig selector
    selector=dantzig(corrMat,t(corrVec),lambda=sparsePar,nlambda=200)
    beta0=selector$BETA0
    # print(beta0)
    #print(proc.time()-start)
    estimateTheta=beta0[,length(beta0[1,])]
    # rm(selector)
    # rm(beta0)
    # Block-wise coordinate descent algorithm
    estimeteF=matrix(0,sampleNum,varNum)
    curEpoch=1
    while(curEpoch<=iter){
        for (var in 1: varNum){
            if(estimateTheta[var]!=0){
                # Calculate the residue
                residue=Y
                for(anotherVar in 1:varNum){
                    if(anotherVar!=var){
                        residue=residue-estimateTheta[anotherVar]*estimeteF[,anotherVar]
                    }
                }
                # PAV step
                reorder=order(X[,var])
                # isoModel=isoreg(Y[reorder])
                isoModel=isoreg((residue/estimateTheta[var])[reorder])
                # plot(isoModel)
                estimateY=isoModel$yf
                # reorderY=matrix(0,1,sampleNum)
                # for (i in 1:sampleNum){
                #     reorderY[i]=estimateY[reorder[i]]
                # }
                estimateY=uniqueProject(estimateY)
                estimeteF[,var][reorder]=estimateY

                # plot(X[,var][reorder], reorderY)
            }
        }
        curEpoch=curEpoch+1
    }
    return (list(theta=estimateTheta,F=estimeteF,order=reorder, selector=selector))
}

SLIM_testt<-function(X,Y,selector,sparsePar=0.1,iter=100){
    "
    Main function for solving SLIM.
    ================================
    Reference:
        ``Sparse Linear Isotonic Models''
        -- http://proceedings.mlr.press/v84/chen18d/chen18d.pdf
    "
    varNum=length(X[1,])
    sampleNum=length(Y)
    # Kendall's tau correlation
    #start=proc.time()
    kendall=sampleCorr(X,Y)
    corrMat=kendall$mat
    corrVec=kendall$vec
    # Dantzig selector
    # selector=dantzig(corrMat,t(corrVec),lambda=sparsePar,nlambda=200)
    beta0=selector$BETA0
    # print(beta0)
    #print(proc.time()-start)
    estimateTheta=beta0[,length(beta0[1,])]
    # rm(selector)
    # rm(beta0)
    # Block-wise coordinate descent algorithm
    estimeteF=matrix(0,sampleNum,varNum)
    curEpoch=1
    # print('======= stop1 ========')
    while(curEpoch<=iter){
        for (var in 1: varNum){
            if(estimateTheta[var]!=0){
                # Calculate the residue
                residue=Y
                for(anotherVar in 1:varNum){
                    if(anotherVar!=var){
                        residue=residue-estimateTheta[anotherVar]*estimeteF[,anotherVar]
                    }
                }
                # PAV step
                reorder=order(X[,var])
                # isoModel=isoreg(Y[reorder])
                isoModel=isoreg((residue/estimateTheta[var])[reorder])
                # plot(isoModel)
                estimateY=isoModel$yf
                # reorderY=matrix(0,1,sampleNum)
                # for (i in 1:sampleNum){
                #     reorderY[i]=estimateY[reorder[i]]
                # }
                estimateY=uniqueProject(estimateY)
                estimeteF[,var][reorder]=estimateY

                # plot(X[,var][reorder], reorderY)
            }
        }
        curEpoch=curEpoch+1
    }
    # print('========stop 2=========')
    return (list(theta=estimateTheta,F=estimeteF,order=reorder))
}

ParSLIM<-function(X,Y,sparsePar=0.1,iter=100){
    "
    Main function for solving SLIM.
    ================================
    Reference:
        ``Sparse Linear Isotonic Models''
        -- http://proceedings.mlr.press/v84/chen18d/chen18d.pdf
    "
    varNum=length(X[1,])
    sampleNum=length(Y)
    # Kendall's tau correlation
    #start=proc.time()
    kendall=sampleCorr(X,Y)
    corrMat=kendall$mat
    corrVec=kendall$vec
    # Dantzig selector
    selector=dantzig(corrMat,t(corrVec),lambda=sparsePar,nlambda=200)
    beta0=selector$BETA0
    # print(beta0)
    #print(proc.time()-start)
    estimateTheta=beta0[,length(beta0[1,])]
    rm(selector)
    rm(beta0)
    # Block-wise coordinate descent algorithm
    estimeteF=matrix(0,sampleNum,varNum)
    varSet = c()
    for(i in 1:varNum){
        if(estimateTheta[i]!=0){
            varSet = c(varSet, i)
        }
    }
    foreach(var=varSet)%dopar%
        eachVar(var, Y, varSet,estimateTheta,estimeteF,X, iter)

    return (list(theta=estimateTheta,F=estimeteF,order=reorder))
}

eachVar<-function(var, Y, varSet,estimateTheta,estimeteF,X,iter){
    # Calculate the residue
    curEpoch=1
    while(curEpoch<=iter){
        residue=Y
        for(anotherVar in varSet){
            if(anotherVar!=var){
                residue=residue-estimateTheta[anotherVar]*estimeteF[,anotherVar]
            }
        }
        # PAV step
        reorder=order(X[,var])
        # isoModel=isoreg(Y[reorder])
        isoModel=isoreg((residue/estimateTheta[var])[reorder])
        # plot(isoModel)
        estimateY=isoModel$yf
        # reorderY=matrix(0,1,sampleNum)
        # for (i in 1:sampleNum){
        #     reorderY[i]=estimateY[reorder[i]]
        # }
        estimateY=uniqueProject(estimateY)
        estimeteF[,var][reorder]=estimateY
        curEpoch=curEpoch+1
    }
}

CyclicSPAV<-function(X,Y,smoothPar=0.1,iter=1,pav_iter=100){
    "
    Main function for Cyclic SPAV.
    ==========================
    Reference:
        ``A smoothed monotonic regression via L2 regularization''
        -- https://link.springer.com/article/10.1007/s10115-018-1201-2
    "
    varNum=length(X[1,])
    sampleNum=length(Y)
    # Block-wise coordinate descent algorithm
    estimeteF=matrix(0,sampleNum,varNum)
    curEpoch=1
    while(curEpoch<=iter){
        for (var in 1: varNum){
            # Calculate the residue
            residue=Y
            for(anotherVar in 1:varNum){
                if(anotherVar!=var){
                    residue=residue-1*estimeteF[,anotherVar]
                }
            }
            # PAV step
            reorder=sort(X[,var],index.return=TRUE)$ix
            # TODO: Change to original SPAV
            sample_indices = order(X[,var])
            estimateY=SPAV(X[reorder][sample_indices],residue[sample_indices],smoothPar,pav_iter)#$Y
            #estimateY=originalSPAV(X[reorder],residue,lambda=smoothPar)$Y
            reorderY=matrix(0,1,sampleNum)
            for (i in 1:sampleNum){
                reorderY[i]=estimateY[reorder[i]]
            }
            reorderY=uniqueProject(reorderY) # make the solution unique
            estimeteF[,var]=reorderY
        }
        curEpoch=curEpoch+1
    }
    return (list(F=estimeteF,order=reorder))
}

"
x=as.matrix(read.csv(sprintf('data/X/X_200p_200n.csv'),header=FALSE))
y=as.matrix(read.csv(sprintf('data/Y/Y_200p_200n.csv'),header=FALSE))
#start=proc.time()
print('LAsso:')
model=IsoLasso(x,y)
print('Dantzig:')
model=SLIM(x,y)
print('EE:')
model=SSAIM(x,y)
#print(proc.time()-start)
#if(length(model)==1){
#    print('High-dim rank-deficient.')
#} else {
    #print(length(model))
#    print(integralPreError(t(y),model$estimatedY))
#}
#model=SLIM(x,y)
#start=proc.time()
#model=CyclicSPAV(x,y)
#print(proc.time()-start)
"