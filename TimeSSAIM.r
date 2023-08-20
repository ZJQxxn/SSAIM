"
Description:  
    Implementation of Smooth and Sparse Additive Isotonic Modeel (SSAIM).

Programmer: 
    Jiaqi Zhang

Reference:

Functions:
    uniqueProject: Project the estimated F to make it unique.

    **SSAIM**: Main function for solving SSAIM.

Update history: 
    --2019/04/07: Transform from Python to R.  
"
#TODO: Need reference

source('VarSelect.r')
source('SPAV.r')
source('originalSPAV.r')

library(foreach)


uniqueProject<-function(vec){
    "
    Project estimated F.
    "
    vec=vec-sum(vec)/length(vec)
    vec=min(c(sqrt(length(vec))/norm(vec,type="2"),1))*vec
    return (vec)
}

SSAIM<-function(X,Y,epsilon=0.1,sparsePar=0.1,smoothPar=0.1,iter=1,pavIter=100){
    "
    Main function for solving SSAIM.
    :param X: Predictors matrix.
    :param Y: Responses vector.
    :param sparsePar: Parameter for sparsity (default=0.1).
    :param smoothPar: Parameter for smoothness (default=0.1).
    :param iter: Number of iterations (default=10).
    :return: Estimated theta, estimated F.
    "
    #TODO: Some checks to determin excecute for additive mdoel

    # Initialization
    sampNum = length(X[,1])
    varNum = length(X[1,])
    varSet = c()
    estimatedTheta = NULL
    F = matrix(0,sampNum,varNum)
    # Variable selection
    #start=proc.time()
    selector=varSelect(X,Y,epsilon,sparsePar)
    estimatedTheta=selector$theta
    varSet=selector$vars
    #print(proc.time()-start)
    #rm(selector)
    printTheta(estimatedTheta)
    print('Finish variable selection!')
    # Block-wise coordinate descent algorithm
    #TODO: Stop criteria
    for (t in 1:iter){ 
        var_sample_indices = list()
        foreach(index=1:length(varSet)) %dopar%
            everyVar(varSet, index, Y, estimatedTheta, X, F, smoothPar, pavIter, var_sample_indices)
    }
    print('Finish estimating F!')
    return (list(theta=estimatedTheta,F=F, indices=var_sample_indices))
}

everyVar <- function(varSet, index, Y, estimatedTheta, X, F, smoothPar, pavIter, var_sample_indices){
    eachVar=varSet[index]
    # Compute the residue
    residue = Y
    otherVar=varSet
    otherVar=otherVar[otherVar!=eachVar]
    for (j in 1:length(otherVar)){
        jVar=otherVar[j]
        residue=residue-estimatedTheta[jVar]*F[,jVar]
    }
    # TODO: Modify the SPAV class
    # SPAV algorithm for each variable

    # TODO: reorder 
    sample_indices = order(X[,eachVar])
    var_sample_indices[[index]] = sample_indices
    F[,eachVar][sample_indices] = SPAV(X[,eachVar][sample_indices], (residue/estimatedTheta[eachVar])[sample_indices], smoothPar, pavIter)
    #F[,eachVar] = originalSPAV(X[,eachVar], residue/estimatedTheta[eachVar], lambda=smoothPar)$Y
    F[,eachVar] = uniqueProject(F[,eachVar])
}
