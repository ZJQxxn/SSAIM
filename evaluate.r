"
Description:  
    Evaluation for the performance of SSAIM.

Programmer: 
    Jiaqi Zhang

Functions:
    thetaError: Error of estimated theta.
    predictError: Error of estimated response.
    smoothness: Evaluate the smoothness.

Update history: 
    -- 2019/04/07: Transform from Python to R.
"
    
Ferror<-function(F, estimate_F){
    return (norm(as.matrix(estimate_F)-as.matrix(F),type="F") / norm(as.matrix(F),type="F"))
}  

thetaError<-function(theta,estimateTheta){
    "
    Evaluate error of estimated theta.
    "
    # return (norm(as.matrix(theta)-as.matrix(estimateTheta),type="F") / norm(as.matrix(theta),type="F"))
    return (norm(as.matrix(theta)-as.matrix(estimateTheta),type="I"))
}

predictError<-function(Y,estimatedTheta,estimatedF){
    "
    Mean l2-norm error of estimated responses.
    "
    return (norm(as.matrix(estimatedF %*% t(estimatedTheta) - Y),type="F")^2/length(Y))
}

yError<-function(Y,estimatedY){
    "
    Mean l2-norm error of estimated responses.
    "
    return (norm(as.matrix(estimatedY - Y),type="F")^2/length(Y))
}

integralPreError<-function(Y,estimatedY){
    "
    Mean l2-norm error of estimated responses.
    "
    return (norm(Y-estimatedY,type="F")^2/length(Y))
}

smoothness<-function(estimateTheta,estimateF){
    "
    Evaluated smoothness of the estimated response 'Y'.
    "
    sampNum=length(estimateF[,1])
    #varSet=c()
    smooth=0
    count=0
    Y=matrix(0,1,sampNum)
    for (s in 1:sampNum){ 
        for (var in 1:length(estimateTheta)){
            if(estimateTheta[var]!=0){
                count=count+1
                Y[s]=Y[s]+(estimateTheta[var]*estimateF[s,var])
            }
        }
    }
    for (s in 1:(sampNum-1)){
        smooth=smooth+(Y[s+1]-Y[s])^2
    }
    # Avoid 0/0.
    if(count==0){
        count=1
    }
    return (smooth/sampNum)
}

Fsmoothness<-function(estimateTheta,estimateF){
    "
    Evaluated smoothness of the estimated latent variable 'F'.
    "
    sampNum=length(estimateF[,1])
    varNum=length(estimateTheta)
    #varSet=c()
    smooth=0
    count=0
    for (var in 1:varNum){
        if(estimateTheta[var]!=0){
            count=count+1
            for (i in 1:(sampNum-1)){
                smooth=smooth+(estimateF[i+1,var]-estimateF[i,var])^2
            }
        }
    }
    # Avoid 0/0.
    if(count==0){
        count=1
    }
    # return (smooth/(count*sampNum))
    return (smooth/(sampNum))
}

ySmoothness<-function(estimateY){
    "
    Evaluated smoothness of the estimated latent variable 'F'.
    "
    sampNum=length(estimateY)
    smoothness=0
    for (i in 1:(sampNum-1)){
        smoothness=smoothness+(estimateY[i+1]-estimateY[i])^2
    }
    return (smoothness/sampNum)
}
    