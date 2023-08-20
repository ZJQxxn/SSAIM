"
Description: 
    EE-Ridge variable selector utilizing Kendall's tau correlations.
    
Programmer:
    Jiaqi Zhang

Reference:
    ``Elementary Estimators for High-Dimensional Linear Regression''
    -- http://proceedings.mlr.press/v32/yangc14.pdf

Functions:
    sampleCorr: Compute sample Kendall's tau correlation;
    softThreshold: Soft-threshold;
    printTheta: Print estimated theta and the sparsity.
    
    **varSelect**: EE-Ridge selector;

"
#TODO: Need reference

library('matlab')
library('VGAM')

kendallCorr<-function(X,Y){
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

softThreshold<-function(estimateTheta,threshold){
    "
    Soft-thresholding the estimated theta.
    "
    for (i in 1:length(estimateTheta)){
        if (abs(estimateTheta[i])<threshold){
            estimateTheta[i]=0
        }
        else{
            estimateTheta[i]=estimateTheta[i]-(sign(estimateTheta[i])*threshold)
        }
    }
    return (estimateTheta)
}

varSelect<-function(X,Y,epsilon=1,threshold=0.1){
    "
    Variable selection.
    "
    varNum=length(X[1,])
    sampleNum=length(Y)
    kendall=kendallCorr(X,Y)
    #print(kendall)
    corr_matrix=kendall$mat
    corr_vector=kendall$vec
    estimateTheta=solve(corr_matrix + epsilon*eye(varNum))
    estimateTheta=t(estimateTheta %*% t(corr_vector))
    #print(estimateTheta)
    responseVar=sqrt(var(Y))
    print(sprintf('Square root of the sample variance is %f',responseVar[1]))
    estimateTheta = softThreshold(estimateTheta,threshold)
    #print(estimateTheta)
    # estimateTheta = responseVar[1]*estimateTheta #TODO: check this
    varSet=c()
    for (index in 1:varNum){
        if (estimateTheta[index]!=0){
            varSet=union(varSet,index)
        }
    }
    return (list(theta=estimateTheta,vars=varSet))
}

printTheta<-function(estimatedTheta){
    "
    Print the estimated theta.
    :return: None.
    "
    nonzero=c()
    for (i in 1:length(estimatedTheta)){
        if (estimatedTheta[i]!=0){
            nonzero<-c(nonzero,i)
        }
    }
    print(sprintf('There are %d useful variables. Indices are',length(nonzero)))
    print(nonzero)
}

"
p=500
path='D:/programming/Python/SSAIM/'
thr=0.07 # 10 vars are selected 
n=500
X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn50s.csv',path,p,n),header=FALSE))
Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn50s.csv',path,p,n),header=FALSE))
res=varSelect(X,Y,1,thr)
write.table(as.matrix(res$theta),file=sprintf('data/estimateTheta/estimateTheta_%dp_%dn50s.csv',p,n),row.names=FALSE,col.names=FALSE,sep=',')
write.table(res$vars,file=sprintf('data/estimateTheta/varSet_%dp_%dn50s.csv',p,n),row.names=FALSE,col.names=FALSE,sep=',')
printTheta(res$theta)
"
