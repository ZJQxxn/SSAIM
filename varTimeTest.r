"
Description:  
    Compare SSAIM and SLIM on the time cost for variable selection.

Programmer: 
    Jiaqi Zhang

Update history: 
    --2019/04/30: Create test codes and do experiments.
"
source('VarSelect.r')
source('OtherMethods.r')
source('evaluate.r')

SLIMVarSelect<-function(X,Y,sparsePar=0.1){
    varNum=length(X[1,])
    sampleNum=length(Y)
    # Kendall's tau correlation
    #start=proc.time()
    kendall=sampleCorr(X,Y)
    corrMat=kendall$mat
    corrVec=kendall$vec
    # Dantzig selector
    selector=dantzig(corrMat,t(corrVec),lambda=sparsePar)
    beta0=selector$BETA0
}

# Initialization
#path='D:/programming/Python/SSAIM/'
path=''
count=1
epsilon=1
lambda=0.1
pList=c(200,400,600,800,1000)
#pList=c()

"
    =========================================
          Varying variable number (n=p)
    =========================================
"
for (p in pList){
    print(sprintf('========= %d Variables when n=p ============',p))
    # Reading data
    n=p
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn.csv',path,p,n),header=FALSE))
    X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn.csv',path,p,n),header=FALSE))
    Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn.csv',path,p,n),header=FALSE))
    print('Finish reading data!')
    
    #========= SSAIM ==============
    print('-------- SSAIM ---------')
    totalTime=proc.time()-proc.time()
    for(i in 1:count){
        start=proc.time()
        ssaim=varSelect(X,Y,epsilon,lambda)
        totalTime=totalTime+(proc.time()-start)
        #printTheta(ssaim$estimateTheta)
    }
    print(totalTime/count)
    
    #========= SLIM ==============
    print('-------- SLIM ---------')
    totalTime=proc.time()-proc.time()
    for(i in 1:count){
        start=proc.time()
        slim=SLIMVarSelect(X,Y,lambda)
        totalTime=totalTime+(proc.time()-start)
    }
    print(totalTime/count)
}

