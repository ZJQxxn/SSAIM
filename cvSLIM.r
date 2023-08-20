"
Description:
    Cross-Validation for SLIM.

Author:
    Jiaqi Zhang
"

source('OtherMethods.r')
source('evaluate.r')
source("ParameterGroup.r")


cvSLIM<-function(X, Y, lambda_list){
    # parameter group
    pars_group = lambda_list
    pars_group_num = length(pars_group)
    print(sprintf("Num of parameter groups : %d", pars_group_num))
    mse_records = c()
    for(i in 1:pars_group_num){
        print("--------------------")
        print(sprintf("|%d| Parameters : ", i))
        pars = pars_group[[i]]
        lambda = pars
        print(pars)
        #estimation
        slim=SLIM(X,Y,lambda,iter=20)
        estimate_F=slim$F
        estimatedTheta=slim$theta
        pred_error = predictError(Y, t(estimatedTheta), estimate_F)
        mse_records = c(mse_records, pred_error)
        print(sprintf("MSE : %f", pred_error))
    }
    best_index = which.min(mse_records)
    best_pars = pars_group[[best_index]]
    print("===================")
    print(sprintf("Best MSE : %f", mse_records[best_index]))
    print(sprintf("Lambda : %f", best_pars))
    return (best_pars)
}


# ================================================================================
# ================================================================================

# Parameter sets
lambda_list = c(1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1) # for sparsity

path = ''
pList = c(200, 400, 600, 800, 1000)
pList = c()
for (p in pList) {
    n = p
    # s = as.integer(p*0.5)
    s = 10
    print(paste(replicate(40, "="), collapse = ""))
    print(sprintf("|%d vars| %d samples ; %d sparsity", p, n, s))
    X = as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    Y = as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    F = as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    # Cross-Validation
    res = cvSLIM(X, Y, lambda_list)
    print(res)
}


p = 200
s = 100
s = 10
nList = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200)
nList = c()
for (n in nList) {
    print(paste(replicate(40, "="), collapse = ""))
    print(sprintf("|%d samples| %d vars ; %d sparsity", n, p, s))
    X = as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    Y = as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    F = as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    # Cross-Validation
    res = cvSLIM(X, Y, lambda_list)
    print(res)
}

p = 200
n = 200
sList = c(20, 40, 60, 80, 100, 120, 140, 160, 180, 200)
# sList = c()
for (s in sList) {
    print(paste(replicate(40, "="), collapse = ""))
    print(sprintf("|%d sparsity| %d vars ; %d samples", s, p, n))
    X = as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    Y = as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    F = as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    # Cross-Validation
    res = cvSLIM(X, Y, lambda_list)
    print(res)
}
