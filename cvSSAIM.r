"
Description:
    Cross-Validation for SSAIM.

Author:
    Jiaqi Zhang
"

source('SSAIM.r')
source('evaluate.r')
source("ParameterGroup.r")
source('OtherMethods.r')


cvSSAIM<-function(X, Y, epsilon_list, lambda_list, delta_list){
    # parameter group
    pars_group = parameterGroup(epsilon_list, lambda_list, delta_list)
    pars_group_num = length(pars_group)
    print(sprintf("Num of parameter groups : %d", pars_group_num))
    mse_records = c()
    for(i in 1:pars_group_num){
        print("--------------------")
        print(sprintf("|%d| Parameters : ", i))
        pars = pars_group[[i]]
        epsilon = pars[1]
        lambda = pars[2]
        delta = pars[3]
        print(pars)
        #estimation
        ssaim=SSAIM(X,Y,epsilon,lambda,delta,iter=20)
        estimate_F=ssaim$F
        estimatedTheta=ssaim$theta
        pred_error = predictError(Y, estimatedTheta, estimate_F)
        mse_records = c(mse_records, pred_error)
        print(sprintf("MSE : %.3f", pred_error))
    }
    best_index = which.min(mse_records)
    best_pars = pars_group[[best_index]]
    print("===================")
    print(sprintf("Best MSE : %.3f", mse_records[best_index]))
    print(sprintf("Epsilon : %f | Lambda : %f | Delta : %f", best_pars[1], best_pars[2], best_pars[3]))
    return (best_pars)
}


# ================================================================================
# ================================================================================


 # Parameter sets
epsilon_list = c(0.1, 0.2, 0.4, 0.5) # for non-deficient rank
lambda_list = c(0.1, 0.11, 0.12, 0.13, 0.14, 0.15) # for sparsity
delta_list = c(1) # for smoothness

# path = ''

# pList = c(400, 600, 800, 1000)
# for (p in pList) {
#     n = p
#     # s = as.integer(p*0.5)
#     s = 10
#     print(paste(replicate(40, "="), collapse = ""))
#     print(sprintf("|%d vars| %d samples ; %d sparsity", p, n, s))
#     X = as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
#     Y = as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
#     F = as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
#     # Cross-Validation
#     res = cvSSAIM(X, Y, epsilon_list, lambda_list, delta_list)
#     print(res)
# }

# path = ''
# p = 200
# # s = 100
# s = 10
# nList = c(50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200)
# for (n in nList) {
#     print(paste(replicate(40, "="), collapse = ""))
#     print(sprintf("|%d samples| %d vars ; %d sparsity", n, p, s))
#     X = as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
#     Y = as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
#     F = as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
#     # Cross-Validation
#     res = cvSSAIM(X, Y, epsilon_list, lambda_list, delta_list)
#     print(res)
# }

path = ''
p = 200
n = 200
sList = c(20, 40, 60, 80, 100, 120, 140, 160, 180, 200)
for (s in sList) {
    print(paste(replicate(40, "="), collapse = ""))
    print(sprintf("|%d sparsity| %d vars ; %d samples", s, p, n))
    X = as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    Y = as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    F = as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    # Cross-Validation
    res = cvSSAIM(X, Y, epsilon_list, lambda_list, delta_list)
    print(res)
}

# path = ''
# p = 200
# n = 200
# s = 10
# X = as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
# Y = as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
# F = as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
# theta = as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
# true_theta_index = c()
# for (each in 1:length(theta)){
#     if(theta[each] != 0){
#         true_theta_index = c(true_theta_index, each)
#     }
# }
# print(true_theta_index)

# epsilon = 0.5
# lambda = 0.1
# delta = 1
# ssaim=SSAIM(X,Y,epsilon,lambda,delta,iter=200, pavIter = 100)
# estimate_F=ssaim$F
# estimatedTheta=ssaim$theta
# pred_error = predictError(Y, estimatedTheta, estimate_F)
# print(sprintf("MSE : %f", pred_error))

# slim=SLIM(X,Y,lambda,iter=50)

# # estimation
# estimate_F=slim$F
# estimatedTheta=slim$theta