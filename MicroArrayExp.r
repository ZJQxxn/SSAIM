"
Description:
    Cross-Validation on real-world dataset.

Author:
    Jiaqi Zhang
"

source('SSAIM.r')
source("OtherMethods.r")
source('evaluate.r')
source("ParameterGroup.r")

library(foreign)


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

library(GEOquery)

gds80606 <- getGEO(filename='./dataset/GSE42866_family.soft.gz')


# for (j in 1:3){
#     X[,j] = X[,j] / norm(X[,j], type = "2")
# }
# Y = Y / norm(Y, type = "2")

# ssaim=SSAIM(X,Y,0.5,0.01,0.0,iter=20)
# estimate_F=ssaim$F
# estimatedTheta=ssaim$theta
# pred_error = predictError(Y, estimatedTheta, estimate_F)
# print(sprintf("MSE : %.3f", pred_error))

# slim=SLIM(X,Y,0.01,iter=20)
# # estimation
# estimate_F=slim$F
# estimatedTheta=slim$theta
# estimatedTheta = t(estimatedTheta) #TODO: check this
# pred_error = predictError(Y,estimatedTheta,estimate_F)
# print(sprintf("SLIM MSE : %f", pred_error))

# scar_model= scar(X, Y, shape = c("cvxin", "cvxin", "cvxin"), family = gaussian())
# # Predicted Y
# scar_predicted_Y = predict(scar_model)
# mse_val = yError(Y, scar_predicted_Y)
# print(sprintf("SCAR MSE : %f", mse_val))


# # Parameter sets
# epsilon_list = c(1e-2, 5e-2, 1e-1, 5e-1, 1) # for non-deficient rank
# lambda_list = c(1e-2, 5e-2, 1e-1, 5e-1, 1) # for sparsity
# delta_list = c(1e-2, 5e-2, 1e-1, 5e-1, 1) # for smoothness

# res = cvSSAIM(X, Y, epsilon_list, lambda_list, delta_list)
# print(res)
