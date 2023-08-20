"
Description:
    Test SSAIM methods.

Author:
    Jiaqi Zhang
"
source('SSAIM_jiaqi/SSAIM.r')
source('SSAIM_jiaqi/OtherMethods.r')
source('SSAIM_jiaqi/evaluate.r')

library(mboost)
library(scar)

# Initialization
#path='D:/programming/Python/SSAIM/'
path='SSAIM_jiaqi/'
count = 1
epsilon = 1
lambda = 0.06
delta_list = c(1e-8, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10)

p = 200
n = 100
s = 10

# 200: c(1, 0.06, 10), 
# pList = c()
# parameter_list = list(c(1, 0.06, 10))

# Reading data
print('Begin reading data....')
theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
X_full=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
Y_full=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
F_full=as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
print('Finish reading data!')

print('Begin splitting data....')
varNum=length(X_full[1,])
sampNum=length(X_full[,1])

drop_ratio = 0.1
drop_index = sample(2:(sampNum-1), ceiling(sampNum * drop_ratio), replace=FALSE)
drop_index = c(29, 57, 42, 21, 56, 33, 69, 83, 3, 47)
print(drop_index)
X_drop = X_full[drop_index, ]
Y_drop = Y_full[drop_index, ]
F_drop = F_full[drop_index, ]

X = X_full[-drop_index, ]
Y = Y_full[-drop_index, ]
F = F_full[-drop_index, ]

save_tag = "0002"

write.csv(X_drop, paste(sprintf("smoothtest_result_synthetic/%s_%dp_%dn_%ds_X_drop.csv", save_tag, p, n, s)))
write.csv(X, paste(sprintf("smoothtest_result_synthetic/%s_%dp_%dn_%ds_X.csv", save_tag, p, n, s)))
write.csv(Y_drop, paste(sprintf("smoothtest_result_synthetic/%s_%dp_%dn_%ds_Y_drop.csv", save_tag, p, n, s)))
write.csv(Y, paste(sprintf("smoothtest_result_synthetic/%s_%dp_%dn_%ds_Y.csv", save_tag, p, n, s)))
write.csv(F_drop, paste(sprintf("smoothtest_result_synthetic/%s_%dp_%dn_%ds_F_drop.csv", save_tag, p, n, s)))
write.csv(F, paste(sprintf("smoothtest_result_synthetic/%s_%dp_%dn_%ds_F.csv", save_tag, p, n, s)))

n_drop = length(X[1, ])
p_drop = length(X[ , 1])
# print(X)
print('Finish splitting data!')

"
    =========================================
          Varying variable number (n=p)
    =========================================
"
mse_history = c()
smoothness_history = c()
for (delta_index in 1:length(delta_list)){
    delta = delta_list[delta_index]
    print(sprintf('========= delta = %f ============', delta))

    #========= SSAIM ==============
    print('-------- SSAIM ---------')
    totalTime=proc.time()-proc.time()
    totalMSe = 0.0
    totalSmoothness = 0.0
    total_theta_error = 0.0
    total_F_error = 0.0
    for(i in 1:count){
        start=proc.time()
        ssaim=SSAIM(X,Y,epsilon,lambda,delta)
        totalTime=totalTime+(proc.time()-start)
        estimate_F=ssaim$F
        estimatedTheta=ssaim$theta
        selector=ssaim$selector
        estimate_Y=estimate_F %*% t(estimatedTheta)
        # print('Theta error:')
        # print(thetaError(t(theta),estimatedTheta))
        # print('Prediction error:')
        # print(predictError(Y,estimatedTheta,estimate_F))
        # print('Smoothness:')
        # print(Fsmoothness(estimatedTheta,estimate_F))
        ssaim_y = estimate_F %*% t(estimatedTheta)
        ssaim_F = estimate_F
        print(sprintf('======= testing result (delta = %f)=======', delta))
        # smoothness test
        ssaim_drop = SSAIM_testt(X_drop, Y_drop, selector, epsilon, lambda, delta)
        estimate_F_drop=ssaim_drop$F
        estimatedTheta_drop=ssaim_drop$theta
        estimate_Y_drop=estimate_F_drop %*% t(estimatedTheta_drop)
        print('Theta error:')
        print(thetaError(t(theta),estimatedTheta_drop))
        print('Prediction error:')
        print(predictError(Y_drop,estimatedTheta_drop,estimate_F_drop))
        print('Smoothness:')
        print(Fsmoothness(estimatedTheta_drop,estimate_F_drop))
        mse_history = c(mse_history, predictError(Y_drop,estimatedTheta_drop,estimate_F_drop))
        smoothness_history = c(smoothness_history, Fsmoothness(estimatedTheta_drop,estimate_F_drop))
        write.csv(estimate_F, paste(sprintf("smoothtest_result_synthetic/%s_%dp_%dn_%ds_%fdelta_EstimateF.csv", save_tag, p, n, s, delta)))
        write.csv(estimate_F_drop, paste(sprintf("smoothtest_result_synthetic/%s_%dp_%dn_%ds_%fdelta_EstimateF_drop.csv", save_tag, p, n, s, delta)))
        write.csv(estimate_Y, paste(sprintf("smoothtest_result_synthetic/%s_%dp_%dn_%ds_%fdelta_EstimateY.csv", save_tag, p, n, s, delta)))
        write.csv(estimate_Y_drop, paste(sprintf("smoothtest_result_synthetic/%s_%dp_%dn_%ds_%fdelta_EstimateY_drop.csv", save_tag, p, n, s, delta)))
        write.csv(estimatedTheta, paste(sprintf("smoothtest_result_synthetic/%s_%dp_%dn_%ds_%fdelta_EstimateTheta.csv", save_tag, p, n, s, delta)))
    }
}
print("Finished saving low-dimensaionl data.")
print("====================================\n")
print(mse_history)
print(smoothness_history)
"
drop_index = c(29, 57, 42, 21, 56, 33, 69, 83, 3, 47)
[1] 0.004162091 0.004162089 0.004162070 0.004161874 0.004160139 0.004156289 0.004340529
[1] 25.18536 25.18532 25.18500 25.18177 25.15029 24.89002 24.11673
"
# 0.001, 0.01, 0.1, 1.0
# 0.2437, 0.2410, 0.2407, 0.2513