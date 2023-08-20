"
Description:  
    Utilize SSAIM on Boston housing dataset.

Programmer: 
    Jiaqi Zhang

Update history: 
    --2019/04/16: Create codes and finish the experiment.
"
source('SSAIM_jiaqi/SSAIM.r')
source('SSAIM_jiaqi/OtherMethods.r')
source('SSAIM_jiaqi/evaluate.r')

library(mboost)
library(scar)

print('Begin reading data....')
theta=as.matrix(read.csv(sprintf('dataset/highdim-housing_theta.csv'),header=FALSE))
X_full=as.matrix(read.csv(sprintf('dataset/highdim-housing_X.csv'),header=FALSE))
Y_full=as.matrix(read.csv(sprintf('dataset/housing_Y.csv'),header=FALSE))
print('Finish reading data!')
print('================ Wine Data ==================')
varNum=length(X_full[1,])
sampNum=length(X_full[,1])

drop_ratio = 0.2
drop_index = sample(2:(sampNum-1), ceiling(sampNum * drop_ratio), replace=FALSE)
print(drop_index)
# drop_index = c(77, 9, 34, 90, 29, 61, 86, 17, 99, 11, 36, 51, 10, 94, 83, 88, 47, 50, 15, 37)
X_drop = X_full[drop_index, ]
Y_drop = Y_full[drop_index, ]

X = X_full[-drop_index, ]
Y = Y_full[-drop_index, ]

save_tag = "0006"

write.csv(X_drop, paste("smoothtest_result_boston/X_drop", "_", save_tag, ".csv", sep = "", collapse = ""))
write.csv(X, paste("smoothtest_result_boston/X", "_", save_tag, ".csv", sep = "", collapse = ""))
write.csv(Y_drop, paste("smoothtest_result_boston/Y_drop", "_", save_tag, ".csv", sep = "", collapse = ""))
write.csv(Y, paste("smoothtest_result_boston/Y", "_", save_tag, ".csv", sep = "", collapse = ""))

n = length(X[1, ])
p = length(X[ , 1])


save_f = 1
save_theta = 0
# save_path = "smoothtest_result_boston/"

count=1
epsilon=0.6
lambda=0.03
delta=2

#========= SSAIM ==============
print('-------- SSAIM ---------')
totalTime=proc.time()-proc.time()
for(i in 1:count){
    start=proc.time()
    ssaim=SSAIM(X,Y,epsilon,lambda,delta)
    totalTime=totalTime+(proc.time()-start)
}
print(totalTime/count)
# Evaluation
estimate_F=ssaim$F
estimatedTheta=ssaim$theta
selector=ssaim$selector
print('Theta error:')
print(thetaError(t(theta),estimatedTheta))
print('Prediction error:')
print(predictError(Y,estimatedTheta,estimate_F))
print('Smoothness:')
print(Fsmoothness(estimatedTheta,estimate_F))
ssaim_y = estimate_F %*% t(estimatedTheta)
ssaim_F = estimate_F
print('======= testing result =======')
# smoothness test
ssaim_drop = SSAIM_testt(X_drop, Y_drop, selector, epsilon, lambda, delta)
estimate_F_drop=ssaim_drop$F
estimatedTheta_drop=ssaim_drop$theta
# print('Theta error:')
# print(thetaError(t(theta),estimatedTheta_drop))
print('Prediction error:')
print(predictError(Y_drop,estimatedTheta_drop,estimate_F_drop))
print('Smoothness:')
print(Fsmoothness(estimatedTheta_drop,estimate_F_drop))
estimate_y_drop = estimate_F_drop %*% t(estimatedTheta_drop)
if(save_f == 1) {
    write.csv(ssaim_y, paste("smoothtest_result_boston/ssaim_y", "_", save_tag, ".csv", sep = "", collapse = ""))
    write.csv(ssaim_F, paste("smoothtest_result_boston/ssaim_f", "_", save_tag, ".csv", sep = "", collapse = ""))
    write.csv(estimate_y_drop, paste("smoothtest_result_boston/ssaim_test_y", "_", save_tag, ".csv", sep = "", collapse = ""))
    write.csv(estimate_F_drop, paste("smoothtest_result_boston/ssaim_test_f", "_", save_tag, ".csv", sep = "", collapse = ""))
    print("F saved!")
}
if(save_theta == 1) {
    write.csv(estimatedTheta, paste("smoothtest_result_boston/ssaim_theta", "_", save_tag, ".csv", sep = "", collapse = ""))
    write.csv(estimatedTheta_drop, paste("smoothtest_result_boston/ssaim_test_theta", "_", save_tag, ".csv", sep = "", collapse = ""))
    print("Theta saved!")
}



#========= SLIM ==============
print('-------- SLIM ---------')
totalTime=proc.time()-proc.time()
for(i in 1:count){
    start=proc.time()
    slim=SLIM(X,Y,2)
    totalTime=totalTime+(proc.time()-start)
}
print(totalTime/count)
# Evaluation
estimate_F=slim$F
estimatedTheta=slim$theta
selector=slim$selector
rm(slim)
print('Theta error:')
print(thetaError(theta,estimatedTheta))
print('Prediction error:')
print(predictError(Y,t(estimatedTheta),estimate_F))
print('Smoothness:')
print(Fsmoothness(estimatedTheta,estimate_F))
slim_y = estimate_F %*% estimatedTheta
slim_F = estimate_F
print('======= testing result =======')
# smoothness test
slim_drop = SLIM_testt(X_drop, Y_drop, selector, 2)
estimate_F_drop=slim_drop$F
estimatedTheta_drop=slim_drop$theta
# print('Theta error:')
# print(thetaError(theta,estimatedTheta_drop))
print('Prediction error:')
print(predictError(Y_drop,t(estimatedTheta_drop),estimate_F_drop))
print('Smoothness:')
print(Fsmoothness(estimatedTheta_drop,estimate_F_drop))
estimate_y_drop = estimate_F_drop %*% estimatedTheta_drop
if(save_f == 1) {
    write.csv(slim_y, paste("smoothtest_result_boston/slim_y", "_", save_tag, ".csv", sep = "", collapse = ""))
    write.csv(slim_F, paste("smoothtest_result_boston/slim_f", "_", save_tag, ".csv", sep = "", collapse = ""))
    write.csv(estimate_y_drop, paste("smoothtest_result_boston/slim_test_y", "_", save_tag, ".csv", sep = "", collapse = ""))
    write.csv(estimate_F_drop, paste("smoothtest_result_boston/slim_test_f", "_", save_tag, ".csv", sep = "", collapse = ""))
    print("F saved!")
}
if(save_theta == 1) {
    write.csv(estimatedTheta, paste("smoothtest_result_boston/slim_theta", "_", save_tag, ".csv", sep = "", collapse = ""))
    write.csv(estimatedTheta_drop, paste("smoothtest_result_boston/slim_test_theta", "_", save_tag, ".csv", sep = "", collapse = ""))
    print("Theta saved!")
}


"
#========= Isotonic Lasso ==============
print('-------- Isotonic Lasso ---------')
totalTime=proc.time()-proc.time()
for(i in 1:count){
    start=proc.time()
    isolasso=IsoLasso(X,Y,lambda)
    if(length(isolasso)==1){
        break
    }
    totalTime=totalTime+(proc.time()-start)
}
if(length(isolasso)>1){
    print(totalTime/count)
    # Evaluation
    estimateY=isolasso$estimatedY
    estimatedTheta=isolasso$theta
    rm(isolasso)
    print('Prediction error:')
    print(integralPreError(t(Y),estimateY))
    print('Smoothness:')
    print(ySmoothness(estimateY))
}
"
#
# print('----------- CPAV ----------')
# print(date())
# start=proc.time()
# model=CyclicSPAV(X,Y,delta)
# print(proc.time()-start)
# estimate_F=model$F
# print('Prediction error:')
# print(predictError(Y,matrix(1,1,varNum),estimate_F))
# print('Smoothness:')
# print(Fsmoothness(matrix(1,1,varNum),estimate_F))
#
#
# #========= MBoost ==============
# print('-------- MBoost ---------')
# mboost_formula <- "Y[,1]~"
# for (j in 2:p-1){
#     mboost_formula <- paste(mboost_formula, sprintf("X[,%d]", j), "+")
# }
# mboost_formula <- paste(mboost_formula, sprintf("X[,%d]", p))
# start=proc.time()
# mboost_model=gamboost(as.formula(mboost_formula), control = boost_control(mstop = 50))
# print(proc.time()-start)
# # Predicted Y
# mboost_predicted_Y = predict(mboost_model)
# mse_val = yError(Y, mboost_predicted_Y)
# smooth_val = ySmoothness(mboost_predicted_Y)
# print('Prediction error:')
# print(mse_val)
# print('Smoothness:')
# print(smooth_val)
#
#
# #========= SCAR ==============
# print('-------- SCAR ---------')
# scar_shape_formula = c()
# for (j in 1:p){
#     scar_shape_formula = c(scar_shape_formula, "cvxin")
# }
# start=proc.time()
# scar_model= scar(X, Y[,1], shape = scar_shape_formula, family = gaussian())
# print(proc.time()-start)
# # Predicted Y
# scar_predicted_Y = predict(scar_model)
# mse_val = yError(Y, scar_predicted_Y)
# smooth_val = ySmoothness(scar_predicted_Y)
# print('Prediction error:')
# print(mse_val)
# print('Smoothness:')
# print(smooth_val)
#
#
# #========= Linear Model ==============
# print('-------- Linear Fit ---------')
# start=proc.time()
# linear_model=lm(as.formula(mboost_formula))
# print(proc.time()-start)
# # Predicted Y
# linear_predicted_Y = predict(linear_model)
# mse_val = yError(Y, linear_predicted_Y)
# smooth_val = ySmoothness(linear_predicted_Y)
# print('Prediction error:')
# print(mse_val)
# print('Smoothness:')
# print(smooth_val)