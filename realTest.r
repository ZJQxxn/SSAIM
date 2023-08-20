"
Description:
    Utilize SSAIM on World Bank dataset.

Programmer:
    Jiaqi Zhang

Update history:
    --2019/04/11: Create codes and finish the experiment.
"
source('SSAIM.r')
source('OtherMethods.r')
source('evaluate.r')


print('Begin reading data....')
theta=as.matrix(read.csv(sprintf('dataset/highdim-env_extract_theta.csv'),header=FALSE))
X=as.matrix(read.csv(sprintf('dataset/highdim-env_extract_X.csv'),header=FALSE))
Y=as.matrix(read.csv(sprintf('dataset/env_extract_Y.csv'),header=FALSE))
for (i in 1:length(X[1,])){
    X[,i]=X[,i]/norm(X[,i], type='2')
}
Y = Y/norm(Y, type='2')
print('Finish reading data!')
print('================ World Bank Data ==================')
varNum=length(X[1,])
sampNum=length(X[,1])

n = sampNum
p = varNum


count=1
epsilon=0.1
lambda=0.1
delta=1

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
print('Theta error:')
print(thetaError(t(theta),estimatedTheta))
print('Prediction error:')
print(predictError(Y,estimatedTheta,estimate_F))
print('Smoothness:')
print(Fsmoothness(estimatedTheta,estimate_F))

# #========= SLIM ==============
# print('-------- SLIM ---------')
# totalTime=proc.time()-proc.time()
# for(i in 1:count){
#     start=proc.time()
#     slim=SLIM(X,Y,lambda)
#     totalTime=totalTime+(proc.time()-start)
# }
# print(totalTime/count)
# # Evaluation
# estimate_F=slim$F
# estimatedTheta=slim$theta
# print('Theta error:')
# print(thetaError(theta,estimatedTheta))
# print('Prediction error:')
# print(predictError(Y,t(estimatedTheta),estimate_F))
# print('Smoothness:')
# print(Fsmoothness(estimatedTheta,estimate_F))

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

print('----------- CPAV ----------')
print(date())
start=proc.time()
model=CyclicSPAV(X,Y,delta)
print(proc.time()-start)
estimate_F=model$F
print('Prediction error:')
print(predictError(Y,matrix(1,1,varNum),estimate_F))
print('Smoothness:')
print(Fsmoothness(matrix(1,1,varNum),estimate_F))

#========= MBoost ==============
print('-------- MBoost ---------')
mboost_formula <- "Y[,1]~"
for (j in 2:p-1){
    mboost_formula <- paste(mboost_formula, sprintf("X[,%d]", j), "+")
}
mboost_formula <- paste(mboost_formula, sprintf("X[,%d]", p))
start=proc.time()
mboost_model=gamboost(as.formula(mboost_formula), control = boost_control(mstop = 50))
print(proc.time()-start)
# Predicted Y
mboost_predicted_Y = predict(mboost_model)
mse_val = yError(Y, mboost_predicted_Y)
smooth_val = ySmoothness(mboost_predicted_Y)
print('Prediction error:')
print(mse_val)
print('Smoothness:')
print(smooth_val)
    
    
#========= SCAR ==============
print('-------- SCAR ---------')
scar_shape_formula = c()
for (j in 1:p){
    scar_shape_formula = c(scar_shape_formula, "cvxin")
}
start=proc.time()
scar_model= scar(X, Y[,1], shape = scar_shape_formula, family = gaussian())
print(proc.time()-start)
# Predicted Y
scar_predicted_Y = predict(scar_model)
mse_val = yError(Y, scar_predicted_Y)
smooth_val = ySmoothness(scar_predicted_Y)
print('Prediction error:')
print(mse_val)
print('Smoothness:')
print(smooth_val)


#========= Linear Model ==============
print('-------- Linear Fit ---------')
start=proc.time()
linear_model=lm(as.formula(mboost_formula))
print(proc.time()-start)
# Predicted Y
linear_predicted_Y = predict(linear_model)
mse_val = yError(Y, linear_predicted_Y)
smooth_val = ySmoothness(linear_predicted_Y)
print('Prediction error:')
print(mse_val)
print('Smoothness:')
print(smooth_val)