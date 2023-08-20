"
Description:
    Compare SSAIM with other related methods.

Programmer:
    Jiaqi Zhang

Update history:
    --2019/04/07: Create test codes and do experiments.
"
source('SSAIM.r')
source('OtherMethods.r')
source('evaluate.r')
# source("mslasso2.r")

library(mboost)
library(scar)


# Initialization
#path='D:/programming/Python/SSAIM/'
path=''
count=5
epsilon=1
lambda=0.1
delta=0.1
pList=c(200,400,600,800,1000)
pList=c(200)

"
    =========================================
          Varying variable number (n=p)
    =========================================
"
p_record = c()
n_record = c()
name_record = c()
time_record = c()
MSE_record = c()
smoothness_record = c()
theta_error_recrd = c()

for (p in pList){
    print(sprintf('========= %d Variables when n=p ============',p))
    # Reading data
    n=p
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn.csv',path,p,n),header=FALSE))
    X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn.csv',path,p,n),header=FALSE))
    Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn.csv',path,p,n),header=FALSE))
    print('Finish reading data!')

    #========= MBoost ==============
    print('-------- MBoost ---------')
    mboost_formula <- "Y[,1]~"
    for (j in 2:p-1){
        mboost_formula <- paste(mboost_formula, sprintf("X[,%d]", j), "+")
    }
    mboost_formula <- paste(mboost_formula, sprintf("X[,%d]", p))
    totalMSe <- 0.0
    totalSmoothness <- 0.0
    totalTime=proc.time()-proc.time()
    for(i in 1:count){
        start=proc.time()
        mboost_model=gamboost(as.formula(mboost_formula))
        totalTime=totalTime+(proc.time()-start)
        # Predicted Y
        mboost_predicted_Y = predict(mboost_model)
        mse_val = yError(Y, mboost_predicted_Y)
        totalMSe = totalMSe + mse_val
        smooth_val = ySmoothness(mboost_predicted_Y)
        totalSmoothness = totalSmoothness + smooth_val
        # reord
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'MBoost')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, mse_val)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, -1.0)
    }
    print("Avg Time Cost:")
    print(totalTime/count)
    print('Avg Prediction error:')
    print(totalMSe/count)
    print('Avg Smoothness:')
    print(totalSmoothness/count)
    # Save estimation result
    if (!('./estimation' %in% list.dirs())){
        dir.create('./estimation', showWarnings = FALSE)
    }
    write.csv(mboost_predicted_Y, sprintf('./estimation/MBoost_%dp_%dn_Y.csv', p, n), row.names = FALSE)
    
    
    #========= SCAR ==============
    print('-------- SCAR ---------')
    scar_shape_formula = c()
    for (j in 1:p){
        scar_shape_formula = c(scar_shape_formula, "cvxin")
    }
    totalMSe <- 0.0
    totalSmoothness <- 0.0
    totalTime=proc.time()-proc.time()
    for(i in 1:count){
        start=proc.time()
        scar_model= scar(X, Y[,1], shape = scar_shape_formula, family = gaussian())
        totalTime=totalTime+(proc.time()-start)
        # Predicted Y
        scar_predicted_Y = predict(scar_model)
        mse_val = yError(Y, scar_predicted_Y)
        totalMSe = totalMSe + mse_val
        smooth_val = ySmoothness(scar_predicted_Y)
        totalSmoothness = totalSmoothness + smooth_val
        # reord
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'SCAR')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, mse_val)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, -1.0)
    }
    print("Avg Time Cost:")
    print(totalTime/count)
    print('Avg Prediction error:')
    print(totalMSe/count)
    print('Avg Smoothness:')
    print(totalSmoothness/count)
    # Save estimation result
    if (!('./estimation' %in% list.dirs())){
        dir.create('./estimation', showWarnings = FALSE)
    }
    write.csv(scar_predicted_Y, sprintf('./estimation/SCAR_%dp_%dn_Y.csv', p, n), row.names = FALSE)
    
    #========= Linear Model ==============
    print('-------- Linear Fit ---------')
    totalMSe <- 0.0
    totalSmoothness <- 0.0
    totalTime=proc.time()-proc.time()
    for(i in 1:count){
        start=proc.time()
        linear_model=lm(as.formula(mboost_formula))
        totalTime=totalTime+(proc.time()-start)
        # Predicted Y
        linear_predicted_Y = predict(linear_model)
        mse_val = yError(Y, linear_predicted_Y)
        totalMSe = totalMSe + mse_val
        smooth_val = ySmoothness(linear_predicted_Y)
        totalSmoothness = totalSmoothness + smooth_val
        # reord
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'linear')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, mse_val)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, -1.0)
    }
    print("Avg Time Cost:")
    print(totalTime/count)
    print('Avg Prediction error:')
    print(totalMSe/count)
    print('Avg Smoothness:')
    print(totalSmoothness/count)
    # Save estimation result
    if (!('./estimation' %in% list.dirs())){
        dir.create('./estimation', showWarnings = FALSE)
    }
    write.csv(linear_predicted_Y, sprintf('./estimation/Linear_%dp_%dn_Y.csv', p, n), row.names = FALSE)
}
# Save to csv file
records = data.frame(p=p_record, n = n_record, name=name_record, time=time_record, mse = MSE_record, smoothness = smoothness_record, thta_error = theta_error_recrd)
write.csv(records, file = "low_dimensional_records_extra.csv", append = FALSE, quote = TRUE, sep = " ")
print("Finished saving low-dimensaionl data.")
print("====================================\n")


"
    =========================================
          Varying sample number (p=100, s=10)
    =========================================
"
path='D:/programming/Python/SSAIM/'
delta=0.4
lambda=0.1
nList=c(10,20,30,40,50,60,70,80,90,100)
nList=c(10)
p=100
p_record = c()
n_record = c()
name_record = c()
time_record = c()
MSE_record = c()
smoothness_record = c()
theta_error_recrd = c()

for (n in nList){
    print(sprintf('========= %d Samples when p=100 s=10 ============',n))
    # Reading data
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn_10s.csv',path,p,n),header=FALSE))
    X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn10s.csv',path,p,n),header=FALSE))
    Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn10s.csv',path,p,n),header=FALSE))
    print('Finish reading data!')
    
        #========= MBoost ==============
    print('-------- MBoost ---------')
    mboost_formula <- "Y[,1]~"
    for (j in 2:p-1){
        mboost_formula <- paste(mboost_formula, sprintf("X[,%d]", j), "+")
    }
    mboost_formula <- paste(mboost_formula, sprintf("X[,%d]", p))
    totalMSe <- 0.0
    totalSmoothness <- 0.0
    totalTime=proc.time()-proc.time()
    for(i in 1:count){
        start=proc.time()
        mboost_model=gamboost(as.formula(mboost_formula))
        totalTime=totalTime+(proc.time()-start)
        # Predicted Y
        mboost_predicted_Y = predict(mboost_model)
        mse_val = yError(Y, mboost_predicted_Y)
        totalMSe = totalMSe + mse_val
        smooth_val = ySmoothness(mboost_predicted_Y)
        totalSmoothness = totalSmoothness + smooth_val
        # reord
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'MBoost')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, mse_val)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, -1.0)
    }
    print("Avg Time Cost:")
    print(totalTime/count)
    print('Avg Prediction error:')
    print(totalMSe/count)
    print('Avg Smoothness:')
    print(totalSmoothness/count)
    # Save estimation result
    if (!('./estimation' %in% list.dirs())){
        dir.create('./estimation', showWarnings = FALSE)
    }
    write.csv(mboost_predicted_Y, sprintf('./estimation/MBost_%dp_%dn_Y.csv', p, n), row.names = FALSE)
    
    #========= Linear Model ==============
    print('-------- Linear Fit ---------')
    totalMSe <- 0.0
    totalSmoothness <- 0.0
    totalTime=proc.time()-proc.time()
    for(i in 1:count){
        start=proc.time()
        linear_model=lm(as.formula(mboost_formula))
        totalTime=totalTime+(proc.time()-start)
        # Predicted Y
        linear_predicted_Y = predict(linear_model)
        mse_val = yError(Y, linear_predicted_Y)
        totalMSe = totalMSe + mse_val
        smooth_val = ySmoothness(linear_predicted_Y)
        totalSmoothness = totalSmoothness + smooth_val
        # reord
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'linear')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, mse_val)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, -1.0)
    }
    print("Avg Time Cost:")
    print(totalTime/count)
    print('Avg Prediction error:')
    print(totalMSe/count)
    print('Avg Smoothness:')
    print(totalSmoothness/count)
    # Save estimation result
    if (!('./estimation' %in% list.dirs())){
        dir.create('./estimation', showWarnings = FALSE)
    }
    write.csv(linear_predicted_Y, sprintf('./estimation/Linear_%dp_%dn_Y.csv', p, n), row.names = FALSE)
}
# Save to csv file
records = data.frame(p=p_record, n = n_record, name=name_record, time=time_record, mse = MSE_record, smoothness = smoothness_record, thta_error = theta_error_recrd)
write.csv(records, file = "high_dimensional_records_extra.csv", append = FALSE, quote = TRUE, sep = " ")
print("Finished saving high-dimensaionl data.")
print("====================================\n")