"
Description:
    Other baseline methods.

Author:
    Jiaqi Zhang
"
source('SSAIM.r')
source('OtherMethods.r')
source('evaluate.r')

library(mboost)
library(scar)

# Initialization
#path='D:/programming/Python/SSAIM/'
path=''
count = 10
pList = c(200, 400, 600, 800, 1000)
pList = c()

p_record = c()
n_record = c()
name_record = c()
time_record = c()
MSE_record = c()
smoothness_record = c()
theta_error_recrd = c()
F_error_record = c()

"
    =========================================
          Varying variable number (n=p)
    =========================================
"
for (p in pList){
    print(sprintf('========= %d Variables when n=p ============',p))
    # Reading data
    n=p
    # s = as.integer(p*0.5)
    s = 10
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    F=as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
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
        mboost_model=gamboost(as.formula(mboost_formula), control = boost_control(mstop = 50))
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
    write.csv(mboost_predicted_Y, sprintf('./estimation/MBoost_%dp_%dn_%ds_Y.csv', p, n, s), row.names = FALSE)
    
    
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
    write.csv(scar_predicted_Y, sprintf('./estimation/SCAR_%dp_%dn_%ds_Y.csv', p, n, s), row.names = FALSE)
    
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
    write.csv(linear_predicted_Y, sprintf('./estimation/Linear_%dp_%dn_%ds_Y.csv', p, n, s), row.names = FALSE)


    #========= Isotonic Lasso ==============
    print('-------- Isotonic Lasso ---------')
    totalTime=proc.time()-proc.time()
    totalMSe = 0.0
    totalSmoothness = 0.0
    total_theta_error = 0.0
    interupt <- 0
    for(i in 1:count){
        start=proc.time()
        isolasso=IsoLasso(X,Y,0.1)
        if(length(isolasso)==1){
            iso_lasso_count <- i
            interupt <- 1
            break
        }
        totalTime=totalTime+(proc.time()-start)
    }
    if (0 == interupt){
        iso_lasso_count <- count
    }
    if(length(isolasso)>1){
        print(totalTime/iso_lasso_count)
        # estimation
        estimateY=isolasso$estimatedY
        estimatedTheta=isolasso$theta
        theta_error = thetaError(t(theta),t(estimatedTheta))
        total_theta_error = total_theta_error + theta_error
        pred_error = yError(Y,t(estimateY))
        totalMSe = totalMSe + pred_error
        smooth_val = ySmoothness(estimateY)
        totalSmoothness = totalSmoothness + smooth_val
        # record
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'IsoLasso')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, pred_error)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, theta_error)
        F_error_record = c(F_error_record, -1)
        # print
        print('Theta error:')
        print(total_theta_error/count)
        print('Prediction error:')
        print(totalMSe / count)
        print('Smoothness:')
        print(totalSmoothness/count)

        if (!('./estimation' %in% list.dirs())){
            dir.create('./estimation', showWarnings = FALSE)
        }
        write.csv(estimatedTheta, sprintf('./estimation/IsoLasso_%dp_%dn_%ds_theta.csv', p, n, s), row.names = FALSE)
        write.csv(estimateY, sprintf('./estimation/IsoLasso_%dp_%dn_%ds_Y.csv', p, n, s), row.names = FALSE)
    }
}
# Save to csv file
if (length(p_record) > 0){
    records = data.frame(p=p_record, n = n_record, name=name_record, time=time_record, mse = MSE_record, smoothness = smoothness_record, thta_error = theta_error_recrd, F_error = F_error_record)
    write.csv(records, file = "low_dimensional_records_other.csv", append = FALSE, quote = TRUE, sep = " ")
    print("Finished saving low-dimensaionl data.")
    print("====================================\n")
}





"
    =========================================
          Varying sample number (p=200, s=100)
    =========================================
"
path=''
p = 200
# s = 100
s = 10
nList = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200)
for (n in nList){
    print(sprintf('========= %d Samples when p=100 s=10 ============',n))
    # Reading data
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    F=as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
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
        mboost_model=gamboost(as.formula(mboost_formula), control = boost_control(mstop = 50))
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
    write.csv(mboost_predicted_Y, sprintf('./estimation/MBoost_%dp_%dn_%ds_Y.csv', p, n, s), row.names = FALSE)
    
    
    # #========= SCAR ==============
    # print('-------- SCAR ---------')
    # scar_shape_formula = c()
    # for (j in 1:p){
    #     scar_shape_formula = c(scar_shape_formula, "cvxin")
    # }
    # totalMSe <- 0.0
    # totalSmoothness <- 0.0
    # totalTime=proc.time()-proc.time()
    # for(i in 1:count){
    #     start=proc.time()
    #     scar_model= scar(X, Y[,1], shape = scar_shape_formula, family = gaussian())
    #     totalTime=totalTime+(proc.time()-start)
    #     # Predicted Y
    #     scar_predicted_Y = predict(scar_model)
    #     mse_val = yError(Y, scar_predicted_Y)
    #     totalMSe = totalMSe + mse_val
    #     smooth_val = ySmoothness(scar_predicted_Y)
    #     totalSmoothness = totalSmoothness + smooth_val
    #     # reord
    #     p_record = c(p_record, p)
    #     n_record = c(n_record, n)
    #     name_record = c(name_record, 'SCAR')
    #     time_record = c(time_record, (proc.time()-start)[3])
    #     MSE_record = c(MSE_record, mse_val)
    #     smoothness_record = c(smoothness_record, smooth_val)
    #     theta_error_recrd = c(theta_error_recrd, -1.0)
    # }
    # print("Avg Time Cost:")
    # print(totalTime/count)
    # print('Avg Prediction error:')
    # print(totalMSe/count)
    # print('Avg Smoothness:')
    # print(totalSmoothness/count)
    # # Save estimation result
    # if (!('./estimation' %in% list.dirs())){
    #     dir.create('./estimation', showWarnings = FALSE)
    # }
    # write.csv(scar_predicted_Y, sprintf('./estimation/SCAR_%dp_%dn_%ds_Y.csv', p, n, s), row.names = FALSE)
    
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
    write.csv(linear_predicted_Y, sprintf('./estimation/Linear_%dp_%dn_%ds_Y.csv', p, n, s), row.names = FALSE)


    #========= Isotonic Lasso ==============
    print('-------- Isotonic Lasso ---------')
    totalTime=proc.time()-proc.time()
    totalMSe = 0.0
    totalSmoothness = 0.0
    total_theta_error = 0.0
    interupt <- 0
    for(i in 1:count){
        start=proc.time()
        isolasso=IsoLasso(X,Y,0.1)
        if(length(isolasso)==1){
            iso_lasso_count <- i
            interupt <- 1
            break
        }
        totalTime=totalTime+(proc.time()-start)
    }
    if (0 == interupt){
        iso_lasso_count <- count
    }
    if(length(isolasso)>1){
        print(totalTime/iso_lasso_count)
        # estimation
        estimateY=isolasso$estimatedY
        estimatedTheta=isolasso$theta
        theta_error = thetaError(t(theta),t(estimatedTheta))
        total_theta_error = total_theta_error + theta_error
        pred_error = yError(Y,t(estimateY))
        totalMSe = totalMSe + pred_error
        smooth_val = ySmoothness(estimateY)
        totalSmoothness = totalSmoothness + smooth_val
        # record
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'IsoLasso')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, pred_error)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, theta_error)
        F_error_record = c(F_error_record, -1)
        # print
        print('Theta error:')
        print(total_theta_error/count)
        print('Prediction error:')
        print(totalMSe / count)
        print('Smoothness:')
        print(totalSmoothness/count)

        if (!('./estimation' %in% list.dirs())){
            dir.create('./estimation', showWarnings = FALSE)
        }
        write.csv(estimatedTheta, sprintf('./estimation/IsoLasso_%dp_%dn_%ds_theta.csv', p, n, s), row.names = FALSE)
        write.csv(estimateY, sprintf('./estimation/IsoLasso_%dp_%dn_%ds_Y.csv', p, n, s), row.names = FALSE)
    } 
}
# Save to csv file
records = data.frame(p=p_record, n = n_record, name=name_record, time=time_record, mse = MSE_record, smoothness = smoothness_record, thta_error = theta_error_recrd, F_error = F_error_record)
write.csv(records, file = "high_dimensional_records_other.csv", append = FALSE, quote = TRUE, sep = " ")
print("Finished saving high-dimensaionl data.")
print("====================================\n")


"
    =========================================
          Varying sparsity level (p=200, n=200)
    =========================================
"
path=''
p = 200
n = 200
sList = c(20, 40, 60, 80, 100, 120, 140, 160, 180, 200)
sList = c()
for (s in sList){
    print(sprintf('========= %d Sparsity when p=200 n=200 ============',s))
    # Reading data
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    F=as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
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
        mboost_model=gamboost(as.formula(mboost_formula), control = boost_control(mstop = 50))
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
    write.csv(mboost_predicted_Y, sprintf('./estimation/MBoost_%dp_%dn_%ds_Y.csv', p, n, s), row.names = FALSE)
    
    
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
    write.csv(scar_predicted_Y, sprintf('./estimation/SCAR_%dp_%dn_%ds_Y.csv', p, n, s), row.names = FALSE)
    
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
    write.csv(linear_predicted_Y, sprintf('./estimation/Linear_%dp_%dn_%ds_Y.csv', p, n, s), row.names = FALSE)


    #========= Isotonic Lasso ==============
    print('-------- Isotonic Lasso ---------')
    totalTime=proc.time()-proc.time()
    totalMSe = 0.0
    totalSmoothness = 0.0
    total_theta_error = 0.0
    interupt <- 0
    for(i in 1:count){
        start=proc.time()
        isolasso=IsoLasso(X,Y,0.1)
        if(length(isolasso)==1){
            iso_lasso_count <- i
            interupt <- 1
            break
        }
        totalTime=totalTime+(proc.time()-start)
    }
    if (0 == interupt){
        iso_lasso_count <- count
    }
    if(length(isolasso)>1){
        print(totalTime/iso_lasso_count)
        # estimation
        estimateY=isolasso$estimatedY
        estimatedTheta=isolasso$theta
        theta_error = thetaError(t(theta),t(estimatedTheta))
        total_theta_error = total_theta_error + theta_error
        pred_error = yError(Y,t(estimateY))
        totalMSe = totalMSe + pred_error
        smooth_val = ySmoothness(estimateY)
        totalSmoothness = totalSmoothness + smooth_val
        # record
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'IsoLasso')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, pred_error)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, theta_error)
        F_error_record = c(F_error_record, -1)
        # print
        print('Theta error:')
        print(total_theta_error/count)
        print('Prediction error:')
        print(totalMSe / count)
        print('Smoothness:')
        print(totalSmoothness/count)

        if (!('./estimation' %in% list.dirs())){
            dir.create('./estimation', showWarnings = FALSE)
        }
        write.csv(estimatedTheta, sprintf('./estimation/IsoLasso_%dp_%dn_%ds_theta.csv', p, n, s), row.names = FALSE)
        write.csv(estimateY, sprintf('./estimation/IsoLasso_%dp_%dn_%ds_Y.csv', p, n, s), row.names = FALSE)
    } 
}
# Save to csv file
records = data.frame(p=p_record, n = n_record, name=name_record, time=time_record, mse = MSE_record, smoothness = smoothness_record, thta_error = theta_error_recrd, F_error = F_error_record)
write.csv(records, file = "vary_sparsity_records_other.csv", append = FALSE, quote = TRUE, sep = " ")
print("Finished saving sparsity level data.")
print("====================================\n")