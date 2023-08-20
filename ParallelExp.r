"
Description:
    Test SSAIM methods.

Author:
    Jiaqi Zhang
"
source('TimeSSAIM.r')
source('OtherMethods.r')
source('evaluate.r')

library(mboost)
library(scar)

# Initialization
#path='D:/programming/Python/SSAIM/'
path=''
count = 1
pList = c(200, 400, 600, 800, 1000)
# pList = c()
parameter_list = list(c(1, 0.06, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1))
lambda_list = c(0.1, 0.1, 0.05, 0.1, 0.1)
# 200: c(1, 0.06, 10), 
# pList = c()
# parameter_list = list(c(1, 0.06, 10))

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
if(length(pList) > 0){
for (p_index in 1:length(pList)){
    p = pList[p_index]
    print(sprintf('========= %d Variables when n=p ============',p))
    # Reading data
    n=p
    # s = as.integer(p*0.5)
    s = 10
    pars = parameter_list[[p_index]]
    epsilon = pars[1]
    lambda = pars[2]
    delta = pars[3]
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    F=as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    print('Finish reading data!')

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
        # estimation
        estimate_F=ssaim$F
        estimatedTheta=ssaim$theta
        theta_error = thetaError(t(theta),estimatedTheta)
        total_theta_error = total_theta_error + theta_error
        pred_error = predictError(Y,estimatedTheta,estimate_F)
        totalMSe = totalMSe + pred_error
        smooth_val = Fsmoothness(estimatedTheta,estimate_F)
        totalSmoothness = totalSmoothness + smooth_val
        F_error = Ferror(F, estimate_F)
        total_F_error = total_F_error + F_error
        # record
         # reord
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'SSAIM')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, pred_error)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, theta_error)
        F_error_record = c(F_error_record, F_error)
    }
    print(totalTime/count)
    # Evaluation
    print('Theta error:')
    print(total_theta_error/count)
    print('Prediction error:')
    print(totalMSe / count)
    print('Smoothness:')
    print(totalSmoothness/count)
    print('F error:')
    print(total_F_error/count)
    
    #========= SLIM ==============
    print('-------- SLIM ---------')
    lambda = lambda_list[p_index]
    totalTime=proc.time()-proc.time()
    totalMSe = 0.0
    totalSmoothness = 0.0
    total_theta_error = 0.0
    total_F_error = 0.0
    for(i in 1:count){
        start=proc.time()
        slim=ParSLIM(X,Y,lambda)
        totalTime=totalTime+(proc.time()-start)
        # estimation
        estimate_F=slim$F
        estimatedTheta=slim$theta
        estimatedTheta = t(estimatedTheta) #TODO: check this
        theta_error = thetaError(theta,t(estimatedTheta))
        total_theta_error = total_theta_error + theta_error
        pred_error = predictError(Y,estimatedTheta,estimate_F)
        totalMSe = totalMSe + pred_error
        smooth_val = Fsmoothness(estimatedTheta,estimate_F)
        totalSmoothness = totalSmoothness + smooth_val
        F_error = Ferror(F, estimate_F)
        total_F_error = total_F_error + F_error
        # record
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'SLIM')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, pred_error)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, theta_error)
        F_error_record = c(F_error_record, F_error)
    }
    print(totalTime/count)
    # Evaluation
    print('Theta error:')
    print(total_theta_error/count)
    print('Prediction error:')
    print(totalMSe / count)
    print('Smoothness:')
    print(totalSmoothness/count)
    print('F error:')
    print(total_F_error/count)

    # # Save estimation result
    # if (!('./estimation' %in% list.dirs())){
    #     dir.create('./estimation', showWarnings = FALSE)
    # }
    # write.csv(t(estimatedTheta), sprintf('./estimation/SSAIM_%dp_%dn_%ds_theta.csv', p, n,s), row.names = FALSE)
    # write.csv(estimate_F, sprintf('./estimation/SSAIM_%dp_%dn_%ds_F.csv', p, n,s), row.names = FALSE)
    # write.csv(as.matrix(estimate_F %*% t(estimatedTheta)), sprintf('./estimation/SSAIM_%dp_%dn_%ds_Y.csv', p, n,s), row.names = FALSE)
}
# Save to csv file
records = data.frame(p=p_record, n = n_record, name=name_record, time=time_record, mse = MSE_record, smoothness = smoothness_record, thta_error = theta_error_recrd, F_error = F_error_record)
# write.csv(records, file = "low_dimensional_records_SSAIM.csv", append = FALSE, quote = TRUE, sep = " ")
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
# nList = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200)
# parameter_list = list(c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1),
# c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1),
# c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1), c(0.1, 0.1, 1))
# parameter_list = list(c(0.05, 0.1, 10), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1), c(0.5, 0.11, 1))

# nList = c(10)
# parameter_list = list(c(0.05, 0.13, 10))

nList = c(30, 50, 70, 100, 150, 200)
parameter_list = list(c(0.5, 0.3, 1), c(0.5, 0.1, 1), c(0.5, 0.08, 1), c(1, 0.06, 1), c(1, 0.06, 1), c(1, 0.06, 1))
nList = c()
if(length(nList) > 0){
for (n_index in 1:length(nList)){
    n = nList[n_index]
    print(sprintf('========= %d Samples when p=100 s=10 ============',n))
    pars = parameter_list[[n_index]]
    epsilon = pars[1]
    lambda = pars[2]
    delta = pars[3]
    # Reading data
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    F=as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    print('Finish reading data!')

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
        # estimation
        estimate_F=ssaim$F
        estimatedTheta=ssaim$theta
        theta_error = thetaError(t(theta),estimatedTheta)
        total_theta_error = total_theta_error + theta_error
        pred_error = predictError(Y,estimatedTheta,estimate_F)
        totalMSe = totalMSe + pred_error
        smooth_val = Fsmoothness(estimatedTheta,estimate_F)
        totalSmoothness = totalSmoothness + smooth_val
        F_error = Ferror(F, estimate_F)
        total_F_error = total_F_error + F_error
        # record
         # reord
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'SSAIM')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, pred_error)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, theta_error)
        F_error_record = c(F_error_record, F_error)
    }
    print(totalTime/count)
    # Evaluation
    print('Theta error:')
    print(total_theta_error/count)
    print('Prediction error:')
    print(totalMSe / count)
    print('Smoothness:')
    print(totalSmoothness/count)
    print('F error:')
    print(total_F_error/count)
    
    #========= SLIM ==============
    print('-------- SLIM ---------')
    totalTime=proc.time()-proc.time()
    totalMSe = 0.0
    totalSmoothness = 0.0
    total_theta_error = 0.0
    total_F_error = 0.0
    for(i in 1:count){
        start=proc.time()
        slim=ParSLIM(X,Y,lambda)
        totalTime=totalTime+(proc.time()-start)
        # estimation
        estimate_F=slim$F
        estimatedTheta=slim$theta
        estimatedTheta = t(estimatedTheta) #TODO: check this
        theta_error = thetaError(theta,t(estimatedTheta))
        total_theta_error = total_theta_error + theta_error
        pred_error = predictError(Y,estimatedTheta,estimate_F)
        totalMSe = totalMSe + pred_error
        smooth_val = Fsmoothness(estimatedTheta,estimate_F)
        totalSmoothness = totalSmoothness + smooth_val
        F_error = Ferror(F, estimate_F)
        total_F_error = total_F_error + F_error
        # record
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'SLIM')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, pred_error)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, theta_error)
        F_error_record = c(F_error_record, F_error)
    }
    print(totalTime/count)
    # Evaluation
    print('Theta error:')
    print(total_theta_error/count)
    print('Prediction error:')
    print(totalMSe / count)
    print('Smoothness:')
    print(totalSmoothness/count)
    print('F error:')
    print(total_F_error/count)

    # # Save estimation result
    # if (!('./estimation' %in% list.dirs())){
    #     dir.create('./estimation', showWarnings = FALSE)
    # }
    # write.csv(t(estimatedTheta), sprintf('./estimation/SSAIM_%dp_%dn_%ds_theta.csv', p, n,s), row.names = FALSE)
    # write.csv(estimate_F, sprintf('./estimation/SSAIM_%dp_%dn_%ds_F.csv', p, n,s), row.names = FALSE)
    # write.csv(as.matrix(estimate_F %*% t(estimatedTheta)), sprintf('./estimation/SSAIM_%dp_%dn_%ds_Y.csv', p, n,s), row.names = FALSE)
}
# Save to csv file
records = data.frame(p=p_record, n = n_record, name=name_record, time=time_record, mse = MSE_record, smoothness = smoothness_record, thta_error = theta_error_recrd, F_error = F_error_record)
# write.csv(records, file = "high_dimensional_records_SSAIM.csv", append = FALSE, quote = TRUE, sep = " ")
print("Finished saving high-dimensaionl data.")
print("====================================\n")
}


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
parameter_list = list()
if(length(sList) > 0){
for (s_index in 1:length(sList)){
    s = sList[s_index]
    print(sprintf('========= %d Sparsity when p=200 n=200 ============',s))
    pars = parameter_list[[s_index]]
    epsilon = pars[1]
    lambda = pars[2]
    delta = pars[3]
    # Reading data
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    F=as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
    print('Finish reading data!')

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
        # estimation
        estimate_F=ssaim$F
        estimatedTheta=ssaim$theta
        theta_error = thetaError(t(theta),estimatedTheta)
        total_theta_error = total_theta_error + theta_error
        pred_error = predictError(Y,estimatedTheta,estimate_F)
        totalMSe = totalMSe + pred_error
        smooth_val = Fsmoothness(estimatedTheta,estimate_F)
        totalSmoothness = totalSmoothness + smooth_val
        F_error = Ferror(F, estimate_F)
        total_F_error = total_F_error + F_error
        # record
         # reord
        p_record = c(p_record, p)
        n_record = c(n_record, n)
        name_record = c(name_record, 'SSAIM')
        time_record = c(time_record, (proc.time()-start)[3])
        MSE_record = c(MSE_record, pred_error)
        smoothness_record = c(smoothness_record, smooth_val)
        theta_error_recrd = c(theta_error_recrd, theta_error)
        F_error_record = c(F_error_record, F_error)
    }
    print(totalTime/count)
    # Evaluation
    print('Theta error:')
    print(total_theta_error/count)
    print('Prediction error:')
    print(totalMSe / count)
    print('Smoothness:')
    print(totalSmoothness/count)
    print('F error:')
    print(total_F_error/count)
    
    # Save estimation result
    if (!('./estimation' %in% list.dirs())){
        dir.create('./estimation', showWarnings = FALSE)
    }
    write.csv(t(estimatedTheta), sprintf('./estimation/SSAIM_%dp_%dn_%ds_theta.csv', p, n,s), row.names = FALSE)
    write.csv(estimate_F, sprintf('./estimation/SSAIM_%dp_%dn_%ds_F.csv', p, n,s), row.names = FALSE)
    write.csv(as.matrix(estimate_F %*% t(estimatedTheta)), sprintf('./estimation/SSAIM_%dp_%dn_%ds_Y.csv', p, n,s), row.names = FALSE)
}
# Save to csv file
records = data.frame(p=p_record, n = n_record, name=name_record, time=time_record, mse = MSE_record, smoothness = smoothness_record, thta_error = theta_error_recrd, F_error = F_error_record)
write.csv(records, file = "vary_sparsity_records_SSAIM.csv", append = FALSE, quote = TRUE, sep = " ")
print("Finished saving sparsity level data.")
print("====================================\n")
}