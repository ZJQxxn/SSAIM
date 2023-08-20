"
Description:
    Compare SSAIM with other related methods.

Programmer:
    Jiaqi Zhang
"
source('SSAIM.r')
source('OtherMethods.r')
source('evaluate.r')


# Initialization
#path='D:/programming/Python/SSAIM/'
path=''
count=1
epsilon=1
lambda=0.05
delta=0.1
# pList=c(200,400,600,800,1000)
pList=c(200)

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
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn.csv',path,p,n),header=FALSE))
    X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn.csv',path,p,n),header=FALSE))
    Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn.csv',path,p,n),header=FALSE))
    F=as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn.csv',path,p,n),header=FALSE))
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
        ssaim=SSAIM(X,Y,epsilon,lambda,delta,iter=20)
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
    write.csv(t(estimatedTheta), sprintf('./estimation/SSAIM_%dp_%dn_theta.csv', p, n), row.names = FALSE)
    write.csv(estimate_F, sprintf('./estimation/SSAIM_%dp_%dn_F.csv', p, n), row.names = FALSE)
    write.csv(estimated_Y, sprintf('./estimation/SSAIM_%dp_%dn_Y.csv', p, n), row.names = FALSE)

    #========= SLIM ==============
    print('-------- SLIM ---------')
    totalTime=proc.time()-proc.time()
    totalMSe = 0.0
    totalSmoothness = 0.0
    total_theta_error = 0.0
    total_F_error = 0.0
    for(i in 1:count){
        start=proc.time()
        slim=SLIM(X,Y,lambda,iter=20)
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
    # Save estimation result
    if (!('./estimation' %in% list.dirs())){
        dir.create('./estimation', showWarnings = FALSE)
    }
    write.csv(estimatedTheta, sprintf('./estimation/SLIM_%dp_%dn_theta.csv', p, n), row.names = FALSE)
    write.csv(estimate_F, sprintf('./estimation/SLIM_%dp_%dn_F.csv', p, n), row.names = FALSE)
    write.csv(estimated_Y, sprintf('./estimation/SLIM_%dp_%dn_Y.csv', p, n), row.names = FALSE)


    #========= Isotonic Lasso ==============
    print('-------- Isotonic Lasso ---------')
    totalTime=proc.time()-proc.time()
    totalMSe = 0.0
    totalSmoothness = 0.0
    total_theta_error = 0.0
    interupt <- 0
    for(i in 1:count){
        start=proc.time()
        isolasso=IsoLasso(X,Y,lambda)
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
        write.csv(estimatedTheta, sprintf('./estimation/IsoLasso_%dp_%dn_theta.csv', p, n), row.names = FALSE)
        write.csv(estimateY, sprintf('./estimation/IsoLasso_%dp_%dn_Y.csv', p, n), row.names = FALSE)
    }
    

    #========= CPAV ==============
    print('-------- CPAV ---------')
    totalTime=proc.time()-proc.time()
    totalMSe = 0.0
    totalSmoothness = 0.0
    total_F_error = 0.0
    for(i in 1:count){
        print(date())
        start=proc.time()
        model=CyclicSPAV(X,Y,delta)
        totalTime=totalTime+(proc.time()-start)
        # estimation
        estimate_F=model$F
        pred_error = predictError(Y,matrix(1,1,p),estimate_F)
        totalMSe = totalMSe + pred_error
        smooth_val = Fsmoothness(matrix(1,1,p),estimate_F)
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
        theta_error_recrd = c(theta_error_recrd, -1)
        F_error_record = c(F_error_record, F_error)
    }
    print(totalTime/count)
    # Evaluation
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
    write.csv(matrix(1,1,p), sprintf('./estimation/CSPAV_%dp_%dn_theta.csv', p, n), row.names = FALSE)
    write.csv(estimate_F, sprintf('./estimation/CSPAV_%dp_%dn_F.csv', p, n), row.names = FALSE)
    write.csv(as.matrix(estimate_F %*% t(matrix(1,1,p)), sprintf('./estimation/CSPAV_%dp_%dn_Y.csv', p, n), row.names = FALSE)
    
}
# Save to csv file
records = data.frame(p=p_record, n = n_record, name=name_record, time=time_record, mse = MSE_record, smoothness = smoothness_record, thta_error = theta_error_recrd, F_error = F_error_record)
write.csv(records, file = "low_dimensional_records.csv", append = FALSE, quote = TRUE, sep = " ")
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
for (n in nList){
    print(sprintf('========= %d Samples when p=100 s=10 ============',n))
    # Reading data
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn_10s.csv',path,p,n),header=FALSE))
    X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn10s.csv',path,p,n),header=FALSE))
    Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn10s.csv',path,p,n),header=FALSE))
    F=as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn10s.csv',path,p,n),header=FALSE))
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
        ssaim=SSAIM(X,Y,epsilon,lambda,delta,iter=20)
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
    write.csv(t(estimatedTheta), sprintf('./estimation/SSAIM_%dp_%dn_theta.csv', p, n), row.names = FALSE)
    write.csv(estimate_F, sprintf('./estimation/SSAIM_%dp_%dn_F.csv', p, n), row.names = FALSE)
    write.csv(estimated_Y, sprintf('./estimation/SSAIM_%dp_%dn_Y.csv', p, n), row.names = FALSE)

    #========= SLIM ==============
    print('-------- SLIM ---------')
    totalTime=proc.time()-proc.time()
    totalMSe = 0.0
    totalSmoothness = 0.0
    total_theta_error = 0.0
    total_F_error = 0.0
    for(i in 1:count){
        start=proc.time()
        slim=SLIM(X,Y,lambda,iter=20)
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
    # Save estimation result
    if (!('./estimation' %in% list.dirs())){
        dir.create('./estimation', showWarnings = FALSE)
    }
    write.csv(estimatedTheta, sprintf('./estimation/SLIM_%dp_%dn_theta.csv', p, n), row.names = FALSE)
    write.csv(estimate_F, sprintf('./estimation/SLIM_%dp_%dn_F.csv', p, n), row.names = FALSE)
    write.csv(estimated_Y, sprintf('./estimation/SLIM_%dp_%dn_Y.csv', p, n), row.names = FALSE)


    # #========= Isotonic Lasso ==============
    # print('-------- Isotonic Lasso ---------')
    # totalTime=proc.time()-proc.time()
    # totalMSe = 0.0
    # totalSmoothness = 0.0
    # total_theta_error = 0.0
    # interupt <- 0
    # for(i in 1:count){
    #     start=proc.time()
    #     isolasso=IsoLasso(X,Y,lambda)
    #     if(length(isolasso)==1){
    #         iso_lasso_count <- i
    #         interupt <- 1
    #         break
    #     }
    #     totalTime=totalTime+(proc.time()-start)
    # }
    # if (0 == interupt){
    #     iso_lasso_count <- count
    # }
    # if(length(isolasso)>1){
    #     print(totalTime/iso_lasso_count)
    #     # estimation
    #     estimateY=isolasso$estimatedY
    #     estimatedTheta=isolasso$theta
    #     theta_error = thetaError(t(theta),t(estimatedTheta))
    #     total_theta_error = total_theta_error + theta_error
    #     pred_error = yError(Y,t(estimateY))
    #     totalMSe = totalMSe + pred_error
    #     smooth_val = ySmoothness(estimateY)
    #     totalSmoothness = totalSmoothness + smooth_val
    #     # record
    #     p_record = c(p_record, p)
    #     n_record = c(n_record, n)
    #     name_record = c(name_record, 'IsoLasso')
    #     time_record = c(time_record, (proc.time()-start)[3])
    #     MSE_record = c(MSE_record, pred_error)
    #     smoothness_record = c(smoothness_record, smooth_val)
    #     theta_error_recrd = c(theta_error_recrd, theta_error)
    #     F_error_record = c(F_error_record, -1)
    #     # print
    #     print('Theta error:')
    #     print(total_theta_error/count)
    #     print('Prediction error:')
    #     print(totalMSe / count)
    #     print('Smoothness:')
    #     print(totalSmoothness/count)

    #     if (!('./estimation' %in% list.dirs())){
    #         dir.create('./estimation', showWarnings = FALSE)
    #     }
    #     write.csv(estimatedTheta, sprintf('./estimation/IsoLasso_%dp_%dn_theta.csv', p, n), row.names = FALSE)
    #     write.csv(estimateY, sprintf('./estimation/IsoLasso_%dp_%dn_Y.csv', p, n), row.names = FALSE)
    # }
    

    #========= CPAV ==============
    print('-------- CPAV ---------')
    totalTime=proc.time()-proc.time()
    totalMSe = 0.0
    totalSmoothness = 0.0
    total_F_error = 0.0
    for(i in 1:count){
        print(date())
        start=proc.time()
        model=CyclicSPAV(X,Y,delta)
        totalTime=totalTime+(proc.time()-start)
        # estimation
        estimate_F=model$F
        pred_error = predictError(Y,matrix(1,1,p),estimate_F)
        totalMSe = totalMSe + pred_error
        smooth_val = Fsmoothness(matrix(1,1,p),estimate_F)
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
        theta_error_recrd = c(theta_error_recrd, -1)
        F_error_record = c(F_error_record, F_error)
    }
    print(totalTime/count)
    # Evaluation
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
    write.csv(matrix(1,1,p), sprintf('./estimation/CSPAV_%dp_%dn_theta.csv', p, n), row.names = FALSE)
    write.csv(estimate_F, sprintf('./estimation/CSPAV_%dp_%dn_F.csv', p, n), row.names = FALSE)
    write.csv(as.matrix(estimate_F %*% t(matrix(1,1,p)), sprintf('./estimation/CSPAV_%dp_%dn_Y.csv', p, n), row.names = FALSE)
    
}
# Save to csv file
records = data.frame(p=p_record, n = n_record, name=name_record, time=time_record, mse = MSE_record, smoothness = smoothness_record, thta_error = theta_error_recrd, F_error = F_error_record)
write.csv(records, file = "high_dimensional_records.csv", append = FALSE, quote = TRUE, sep = " ")
print("Finished saving high-dimensaionl data.")
print("====================================\n")


"
    =========================================
          Varying lambda (p=100, n=100, s=10)
    =========================================
"
#lambdaList=c(0.05,0.06,0.07,0.08,0.09,0.1)
lambdaList=c()
delta=0.1
p=100
n=100
# Reading data
if(length(lambdaList)>0){
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('data/theta/Theta_%dp_%dn_10s.csv',p,n),header=FALSE))
    X=as.matrix(read.csv(sprintf('data/X/X_%dp_%dn10s.csv',p,n),header=FALSE))
    Y=as.matrix(read.csv(sprintf('data/Y/Y_%dp_%dn10s.csv',p,n),header=FALSE))
    print('Finish reading data!')
}
for (lambda in lambdaList){
    print(sprintf('========= %.2f Lambda when p=100 n=100 s=10 ============',lambda))

    #========= SSAIM ==============
    print('-------- SSAIM ---------')
    totalTime=proc.time()-proc.time()
    for(i in 1:count){
        start=proc.time()
        ssaim=SSAIM(X,Y,epsilon,lambda,delta,iter=20)
        totalTime=totalTime+(proc.time()-start)
    }
    print(totalTime/count)
    # Evaluation
    estimate_F=ssaim$F
    estimatedTheta=ssaim$theta
    rm(ssaim)
    print('Theta error:')
    print(thetaError(t(theta),estimatedTheta))
    print('Prediction error:')
    print(predictError(Y,estimatedTheta,estimate_F))
    print('Smoothness:')
    print(Fsmoothness(estimatedTheta,estimate_F))

    #========= SLIM ==============
    print('-------- SLIM ---------')
    totalTime=proc.time()-proc.time()
    for(i in 1:count){
        start=proc.time()
        slim=SLIM(X,Y,lambda,iter=20)
        totalTime=totalTime+(proc.time()-start)
    }
    print(totalTime/count)
    # Evaluation
    estimate_F=slim$F
    estimatedTheta=slim$theta
    rm(slim)
    print('Theta error:')
    print(thetaError(t(theta),estimatedTheta))
    print('Prediction error:')
    print(predictError(Y,t(estimatedTheta),estimate_F))
    print('Smoothness:')
    print(Fsmoothness(estimatedTheta,estimate_F))
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
        #TODO: No need to compare theta error for Isotonic Lasso
        #print('Theta error:')
        #print(thetaError(theta,estimatedTheta))
        print('Prediction error:')
        print(integralPreError(t(Y),estimateY))
        print('Smoothness:')
        print(ySmoothness(estimateY))
    }
    "
}


"
    =========================================
          Varying delta (p=100, n=100, s=10)
    =========================================
"
#deltaList=c(0.1,0.2,0.3,0.4,0.5)
deltaList=c()
p=100
n=100
epsilon=1
lambda=0.05
# Reading data
if(length(deltaList)>0){
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('data/theta/Theta_%dp_%dn_10s.csv',p,n),header=FALSE))
    X=as.matrix(read.csv(sprintf('data/X/X_%dp_%dn10s.csv',p,n),header=FALSE))
    Y=as.matrix(read.csv(sprintf('data/Y/Y_%dp_%dn10s.csv',p,n),header=FALSE))
    print('Finish reading data!')
}
for (delta in deltaList){
    print(sprintf('========= %.2f Delta when p=100 n=100 s=10 ============',delta))

    #========= SSAIM ==============
    print('-------- SSAIM ---------')
    totalTime=proc.time()-proc.time()
    for(i in 1:count){
        start=proc.time()
        ssaim=SSAIM(X,Y,epsilon,lambda,delta,iter=20)
        totalTime=totalTime+(proc.time()-start)
    }
    print(totalTime/count)
    # Evaluation
    estimate_F=ssaim$F
    estimatedTheta=ssaim$theta
    rm(ssaim)
    print('Theta error:')
    print(thetaError(t(theta),estimatedTheta))
    print('Prediction error:')
    print(predictError(Y,estimatedTheta,estimate_F))
    print('Smoothness:')
    print(Fsmoothness(estimatedTheta,estimate_F))

    #============ CPAV ===============
    print('----------- CPAV ----------')
    print(date())
    start=proc.time()
    model=CyclicSPAV(X,Y,delta)
    print(proc.time()-start)
    estimate_F=model$F
    print('Prediction error:')
    print(predictError(Y,matrix(1,1,p),estimate_F))
    print('Smoothness:')
    print(Fsmoothness(matrix(1,1,p),estimate_F))
}


"
    =========================================
          Varying epsilon (p=100, n=100, s=10)
    =========================================
"
#epsilonList=c(0.4,0.6,0.8,1.0,1.2,1.4)
epsilonList=c()
p=100
n=100
delta=0.1
lambda=0.05
# Reading data
if(length(epsilonList)>0){
    print('Begin reading data....')
    theta=as.matrix(read.csv(sprintf('data/theta/Theta_%dp_%dn_10s.csv',p,n),header=FALSE))
    X=as.matrix(read.csv(sprintf('data/X/X_%dp_%dn10s.csv',p,n),header=FALSE))
    Y=as.matrix(read.csv(sprintf('data/Y/Y_%dp_%dn10s.csv',p,n),header=FALSE))
    print('Finish reading data!')
}
for (epsilon in epsilonList){
    print(sprintf('========= %.2f Epsilon when p=100 n=100 s=10 ============',epsilon))

    #========= SSAIM ==============
    print('-------- SSAIM ---------')
    totalTime=proc.time()-proc.time()
    for(i in 1:count){
        start=proc.time()
        ssaim=SSAIM(X,Y,epsilon,lambda,delta,iter=20)
        totalTime=totalTime+(proc.time()-start)
    }
    print(totalTime/count)
    # Evaluation
    estimate_F=ssaim$F
    estimatedTheta=ssaim$theta
    rm(ssaim)
    print('Theta error:')
    print(thetaError(t(theta),estimatedTheta))
    print('Prediction error:')
    print(predictError(Y,estimatedTheta,estimate_F))
    print('Smoothness:')
    print(Fsmoothness(estimatedTheta,estimate_F))

    #========= SLIM ==============
    print('-------- SLIM ---------')
    totalTime=proc.time()-proc.time()
    for(i in 1:count){
        start=proc.time()
        slim=SLIM(X,Y,lambda,iter=20)
        totalTime=totalTime+(proc.time()-start)
    }
    print(totalTime/count)
    # Evaluation
    estimate_F=slim$F
    estimatedTheta=slim$theta
    rm(slim)
    print('Theta error:')
    print(thetaError(t(theta),estimatedTheta))
    print('Prediction error:')
    print(predictError(Y,t(estimatedTheta),estimate_F))
    print('Smoothness:')
    print(Fsmoothness(estimatedTheta,estimate_F))
}
