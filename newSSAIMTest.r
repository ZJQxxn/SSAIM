source('SSAIM.r')
source('OtherMethods.r')
source('evaluate.r')

path=''
count = 1

pList = c(200)
parameter_list = list(c(1, 0.06, 10))

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
        ssaim=SSAIM(X,Y,epsilon,lambda,delta, iter = 1)
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
    
    indices = ssaim$indices
}
}
    