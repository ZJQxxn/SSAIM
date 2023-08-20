"
Description:
    Test SSAIM methods.

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
count = 1

epsilon = 0.5
lambda_list = c(0.02, 0.04, 0.06, 0.08, 0.1)
delta = 0.1

p = 200
n = 100
s = 10

p_record = c()
n_record = c()
name_record = c()
time_record = c()
MSE_record = c()
smoothness_record = c()
theta_error_recrd = c()
F_error_record = c()

print('Begin reading data....')
theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
F=as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn_%ds.csv',path,p,n,s),header=FALSE))
print('Finish reading data!')


for (lambda in lambda_list){
    print(sprintf('========= %f Lambda ============', lambda))

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
    write.csv(t(estimatedTheta), sprintf('./estimation/SSAIM_%dp_%dn_%.3flambda_theta.csv', p, n, lambda), row.names = FALSE)
    write.csv(estimate_F, sprintf('./estimation/SSAIM_%dp_%dn_%.3flambda_F.csv', p, n ,lambda), row.names = FALSE)
    write.csv(as.matrix(estimate_F %*% t(estimatedTheta)), sprintf('./estimation/SSAIM_%dp_%dn_%.3flambda_Y.csv', p, n, lambda), row.names = FALSE)
}
# Save to csv file
records = data.frame(p=p_record, n = n_record, name=name_record, time=time_record, mse = MSE_record, smoothness = smoothness_record, thta_error = theta_error_recrd, F_error = F_error_record)
write.csv(records, file = "diff_lambda_records_SSAIM.csv", append = FALSE, quote = TRUE, sep = " ")
print("Finished saving diff lambda data.")
print("====================================\n")