"
Description:  
    Utilize SSAIM on Noise dataset.

Programmer: 
    Jiaqi Zhang

Update history: 
    --2019/04/16: Create codes and finish the experiment.
"
source('SSAIM.r')
source('OtherMethods.r')
source('evaluate.r')


print('Begin reading data....')
theta=as.matrix(read.csv(sprintf('dataset/noise_theta.csv'),header=FALSE))
X=as.matrix(read.csv(sprintf('dataset/noise_X.csv'),header=FALSE))
Y=as.matrix(read.csv(sprintf('dataset/noise_Y.csv'),header=FALSE))
print('Finish reading data!')
for (i in 1:length(X[1,])){
    X[,i]=X[,i]/norm(X[,i],type='2')
}
Y = Y/norm(Y, type='2')
print('================ Noise Data ==================')
varNum=length(X[1,])
sampNum=length(X[,1])
count=1
epsilon=1
lambda=0.01
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

#========= SLIM ==============
print('-------- SLIM ---------')
totalTime=proc.time()-proc.time()
for(i in 1:count){
    start=proc.time()
    slim=SLIM(X,Y,lambda)
    totalTime=totalTime+(proc.time()-start)
}
print(totalTime/count)
# Evaluation
estimate_F=slim$F
estimatedTheta=slim$theta
print('Theta error:')
print(thetaError(theta,estimatedTheta))
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
