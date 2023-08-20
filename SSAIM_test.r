"
Description:  
    Test SSAIM implementation.

Programmer: 
    Jiaqi Zhang

"

source('SSAIM.r')
source('evaluate.r')

library(scar)


# Generate a 2-D dataset
X=matrix(0,10,2)
X[,1]=1:10
X[,2]=1:10
theta=c(1,1)
Y=X[,1]+X[,2]+rnorm(10)
smoothPar=0.1
sparsePar=0.1
epsilon=0.1

# Solve SSAIM
model=SSAIM(X,Y,epsilon,sparsePar,smoothPar)
estimate_F=model$F
estimatedTheta=model$theta

# Evaluation
print('Theta error:')
print(thetaError(theta,estimatedTheta))
print('Prediction error:')
print(predictError(Y,estimatedTheta,estimate_F))
print('Smoothness:')
print(smoothness(estimatedTheta,estimate_F))

plot(0,0,xlim = c(0,10),ylim = c(0,5),type = "n")
lines(X[,1],estimate_F[,1],type='b', col = "blue")


scar_model= scar(X, Y, shape = c("cvxin", "cvxin"), family = gaussian())
scar_esimated_F = scar_model$componentfit
lines(X[,1],scar_esimated_F[,1],type='b', col = "red")


