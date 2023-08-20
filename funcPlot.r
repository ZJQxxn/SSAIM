source('SSAIM.r')

# settings
n = 200
p = 200
path = ''

theta=as.matrix(read.csv(sprintf('%sdata/theta/Theta_%dp_%dn.csv',path,p,n),header=FALSE))
X=as.matrix(read.csv(sprintf('%sdata/X/X_%dp_%dn.csv',path,p,n),header=FALSE))
Y=as.matrix(read.csv(sprintf('%sdata/Y/Y_%dp_%dn.csv',path,p,n),header=FALSE))
F=as.matrix(read.csv(sprintf('%sdata/F/F_%dp_%dn.csv',path,p,n),header=FALSE))

plt_x = X[, 1]
plt_y = F[, 1]

plot(plt_x, plt_y, type = "p")

# fit SSAIM model
epsilon=1
lambda=0.01
delta=0.01
ssaim=SSAIM(X,Y,epsilon,lambda,delta,iter=20)
estimate_F=ssaim$F
estimatedTheta=ssaim$theta

plot(plt_x, estimate_F[,1], col = "red", type = "p")