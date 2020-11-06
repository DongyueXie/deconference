# test mvmeta package
library(mvmeta)
library(MASS)
set.seed(12345)
k = 4
n = 100
mu = c(1,1,2,2)

Sigma = (rWishart(1,k,diag(k)/k))[,,1]
S = matrix(runif(n*k),nrow=n,ncol=k)

reps = 50

mse = matrix(nrow=reps,ncol=k)

for(rr in 1:reps){

  y = matrix(nrow=n,ncol=k)
  for(i in 1:n){
    y[i,] = MASS::mvrnorm(1,mu,Sigma+diag(S[i,]))
  }

  fit1 = mvmeta(y~1,S=S,data=data.frame(y=y),method = 'reml')
  mse[rr,1] = norm(fit1$Psi-Sigma,'2')^2

  fit2 = mvmeta(y~1,S=S,data=data.frame(y=y),method = 'ml')
  mse[rr,2] = norm(fit2$Psi-Sigma,'2')^2

  fit3 = mvmeta(y~1,S=S,data=data.frame(y=y),method = 'mm')
  mse[rr,3] = norm(fit3$Psi-Sigma,'2')^2

  fit4 = mvmeta(y~1,S=S,data=data.frame(y=y),method = 'vc')
  mse[rr,4] = norm(fit4$Psi-Sigma,'2')^2


}

sqrt(colMeans(mse))

