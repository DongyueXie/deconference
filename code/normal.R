
## check everthing using normal distribution

n = 100
p = 3
Q = matrix(0.5,nrow=p,ncol=p) + diag(0.5,p)*3
V = matrix(0.5,nrow=p,ncol=p) + diag(0.5,p)
b = c(1,-2,3)
b2 = p:1
b3 = rep(1,p)

sigma_y = c(2,3,4)
set.seed(12345)
x = MASS::mvrnorm(n,rep(2,p),Q)
X = x + MASS::mvrnorm(n,rep(0,p),V)

y = x%*%b + rnorm(n,0,sigma_y[1])
y2 = x%*%b2 + rnorm(n,0,sigma_y[2])
y3 = x%*%b3 + rnorm(n,0,sigma_y[3])

Y = cbind(y,y2,y3)

A = t(X)%*%X/n
kappa(A-V)
kappa(A-diag(diag(V)))

eigen(A-V)$values

eigen(solve(V)%*%A)$values

Z = cbind(Y,X)
S = matrix(0,nrow=ncol(Y)+p,ncol=ncol(Y)+p)
S[((ncol(Y)+1):nrow(S)),((ncol(Y)+1):ncol(S))] = V
S[1:ncol(Y),1:ncol(Y)] = diag(sigma_y^2)

Mzz = t(Z)%*%Z/n
(eigen(solve(S)%*%Mzz)$values)



S = V
S = cbind(rep(0,p),S)
S = rbind(c(sigma_y[1]^2,rep(0,p)),S)
min(eigen(solve(S)%*%(crossprod(cbind(y,X))/n))$values)

S[1,1] = sigma_y[2]^2
min(eigen(solve(S)%*%(crossprod(cbind(y2,X))/n))$values)

S[1,1] = sigma_y[3]^2
min(eigen(solve(S)%*%(crossprod(cbind(y3,X))/n))$values)


S2 = S
S2[1,1] = 1e-5
min(eigen(solve(S2)%*%Mzz)$values)

S3 = diag(diag(V))
S3 = cbind(rep(0,p),S3)
S3 = rbind(c(sigma_y^2,rep(0,p)),S3)
min(eigen(solve(S3)%*%Mzz)$values)

solve(A-V)%*%t(X)%*%y/n
solve(A)%*%t(X)%*%y/n




##
lambda = 10
n=10
reps = 100
out = c()
for(i in  1:reps){
  x = rpois(n,lambda)
  out = rbind(out,c(mean(x),var(x)))
}
apply(out, 2, mean)
apply(out,2,sd)




##

fit1 = (lm(y~.-1,data.frame(y=y,x=x)))
fit2 = (lm(y~.-1,data.frame(y=y/abs(sum(y))*n,x=x)))
