source('code/deconference_main.R')
source('code/utils.R')
source('code/wols.R')
source('code/simulation/simu_correlation.R')

n = 500
K = 4
set.seed(12345)
ref = matrix(rnorm(n*K),nrow=n)
ref = abs(ref)
ref = apply(ref, 2, function(z){z/sum(z)})*n

sigma2 = ref/2
#tau2 = sigma2/100
rownames(ref) = 1:n

# library(Matrix)
# set.seed(12345)
# L = rsparsematrix(n,n,0.01,symmetric = T)
# L[upper.tri(L)] = 0
# diag(L) = 0.5
# non0.idx = which(L<0)
# temp.idx = sample(non0.idx,0.8*length(non0.idx))
# L[temp.idx] = abs(L[temp.idx])
# R = tcrossprod(L)
# R = R + 0.1*diag(n)
# R = cov2cor(R)

R = as.matrix(bandSparse(n,k=c(0,1),diagonals = list(rep(1,n),rep(0.45,n)),symmetric = T))


b1 = c(0.1,0.1,0.3,0.5)
b2 = c(0.1,0.2,0.5,0.2)

nb = 4
b.m = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))

set.seed(12345)
sim0 = simu_corr_simple(ref,b.m,R=diag(n),sigma2 = sigma2,
                        nreps = 500,n_indi = 10,verbose = F,printevery = 50)

mean(sim0$coverage_unadj_hc0)


set.seed(12345)
sim1 = simu_corr_simple(ref,b.m,R=R,sigma2 = sigma2,
                        nreps = 500,n_indi = 10,verbose = F,printevery = 50)

#tt = simu_correlation(ref,b.m,R,nreps = 10,tau2=tau2,sigma2=sigma2,verbose = FALSE,n_indi = 8)
mean(sim1$coverage_unadj_hc0)
mean(sim1$coverage_adj_hc0)



p = 500
A = matrix(0,nrow=p,ncol=p)

for(i in 1:p){
  for(j in 1:p){
    A[i,j] = max(1-abs(i-j)/2,0)
  }
}


set.seed(12345)
sim2 = simu_corr_simple(ref,b.m,R=A,sigma2 = sigma2,
                        nreps = 500,n_indi = 10,verbose = F,printevery = 50)

#tt = simu_correlation(ref,b.m,R,nreps = 10,tau2=tau2,sigma2=sigma2,verbose = FALSE,n_indi = 8)
mean(sim2$coverage_unadj_hc0)
mean(sim2$coverage_adj_hc0)



