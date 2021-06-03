
source('code/deconference_main.R')
source('code/utils.R')
source('code/wols.R')
source('code/simulation/simu_correlation.R')

genecorr = readRDS("output/geneCor_gtexpancreas_tpm.rds")

xin_raw <- readRDS("data/pancreas/xin_raw.rds")
cell_types = c('alpha', 'beta', 'delta', 'gamma')
K = length(cell_types)
rm.indi = c("Non T2D 4","Non T2D 7","Non T2D 10","Non T2D 12")
rm.indi.idx = which(xin_raw$individual%in%rm.indi)

datax.xin = set_data_decon(Y = xin_raw[,-rm.indi.idx],cell_types = cell_types, gene_thresh = 0.05,max_count_quantile = 0.95,w=1)
design.mat.xin = scRef_multi_proc(datax.xin$Y,datax.xin$cell_type_idx,datax.xin$indi_idx,estimator="separate",est_sigma2 = TRUE,meta_mode = 'local',smooth.sigma = FALSE)

ref = design.mat.xin$X

inter.gene = intersect(genecorr$genes,rownames(ref))
ref = ref[match(inter.gene,rownames(ref)),]
sigma2 = design.mat.xin$Sigma
sigma2 = sigma2[match(inter.gene,rownames(ref)),]

rm.idx = which(rowSums(sigma2)==0)
if(length(rm.idx)>0){
  ref=ref[-rm.idx,]
  sigma2=sigma2[-rm.idx,]
}
ref = ref+1/nrow(ref)
sigma2 = sigma2 + 1/nrow(ref)


rm(datax.xin)
rm(design.mat.xin)

gene.Cov = genecorr$S
idx = match(inter.gene,genecorr$genes)
gene.R = cov2cor(gene.Cov[idx,idx])
rm(genecorr)

# make gene.R sparser
#'@param n keep #top pairs correlations and set all others to 0
#'@param R
#'@param random If true, random draw n pairs to keep; otherwise, top n pairs to keep
pruneR = function(R,n=nrow(R),random=F){

  R.up = R
  R.up[lower.tri(R.up)] = 0
  diag(R.up) = 0
  non0idx = which(R.up!=0)

  if(random){
    idx.temp = sample(non0idx,n)
  }else{
    r.non0 = R.up[non0idx]
    t = quantile(abs(r.non0),1-n/length(non0idx))
    idx.temp = non0idx[which(abs(r.non0)>t)]
  }

  R.up[-idx.temp] = 0
  R = R.up+t(R.up)
  diag(R) = 1
  # make sure R is positive definite
  min.ev = RSpectra::eigs(Matrix::Matrix(R,sparse = T),1,which = 'SA')$values
  if(min.ev<0.1){
    diag(R) = 1+abs(min.ev)+0.1
    adj = sqrt(1+abs(min.ev)+0.1)
    R = t(R/(rep(adj,nrow(R))))/rep(adj,nrow(R))
  }

  R

}


## 8724 is too large, sub-sample 1000 genes

set.seed(12345)
ii = sample(1:8724,1000)
ii = sort(ii)
ref.sub = ref[ii,]
sigma2.sub=sigma2[ii,]
gene.R.sub = gene.R[ii,ii]


#hist(R1[upper.tri(R1)&R1!=0],breaks = 100)
#sum(R1[upper.tri(R1)]!=0)/(nrow(R1)*(nrow(R1)+1)/2-nrow(R1))
#RSpectra::eigs(Matrix::Matrix(R1,sparse = T),1,which = 'SA')$values

K = 4
#set.seed(12345)
#b = matrix(runif(n*K),nrow=4)
#b = apply(b, 2, function(z){z/sum(z)})
b1 = c(0.1,0.1,0.3,0.5)
b2 = c(0.1,0.2,0.5,0.2)

#dim(ref)
#ref_test = ref[1:300,]
#tau2_test = tau2[1:300,]
#sigma2_test = sigma2[1:300,]
#simu_test = simu_twosample(ref_test,b.m,nreps=2,tau2=tau2_test,sigma2=sigma2_test,sc_lib_size = 0.1,n_indi = 6,nk=50,verbose = F,printevery = 1)
#get_summary(simu1,b=b.m,trans = T)
nb = 10
b.m = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))
set.seed(12345)
R1 = pruneR(gene.R.sub,random = T,n=1000)
simu = simu_corr_simple(ref.sub,b.m,nreps=100,
                        sigma2=sigma2.sub,R=R1,n_indi = 10,verbose = F,printevery = 1)



saveRDS(simu,file='output/simu_corr_xin_G1000_corpair1000.rds')
rm(simu)


set.seed(12345)
R1 = pruneR(gene.R.sub,random = T,n=5000)
simu = simu_corr_simple(ref.sub,b.m,nreps=100,
                        sigma2=sigma2.sub,R=R1,n_indi = 10,verbose = F,printevery = 1)



saveRDS(simu,file='output/simu_corr_xin_G1000_corpair5000.rds')
rm(simu)

set.seed(12345)
R1 = pruneR(gene.R.sub,random = T,n=10000)
simu = simu_corr_simple(ref.sub,b.m,nreps=100,
                        sigma2=sigma2.sub,R=R1,n_indi = 10,verbose = F,printevery = 1)



saveRDS(simu,file='output/simu_corr_xin_G1000_corpair10000.rds')
rm(simu)


#######
#simu.check = simu_twosample(ref.sub,b.m,nreps=10,tau2=sigma2.sub/100,sigma2=sigma2.sub,sc_lib_size = 0.1,n_indi = 10,nk=100,verbose = F,printevery = 1)
######



#
#
# n = 500
# K = 4
# set.seed(12345)
# ref = matrix(rnorm(n*K),nrow=n)
# ref = abs(ref)
# ref = apply(ref, 2, function(z){z/sum(z)})*n
#
# sigma2 = ref/2
# tau2 = sigma2/100
# rownames(ref) = 1:n
#
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
#
# hist(R[upper.tri(R)&R!=0],breaks = 100)
#
# sum(R[upper.tri(R)]!=0)/(n*(n+1)/2-n)
#
# b1 = c(0.1,0.1,0.3,0.5)
# b2 = c(0.1,0.2,0.5,0.2)
#
# nb = 4
# b.m = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))
#
# tt = simu_corr_simple(ref,b.m,R,nreps = 10,printevery = 1)
#
# #tt = simu_correlation(ref,b.m,R,nreps = 10,tau2=tau2,sigma2=sigma2,verbose = FALSE,n_indi = 8)
# tt$coverage_unadj_hc0
# tt$coverage_adj_hc0
