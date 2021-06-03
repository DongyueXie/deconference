source('code/deconference_main.R')
source('code/utils.R')
source('code/simulation/simu_twoSample.R')


xin_raw <- readRDS("data/pancreas/xin_raw.rds")
cell_types = c('alpha', 'beta', 'delta', 'gamma')
K = length(cell_types)
rm.indi = c("Non T2D 4","Non T2D 7","Non T2D 10","Non T2D 12")
rm.indi.idx = which(xin_raw$individual%in%rm.indi)

datax.xin = set_data_decon(Y = xin_raw[,-rm.indi.idx],cell_types = cell_types, gene_thresh = 0.05,max_count_quantile = 0.95,w=1)
design.mat.xin = scRef_multi_proc(datax.xin$Y,datax.xin$cell_type_idx,datax.xin$indi_idx,estimator="separate",est_sigma2 = TRUE,meta_mode = 'local',smooth.sigma = FALSE)

ref = design.mat.xin$X
sigma2 = design.mat.xin$Sigma
rm.idx = which(rowSums(sigma2)==0)
if(length(rm.idx)>0){
  ref=ref[-rm.idx,]
  sigma2=sigma2[-rm.idx,]
}
ref = ref+1/nrow(ref)
sigma2 = sigma2 + 1/nrow(ref)
tau2 = sigma2/100

rm(datax.xin)
rm(design.mat.xin)



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
simu = simu_twosample(ref,b.m,nreps=100,tau2=tau2,sigma2=sigma2,sc_lib_size = 0.1,n_indi = 10,nk=100,verbose = F,printevery = 1)
saveRDS(simu,file='output/simu_twosample_xin_nb10.rds')
rm(simu)

nb = 50
b.m = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))
set.seed(12345)
simu = simu_twosample(ref,b.m,nreps=100,tau2=tau2,sigma2=sigma2,sc_lib_size = 0.1,n_indi = 10,nk=100,verbose = F,printevery = 1)
saveRDS(simu,file='output/simu_twosample_xin_nb50.rds')
rm(simu)




seger <- readRDS("data/pancreas/segerstolpe_raw.rds")
datax.seger = set_data_decon(Y = seger,cell_types = cell_types, gene_thresh = 0.05,max_count_quantile = 0.95,w=1)
design.mat.seger = scRef_multi_proc(datax.seger$Y,datax.seger$cell_type_idx,datax.seger$indi_idx,estimator="separate",est_sigma2 = TRUE,meta_mode = 'local',smooth.sigma = FALSE)

ref = design.mat.seger$X
sigma2 = design.mat.seger$Sigma
rm.idx = which(rowSums(sigma2)==0)
if(length(rm.idx)>0){
  ref=ref[-rm.idx,]
  sigma2=sigma2[-rm.idx,]
}
ref = ref+1/nrow(ref)
sigma2 = sigma2 + 1/nrow(ref)
tau2 = sigma2/100

rm(datax.seger)
rm(design.mat.seger)




nb = 10
b.m = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))
set.seed(12345)
simu = simu_twosample(ref,b.m,nreps=100,tau2=tau2,sigma2=sigma2,sc_lib_size = 0.1,n_indi = 10,nk=100,verbose = F,printevery = 1)
saveRDS(simu,file='output/simu_twosample_seger_nb10.rds')
rm(simu)


nb = 50
b.m = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))
set.seed(12345)
simu = simu_twosample(ref,b.m,nreps=100,tau2=tau2,sigma2=sigma2,sc_lib_size = 0.1,n_indi = 10,nk=100,verbose = F,printevery = 1)
saveRDS(simu,file='output/simu_twosample_seger_nb50.rds')
rm(simu)




baron = readRDS("data/pancreas/baron.rds")
datax = set_data_decon(Y = baron,cell_types = cell_types, gene_thresh = 0.05,max_count_quantile = 0.95,w=1)
design.mat = scRef_multi_proc(datax$Y,datax$cell_type_idx,datax$indi_idx,estimator="separate",est_sigma2 = TRUE,meta_mode = 'local',smooth.sigma = FALSE)

ref = design.mat$X
sigma2 = design.mat$Sigma
rm.idx = which(rowSums(sigma2)==0)
if(length(rm.idx)>0){
  ref=ref[-rm.idx,]
  sigma2=sigma2[-rm.idx,]
}
ref = ref+1/nrow(ref)
sigma2 = sigma2 + 1/nrow(sigma2)
tau2 = sigma2/100

rm(datax)
rm(design.mat)


nb = 10
b.m = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))
set.seed(12345)
simu = simu_twosample(ref,b.m,nreps=100,tau2=tau2,sigma2=sigma2,sc_lib_size = 0.1,n_indi = 10,nk=100,verbose = F,printevery = 1)
saveRDS(simu,file='output/simu_twosample_baron_nb10.rds')
rm(simu)


nb = 50
b.m = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))
set.seed(12345)
simu = simu_twosample(ref,b.m,nreps=100,tau2=tau2,sigma2=sigma2,sc_lib_size = 0.1,n_indi = 10,nk=100,verbose = F,printevery = 1)
saveRDS(simu,file='output/simu_twosample_baron_nb50.rds')
rm(simu)
