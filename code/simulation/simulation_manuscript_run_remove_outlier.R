xin_raw <- readRDS("data/pancreas/xin_raw.rds")
cell_types = c('alpha', 'beta', 'delta', 'gamma')
K = length(cell_types)
rm.indi = c("Non T2D 4","Non T2D 7","Non T2D 10","Non T2D 12")
rm.indi.idx = which(xin_raw$individual%in%rm.indi)

datax.xin = set_data_decon(Y = xin_raw[,-rm.indi.idx],cell_types = cell_types,
                           gene_thresh = 0.05,max_count_quantile_celltype = 0.95,
                           max_count_quantile_indi = 0.95,
                           w=1)
design.mat.xin = scRef_multi_proc(datax.xin$Y,datax.xin$cell_type_idx,
                                  datax.xin$indi_idx,estimator="separate",
                                  est_sigma2 = TRUE,meta_mode = 'local',smooth.sigma = F)

b1 = c(0.1,0.1,0.3,0.5)
b2 = c(0.1,0.2,0.5,0.2)
nb = 10
b = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))

ref = design.mat.xin$X
sigma2 = design.mat.xin$Sigma

# h = diag(ref%*%solve(t(ref)%*%ref - matrix(colSums(design.mat.xin$Vg),nrow=4))%*%t(ref))
# plot(h)
# idx = which(h>0.01)

h2 = diag(sigma2%*%solve(t(sigma2)%*%sigma2)%*%t(sigma2))
plot(h2)
idx = which(h2>0.005)

round(ref[idx,],2)
round(sigma2[idx,],2)

ref = ref[-idx,]
sigma2 = sigma2[-idx,]

ref = ref+1/nrow(ref)
sigma2 = sigma2 + 1/nrow(ref)



# G = nrow(ref)
# K = 4
# d = 500
# A = matrix(0,nrow=G,ncol=G)
#
# for(i in 1:G){
#   for(j in i:min(i+d,G)){
#     A[i,j] = max(1-abs(i-j)/d,0)
#   }
# }
# A = A+t(A) - diag(G)
# library(Matrix)
# A = Matrix(A,sparse = TRUE)
#
# set.seed(12345)
# simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1)
# saveRDS(simu_out,file = 'output/manuscript/simulation_10bulk_500genecor_fdr05.rds')
#
#
# d = 300
# A = matrix(0,nrow=G,ncol=G)
#
# for(i in 1:G){
#   for(j in i:min(i+d,G)){
#     A[i,j] = max(1-abs(i-j)/d,0)
#   }
# }
# A = A+t(A) - diag(G)
# library(Matrix)
# A = Matrix(A,sparse = TRUE)
# set.seed(12345)
# simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1)
# saveRDS(simu_out,file = 'output/manuscript/simulation_10bulk_300genecor_fdr05.rds')

set.seed(12345)
simu_out = simu_study(ref,b,R=NULL,sigma2,printevery = 1)
saveRDS(simu_out,file = 'output/manuscript/simulation_10bulk_0genecor_fdr05_remove_large_var.rds')

##################################################
ref = design.mat.xin$X
sigma2 = design.mat.xin$Sigma


h2 = diag(sigma2%*%solve(t(sigma2)%*%sigma2)%*%t(sigma2))
idx = which(h2>0.001)


ref = ref[-idx,]
sigma2 = sigma2[-idx,]

ref = ref+1/nrow(ref)
sigma2 = sigma2 + 1/nrow(ref)
set.seed(12345)
simu_out = simu_study(ref,b,R=NULL,sigma2,printevery = 1)
saveRDS(simu_out,file = 'output/manuscript/simulation_10bulk_0genecor_fdr05_remove_large_var2.rds')
##################################################
ref = design.mat.xin$X
sigma2 = design.mat.xin$Sigma
h2 = diag(sigma2%*%solve(t(sigma2)%*%sigma2)%*%t(sigma2))
idx = which(h2>0.0001)
ref = ref[-idx,]
sigma2 = sigma2[-idx,]

ref = ref+1/nrow(ref)
sigma2 = sigma2 + 1/nrow(ref)
set.seed(12345)
simu_out = simu_study(ref,b,R=NULL,sigma2,printevery = 1)
saveRDS(simu_out,file = 'output/manuscript/simulation_10bulk_0genecor_fdr05_remove_large_var3.rds')
