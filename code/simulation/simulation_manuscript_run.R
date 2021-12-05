#
# xin_raw <- readRDS("data/pancreas/xin_raw.rds")
# cell_types = c('alpha', 'beta', 'delta', 'gamma')
# K = length(cell_types)
# rm.indi = c("Non T2D 4","Non T2D 7","Non T2D 10","Non T2D 12")
# rm.indi.idx = which(xin_raw$individual%in%rm.indi)
#
# datax.xin = set_data_decon(Y = xin_raw[,-rm.indi.idx],cell_types = cell_types,
#                            gene_thresh = 0.05,max_count_quantile_celltype = 0.95,
#                            max_count_quantile_indi = 0.95,
#                            w=1)
# design.mat.xin = scRef_multi_proc(datax.xin$Y,datax.xin$cell_type_idx,
#                                   datax.xin$indi_idx,estimator="separate",
#                                   est_sigma2 = TRUE,meta_mode = 'local',smooth.sigma = F)
#
# ref = design.mat.xin$X
# sigma2 = design.mat.xin$Sigma
#
# ref = ref+1/nrow(ref)
# sigma2 = sigma2 + 1/nrow(ref)
#
# saveRDS(list(ref=ref,sigma2=sigma2),file='data/pancreas/xin_ref_sigma9496.rds')

########################################################
########### use neuron data for simulation #############
source('code/simulation/simulation_manuscript.R')
indis_ref = readRDS('data/neuron/indis_ref_12400by6by97.rds')
ref = apply(indis_ref,c(1,2),mean,na.rm=TRUE)
sigma2 = apply(indis_ref,c(1,2),var,na.rm=TRUE)


G = nrow(ref)
K = ncol(ref)
d = 500
A = matrix(0,nrow=G,ncol=G)
for(i in 1:G){
  for(j in i:min(i+d,G)){
    A[i,j] = max(1-abs(i-j)/d,0)
  }
}
A = A+t(A) - diag(G)
library(Matrix)
A = Matrix(A,sparse = TRUE)
alpha.cors = c(0)
cases = c("null")
nbs = c(4)
dirichlet.scale = c(5)

# test
# temp = simu_study(ref[1:100,],sigma2[1:100,],c(0.5,0.3,0.1,0.1),c(0.5,0.3,0.1,0.1),
#                   n_bulk = 50,dirichlet.scale=5,
#                   R=A[1:100,1:100],printevery = 1,est_cor=FALSE,nreps = 5)

set.seed(12345)
for(case in cases){

  if(case=='null'){
    #p1 = c(0.3,0.2,0.15,0.15,0.1,0.1)
    #p2 = c(0.3,0.2,0.15,0.15,0.1,0.1)
    p1 = c(0.15,0.15,0.1,0.1,0.2,0.3)
    p2 = c(0.15,0.15,0.1,0.1,0.2,0.3)
  }else if(case=='all_diff'){
    p1 = c(0.15,0.15,0.1,0.1,0.2,0.3)
    p2 = c(0.1,0.1,0.2,0.3,0.15,0.15)
  }

  for(nb in nbs){
    for(aa in dirichlet.scale){
      for(alpha.cor in alpha.cors){

        if(alpha.cor==0){
          est_cor = FALSE
          cor.status = 'trueR'
        }else{
          est_cor=TRUE
          cor.status = paste('cor0',alpha.cor*10,sep = '')
        }

        print(paste('Running:',case,'nb=',nb,'cor:',cor.status,'aa=',aa))

        simu_out2 = simu_study(ref,sigma2,p1,p2,n_bulk = nb,dirichlet.scale=aa,nreps = 100,
                              R=A,printevery = 1,alpha.cor = alpha.cor,est_cor=est_cor)
        #saveRDS(simu_out,file = paste('output/manuscript/simulation/simulation_',nb,'bulk_500genecor_',cor.status,'_',case,'_dirichlet',aa,'no_pd.rds',sep=''))

      }
    }
  }
}

########10/20evening/2021#####################

# do not use real reference and sigma2
source('code/simulation/simulation_manuscript.R')
G = 500
K = 4
set.seed(12345)
ref = matrix(rnorm(G*K),nrow=G)
ref = abs(ref)
ref = apply(ref, 2, function(z){z/sum(z)})*G
sigma2 = ref/2


d = 25
A = matrix(0,nrow=G,ncol=G)
for(i in 1:G){
  for(j in i:min(i+d,G)){
    A[i,j] = max(1-abs(i-j)/d,0)
  }
}
A = A+t(A) - diag(G)
library(Matrix)
A = Matrix(A,sparse = TRUE)
alpha.cors = c(0,0.1,0.5)
cases = c("null","all_diff")
nbs = c(100)
dirichlet.scale = c(10,5)

set.seed(12345)
for(case in cases){

  if(case=='null'){
    p1 = c(0.5,0.3,0.1,0.1)
    p2 = c(0.5,0.3,0.1,0.1)
  }else if(case=='all_diff'){
    p1 = c(0.15,0.2,0.45,0.2)
    p2 = c(0.1,0.1,0.3,0.5)
  }

  for(nb in nbs){
    for(aa in dirichlet.scale){
      for(alpha.cor in alpha.cors){

        if(alpha.cor==0){
          est_cor = FALSE
          cor.status = 'trueR'
        }else{
          est_cor=TRUE
          cor.status = paste('cor0',alpha.cor*10,sep = '')
        }

        print(paste('Running:',case,'nb=',nb,'cor:',cor.status,'aa=',aa))

        simu_out = simu_study(ref,sigma2,p1,p2,n_bulk = nb,dirichlet.scale=aa,
                              R=A,printevery = 1,alpha.cor = alpha.cor,est_cor=est_cor,nreps = 100)
        saveRDS(simu_out,file = paste('output/manuscript/simulation/test/',nb,'bulk_500genecor_',cor.status,'_',case,'_dirichlet',aa,'no_pd.rds',sep=''))

      }
    }
  }
}


#######################################

########10/20/2021#####################

# random p, do not make.pos

library(gtools)
source('code/simulation/simulation_manuscript.R')
xin = readRDS('data/pancreas/xin_ref_sigma9496.rds')
ref = xin$ref
sigma2 = xin$sigma2
G = nrow(ref)
K = ncol(ref)
d = 500
A = matrix(0,nrow=G,ncol=G)
for(i in 1:G){
  for(j in i:min(i+d,G)){
    A[i,j] = max(1-abs(i-j)/d,0)
  }
}
A = A+t(A) - diag(G)
library(Matrix)
A = Matrix(A,sparse = TRUE)
alpha.cors = c(0,0.1,0.5,0.8)
cases = c("null","all_diff")
nbs = c(100)
dirichlet.scale = c(5,10,100)

# test
# temp = simu_study(ref[1:100,],sigma2[1:100,],c(0.5,0.3,0.1,0.1),c(0.5,0.3,0.1,0.1),
#                   n_bulk = 50,dirichlet.scale=5,
#                   R=A[1:100,1:100],printevery = 1,est_cor=FALSE,nreps = 5)

set.seed(12345)
for(case in cases){

  if(case=='null'){
    p1 = c(0.5,0.3,0.1,0.1)
    p2 = c(0.5,0.3,0.1,0.1)
  }else if(case=='all_diff'){
    p1 = c(0.15,0.2,0.45,0.2)
    p2 = c(0.1,0.1,0.3,0.5)
  }

  for(nb in nbs){
    for(aa in dirichlet.scale){
      for(alpha.cor in alpha.cors){

        if(alpha.cor==0){
          est_cor = FALSE
          cor.status = 'trueR'
        }else{
          est_cor=TRUE
          cor.status = paste('cor0',alpha.cor*10,sep = '')
        }

        print(paste('Running:',case,'nb=',nb,'cor:',cor.status,'aa=',aa))

        simu_out = simu_study(ref,sigma2,p1,p2,n_bulk = nb,dirichlet.scale=aa,
                              R=A,printevery = 1,alpha.cor = alpha.cor,est_cor=est_cor)
        saveRDS(simu_out,file = paste('output/manuscript/simulation/simulation_',nb,'bulk_500genecor_',cor.status,'_',case,'_dirichlet',aa,'no_pd.rds',sep=''))

      }
    }
  }
}

#######################################

########10/19/2021#####################

# try to fix p while keep all other things unchanged

source('code/simulation/simulation_manuscript.R')
xin = readRDS('data/pancreas/xin_ref_sigma9496.rds')
ref = xin$ref
sigma2 = xin$sigma2
G = nrow(ref)
K = 4
d = 500
A = matrix(0,nrow=G,ncol=G)
for(i in 1:G){
  for(j in i:min(i+d,G)){
    A[i,j] = max(1-abs(i-j)/d,0)
  }
}
A = A+t(A) - diag(G)
library(Matrix)
A = Matrix(A,sparse = TRUE)
alpha.cors = c(0,0.1,0.5)
cases = c("null","all_diff")
nbs = c(100)


# test
# temp = simu_study(ref[1:100,],sigma2[1:100,],c(0.5,0.3,0.1,0.1),c(0.5,0.3,0.1,0.1),
#                   n_bulk = 50,dirichlet.scale=5,
#                   R=A[1:100,1:100],printevery = 1,est_cor=FALSE,nreps = 5)

set.seed(12345)
for(case in cases){

  if(case=='null'){
    p1 = c(0.5,0.3,0.1,0.1)
    p2 = c(0.5,0.3,0.1,0.1)
  }else if(case=='all_diff'){
    p1 = c(0.15,0.2,0.45,0.2)
    p2 = c(0.1,0.1,0.3,0.5)
  }

  for(nb in nbs){

      for(alpha.cor in alpha.cors){

        if(alpha.cor==0){
          est_cor = FALSE
          cor.status = 'trueR'
        }else{
          est_cor=TRUE
          cor.status = paste('cor0',alpha.cor*10,sep = '')
        }

        print(paste('Running:',case,'nb=',nb,'cor:',cor.status))

        simu_out = simu_study(ref,sigma2,p1,p2,n_bulk = nb,dirichlet=FALSE,
                              R=A,printevery = 1,alpha.cor = alpha.cor,est_cor=est_cor)
        saveRDS(simu_out,file = paste('output/manuscript/simulation/fixp/simulation_',nb,'bulk_500genecor_',cor.status,'_',case,'_fixp','.rds',sep=''))

      }

  }
}
##################################

########10/12/2021#####################

# do not use real reference and sigma2

G = 500
set.seed(12345)
ref = matrix(rnorm(G*K),nrow=G)
ref = abs(ref)
ref = apply(ref, 2, function(z){z/sum(z)})*G
sigma2 = ref/2

K = 4
d = 25
A = matrix(0,nrow=G,ncol=G)
for(i in 1:G){
  for(j in i:min(i+d,G)){
    A[i,j] = max(1-abs(i-j)/d,0)
  }
}
A = A+t(A) - diag(G)
library(Matrix)
A = Matrix(A,sparse = TRUE)
alpha.cors = c(0)
cases = c("null","all_diff")
nbs = c(50,100)
dirichlet.scale = c(5,10)

set.seed(12345)
for(case in cases){

  if(case=='null'){
    p1 = c(0.5,0.3,0.1,0.1)
    p2 = c(0.5,0.3,0.1,0.1)
  }else if(case=='all_diff'){
    p1 = c(0.15,0.2,0.45,0.2)
    p2 = c(0.1,0.1,0.3,0.5)
  }

  for(nb in nbs){
    for(aa in dirichlet.scale){
      for(alpha.cor in alpha.cors){

        if(alpha.cor==0){
          est_cor = FALSE
          cor.status = 'trueR'
        }else{
          est_cor=TRUE
          cor.status = paste('cor0',alpha.cor*10,sep = '')
        }

        print(paste('Running:',case,'nb=',nb,'cor:',cor.status,'aa=',aa))

        simu_out = simu_study(ref,sigma2,p1,p2,n_bulk = nb,dirichlet.scale=aa,
                              R=A,printevery = 1,alpha.cor = alpha.cor,est_cor=est_cor,nreps = 100)
        saveRDS(simu_out,file = paste('output/manuscript/simulation/test/',nb,'bulk_500genecor_',cor.status,'_',case,'_dirichlet',aa,'.rds',sep=''))

      }
    }
  }
}


#######################################

########10/06/2021#####################
library(gtools)
source('code/simulation/simulation_manuscript.R')
xin = readRDS('data/pancreas/xin_ref_sigma9496.rds')
ref = xin$ref
sigma2 = xin$sigma2
G = nrow(ref)
K = 4
d = 500
A = matrix(0,nrow=G,ncol=G)
for(i in 1:G){
  for(j in i:min(i+d,G)){
    A[i,j] = max(1-abs(i-j)/d,0)
  }
}
A = A+t(A) - diag(G)
library(Matrix)
A = Matrix(A,sparse = TRUE)
alpha.cors = c(0,0.1,0.5)
cases = c("null","all_diff")
nbs = c(10,50,100)
dirichlet.scale = c(5,10)

# test
# temp = simu_study(ref[1:100,],sigma2[1:100,],c(0.5,0.3,0.1,0.1),c(0.5,0.3,0.1,0.1),
#                   n_bulk = 50,dirichlet.scale=5,
#                   R=A[1:100,1:100],printevery = 1,est_cor=FALSE,nreps = 5)

set.seed(12345)
for(case in cases){

  if(case=='null'){
    p1 = c(0.5,0.3,0.1,0.1)
    p2 = c(0.5,0.3,0.1,0.1)
  }else if(case=='all_diff'){
    p1 = c(0.15,0.2,0.45,0.2)
    p2 = c(0.1,0.1,0.3,0.5)
  }

  for(nb in nbs){
    for(aa in dirichlet.scale){
      for(alpha.cor in alpha.cors){

        if(alpha.cor==0){
          est_cor = FALSE
          cor.status = 'trueR'
        }else{
          est_cor=TRUE
          cor.status = paste('cor0',alpha.cor*10,sep = '')
        }

        print(paste('Running:',case,'nb=',nb,'cor:',cor.status,'aa=',aa))

        simu_out = simu_study(ref,sigma2,p1,p2,n_bulk = nb,dirichlet.scale=aa,
                              R=A,printevery = 1,alpha.cor = alpha.cor,est_cor=est_cor)
        saveRDS(simu_out,file = paste('output/manuscript/simulation/simulation_',nb,'bulk_500genecor_',cor.status,'_',case,'_dirichlet',aa,'.rds',sep=''))

      }
    }
  }
}

#######################################


########08/18/2021###################
source('code/simulation/simulation_manuscript.R')
xin = readRDS('data/pancreas/xin_ref_sigma9496.rds')
ref = xin$ref
sigma2 = xin$sigma2
b1 = c(0.1,0.1,0.3,0.5)
b2 = c(0.1,0.2,0.5,0.2)
nb = 10
b = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))

G = nrow(ref)
K = 4
d = 500
A = matrix(0,nrow=G,ncol=G)


for(i in 1:G){
  for(j in i:min(i+d,G)){
    A[i,j] = max(1-abs(i-j)/d,0)
  }
}
A = A+t(A) - diag(G)
library(Matrix)
A = Matrix(A,sparse = TRUE)

set.seed(12345)
simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1)
saveRDS(simu_out,file = 'output/manuscript/simulation_10bulk_500genecor_fdr05.rds')




##############09/28/2021###############
b1 = c(0.5,0.3,0.1,0.1)
b2 = c(0.5,0.3,0.1,0.1)
nb = 50
b = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))
set.seed(12345)
simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1)
saveRDS(simu_out,file = 'output/manuscript/simulation_50bulk_500genecor_fdr05_null.rds')

b1 = c(0.5,0.3,0.1,0.1)
b2 = c(0.1,0.1,0.3,0.5)
nb = 30
b = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))
set.seed(12345)
simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1)
saveRDS(simu_out,file = 'output/manuscript/simulation_50bulk_500genecor_fdr05_all_diff.rds')

b1 = c(0.1,0.1,0.3,0.5)
b2 = c(0.1,0.15,0.4,0.35)
nb = 30
b = cbind(b1%*%t(rep(1,nb/2)),b2%*%t(rep(1,nb/2)))
set.seed(12345)
simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1)
saveRDS(simu_out,file = 'output/manuscript/simulation_50bulk_500genecor_fdr05_one_diff.rds')
#######################################


##############10/05/2021###############
library(gtools)
source('code/simulation/simulation_manuscript.R')
xin = readRDS('data/pancreas/xin_ref_sigma9496.rds')
ref = xin$ref
sigma2 = xin$sigma2
G = nrow(ref)
K = 4
d = 500
A = matrix(0,nrow=G,ncol=G)
for(i in 1:G){
  for(j in i:min(i+d,G)){
    A[i,j] = max(1-abs(i-j)/d,0)
  }
}
A = A+t(A) - diag(G)
library(Matrix)
A = Matrix(A,sparse = TRUE)


alpha.cors = 0.5


for(nb in c(50,100)){
  set.seed(12345)
  b1 = t(rdirichlet(nb/2,p1*10))
  b2 = t(rdirichlet(nb/2,p2*10))
  b = cbind(b1,b2)
  set.seed(12345)
  simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1,alpha.cor = alpha.cor)
  saveRDS(simu_out,file = paste('output/manuscript/simulation/simulation_',nb,'bulk_500genecor_fdr0',alpha.cor*10,'_null_dirichlet.rds',sep=''))
}


## all diff
p1 = c(0.5,0.3,0.1,0.1)
p2 = c(0.1,0.1,0.3,0.5)
nb = 100
# generate group proportions using dirichlet distribution
set.seed(12345)
b1 = t(rdirichlet(nb/2,p1*10))
b2 = t(rdirichlet(nb/2,p2*10))
b = cbind(b1,b2)


set.seed(12345)
simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1,alpha.cor = alpha.cor)
saveRDS(simu_out,file = paste('output/manuscript/simulation/simulation_',nb,'bulk_500genecor_fdr0',alpha.cor*10,'_all_diff_dirichlet.rds',sep=''))

nb = 50
# generate group proportions using dirichlet distribution
set.seed(12345)
b1 = t(rdirichlet(nb/2,p1*10))
b2 = t(rdirichlet(nb/2,p2*10))
b = cbind(b1,b2)


set.seed(12345)
simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1,alpha.cor = alpha.cor)
saveRDS(simu_out,file = paste('output/manuscript/simulation/simulation_',nb,'bulk_500genecor_fdr0',alpha.cor*10,'_all_diff_dirichlet.rds',sep=''))

## one diff
p1 = c(0.1,0.1,0.3,0.5)
p2 = c(0.1,0.15,0.4,0.35)
nb = 100
# generate group proportions using dirichlet distribution
set.seed(12345)
b1 = t(rdirichlet(nb/2,p1*10))
b2 = t(rdirichlet(nb/2,p2*10))
b = cbind(b1,b2)


set.seed(12345)
simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1,alpha.cor = alpha.cor)
saveRDS(simu_out,file = paste('output/manuscript/simulation/simulation_',nb,'bulk_500genecor_fdr0',alpha.cor*10,'_one_diff_dirichlet.rds',sep=''))

nb = 50
# generate group proportions using dirichlet distribution
set.seed(12345)
b1 = t(rdirichlet(nb/2,p1*10))
b2 = t(rdirichlet(nb/2,p2*10))
b = cbind(b1,b2)


set.seed(12345)
simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1,alpha.cor = alpha.cor)
saveRDS(simu_out,file = paste('output/manuscript/simulation/simulation_',nb,'bulk_500genecor_fdr0',alpha.cor*10,'_one_diff_dirichlet.rds',sep=''))

#######################################

##############10/06/2021###############
library(gtools)
source('code/simulation/simulation_manuscript.R')
xin = readRDS('data/pancreas/xin_ref_sigma9496.rds')
ref = xin$ref
sigma2 = xin$sigma2
G = nrow(ref)
K = 4
d = 500
A = matrix(0,nrow=G,ncol=G)
for(i in 1:G){
  for(j in i:min(i+d,G)){
    A[i,j] = max(1-abs(i-j)/d,0)
  }
}
A = A+t(A) - diag(G)
library(Matrix)
A = Matrix(A,sparse = TRUE)




alpha.cors = c(0,0.5)
cases = c("null","all_diff")
nbs = c(50,100)
dirichlet.scale = c(5,10)

set.seed(12345)
for(case in cases){

  if(case=='null'){
    p1 = c(0.5,0.3,0.1,0.1)
    p2 = c(0.5,0.3,0.1,0.1)
  }else if(case=='all_diff'){
    p1 = c(0.4,0.3,0.2,0.1)
    p2 = c(0.1,0.1,0.3,0.5)
  }

  for(nb in nbs){
    for(aa in dirichlet.scale){


      b1 = t(rdirichlet(nb/2,p1*aa))
      b2 = t(rdirichlet(nb/2,p2*aa))
      b = cbind(b1,b2)

      for(alpha.cor in alpha.cors){

        if(alpha.cor==0){
          est_cor = FALSE
          cor.status = 'trueR'
        }else{
          est_cor=TRUE
          cor.status = paste('cor0',alpha.cor*10,sep = '')
        }

        print(paste('Running:',case,'nb=',nb,'cor:',cor.status,'aa=',aa))

        simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1,alpha.cor = alpha.cor,est_cor=est_cor)
        saveRDS(simu_out,file = paste('output/manuscript/simulation/simulation_',nb,'bulk_500genecor_',cor.status,'_',case,'_dirichlet',aa,'.rds',sep=''))

      }
    }
  }
}
#######################################




d = 300
A = matrix(0,nrow=G,ncol=G)

for(i in 1:G){
  for(j in i:min(i+d,G)){
    A[i,j] = max(1-abs(i-j)/d,0)
  }
}
A = A+t(A) - diag(G)
library(Matrix)
A = Matrix(A,sparse = TRUE)
set.seed(12345)
simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1)
saveRDS(simu_out,file = 'output/manuscript/simulation/simulation_10bulk_300genecor_fdr05.rds')

set.seed(12345)
simu_out = simu_study(ref,b,R=NULL,sigma2,printevery = 1)
saveRDS(simu_out,file = 'output/manuscript/simulation/simulation_10bulk_0genecor_fdr05.rds')


# for(t in 1:7){
#   temp = simu_out[[t]]
#   temp_out = array(dim=c(4,10,100))
#   for(r in 1:100){
#     temp_out[,,r] = matrix(c(temp[,,r]),nrow=4)
#   }
#   simu_out[[t]] = temp_out
# }
# saveRDS(simu_out,file = 'output/manuscript/simulation1.rds')

#

# diff_hat_se = matrix(nrow=100,ncol=4)
# diff_hat_se_cor = matrix(nrow=100,ncol=4)
# diff_hat_weight_se_cor = matrix(nrow=100,ncol=4)
# for(i in 1:100){
#   diff_hat_se[i,] = two_group_test(simu_out$all_fit[[i]]$fit.err,c(1,1,1,1,1,2,2,2,2,2))$diff_se
#   diff_hat_se_cor[i,] = two_group_test(simu_out$all_fit[[i]]$fit.err.cor,c(1,1,1,1,1,2,2,2,2,2))$diff_se
#   diff_hat_weight_se_cor[i,] = two_group_test(simu_out$all_fit[[i]]$fit.err.cor.weight,c(1,1,1,1,1,2,2,2,2,2))$diff_se
# }
# simu_out$diff_hat_se = diff_hat_se
# simu_out$diff_hat_se_cor = diff_hat_se_cor
# simu_out$diff_hat_weight_se_cor = diff_hat_weight_se_cor

## Look at rmse
rmse = function(x,y){sqrt(mean((x-y)^2))}

get_rmse = function(p_hat,b){
  K = dim(p_hat)[1]
  nb = dim(p_hat)[2]
  nreps = dim(p_hat)[3]
  rmses = c()
  for(i in 1:nb){
    err = c()
    for(j in 1:nreps){
      err[j] = sum((p_hat[,i,j]-b[,i])^2)
    }
    rmses[i] = sqrt(mean(err))
  }
  names(rmses) = paste('bulk',1:nb)
  rmses

}

get_rmse(simu_out$p_hat_ols,b)
get_rmse(simu_out$p_hat,b)
get_rmse(simu_out$p_hat_weight,b)


# Look at coverage

get_coverage_p = function(p_hat,p_hat_se,b){

  K = dim(z)[1]
  nb = dim(z)[2]
  z = array(dim = dim(p_hat))
  for(i in 1:dim(z)[3]){
    z[,,i] = (p_hat[,,i]-b)/p_hat_se[,,i]
  }
  crg = apply(z,c(1,2),function(z){round(mean(abs(z)<1.96,na.rm=T),3)})
  rownames(crg) = paste('cell',1:K)
  colnames(crg) = paste('bulk',1:nb)
  crg
}

get_coverage_p(simu_out$p_hat_ols,simu_out$p_hat_ols_se,b)
get_coverage_p(simu_out$p_hat,simu_out$p_hat_se,b)
get_coverage_p(simu_out$p_hat,simu_out$p_hat_se_cor,b)
get_coverage_p(simu_out$p_hat_weight,simu_out$p_hat_weight_se_cor,b)

get_power_diff = function(diff_hat,diff_hat_se){
  colMeans(abs(diff_hat/diff_hat_se)>1.96,na.rm=TRUE)
}

get_power_diff(simu_out$diff_hat_ols,simu_out$diff_hat_ols_se)
get_power_diff(simu_out$diff_hat,simu_out$diff_hat_se)
get_power_diff(simu_out$diff_hat,simu_out$diff_hat_se_cor)
get_power_diff(simu_out$diff_hat_weight,simu_out$diff_hat_weight_se_cor)

