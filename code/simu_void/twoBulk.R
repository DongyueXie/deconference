############## test ################
source('code/deconference.R')
source('code/simu_func1.R')
source('code/simu_func_multi.R')


#### case 1 : bulk reference data ####

G = 1000
K = 4
b = c(0.05,0.05,0.1,0.8)
b2 = c(0.1,0.2,0.3,0.4)
library(gtools)
set.seed(12345)
ref = t(rdirichlet(K,rep(1,G)))
G_list = round(seq(50,nrow(ref),length.out = 10))

results = list()
results_null = list()

for(i in 1:length(G_list)){

  if(i%%2==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}

  results[[i]] = simu_study(ref,G_list[i],b,
                            ref_type='bulk',
                            printevery = 1e5,
                            b2=b2)

  results_null[[i]] = simu_study(ref,G_list[i],b,
                                 ref_type='bulk',
                                 printevery = 1e5,
                                 b2=b)

}

save(results,file=paste('output/twoBulk_bulkref_simu_G',G,'_K',K,'_refls',30,'.RData',sep = ''))
save(results_null,file=paste('output/twoBulk_null_bulkref_simu_G',G,'_K',K,'_refls',30,'.RData',sep = ''))

#### case 2 : single cell reference data ####

G = 1000
K = 4
b = c(0.05,0.05,0.1,0.8)
b2 = c(0.1,0.2,0.3,0.4)
library(gtools)
set.seed(12345)
ref = t(rdirichlet(K,rep(1,G)))
G_list = round(seq(50,nrow(ref),length.out = 10))

results = list()
results_null = list()

for(i in 1:length(G_list)){

  if(i%%2==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}

  results[[i]] = simu_study(ref,G_list[i],b,
                            ref_type='sc',
                            x_estimator = 'aggregate',
                            printevery = 1e5,
                            b2=b2)
  results_null[[i]] = simu_study(ref,G_list[i],b,
                            ref_type='sc',
                            x_estimator = 'aggregate',
                            printevery = 1e5,
                            b2=b)

}

save(results,file=paste('output/twoBulk_scref_simu_G',G,'_K',K,'.RData',sep = ''))
save(results_null,file=paste('output/twoBulk_null_scref_simu_G',G,'_K',K,'.RData',sep = ''))


#### case 3 : multiple single cell reference data ####

G = 1000
K = 4
b = c(0.05,0.05,0.1,0.8)
b2 = c(0.1,0.2,0.3,0.4)
library(gtools)
set.seed(12345)
ref = t(rdirichlet(K,rep(1,G)))
G_list = round(seq(50,nrow(ref),length.out = 10))

results = list()
results_null = list()

for(i in 1:length(G_list)){

  if(i%%2==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}

  results[[i]] = simu_study_multiInd(ref,G_list[i],b,
                            ref_type='sc',
                            x_estimator = 'aggregate',
                            printevery = 1e5,
                            b2=b2)

  results_null[[i]] = simu_study_multiInd(ref,G_list[i],b,
                                     ref_type='sc',
                                     x_estimator = 'aggregate',
                                     printevery = 1e5,
                                     b2=b)

}

save(results,file=paste('output/twoBulk_multiscref_simu_G',G,'_K',K,'.RData',sep = ''))

save(results_null,file=paste('output/twoBulk_null_multiscref_simu_G',G,'_K',K,'.RData',sep = ''))


