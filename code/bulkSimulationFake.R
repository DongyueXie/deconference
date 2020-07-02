source('code/deconference.R')

set.seed(12345)

G = 1000
K = 4
b = 1:K
b = b/sum(b)
library(gtools)
ref = t(rdirichlet(K,rep(1,G)))
s=rep(1,K)
bulk_lib_size = 50

G_list = round(seq(50,nrow(ref),length.out = 50))

ref_lib_size_list = c(1,3,5,10)

for(r in ref_lib_size_list){

  print(paste('Running ref libary size:',r))

  results = list()

  for(i in 1:length(G_list)){

    if(i%%10==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}

    results[[i]] = simu_study(ref,G_list[i],(b/s)/sum(b/s),
                              ref_type='bulk',
                              addw=TRUE,
                              bulk_lib_size=bulk_lib_size,
                              ref_lib_size=r)

  }

  if(r<1){r=paste('0',r*10,sep='')}
  save(results,file=paste('output/bulkref_simu_addw_G',G,'_K',K,'_refls',r,'_bulkls',bulk_lib_size,'.RData',sep = ''))
}


########## correction ###################
set.seed(12345)
ref_lib_size_list = c(1,3)

for(r in ref_lib_size_list){

  print(paste('Running ref libary size:',r))

  results = list()

  for(i in 1:length(G_list)){

    if(i%%10==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}

    results[[i]] = simu_study(ref,G_list[i],(b/s)/sum(b/s),
                              ref_type='bulk',
                              addw=TRUE,
                              bulk_lib_size=bulk_lib_size,
                              ref_lib_size=r,
                              correction = TRUE)

  }

  if(r<1){r=paste('0',r*10,sep='')}
  save(results,file=paste('output/bulkref_simu_correction_addw_G',G,'_K',K,'_refls',r,'_bulkls',bulk_lib_size,'.RData',sep = ''))
}
