source('code/deconference.R')

set.seed(12345)

G = 1000
K = 4
b = 1:K
b = b/sum(b)
library(gtools)
ref = t(rdirichlet(K,rep(1,G)))
bulk_lib_size = 50

G_list = round(seq(50,nrow(ref),length.out = 50))

sc_lib_size_list = c(0.1,0.2,0.5)
sc_number_list = c(10,30,300)

for(n in sc_number_list){
  for(r in sc_lib_size_list){

    print(paste('Running number of cell:',n, ' and sc libary size:',r))

    results = list()

    for(i in 1:length(G_list)){

      if(i%%10==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}

      results[[i]] = simu_study(ref,G_list[i],(b)/sum(b),
                                ref_type='sc',
                                tau2known=FALSE,addw=TRUE,
                                bulk_lib_size=bulk_lib_size,
                                sc_lib_size=r,
                                nk=n)

    }

    if(r<1){r=paste('0',r*10,sep='')}
    save(results,file=paste('output/scref_simu_G',G,'_K',K,'_scls',r,'_scN',n,'_bulkls',bulk_lib_size,'_addw_tauunknown','.RData',sep = ''))
  }
}



########## correction ###################


# set.seed(12345)
# sc_lib_size_list = c(0.1,0.2)
# sc_number_list = c(10)
#
# for(n in sc_number_list){
#   for(r in sc_lib_size_list){
#
#     print(paste('Running number of cell:',n, ' and sc libary size:',r))
#
#     results = list()
#
#     for(i in 1:length(G_list)){
#
#       if(i%%10==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}
#
#       results[[i]] = simu_study_sc(ref,G_list[i],(b)/sum(b),
#                                    bulk_lib_size,r,n,correction = T)
#
#     }
#
#     if(r<1){r=paste('0',r*10,sep='')}
#     save(results,file=paste('output/scref_simu_correction_G',G,'_K',K,'_scls',r,'_scN',n,'_bulkls',bulk_lib_size,'.RData',sep = ''))
#   }
# }
