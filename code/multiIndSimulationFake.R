source('code/deconference.R')


# 1. n_indi = 10; est_sigma2 / do not est_sigma2
# 2. n_indi = 50; est_sigma2 / do not est_sigma2

G = 300
K = 4
b = 1:K
b = b/sum(b)
library(gtools)
set.seed(12345)
ref = t(rdirichlet(K,rep(1,G)))
bulk_lib_size = 50
G_list = round(seq(50,nrow(ref),length.out = 10))



n_indi = 10

# est_sigma2

results = list()
for(i in 1:length(G_list)){

  results[[i]] = simu_study_multiInd(ref,G_list[i],(b)/sum(b),
                                     ref_type='sc',
                                     n_indi = n_indi,
                                     tau2known=FALSE,addw=TRUE,
                                     bulk_lib_size=bulk_lib_size,
                                     sc_lib_size=0.2,
                                     est_sigma2 = TRUE,
                                     nk=100,snr=100)

  save(results,file=paste('scref_multiInd',n_indi,'_simu_G',G,'_K',K,'_scls02','_scN100','_bulkls',bulk_lib_size,'_addw_tauunknown','.RData',sep = ''))
  if(i%%1==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}
}


# est_sigma2

results = list()
for(i in 1:length(G_list)){

  results[[i]] = simu_study_multiInd(ref,G_list[i],(b)/sum(b),
                                     ref_type='sc',
                                     n_indi = n_indi,
                                     tau2known=FALSE,addw=TRUE,
                                     bulk_lib_size=bulk_lib_size,
                                     sc_lib_size=0.2,
                                     est_sigma2 = FALSE,
                                     nk=100,snr=100)

  save(results,file=paste('scref_multiInd',n_indi,'_simu_G',G,'_K',K,'_scls02','_scN100','_bulkls',bulk_lib_size,'_addw_tauunknown_Nosigma2','.RData',sep = ''))
  if(i%%1==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}
}

#########################################
#########################################




n_indi = 100

# est_sigma2

results = list()
for(i in 1:length(G_list)){

  results[[i]] = simu_study_multiInd(ref,G_list[i],(b)/sum(b),
                                     ref_type='sc',
                                     n_indi = n_indi,
                                     tau2known=FALSE,addw=TRUE,
                                     bulk_lib_size=bulk_lib_size,
                                     sc_lib_size=0.2,
                                     est_sigma2 = TRUE,
                                     nk=100,snr=100)

  save(results,file=paste('output/scref_multiInd',n_indi,'_simu_G',G,'_K',K,'_scls02','_scN100','_bulkls',bulk_lib_size,'_addw_tauunknown','.RData',sep = ''))
  if(i%%1==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}
}


# est_sigma2

results = list()
for(i in 1:length(G_list)){

  results[[i]] = simu_study_multiInd(ref,G_list[i],(b)/sum(b),
                                     ref_type='sc',
                                     n_indi = n_indi,
                                     tau2known=FALSE,addw=TRUE,
                                     bulk_lib_size=bulk_lib_size,
                                     sc_lib_size=0.2,
                                     est_sigma2 = FALSE,
                                     nk=100,snr=100)

  save(results,file=paste('scref_multiInd',n_indi,'_simu_G',G,'_K',K,'_scls02','_scN100','_bulkls',bulk_lib_size,'_addw_tauunknown_Nosigma2','.RData',sep = ''))
  if(i%%1==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}
}






