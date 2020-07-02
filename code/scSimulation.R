source('code/deconference.R')

################################################
############### PBMC DATA ######################
###############   UMI     ######################
################################################


load("~/misc/data/scde/scCD14.RData")
load("~/misc/data/scde/scCD4.RData")
load("~/misc/data/scde/scCD8.RData")
load("~/misc/data/scde/scMB.RData")

CD14 = as.matrix(CD14)
CD4 = as.matrix(CD4)
CD8 = as.matrix(CD8)
MB = as.matrix(MB)

ref = cbind(rowSums(CD14),rowSums(CD4),rowSums(CD8),rowSums(MB))
rm.idx = which(rowSums(ref)==0)

CD14 = CD14[-rm.idx,]
CD4 = CD4[-rm.idx,]
CD8 = CD8[-rm.idx,]
MB = MB[-rm.idx,]

b = c(1,2,3,4)

bulk_lib_size = 50




# #################################################
# ## choose gene method 0:
# ## genes expressed at extremely high levels in the reference dataset will dominate the inference results.
# ## Since highly expressed genes have large variance, our inferences become very sensitive to these outliers.
# ## We therefore remove them to acquire more robust estimations.
# ## The total removed genes were the union of the top 1 %.
#
# ## see https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1028-7#Sec9
# ################################################
#

ref = cbind(rowSums(CD14),rowSums(CD4),rowSums(CD8),rowSums(MB))
rm.idx = unique(unlist(apply(ref,2,function(x) which(x>=quantile(x,0.99)))))


s = c(sum(CD14[-rm.idx,])/ncol(CD14[-rm.idx,]),
      sum(CD4[-rm.idx,])/ncol(CD4[-rm.idx,]),
      sum(CD8[-rm.idx,])/ncol(CD8[-rm.idx,]),
      sum(MB[-rm.idx,])/ncol(MB[-rm.idx,]))

####
ref = cbind(rowSums(CD14[-rm.idx,]),rowSums(CD4[-rm.idx,]),rowSums(CD8[-rm.idx,]),rowSums(MB[-rm.idx,]))
####

ref = apply(ref,2,function(z){z/sum(z)})

# set.seed(12345)
#
# G_list = round(seq(50,nrow(ref)/2,length.out = 100))
#
# sc_lib_size_list = c(0.5,1,3,5,10,50)
# sc_number_list = c(10,30,50)
#
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
#                                    bulk_lib_size,r,n)
#
#     }
#
#     if(r<1){r=paste('0',r*10,sep='')}
#     save(results,file=paste('output/scref_pbmc_allgene_K',K,'_scls',r,'_scN',n,'_bulkls',bulk_lib_size,'.RData',sep = ''))
#   }
# }


## correction

set.seed(12345)

G_list = round(seq(50,nrow(ref)/2,length.out = 100))

sc_lib_size_list = c(0.5)
sc_number_list = c(10,30)


for(n in sc_number_list){
  for(r in sc_lib_size_list){

    print(paste('Running number of cell:',n, ' and sc libary size:',r))

    results = list()

    for(i in 1:length(G_list)){

      if(i%%10==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}

      results[[i]] = simu_study_sc(ref,G_list[i],(b)/sum(b),
                                   bulk_lib_size,r,n,correction = T)

    }

    if(r<1){r=paste('0',r*10,sep='')}
    save(results,file=paste('output/scref_pbmc_allgene_correction_K',K,'_scls',r,'_scN',n,'_bulkls',bulk_lib_size,'.RData',sep = ''))
  }
}





##################################################################
##############choose gene method 2: use only 'marker' gene###########

## attention: need to select marker genes with existing deconv method?

###################################################################

ref = cbind(rowSums(CD14),rowSums(CD4),rowSums(CD8),rowSums(MB))

# marker gene: only expresses in one cell type.

marker_gene_idx = which(rowSums(ref!=0)==1)

s = c(sum(CD14[marker_gene_idx,])/ncol(CD14[marker_gene_idx,]),
      sum(CD4[marker_gene_idx,])/ncol(CD4[marker_gene_idx,]),
      sum(CD8[marker_gene_idx,])/ncol(CD8[marker_gene_idx,]),
      sum(MB[marker_gene_idx,])/ncol(MB[marker_gene_idx,]))

ref = ref[marker_gene_idx,]


set.seed(12345)

G_list = round(seq(50,nrow(ref),length.out = 20))

sc_lib_size_list = c(0.5,1,3,5,10,50)
sc_number_list = c(10,30,50)


for(n in sc_number_list){
  for(r in sc_lib_size_list){

    print(paste('Running number of cell:',n, ' and sc libary size:',r))

    results = list()

    for(i in 1:length(G_list)){

      if(i%%10==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}

      results[[i]] = simu_study_sc(ref,G_list[i],(b)/sum(b),
                                   bulk_lib_size,r,n)

    }

    if(r<1){r=paste('0',r*10,sep='')}
    save(results,file=paste('output/scref_pbmc_markergene_K',K,'_scls',r,'_scN',n,'_bulkls',bulk_lib_size,'.RData',sep = ''))
  }
}


## correction

set.seed(12345)

sc_lib_size_list = c(0.5,1)
sc_number_list = c(10,30)


for(n in sc_number_list){
  for(r in sc_lib_size_list){

    print(paste('Running number of cell:',n, ' and sc libary size:',r))

    results = list()

    for(i in 1:length(G_list)){

      if(i%%10==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}

      results[[i]] = simu_study_sc(ref,G_list[i],(b)/sum(b),
                                   bulk_lib_size,r,n,correction = T)

    }

    if(r<1){r=paste('0',r*10,sep='')}
    save(results,file=paste('output/scref_pbmc_markergene_correction_K',K,'_scls',r,'_scN',n,'_bulkls',bulk_lib_size,'.RData',sep = ''))
  }
}

