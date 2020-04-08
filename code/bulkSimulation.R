source('code/gls.R')

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



#################################################################
#################################################################
### choose gene method 1: remove genes, using my original methods
#################################################################
#################################################################

rm.idx1 = which(((rowSums(CD14!=0)<=10)|(rowSums(CD14)>=quantile(rowSums(CD14),0.95))))
rm.idx2 = which(((rowSums(CD4!=0)<=10)|(rowSums(CD4)>=quantile(rowSums(CD4),0.95))))
rm.idx3 = which(((rowSums(CD8!=0)<=10)|(rowSums(CD8)>=quantile(rowSums(CD8),0.95))))
rm.idx4 = which(((rowSums(MB!=0)<=10)|(rowSums(MB)>=quantile(rowSums(MB),0.95))))
rm.idx = unique(c(rm.idx1,rm.idx2,rm.idx3,rm.idx4))

s = c(sum(CD14[-rm.idx,])/ncol(CD14[-rm.idx,]),
      sum(CD4[-rm.idx,])/ncol(CD4[-rm.idx,]),
      sum(CD8[-rm.idx,])/ncol(CD8[-rm.idx,]),
      sum(MB[-rm.idx,])/ncol(MB[-rm.idx,]))

####
ref = cbind(rowSums(CD14[-rm.idx,]),rowSums(CD4[-rm.idx,]),rowSums(CD8[-rm.idx,]),rowSums(MB[-rm.idx,]))
####


set.seed(12345)

G_list = round(seq(50,nrow(ref),length.out = 100))

ref_lib_size_list = c(0.1,0.3,0.5,1,3,5,10,50)

for(r in ref_lib_size_list){

  print(paste('Running ref libary size:',r))

  results = list()

  for(i in 1:length(G_list)){

    if(i%%10==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}

    results[[i]] = simu_study(ref,s,G_list[i],(b/s)/sum(b/s),
                              bulk_lib_size,r)

  }

  if(r<1){r=paste('0',r*10,sep='')}
  save(results,file=paste('output/bulkref_pbmc_fewgene_refls',r,'_bulkls',bulk_lib_size,'.RData',sep = ''))
}




#################################################
## choose gene method 0:
## genes expressed at extremely high levels in the reference dataset will dominate the inference results.
## Since highly expressed genes have large variance, our inferences become very sensitive to these outliers.
## We therefore remove them to acquire more robust estimations.
## The total removed genes were the union of the top 1 %.

## see https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1028-7#Sec9
################################################

ref = cbind(rowSums(CD14),rowSums(CD4),rowSums(CD8),rowSums(MB))
rm.idx = unique(unlist(apply(ref,2,function(x) which(x>=quantile(x,0.99)))))


s = c(sum(CD14[-rm.idx,])/ncol(CD14[-rm.idx,]),
      sum(CD4[-rm.idx,])/ncol(CD4[-rm.idx,]),
      sum(CD8[-rm.idx,])/ncol(CD8[-rm.idx,]),
      sum(MB[-rm.idx,])/ncol(MB[-rm.idx,]))

####
ref = cbind(rowSums(CD14[-rm.idx,]),rowSums(CD4[-rm.idx,]),rowSums(CD8[-rm.idx,]),rowSums(MB[-rm.idx,]))
####


set.seed(12345)

G_list = round(seq(50,nrow(ref)/2,length.out = 100))

ref_lib_size_list = c(0.1,0.3,0.5,1,3,5,10,50)

for(r in ref_lib_size_list){

  print(paste('Running ref libary size:',r))

  results = list()

  for(i in 1:length(G_list)){

    if(i%%10==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}

    results[[i]] = simu_study(ref,s,G_list[i],(b/s)/sum(b/s),
                              bulk_lib_size,r)

  }

  if(r<1){r=paste('0',r*10,sep='')}
  save(results,file=paste('output/bulkref_pbmc_allgene_refls',r,'_bulkls',bulk_lib_size,'.RData',sep = ''))
}





##################################################################
##############choose gene method 2: use only marker gene###########

## attention: need to select marker genes with existing deconv method?

###################################################################

# ref = cbind(rowSums(CD14),rowSums(CD4),rowSums(CD8),rowSums(MB))
#
# # marker gene: only expresses in one cell type.
#
# marker_gene_idx = which(rowSums(ref!=0)==1)
#
# s = c(sum(CD14[marker_gene_idx,])/ncol(CD14[marker_gene_idx,]),
#       sum(CD4[marker_gene_idx,])/ncol(CD4[marker_gene_idx,]),
#       sum(CD8[marker_gene_idx,])/ncol(CD8[marker_gene_idx,]),
#       sum(MB[marker_gene_idx,])/ncol(MB[marker_gene_idx,]))
#
#
# G_list = round(seq(50,nrow(ref),length.out = 50))
#
#
#
# set.seed(12345)
#
# results = list()
#
# for(i in 1:length(G_list)){
#
#   print(i)
#
#   results[[i]] = simu_study(ref,s,G_list[i],(b/s)/sum(b/s),
#                             bulk_lib_size,ref_lib_size)
#
# }
#
#
#
#
#
#
# ###############################################################
# ############### Segerstolpe et al. (2016) #####################
# ###############         non-UMI           #####################
# ###############################################################
#
# EMTAB.eset = readRDS('data/deconv/EMTABesethealthy.rds')
# cell_types = c('acinar','alpha','beta','ductal')
#
# acinar = exprs(EMTAB.eset)[,which(EMTAB.eset$cellType=='acinar')]
# alpha = exprs(EMTAB.eset)[,which(EMTAB.eset$cellType=='alpha')]
# beta = exprs(EMTAB.eset)[,which(EMTAB.eset$cellType=='beta')]
# ductal = exprs(EMTAB.eset)[,which(EMTAB.eset$cellType=='ductal')]
#
# s = c(sum(acinar)/ncol(acinar),
#       sum(alpha)/ncol(alpha),
#       sum(beta)/ncol(beta),
#       sum(ductal)/ncol(ductal))
#
#
# ref = cbind(rowSums(acinar),rowSums(alpha),rowSums(beta),rowSums(ductal))
#
# ref = select_gene(ref)
#
#
# G_list = round(seq(50,nrow(ref),length.out = 50))
#
#
# set.seed(12345)
# b = c(1,2,3,4)
# Nr = rep(30,4)
#
# results_nonUMI = list()
#
# for(i in 1:length(G_list)){
#
#   print(i)
#
#   results_nonUMI[[i]] = simu_study(ref,s,G_list[i],b,Nr)
#
# }
