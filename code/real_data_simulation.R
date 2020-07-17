############## test ################
source('code/deconference_main.R')
source('code/simu_func1.R')



######################################
######################################
# start with LM6 matrix from cibersort

# LM22 = read.table('data/cibersort/LM22.txt',skip = 1,row.names = 1)
LM6 = read.table('data/cibersort/signature_rnaseq_geo60424_LM6.txt',header = TRUE,sep='\t',row.names = 1)
LM6_type = c("B","CD8","CD4","NK","Monocytes","Neutrophils")

ref = apply(LM6,2,function(z){z/sum(z)})

#### case 1 : bulk reference data ####

G = nrow(ref)
K = ncol(ref)
b = c(0.05,0.05,0.05,0.05,0.1,0.7)
b2 = c(0.05,0.06,0.08,0.1,0.2,0.51)


set.seed(12345)

results = list()

# for(i in 1:length(G_list)){

  # if(i%%2==0){print(sprintf("done %d (out of %d)",i,length(G_list)))}

  results[[1]] = simu_study(ref,G,b,
                            ref_type='bulk',
                            printevery = 10,
                            b2=b2)

  results[[2]] = simu_study(ref,G,b,
                       ref_type='sc',
                       printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1/100,
                       b2=b2,weight = 'equal',hc.type = 'hc0')

  results[[3]] = simu_study(ref,G,b,
                            ref_type='sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1/100,
                            b2=b2,weight = 'equal',hc.type = 'hc3')

  results[[4]] = simu_study(ref,G,b,
                            ref_type='sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1/100,
                            b2=b2,weight = 'default',hc.type = 'hc0')

  results[[5]] = simu_study(ref,G,b,
                            ref_type='sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1/100,
                            b2=b2,weight = 'default',hc.type = 'hc3')


  save(results,file = 'output/LM6_simu.RData')

  results[[6]] = simu_study(ref,G,b,
                            ref_type='sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1,
                            b2=b2,weight = 'equal',hc.type = 'hc0')

  results[[7]] = simu_study(ref,G,b,
                            ref_type='sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1,
                            b2=b2,weight = 'equal',hc.type = 'hc3')

  results[[8]] = simu_study(ref,G,b,
                            ref_type='sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1,
                            b2=b2,weight = 'default',hc.type = 'hc0')

  results[[9]] = simu_study(ref,G,b,
                            ref_type='sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1,
                            b2=b2,weight = 'default',hc.type = 'hc3')

  results[[10]] = simu_study(ref,G,b,
                            ref_type='sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1/G,
                            b2=b2,weight = 'equal',hc.type = 'hc0')

  save(results,file = 'output/LM6_simu.RData')

  results[[11]] = simu_study(ref,G,b,
                            ref_type='sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1/G,
                            b2=b2,weight = 'equal',hc.type = 'hc3')

  results[[12]] = simu_study(ref,G,b,
                            ref_type='sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1/G,
                            b2=b2,weight = 'default',hc.type = 'hc0')

  results[[13]] = simu_study(ref,G,b,
                            ref_type='sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1/G,
                            b2=b2,weight = 'default',hc.type = 'hc3')


  save(results,file = 'output/LM6_simu.RData')

  results[[14]] = simu_study(ref,G,b,
                            ref_type='multi_sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1/G,
                            b2=b2,weight = 'equal',hc.type = 'hc3')

  results[[15]] = simu_study(ref,G,b,
                            ref_type='multi_sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 10,mean_to_var_tau = 1/G,
                            b2=b2,weight = 'default',hc.type = 'hc3')

  results[[16]] = simu_study(ref,G,b,
                            ref_type='multi_sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 1,mean_to_var_tau = 1/G,
                            b2=b2,weight = 'equal',hc.type = 'hc3')

  results[[17]] = simu_study(ref,G,b,
                            ref_type='multi_sc',
                            printevery = 10,same_indi = F,mean_to_var_sigma = 1,mean_to_var_tau = 1/G,
                            b2=b2,weight = 'default',hc.type = 'hc3')

  save(results,file = 'output/LM6_simu.RData')

  results[[18]] = simu_study(ref,G,b,
                             ref_type='multi_sc',
                             printevery = 10,same_indi = F,mean_to_var_sigma = 1/100,mean_to_var_tau = 1/G,
                             b2=b2,weight = 'equal',hc.type = 'hc3')

  results[[19]] = simu_study(ref,G,b,
                             ref_type='multi_sc',
                             printevery = 10,same_indi = F,mean_to_var_sigma = 1/100,mean_to_var_tau = 1/G,
                             b2=b2,weight = 'default',hc.type = 'hc3')

  results[[20]] = simu_study(ref,G,b,
                             ref_type='multi_sc',
                             printevery = 10,same_indi = F,mean_to_var_sigma = 1/G,mean_to_var_tau = 1/G,
                             b2=b2,weight = 'equal',hc.type = 'hc3')

  results[[21]] = simu_study(ref,G,b,
                             ref_type='multi_sc',
                             printevery = 10,same_indi = F,mean_to_var_sigma = 1/G,mean_to_var_tau = 1/G,
                             b2=b2,weight = 'default',hc.type = 'hc3')

  save(results,file = 'output/LM6_simu.RData')

# }



######################################
######################################
#

results_nsclc = list()

nsclc = read.table('data/cibersort/Fig2ab-NSCLC_PBMCs/Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt',row.names = 1,skip=1)
cells = read.table('data/cibersort/Fig2ab-NSCLC_PBMCs/Fig2ab-NSCLC_PBMCs_scRNAseq_refsample.txt',row.names = 1,nrows = 1,sep="\t")

cell = c()
for(i in 1:length(cells)){
  cell[i] = (cells[1,i])
}

cell_type = unique(cell)


ref_cibersort = read.table('data/cibersort/Fig2ab-NSCLC_PBMCs/Fig2ab-NSCLC_PBMCs_scRNAseq_sigmatrix.txt',row.names = 1,skip=1)

marker_gene = rownames(ref_cibersort)

# gene_idx = match(marker_gene,rownames(mat))
#
#
# ## remove genes
# # rm.gene = which(rowSums(nsclc>0)<=10)
# # mat = nsclc[-rm.gene,]
#
# mat = nsclc[gene_idx,]
# X = matrix(nrow=nrow(mat),ncol=length(cell_type))
# rownames(X) = rownames(mat)
# colnames(X) = cell_type
#
# for(i in 1:length(cell_type)){
#   cell_idx = which(cell == cell_type[i])
#   cell_mat = mat[,cell_idx] + 1
#   X[,i] = rowSums(cell_mat)/sum(cell_mat)
# }


mat = nsclc
X = matrix(nrow=nrow(mat),ncol=length(cell_type))
rownames(X) = rownames(mat)
colnames(X) = cell_type

for(i in 1:length(cell_type)){
  cell_idx = which(cell == cell_type[i])
  cell_mat = mat[,cell_idx] + 1
  X[,i] = rowSums(cell_mat)/sum(cell_mat)
}

# rm.gene = which(rowSums(X==0)>0)
# X = X[-rm.gene,]



#### case 1 : bulk reference data ####

G = nrow(X)
K = ncol(X)
b = c(0.05,0.05,0.05,0.05,0.1,0.7)
b2 = c(0.05,0.1,0.15,0.2,0.3,0.2)

set.seed(12345)


results = simu_study(X,G,b,
                     ref_type='bulk',
                     printevery = 10,
                     marker_gene=NULL,
                     b2=b2,bulk_lib_size = 100,ref_lib_size = 100,same_indi = TRUE,weight = 'equal')

results_nsclc[[1]] = results
save(results_nsclc,file="output/results_nsclc.RData")

results = simu_study(X,G,b,
                     ref_type='bulk',
                     printevery = 10,
                     marker_gene=NULL,
                     b2=b2,bulk_lib_size = 100,ref_lib_size = 100,same_indi = TRUE,weight = 'default')

results_nsclc[[2]] = results
save(results_nsclc,file="output/results_nsclc.RData")


# use marker genes selected by cibersort
results = simu_study(X,G,b,
                     ref_type='bulk',
                     printevery = 10,
                     marker_gene=marker_gene,
                     b2=b2,bulk_lib_size = 100,ref_lib_size = 100,same_indi = TRUE,weight = 'default')

results_nsclc[[5]] = results
save(results_nsclc,file="output/results_nsclc.RData")

#### case 2 : sc reference data ####

results = simu_study(X,G,b,
                     ref_type='sc',
                     printevery = 10,
                     marker_gene=NULL,
                     b2=b2,same_indi = TRUE,weight = 'default')

results_nsclc[[3]] = results
save(results_nsclc,file="output/results_nsclc.RData")

# use marker genes selected by cibersort
results = simu_study(X,G,b,
                     ref_type='sc',
                     printevery = 10,
                     marker_gene=marker_gene,
                     b2=b2,same_indi = TRUE,weight = 'default')

results_nsclc[[6]] = results
save(results_nsclc,file="output/results_nsclc.RData")

#### case 3 : multi sc reference data ####

results = simu_study(X,G,b,
                     ref_type='multi_sc',
                     printevery = 10,
                     marker_gene=NULL,
                     b2=b2,same_indi = F,weight = 'default')

results_nsclc[[4]] = results
save(results_nsclc,file="output/results_nsclc.RData")

results = simu_study(X,G,b,
                     ref_type='multi_sc',
                     printevery = 10,
                     marker_gene=marker_gene,
                     b2=b2,same_indi = F,weight = 'default')

results_nsclc[[7]] = results
save(results_nsclc,file="output/results_nsclc.RData")

######################################
######################################
#

library(Seurat)
library(Matrix)
pbmc = pbmc3k.SeuratData::pbmc3k.final
cell_types = levels(pbmc@active.ident)[1:4]
G = dim(pbmc)[1]
K = length(cell_types)
X = matrix(nrow=G,ncol=K)
for(k in 1:K){
  cell_data = pbmc@assays$RNA@counts[,which(pbmc@active.ident==cell_types[k])]
  X[,k] = rowSums(cell_data)/sum(cell_data)
}
X = X[-which(rowSums(X)==0),]

G = nrow(X)
K = ncol(X)
b = c(0.1,0.2,0.3,0.4)
b2 = c(0.05,0.05,0.1,0.8)

results_pbmc = list()

set.seed(12345)

results = simu_study(X,G,b,
                     ref_type='bulk',
                     printevery = 10,
                     marker_gene=NULL,
                     weight = 'default',
                     same_indi = T,
                     b2=b2)

results_pbmc[[1]] = results
save(results_pbmc,file="output/results_pbmc.RData")

results = simu_study(X,G,b,
                     ref_type='sc',
                     printevery = 10,
                     marker_gene=NULL,
                     weight = 'default',
                     same_indi = T,
                     b2=b2)

results_pbmc[[2]] = results
save(results_pbmc,file="output/results_pbmc.RData")

results = simu_study(X,G,b,
                     ref_type='multi_sc',
                     printevery = 10,
                     marker_gene=NULL,
                     weight = 'default',
                     same_indi = F,
                     b2=b2)

results_pbmc[[3]] = results
save(results_pbmc,file="output/results_pbmc.RData")



