source('code/deconference_main.R')
source('code/utils.R')
source('code/wols.R')
source('code/simulation/neuron/simu_neuron.R')

# remove genes

gene_name = readRDS('data/neuron/gene_name_20293.rds')
indis_ref = readRDS('data/neuron/indis_ref.rds')
indis_ref = indis_ref[match(gene_name,dimnames(indis_ref)[[1]]),,]
dim(indis_ref)

# remove individuals that has no or too few certain cell types
cell_ann = readRDS('data/neuron/cell_ann.rds')
indi_by_cell = table(cell_ann$donor_id,cell_ann$celltype)
indi_by_cell[,6] = indi_by_cell[,6] + indi_by_cell[,7]
indi_by_cell = indi_by_cell[,1:6]
colnames(indi_by_cell)[6] = 'U_Neur'

indi_to_use = names(which(rowSums(indi_by_cell>10)==6))

##

indis_ref_filter = indis_ref[,,match(indi_to_use,dimnames(indis_ref)[[3]])]
tt = apply(indis_ref_filter,3,rowSums)
ii = rowSums(tt>0)

gene_to_use = which(ii>=97)
gene_name_12406 = names(gene_to_use)
saveRDS(gene_name_12406,file = 'data/neuron/gene_name_12406.rds')

indis_ref_filter = indis_ref_filter[match(gene_name_12406,dimnames(indis_ref_filter)[[1]]),,]


## run
G = dim(indis_ref_filter)[1]
corGene_idx = readRDS("data/neuron/corGene_idx_lower_cpm_alpha005.rds")
library(Matrix)
A.indicator = sparseMatrix(i = corGene_idx[,1],j = corGene_idx[,2],x=1,dims=c(20293,20293))
idx = match(gene_name_12406,gene_name)
A.indicator = A.indicator[idx,idx]
cor.idx = which(A.indicator!=0,arr.ind = TRUE)
cor.idx = rbind(cor.idx,cbind(cor.idx[,2],cor.idx[,1]))

b1 = c(0.1,0.1,0.15,0.15,0.2,0.3)
b2 = c(0.1,0.15,0.25,0.3,0.1,0.1)


n = dim(indis_ref_filter)[3]
n_ref = 11
n_bulk = n-n_ref
b = cbind(b1%*%t(rep(1,n_bulk/2)),b2%*%t(rep(1,n_bulk/2)))


set.seed(12345)
out10 = list()
for(i in 1:10){
  print(i)
  ref.idx = sort(sample(1:n,n_ref))
  out10[[i]] = simu_neuron(indis_ref_filter,ref.idx,b,cor.idx,calc_cov = FALSE,verbose=F,weighted = TRUE)
  saveRDS(out10,file='output/neuron/neuron_simu_ref11_weight.rds')
}




n_ref = 21
n_bulk = n-n_ref
b = cbind(b1%*%t(rep(1,n_bulk/2)),b2%*%t(rep(1,n_bulk/2)))

set.seed(12345)
out20 = list()
for(i in 1:10){
  print(i)
  ref.idx = sort(sample(1:n,n_ref))
  out20[[i]] = simu_neuron(indis_ref_filter,ref.idx,b,cor.idx,calc_cov = FALSE,verbose=F,weighted = TRUE)
  saveRDS(out20,file='output/neuron/neuron_simu_ref21_weight.rds')
}


##############vash V###################
n = dim(indis_ref_filter)[3]
n_ref = 11
n_bulk = n-n_ref
b = cbind(b1%*%t(rep(1,n_bulk/2)),b2%*%t(rep(1,n_bulk/2)))


set.seed(12345)
out10.eb = list()
for(i in 1:10){
  print(i)
  ref.idx = sort(sample(1:n,n_ref))
  out10.eb[[i]] = simu_neuron(indis_ref_filter,ref.idx,b,cor.idx,calc_cov = FALSE,verbose=F,weighted = FALSE,eb.V = TRUE)
  saveRDS(out10.eb,file='output/neuron/neuron_simu_ref11_ebV.rds')
}




n_ref = 21
n_bulk = n-n_ref
b = cbind(b1%*%t(rep(1,n_bulk/2)),b2%*%t(rep(1,n_bulk/2)))

set.seed(12345)
out20.eb = list()
for(i in 1:10){
  print(i)
  ref.idx = sort(sample(1:n,n_ref))
  out20.eb[[i]] = simu_neuron(indis_ref_filter,ref.idx,b,cor.idx,calc_cov = FALSE,verbose=F,weighted = FALSE,eb.V = TRUE)
  saveRDS(out20.eb,file='output/neuron/neuron_simu_ref21_ebV.rds')
}
##############################################


# n_ref = 31
# n_bulk = n-n_ref
# b = cbind(b1%*%t(rep(1,n_bulk/2)),b2%*%t(rep(1,n_bulk/2)))
#
# set.seed(12345)
# out30 = list()
# for(i in 1:10){
#   print(i)
#   ref.idx = sort(sample(1:n,n_ref))
#   out30[[i]] = simu_neuron(indis_ref_filter,ref.idx,b,cor.idx,calc_cov = FALSE,verbose=F,weighted = TRUE)
#   saveRDS(out30,file='output/neuron/neuron_simu_ref31_weight.rds')
# }



## MuSiC

n_ref = 31
n_bulk = n-n_ref
b = cbind(b1%*%t(rep(1,n_bulk/2)),b2%*%t(rep(1,n_bulk/2)))
set.seed(12345)
out.music = list()
for(i in 1:10){
  print(i)
  ref.idx = sort(sample(1:n,n_ref))
  out.music[[i]] = simu_neuron_music(indis_ref_filter,ref.idx,b)
  saveRDS(out.music,file='output/neuron/neuron_simu_ref31_music.rds')
}

##

