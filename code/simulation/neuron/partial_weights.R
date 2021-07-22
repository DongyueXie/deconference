source('code/deconference_main.R')
source('code/simulation/neuron/simu_neuron.R')
#source('code/simulation/simu_correlation_ult.R')

# # remove genes
#
# gene_name = readRDS('data/neuron/gene_name_20293.rds')
# indis_ref = readRDS('data/neuron/indis_ref.rds')
# indis_ref = indis_ref[match(gene_name,dimnames(indis_ref)[[1]]),,]
# dim(indis_ref)
#
# # remove individuals that has no or too few certain cell types
# cell_ann = readRDS('data/neuron/cell_ann.rds')
# indi_by_cell = table(cell_ann$donor_id,cell_ann$celltype)
# indi_by_cell[,6] = indi_by_cell[,6] + indi_by_cell[,7]
# indi_by_cell = indi_by_cell[,1:6]
# colnames(indi_by_cell)[6] = 'U_Neur'
#
# indi_to_use = names(which(rowSums(indi_by_cell>10)==6))
#
# ##
#
# indis_ref_filter = indis_ref[,,match(indi_to_use,dimnames(indis_ref)[[3]])]
# tt = apply(indis_ref_filter,3,rowSums)
# ii = rowSums(tt>0)
#
# gene_to_use = which(ii>=97)
# gene_name_12406 = names(gene_to_use)
# indis_ref_filter = indis_ref_filter[match(gene_name_12406,dimnames(indis_ref_filter)[[1]]),,]
# X = apply(indis_ref_filter,c(1,2),mean,na.rm=TRUE)
# V = t(apply(indis_ref_filter,c(1),function(z){(cov(t(z),use = 'complete.obs'))}))/dim(indis_ref_filter)[3]
# h = rowSums((X%*%solve(crossprod(X)-matrix(colSums(V),ncol=6)))*X)
# plot(h)
#
# # remove genes whose h is larger than 0.2
# rm.idx = which(h>0.2)
#
# indis_ref_filter = indis_ref_filter[-rm.idx,,]
#
# G = dim(indis_ref_filter)[1]
# corGene_idx = readRDS("data/neuron/corGene_idx_lower_cpm_alpha05.rds")
# library(Matrix)
# A.indicator = sparseMatrix(i = corGene_idx[,1],j = corGene_idx[,2],x=1,dims=c(20293,20293))
# idx = match(gene_name_12406[-rm.idx],gene_name)
# A.indicator = A.indicator[idx,idx]
# cor.idx = which(A.indicator!=0,arr.ind = TRUE)
# cor.idx = rbind(cor.idx,cbind(cor.idx[,2],cor.idx[,1]))
# saveRDS(cor.idx,file='data/neuron/gene12400_cor_idx_alpha05.rds')

indis_ref_filter = readRDS('data/neuron/indis_ref_12400by6by97.rds')


b1 = c(0.1,0.1,0.15,0.15,0.2,0.3)
b2 = c(0.1,0.15,0.25,0.3,0.1,0.1)

n = dim(indis_ref_filter)[3]
n_ref = 11
n_bulk = n-n_ref
b = cbind(b1%*%t(rep(1,n_bulk/2)),b2%*%t(rep(1,n_bulk/2)))



#############################################################
#############################################################

# only add positive residuals

## alpha = 0.05
cor.idx = readRDS('data/neuron/gene12400_cor_idx_alpha005.rds')
set.seed(12345)
out10 = list()
for(i in 1:10){
  print(i)
  ref.idx = sort(sample(1:n,n_ref))
  out10[[i]] = simu_neuron(indis_ref_filter,ref.idx,b,cor.idx,
                           calc_cov = FALSE,verbose=F,weighted = TRUE,only.add.pos.res = TRUE)
  saveRDS(out10,file='output/neuron/neuron_simu_ref11_gene12400_weight_alpha005_add_pos_res.rds')
}

## alpha = 0.5
cor.idx = readRDS('data/neuron/gene12400_cor_idx_alpha05.rds')
set.seed(12345)
out05p = list()
for(i in 1:5){
  print(i)
  ref.idx = sort(sample(1:n,n_ref))
  out05p[[i]] = simu_neuron(indis_ref_filter,ref.idx,b,cor.idx,
                           calc_cov = FALSE,verbose=F,weighted = TRUE,only.add.pos.res = TRUE)
  saveRDS(out05p,file='output/neuron/neuron_simu_ref11_gene12400_weight_alpha05_add_pos_res.rds')
}



set.seed(12345)
out05 = list()
for(i in 1:5){
  print(i)
  ref.idx = sort(sample(1:n,n_ref))
  out05[[i]] = simu_neuron(indis_ref_filter,ref.idx,b,cor.idx,
                           calc_cov = FALSE,verbose=F,weighted = TRUE,only.add.pos.res = FALSE)
  saveRDS(out05,file='output/neuron/neuron_simu_ref11_gene12400_weight_alpha05.rds')
}

#############################################################
#############################################################
# use all individual for obtaining weights
V.temp = t(apply(indis_ref_filter,c(1),function(z){(cov(t(z),use = 'complete.obs'))}))
fit.vash = vashr::vash(sqrt(rowSums(V.temp)),df=n-1)
w = 1/(fit.vash$sd.post)^2
w = w/sum(w)*G

set.seed(12345)
out10 = list()
for(i in 1:10){
  print(i)
  ref.idx = sort(sample(1:n,n_ref))
  out10[[i]] = simu_neuron(indis_ref_filter,ref.idx,b,cor.idx,
                           calc_cov = FALSE,verbose=F,weighted = TRUE,w=w)
  saveRDS(out10,file='output/neuron/neuron_simu_ref11_rm_outlier_weight_allcalcweight_alpha005.rds')
}




#############################################################
#############################################################
set.seed(12345)
out10 = list()
for(i in 1:10){
  print(i)
  ref.idx = sort(sample(1:n,n_ref))
  out10[[i]] = simu_neuron(indis_ref_filter,ref.idx,b,cor.idx,
                           calc_cov = FALSE,verbose=F,weighted = TRUE)
  saveRDS(out10,file='output/neuron/neuron_simu_ref11_rm_outlier_weight.rds')
}

#############################################################
#############################################################

set.seed(12345)
out10 = list()
for(i in 1:10){
  print(i)
  ref.idx = sort(sample(1:n,n_ref))
  out10[[i]] = simu_neuron(indis_ref_filter,ref.idx,b,cor.idx,
                           calc_cov = FALSE,verbose=F,weighted = FALSE,
                           eb.V = TRUE)
  saveRDS(out10,file='output/neuron/neuron_simu_ref11_rm_outlier_ebV.rds')
}

