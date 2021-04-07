


refs_raw <- readRDS("data/pancreas/refs_raw.rds")

Xin_raw_filtered = refs_raw$xin
Xin_raw_filtered_es = sce_to_es(Xin_raw_filtered)
Xin_bulk = bulk_construct(Xin_raw_filtered_es,'cell_type','individual')

saveRDS(Xin_raw_filtered,file='data/pancreas/rnasieve/Xin_raw_filtered.rds')


## create data for rna-seieve

Xin_bulk_count = exprs(Xin_bulk$Bulk.counts)

Xin_ref_list = list()
cell_types = c('alpha', 'beta', 'delta', 'gamma')
for(i in 1:4){
  cell_mat = counts(Xin_raw_filtered)[,which(Xin_raw_filtered$cell_type==cell_types[i])]
  saveRDS(cell_mat,file = paste('data/pancreas/rnasieve/',cell_types[i],'.rds',sep = ''))
  Xin_ref_list[[i]] = cell_mat
}

saveRDS(Xin_bulk_count,file='data/pancreas/rnasieve/Xin_bulk_count.rds')



ref_data = create_bulk_abitrary_prop(list(Xin_raw_filtered),cell_types=cell_types,c(0.4,0.3,0.2,0.1),c(0.1,0.1,0.2,0.6),ncells = 5000)
saveRDS(counts(ref_data$bulk),file='data/pancreas/rnasieve/Xin_bulk_count_arbitrary_prop.rds')

saveRDS(Xin_bulk$num.real/rowSums(Xin_bulk$num.real),file='data/pancreas/rnasieve/Xin_bulk_count_true_beta.rds')
saveRDS(ref_data$bulk_prop,file='data/pancreas/rnasieve/Xin_bulk_count_arbitrary_prop_true_beta.rds')




my_out = deconference_multi_ref(list(Xin_raw_filtered),ref_data$bulk)

bhat = t(my_out$beta_hat)
colnames(bhat) = cell_types

Eval_multi(ref_data$bulk_prop,list(music_out$Est.prop.weighted, b = bhat))



my_out2 = deconference_multi_ref(list(Xin_raw_filtered),SingleCellExperiment(assays = list(counts=exprs(Xin_bulk$Bulk.counts))))
summary_real(my_out2,as.matrix(Xin_bulk$num.real/rowSums(Xin_bulk$num.real)),groups = c(rep(1,9),rep(2,9)))
