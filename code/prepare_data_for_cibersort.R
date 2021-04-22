

## seger data

seger <- readRDS("data/pancreas/segerstolpe_raw.rds")
cell_types = c('alpha', 'beta', 'delta', 'gamma')
seger.mat = counts(seger)
colnames(seger.mat) = seger$cell_type
seger.mat = seger.mat[,colnames(seger.mat)%in%cell_types]

# remember to add gene at the very beginning manually in the file
write.table(seger.mat,file='data/pancreas/seger_cibersort.txt',quote = FALSE,sep = "\t")

## baron data

baron.mat = counts(baron)
colnames(baron.mat) = baron$cell_type
baron.mat = baron.mat[,colnames(baron.mat)%in%cell_types]
write.table(baron.mat,file='data/pancreas/baron_cibersort.txt',quote = FALSE,sep = "\t")

xin_raw <- readRDS("data/pancreas/xin_raw.rds")
createBulk_tech.noise = function(X_array, b, gene_names, ls_per_gene = 500,seed=12345){

  G = dim(X_array)[1]
  K = dim(X_array)[2]
  Ni = dim(X_array)[3]

  if(is.null(dim(b))){
    y = apply(X_array,3,function(z){z%*%b})
  }else{
    y = matrix(nrow=G,ncol=Ni)
    for(i in 1:Ni){
      y[,i] = X_array[,,i]%*%b[,i]
    }
  }
  y = y/colSums(y)*ls_per_gene*G
  set.seed(seed)
  y = matrix(rpois(prod(dim(y)),y),nrow=nrow(y),ncol=ncol(y))
  rownames(y) = gene_names
  list(y=y,b=b*ls_per_gene)
}
K = length(cell_types)
rm.indi = c("Non T2D 4","Non T2D 7","Non T2D 10","Non T2D 12")
#rm.indi = levels(XinT2D.eset$SubjectName)[rm.indi]
rm.indi.idx = which(xin_raw$individual%in%rm.indi)


datax.xin = set_data_decon(Y = xin_raw[,-rm.indi.idx],cell_types = cell_types, gene_thresh = 0,max_count_quantile = 1,w=1)

design.mat.xin = scRef_multi_proc(datax.xin$Y,datax.xin$cell_type_idx,datax.xin$indi_idx,
                                  estimator="separate",est_sigma2 = FALSE,
                                  meta_var='plug_in',meta_mode='smooth',
                                  verbose = F)

# create bulk data
b = c(0.1,0.1,0.3,0.5)
bulk.xin = createBulk_tech.noise(design.mat.xin$X_array,b,gene_names = rownames(xin_raw))


bulk.xin = bulk.xin$y
colnames(bulk.xin) = 1:ncol(bulk.xin)
write.table(bulk.xin,file='data/pancreas/bulk_xin_cibersort.txt',quote = FALSE,sep = "\t")
