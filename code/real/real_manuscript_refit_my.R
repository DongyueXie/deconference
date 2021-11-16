library(MuSiC)
library(Matrix)
library(reticulate)
library(gtools)




############ revise Dirichlet, re-generate b each time!!!!!!!!!
############ need to refit my, music, re-store bulk data and b etc.

neuron_simu_study = function(indis_ref,
                             day30sce,
                             out,
                             p1,
                             p2,
                             dirichlet = TRUE,
                             dirichlet.scale=10,
                             p1.prop=0.5,
                             n_rep=1,
                             n_ref = 11,
                             folds = NULL,
                             bulk_lib_size = 500,
                             seed=12345,
                             verbose=FALSE,
                             R01=NULL){
  gene_to_use = dimnames(indis_ref)[[1]]
  individuals = dimnames(indis_ref)[[3]]
  celltypes = dimnames(indis_ref)[[2]]
  n_indi = length(individuals)
  G = length(gene_to_use)
  K = length(p1)
  n_bulk = n_indi - n_ref
  n.p1.temp =  round(n_bulk*p1.prop)


  set.seed(seed)

  my_fit = list()
  music_fit = list()
  rep_info = list()
  for(r in 1:n_rep){
    if(verbose){
      print(paste('running rep', r))
    }

    # randomly draw reference individuals
    ref.idx = out$rep_info[[r]]$ref.idx


    X_array_ref = indis_ref[,,ref.idx]
    X_array_bulk = indis_ref[,,-ref.idx]

    # create bulk data

    if(verbose){
      print('creating bulk data')
    }


    # generate bulk true cell proportions

    groups = c(rep(1,n.p1.temp),rep(2,n_bulk-n.p1.temp))
    if(dirichlet){
      b = cbind(t(rdirichlet(n.p1.temp,p1*dirichlet.scale)),t(rdirichlet(n_bulk-n.p1.temp,p2*dirichlet.scale)))
    }else{
      b = cbind(p1%*%t(rep(1,n.p1.temp)),p2%*%t(rep(1,n_bulk-n.p1.temp)))
    }
    perm.temp = sample(1:n_bulk)
    groups = groups[perm.temp]
    b = b[,perm.temp]

    mb = lapply(1:n_bulk,function(i){X_array_bulk[,,i]%*%b[,i]})
    mb = do.call(cbind,mb)
    thetab = apply(mb,2,function(z){z/sum(z)})
    y = matrix(rpois(G*n_bulk,bulk_lib_size*G*thetab),nrow=G)
    rownames(y) = gene_to_use
    colnames(y) = individuals[-ref.idx]

    # fit our model
    X = apply(X_array_ref,c(1,2),mean,na.rm=TRUE)
    V = t(apply(X_array_ref,c(1),function(z){(cov(t(z),use = 'complete.obs'))}))/n_ref
    w = 1/((vashr::vash(sqrt(rowSums(V)),df=n_ref-1))$sd.post)^2

    if(verbose){
      print('fitting measurement model')
    }

    fit.err.cor.weight.cv = estimation_func2(y=y,
                                             X=X,
                                             Vg=V,
                                             w=w,
                                             hc.type='jackknife_indep',
                                             verbose=verbose,
                                             R01=R01,
                                             folds=folds,
                                             groups=groups)

    my_fit[[r]] = fit.err.cor.weight.cv
    # fit music
    # create bulk expression set object

    if(verbose){
      print('preparing data for music')
    }

    phenoData = data.frame(SubjectName = individuals[-ref.idx])
    rownames(phenoData) = colnames(y)
    bulk_es = ExpressionSet(assayData = y,phenoData = AnnotatedDataFrame(phenoData))

    # create reference expression set object
    ref_sce = day30sce[gene_to_use,day30sce$individual%in%(individuals[ref.idx])]
    ref_sce$individual = droplevels(ref_sce$individual)
    ref_es = sce_to_es(ref_sce)
    if(verbose){
      print('fitting music')
    }
    fit_music = music_prop(bulk_es,ref_es,clusters='cell_type',
                           samples='individual',select.ct = celltypes,
                           verbose=verbose)

    music_fit[[r]] = fit_music

    rep_info[[r]] = list(ref.idx=ref.idx,b=b,groups=groups)
    # save bulk and ref for running RNA-sieve and cibersort
    # cibersort
    if(verbose){
      print('saving data for cibersort')
    }
    #ref_mat = exprs(ref_es)
    rm(ref_es)
    #rownames(ref_mat) = rownames(ref_sce)
    #colnames(ref_mat) = ref_sce$cell_type
    #ref_mat = data.frame('gene' = rownames(ref_mat),ref_mat)
    #write.table(ref_mat,
    #            file=paste('data/neuron/real_manu/cibersort/Ref_nref',n_ref,'_rep',r,'.txt',sep = ''),
    #            row.names = FALSE,col.names = c("gene",as.character(ref_sce$cell_type)),quote = FALSE,sep = "\t")

    write.table(data.frame('gene'=rownames(y),y),
                file=paste('data/neuron/real_manu/cibersort/Bulk_nref',n_ref,'_rep',r,'.txt',sep = ''),
                row.names = FALSE,quote = FALSE,sep = "\t")

    #rm(ref_mat)
    rm(y)

    # # RNA-sieve
    # if(verbose){
    #   print('saving data for rna-sieve')
    # }
    # ref_list = list()
    # for(k in 1:K){
    #   ref_list[[k]] = as.matrix(counts(ref_sce[,ref_sce$cell_type==celltypes[k]]))
    # }
    # names(ref_list) = celltypes
    # ref_list = r_to_py(ref_list)
    # py_save_object(ref_list,filename=paste('data/neuron/real_manu/rnasieve/Ref_nref',n_ref,'_rep',r,sep = ''))
    #
    # rm(ref_list)
    # for bulk data run:
    # temp = pd.read_csv('bulk.txt',sep="\t",index_col=0)
    # temp.values



  }
  return(list(my_fit=my_fit,music_fit=music_fit,rep_info=rep_info,p1=p1,p2=p2))

}



out = readRDS('output/manuscript/real/neuron_ref11_rep10_fdr005.rds')
source('code/real/real_manuscript.R')
source('code/deconference_main.R')
# For neuron dataset,
indis_ref = readRDS('data/neuron/indis_ref_12400by6by97.rds')
day30sce = readRDS('data/neuron/day30sce.rds')
R01 = readRDS("data/neuron/R01_neuron_G12400_alpha005.rds")

R01.dist = as.dist(1-R01)
clusters01 = cluster::pam(R01.dist,10,pamonce=5)
folds01 = clusters01$clustering
table(folds01)

p1 = c(0.1,0.1,0.15,0.15,0.2,0.3)
p2 = c(0.1,0.15,0.25,0.3,0.1,0.1)
real_out = neuron_simu_study(indis_ref,day30sce,out,
                             p1,p2,dirichlet=TRUE,p1.prop=0.5,
                             bulk_lib_size = 500,
                             n_ref=11,n_rep=10,seed=12345,
                             verbose=TRUE,folds=folds01,
                             R01=R01)
saveRDS(real_out,file = 'output/manuscript/real/neuron_ref11_rep10_fdr005_dirichlet10.rds')

