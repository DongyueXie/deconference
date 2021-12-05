

source('code/deconference_main.R')
source('code/ols_hc3.R')
source('code/CIBERSORT.R')

# This function takes generated ref, p, groups as input


library(MuSiC)
library(Matrix)
library(reticulate)
library(gtools)
library(vashr)

sce_to_es = function(sce){
  pheno_Data = data.frame(cell_type = sce$cell_type,individual = sce$individual)
  rownames(pheno_Data) = colnames(sce)
  ref.data = ExpressionSet(assayData = as.matrix(counts(sce)),phenoData = AnnotatedDataFrame(pheno_Data))
  ref.data
}

neuron_simu_study = function(indis_ref,
                             ref_idx_mat,
                             bulk_idx_mat,
                             bulk_p_array,
                             groups_mat,
                             case,
                             R01,
                             folds,
                             add_bulk_bias,
                             cor_fdr = 0.05,
                             output.path = 'output/manuscript/real/',
                             method_list = c('my','music','cibersort','ols'),
                             bulk_bias_sd = 0.1,
                             p_null=c(0.3,0.2,0.15,0.15,0.1,0.1),
                             p1_diff=c(0.15,0.15,0.1,0.1,0.2,0.3),
                             p2_diff=c(0.1,0.1,0.2,0.3,0.15,0.15),
                             bulk_lib_size = 500,
                             seed=12345,
                             verbose=FALSE,
                             calc_var=TRUE){
  gene_to_use = dimnames(indis_ref)[[1]]
  individuals = dimnames(indis_ref)[[3]]
  celltypes = dimnames(indis_ref)[[2]]
  n_indi = length(individuals)
  n_ref = ncol(ref_idx_mat)
  n_rep = nrow(ref_idx_mat)
  G = length(gene_to_use)
  K = length(celltypes)
  n_bulk = ncol(groups_mat)

  if(case=='null'){
    p1 = p_null
    p2 = p_null
  }else if(case=='all_diff'){
    p1 = p1_diff
    p2 = p2_diff
  }


  set.seed(seed)

  my_fit = list()
  music_fit = list()
  cibersort_fit =  list()
  ols_fit = list()
  rep_info = list()
  for(r in 1:n_rep){

    if(verbose){
      print(paste('running rep', r))
    }


    ref.idx = ref_idx_mat[r,]
    X_array_ref = indis_ref[,,ref.idx]
    X_array_bulk = indis_ref[,,bulk_idx_mat[r,]]

    # create bulk data

    if(verbose){
      print('creating bulk data')
    }

    groups = groups_mat[r,]
    b = bulk_p_array[,,r]
    #rownames(b) = celltypes
    #colnames(b) = individuals[-ref.idx]
    #perm.temp = sample(1:n_bulk)
    #groups = groups[perm.temp]
    #b = b[,perm.temp]

    mb = lapply(1:n_bulk,function(i){X_array_bulk[,,i]%*%b[,i]})
    mb = do.call(cbind,mb)
    if(add_bulk_bias){
      mb = mb*rnorm(G,1,bulk_bias_sd)
    }
    thetab = apply(mb,2,function(z){z/sum(z)})
    y = matrix(rpois(G*n_bulk,bulk_lib_size*G*thetab),nrow=G)
    rownames(y) = gene_to_use
    colnames(y) = 1:n_bulk

    rep_info[[r]] = list(ref.idx=ref.idx,b=b,groups=groups)

    #browser()

    ############################################################################
    ############# fit music methods ###############################################
    if('music'%in%method_list){
      if(verbose){
        print('preparing data for music')
      }

      phenoData = data.frame(SubjectName = individuals[bulk_idx_mat[r,]])
      rownames(phenoData) = colnames(y)
      bulk_es = ExpressionSet(assayData = y,phenoData = AnnotatedDataFrame(phenoData))



      # create reference expression set object
      #ref_sce = day30sce[gene_to_use,day30sce$individual%in%(individuals[ref.idx])]
      ref_sce = readRDS(paste('data/neuron/real_manu/sce_by_97_indis/',individuals[ref.idx[1]],'_sce.rds',sep=''))
      ref_samples_cibersort = readRDS(paste('data/neuron/real_manu/sce_by_97_indis/',individuals[ref.idx[1]],'_by_cell.rds',sep=''))
      for(i in ref.idx[-1]){
        ref_sce = cbind(ref_sce,
                        readRDS(paste('data/neuron/real_manu/sce_by_97_indis/',individuals[i],'_sce.rds',sep='')))
        ref_samples_cibersort = cbind(ref_samples_cibersort,
                                      readRDS(paste('data/neuron/real_manu/sce_by_97_indis/',individuals[i],'_by_cell.rds',sep='')))
      }
      #browser()
      ref_sce$individual = droplevels(ref_sce$individual)
      ref_sce = ref_sce[gene_to_use,]
      ref_es = sce_to_es(ref_sce)
      if(verbose){
        print('fitting music')
      }
      fit_music = music_prop(bulk_es,ref_es,clusters='cell_type',
                             samples='individual',select.ct = celltypes,
                             verbose=FALSE)

      music_fit[[r]] = fit_music

      if(r%%10==0){

        output.name = paste("music/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,sep='')
        saveRDS(list(music_fit=music_fit,
                     rep_info=rep_info,
                     p1=p1,
                     p2=p2),file = paste(output.path,output.name,'.rds',sep=''))
      }

    }

    ############################################################################
    ############# fit my methods ###############################################
    if('my'%in%method_list){
      if(verbose){
        print('fitting measurement model')
      }

      # fit our model
      X = apply(X_array_ref,c(1,2),mean,na.rm=TRUE)
      colnames(X) = celltypes
      V = t(apply(X_array_ref,c(1),function(z){(cov(t(z),use = 'complete.obs'))}))/n_ref
      #browser()
      w = 1/((vashr::vash(sqrt(rowSums(V)),df=n_ref-1))$sd.post)^2

      fit.err.cor.weight.cv = estimation_func2(y=y,
                                               X=X,
                                               Vg=V,
                                               w=w,
                                               hc.type='jackknife_indep',
                                               verbose=F,
                                               R01=R01,
                                               folds=folds,
                                               groups=groups,
                                               calc_var=calc_var)
      my_fit[[r]] = fit.err.cor.weight.cv

      if(r%%10==0){

        output.name = paste("my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr",paste(strsplit(as.character(cor_fdr),split='')[[1]][-2],collapse = ''),"_",case,sep='')
        saveRDS(list(my_fit=my_fit,
                     rep_info=rep_info,
                     p1=p1,
                     p2=p2),file = paste(output.path,output.name,'.rds',sep=''))
      }

    }

    ############################################################################
    ############# fit my methods - no weight ###############################################
    if('my_unweighted'%in%method_list){
      if(verbose){
        print('fitting measurement model, unweighted')
      }

      # fit our model
      X = apply(X_array_ref,c(1,2),mean,na.rm=TRUE)
      colnames(X) = celltypes
      V = t(apply(X_array_ref,c(1),function(z){(cov(t(z),use = 'complete.obs'))}))/n_ref
      #browser()
      #w = 1/((vashr::vash(sqrt(rowSums(V)),df=n_ref-1))$sd.post)^2

      fit.err.cor.weight.cv = estimation_func2(y=y,
                                               X=X,
                                               Vg=V,
                                               w=1,
                                               hc.type='jackknife_indep',
                                               verbose=F,
                                               R01=R01,
                                               folds=folds,
                                               groups=groups,
                                               calc_var=calc_var)
      my_fit[[r]] = fit.err.cor.weight.cv

      if(r%%10==0){

        output.name = paste("my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr",paste(strsplit(as.character(cor_fdr),split='')[[1]][-2],collapse = ''),"_unweighted_",case,sep='')
        saveRDS(list(my_fit=my_fit,
                     rep_info=rep_info,
                     p1=p1,
                     p2=p2),file = paste(output.path,output.name,'.rds',sep=''))
      }

    }



    ############################################################################
    ############# fit cibersort methods ###############################################

    if('cibersort'%in%method_list){
      # formulate signature matrix for cibersort
      print('running cibersort')
      ref_samples_cibersort = ref_samples_cibersort[match(gene_to_use,rownames(ref_samples_cibersort)),]
      cibersort_sig_mat = build_signature_matrix_CIBERSORT(ref_samples_cibersort)
      fit_cibersort = CIBERSORT(cibersort_sig_mat,data.frame(GeneSymbol = rownames(y),y))
      cibersort_fit[[r]] = fit_cibersort

      if(r%%10==0){

        output.name = paste("cibersort/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,sep='')
        saveRDS(list(cibersort_fit=cibersort_fit,
                     rep_info=rep_info,
                     p1=p1,
                     p2=p2),file = paste(output.path,output.name,'.rds',sep=''))
      }

    }

    ############################################################################
    ############# fit ols methods ###############################################

    if('ols'%in%method_list){
      # formulate signature matrix for cibersort
      print('running ols')

      X = apply(X_array_ref,c(1,2),mean,na.rm=TRUE)
      colnames(X) = celltypes

      fit.ols = ols_hc3(y,X,groups=groups,calc_var=calc_var)
      ols_fit[[r]] = fit.ols
      if(r%%10==0){

        output.name = paste("ols/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,sep='')
        saveRDS(list(ols_fit=ols_fit,
                     rep_info=rep_info,
                     p1=p1,
                     p2=p2),file = paste(output.path,output.name,'.rds',sep=''))
      }

    }


  }


}
