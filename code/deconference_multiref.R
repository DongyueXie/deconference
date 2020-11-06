

#'@param ref.obj either raw counts or the output by set_data_multiref; if raw, set gene-length_adjust to true and input gene length; otherwise set to false.
deconference_multi_ref = function(ref.obj,bulk.obj,tau2=NULL,cell_types=NULL,sigma2=NULL,
                          est_sigma2=TRUE,meta_var='adjust',meta_mode='smooth',
                          gene_length_adjust = FALSE,
                          gene_length=NULL,protocol=NULL,
                          correction=FALSE,cellsize_est='glm',
                          #marker_gene = NULL,
                          a=length(cell_types)+4,
                          calc_cov=TRUE,
                          hc.type = 'hc3',w = 1,verbose=FALSE){


  if(verbose){
    message('constructing reference matrix')
    message('...going through individuals')
  }

  temp = getXV_array_all(ref.obj=ref.obj,tau2=tau2,cell_types=cell_types,
                         indis = NULL,gene_length_adjust=gene_length_adjust,
                         gene_length=gene_length,protocol=protocol,cellsize_est=cellsize_est)

  if(verbose){
    message('...merging individuals')
  }
  design.mat = getXV(temp$all_X_array,temp$all_Vg_array,
                     S=temp$S,sigma2=sigma2,
                   est_sigma2=est_sigma2,meta_var=meta_var,meta_mode=meta_mode,
                   cell_types=cell_types)
  #browser()

  if(gene_length_adjust){
    if(!is.null(gene_length)){
      gene_length = gene_length/sum(gene_length)*length(gene_length)
      bulk_counts = counts(bulk.obj)/gene_length
    }
  }else{
    bulk_counts = counts(bulk.obj)
  }


  out = estimation_func2(y=bulk_counts,X=design.mat$X,Vg=design.mat$Vg,design.mat$Sigma,
                        w=w,hc.type=hc.type,correction=correction,S=design.mat$S,calc_cov=calc_cov,a=a,verbose=verbose)
  rownames(out$beta_hat) = colnames(design.mat$X)
  return(out)

}


#'@param gene_thresh genes that appear less than ...
#'@param cell_thresh cells that have too few expressed genes removed
preprocess_sc = function(Y,
                         gene_length=NULL,
                         protocol = NULL,
                         gene_thresh=0.05,
                         max_count_quantile=0.99,
                         cell_types=NULL){
  temp_idx = which(Y$cell_type%in%cell_types)
  Y = Y[,temp_idx]
  K = length(cell_types)

  if(protocol=='nonUMI'&(!is.null(gene_length))){
    cm_gene = intersect(rownames(Y),names(gene_length))
    gene_y_idx = match(cm_gene,rownames(Y))
    gene_len_idx = match(cm_gene,names(gene_length))
    Y = Y[gene_y_idx,]
    gene_len = gene_length[gene_len_idx]
    gene_len = gene_len/sum(gene_len)*length(gene_len)
    counts(Y) = counts(Y)/gene_len
  }

  if(gene_thresh<1){
    gene_thresh = round(gene_thresh*ncol(Y))
  }

  rm.gene.low = which(rowSums(counts(Y)>0)<gene_thresh)

  # remove union of genes that expressed more than max_count_quantle each cell type

  rm.gene.high = c()
  for(k in 1:K){
    cell_k_idx = which(Y$cell_type==cell_types[k])
    if(length(cell_k_idx)!=0){
      gene_counts = rowSums(counts(Y)[,cell_k_idx])
      rm.gene.high = c(rm.gene.high,which(gene_counts>quantile(gene_counts,max_count_quantile)))
    }

  }

  rm.gene = unique(c(rm.gene.low,rm.gene.high))

  if(length(rm.gene)!=0){
    Y = Y[-rm.gene,]
  }

  rm.cell = which(colSums(counts(Y))==0)
  if(length(rm.cell)!=0){
    Y = Y[,-rm.cell]
  }

  if(is.factor(Y$cell_type)){
    Y$cell_type = droplevels(Y$cell_type)
  }

  Y

}

# Y is a combination of all indis/studies expressions
filter_gene_all = function(Y,gene_thresh){
  if(gene_thresh<1){
    gene_thresh = round(gene_thresh*ncol(Y))
  }

  rm.gene.low = which(rowSums(Y>0)<gene_thresh)
  rm.gene.low
}

# Y is a sce obj
filter_cell = function(Y,cell_thresh){
  if(cell_thresh<1){
    cell_thresh = round(cell_thresh*nrow(Y))
  }
  rm.cell = which(colSums(counts(Y)>0)<cell_thresh)
  if(length(rm.cell)!=0){
    Y = Y[,-rm.cell]
  }
  Y
}

# input raw ref and bulk dataset, return them with common genes.
# note for full length dataset, adjusted for gene length
set_data_multiref = function(ref.obj,bulk.obj = NULL,
                             gene_length=NULL,protocols=NULL,
                             marker_gene=NULL,
                             gene_thresh=0.05,
                             cell_thresh = 0.01,
                             max_count_quantile=0.99,
                             cell_types = c('acinar','alpha','beta','delta','ductal','gamma')){
  n_study = length(ref.obj)

  # pre-process each data set(remove too high expressed genes), and get common genes among reference dataset
  n_indis_all = 0
  if(!is.null(gene_length)){
    common_genes = names(gene_length)
  }else{
    common_genes = rownames(ref.obj[[1]])
  }

  for(i in 1:n_study){
    temp = preprocess_sc(ref.obj[[i]],
                         gene_length=gene_length,
                         protocol = protocols[i],
                         gene_thresh=0,
                         max_count_quantile=max_count_quantile,
                         cell_types=cell_types)

    ref.obj[[i]] = temp
    n_indis_all = n_indis_all + ncol(temp)
    common_genes = intersect(common_genes,rownames(temp))
  }

  #browser()

  # get common genes among ref and bulk data
  if(!is.null(bulk.obj)){
    common_genes = intersect(common_genes,rownames(bulk.obj))
  }

  # get marker genes
  if(!is.null(marker_gene)){
    common_genes = intersect(common_genes,marker_gene)
  }


  # remove cells with no gene expressed and genes with little expression
  all_counts = c()
  for(i in 1:n_study){
    gene_idx = match(common_genes,rownames(ref.obj[[i]]))
    ref.obj[[i]] = filter_cell((ref.obj[[i]])[gene_idx,],cell_thresh)
    all_counts = cbind(all_counts,counts(ref.obj[[i]]))
  }
  rm.gene = filter_gene_all(all_counts,gene_thresh)
  if(length(rm.gene)!=0){
    for(i in 1:n_study){
      ref.obj[[i]] = (ref.obj[[i]])[-rm.gene,]
    }
    common_genes = common_genes[-rm.gene]
  }

  # now finally we have genes to use and we can output ref and bulk data
  if(!is.null(bulk.obj)){
    gene_idx = match(common_genes,rownames(bulk.obj))
    bulk.obj = (bulk.obj)[gene_idx,]
    #browser()
    if(!is.null(gene_length)){
      gene_len = gene_length[match(common_genes,names(gene_length))]
      gene_len = gene_len/sum(gene_len)*length(gene_len)
      counts(bulk.obj) = counts(bulk.obj)/gene_len
    }

    return(list(ref.obj = ref.obj,bulk.obj = bulk.obj))
  }else{
    return(ref.obj)
  }

}


# input reference data object
# output arrays of all individuals' X and V, and estimate of cell size.
getXV_array_all = function(ref.obj,tau2,cell_types,indis,gene_length_adjust,gene_length,protocol,cellsize_est=NULL){

  library(abind)
  n_ref = length(ref.obj)

  all_X_array = c()
  all_Vg_array = c()

  S = c()



  for(i in 1:n_ref){
    out_array = getXV_array(ref.obj[[i]],tau2=tau2,cell_types = cell_types,indis=NULL,
                            gene_length_adjust=gene_length_adjust,
                            gene_length = gene_length,protocol = protocol[i],cellsize_est=cellsize_est)
    all_X_array = abind(all_X_array,out_array$X_array)
    all_Vg_array = abind(all_Vg_array,out_array$Vg_array)
    S = rbind(S,out_array$S)
    #S_glms = rbind(S_glms,out_array$S_glm)
  }

  return(list(all_X_array=all_X_array,
              all_Vg_array=all_Vg_array,
              S=S))

}

# input one reference dataset
# output X,V,S
getXV_array = function(Y,cell_type_idx=NULL,indi_idx=NULL,
                       estimator='separate',eps=0,
                       tau2=NULL,cell_types = NULL,indis=NULL,
                       gene_length_adjust = FALSE,
                       gene_length=NULL,
                       protocol=NULL,
                       cellsize_est=NULL){

  #browser()

  if(is.null(cell_type_idx)&is.null(indi_idx)&class(Y) == "SingleCellExperiment"){
    cell_type_idx = Y$cell_type
    indi_idx = Y$individual
    Y = counts(Y)
  }
  if(gene_length_adjust){
    gene_length = gene_length/sum(gene_length)*length(gene_length)
    if(!is.null(protocol)&!is.null(gene_length)){
      if(protocol=='nonUMI'){
        Y = Y/gene_length
      }
    }
  }


  G = nrow(Y)

  if(is.null(indis)){
    if(is.factor(indi_idx)){
      indis = levels(indi_idx)
    }else{
      indis = levels(as.factor(indi_idx))
    }
  }


  if(is.null(cell_types)){
    if(is.factor(cell_type_idx)){
      cell_types = levels(cell_type_idx)
    }else{
      cell_types = levels(as.factor(cell_type_idx))
    }
  }


  K = length(cell_types)
  NI = length(indis)


  # first, for each individual, obtain X and Vg
  # then, for each g and k, perform PM method


  X_array = array(dim = c(G,K,NI))
  Vg_array = array(dim = c(G,K,NI))
  S_mat = matrix(nrow=NI,ncol=K)


  #browser()

  for(i in 1:NI){
    indi = indis[i]
    indi_cell_idx = which(indi_idx==indi)
    Yi = Y[,indi_cell_idx]
    cell_type_i = cell_type_idx[indi_cell_idx]
    indi_design.mat = scRef1_proc(Yi,cell_type_i,estimator=estimator,
                                  tau2=tau2,cell_types=cell_types)
    #browser()
    #print(dim(indi_design.mat$X))
    #browser()
    X_array[,,i] = indi_design.mat$X
    Vg_array[,,i] = indi_design.mat$Vg
    if(eps!=0){
      (Vg_array[,,i])[(Vg_array[,,i])==0] = eps
    }
    S_mat[i,] = indi_design.mat$S
  }

  #browser()

  if(!is.null(cellsize_est)){
    S_mat = round(S_mat)
    rownames(S_mat) = indis
    colnames(S_mat) = cell_types
    S = colMeans(S_mat,na.rm = TRUE)
    S = S/S[1]
    if(cellsize_est=='glm'&nrow(S_mat)>1){
      # ## method 2: glm
      #
      S_mat_dataframe = data.frame(y = c(S_mat),
                                   indi = factor(rep(indis,ncol(S_mat))),
                                   type = factor(rep(cell_types,each = nrow(S_mat)),levels = cell_types))
      suppressWarnings({fit = try(MASS::glm.nb(y~.,S_mat_dataframe),silent = TRUE)})

      if(class(fit)[1]=='try-error'){

        fit = glm(y~.,S_mat_dataframe,family = 'poisson')

      }

      S[which(!is.nan(S))] = c(1,exp(fit$coefficients[-c(1:nrow(S_mat))]))
      names(S) = cell_types
    }
  }else{
    S = rep(1,length(cell_types))
    names(S) = cell_types
  }

  return(list(X_array = X_array,Vg_array = Vg_array,S=S))

}

# input X_array, V_array, S
# output X,V,S
getXV = function(X_array,Vg_array,S=NULL,sigma2=NULL,
                 est_sigma2=TRUE,meta_var='adjust',meta_mode='smooth',
                 diag_cov=FALSE,cell_types){

  GKN = dim(X_array)
  G = GKN[1]
  K = GKN[2]
  NI = GKN[3]

  if(!is.null(sigma2)){
    X = matrix(nrow = G,ncol = K)
    Vg = matrix(nrow = G,ncol = K)

    for(g in 1:G){
      for(k in 1:K){
        x = X_array[g,k,]
        v = Vg_array[g,k,]
        #browser()
        pm_out = meta_analysis(x,v,sigma2=sigma2[g,k],meta_var='plug_in')
        X[g,k] = pm_out$mu_hat
        Vg[g,k] = pm_out$var_mu_hat
      }
    }
    Sigma = sigma2
  }else{

    if(est_sigma2){

      Sigma = matrix(nrow = G,ncol = K)

      #browser()

      if(meta_mode=='universal'){

        #browser()

        Sigma = matrix(PMmeta_array(X_array,Vg_array),nrow = G,ncol = K)

      }else if(meta_mode=='by_celltype'){

        for(k in 1:K){
          Sigma[,k] = PMmeta_matrix(X_array[,k,],Vg_array[,k,])
        }

      }else if(meta_mode=='by_gene'){

        for(g in 1:G){
          Sigma[g,] = PMmeta_matrix(X_array[g,,],Vg_array[g,,])
        }

      }else if(meta_mode=='local'){

        for(g in 1:G){
          for(k in 1:K){
            x = X_array[g,k,]
            v = Vg_array[g,k,]
            Sigma[g,k] = PMmeta_vector(x,v)
          }
        }
      }else if(meta_mode=="smooth"){
        X = apply(X_array,c(1,2),mean,na.rm=TRUE)
        Vg = t(apply(X_array,c(1),function(z){diag(cov(t(z),use = 'pairwise.complete.obs'))}))
        Sigma = pmax(Vg - apply(Vg_array,c(1,2),mean,na.rm=TRUE),0)
        for(k in 1:K){
          if(!all(is.na(X[,k]))){
            loess_fit = loess(y~.,
                              data.frame(x=X[,k],y=Sigma[,k]))
            Sigma[,k] = pmax(loess_fit$fitted,0)
          }
        }
      }else{
        stop('unsupported meta analysis method')
      }

      X = matrix(nrow = G,ncol = K)
      Vg = matrix(nrow = G,ncol = K)
      for(g in 1:G){
        for(k in 1:K){
          x = X_array[g,k,]
          v = Vg_array[g,k,]
          pm_out = meta_analysis(x,v,Sigma[g,k],meta_var=meta_var)
          X[g,k] = pm_out$mu_hat
          Vg[g,k] = pm_out$var_mu_hat
        }
      }

    }else{
      X = apply(X_array,c(1,2),mean,na.rm=TRUE)
      #browser()
      if(diag_cov){
        Vg = t(apply(X_array,c(1),function(z){diag(cov(t(z),use = 'pairwise.complete.obs'))}))/NI
        Sigma = pmax(Vg*NI - apply(Vg_array,c(1,2),mean,na.rm=TRUE),0)
      }else{
        Vg = t(apply(X_array,c(1),function(z){(cov(t(z),use = 'pairwise.complete.obs'))}))/NI
        Sigma=NULL
      }


    }

  }

  #browser()

  colnames(X) = cell_types
  #colnames(Vg) = cell_types


  # estimate cell size

  ## method 1: ols

  #browser()
  S = apply(rbind(S),2,mean,na.rm=T)

  return(list(X=X,Vg=Vg,Sigma=Sigma,S=S,
              X_array=X_array,Vg_array=Vg_array))

}



create_sce_from_counts = function(counts, colData, rowData = NULL) {
  if(is.null(rowData)) {
    sceset <- SingleCellExperiment(assays = list(counts = as.matrix(counts)),
                                   colData = colData)
  } else {
    sceset <- SingleCellExperiment(assays = list(counts = as.matrix(counts)),
                                   colData = colData,
                                   rowData = rowData)
  }
  # use gene names as feature symbols
  rowData(sceset)$feature_symbol <- rownames(sceset)
  # remove features with duplicated names
  if(is.null(rowData)) {
    sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
  }
  # QC
  isSpike(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
  sceset <- calculateQCMetrics(sceset, feature_controls = list("ERCC" = isSpike(sceset, "ERCC")))
  return(sceset)
}






#'@title void voild voild
#'@param ref.obj a list of singlecellexperiment object, each with raw counts, and colData cell_type and individual.
#'@param genes external gene annotation file
multi_ref_proc = function(ref.obj,bulk.obj = NULL,gene_length=NULL,protocols=NULL,marker_gene=NULL,
                          cell_types = c('alpha','acinar','beta','delta','ductal','gamma')){

  n_study = length(ref.obj)

  common_genes = set_data_decon(Y = ref.obj[[1]],cell_types=cell_types)$genes

  for(i in 2:n_study){

    datax = set_data_decon(Y = ref.obj[[i]],cell_types=cell_types)
    common_genes = intersect(common_genes,set_data_decon(Y = ref.obj[[i]],cell_types=cell_types)$genes)

  }


  if(!is.null(bulk.obj)){
    common_genes = intersect(common_genes,rownames(bulk.obj))
  }

  if(!is.null(genes)){
    common_genes = intersect(common_genes,genes$featurename)
    gene_length = (genes$featureend-genes$featurestart)[match(common_genes,genes$featurename)]
  }


  for(i in 1:n_study){
    gene_idx = match(common_genes,rownames(ref.obj[[i]]))
    ref.obj[[i]] = (ref.obj[[i]])[gene_idx,]
  }

  if(!is.null(bulk.obj)){
    gene_idx = match(common_genes,rownames(bulk.obj))
    bulk.obj = (bulk.obj)[gene_idx,]
    if(!is.null(genes)){
      counts(bulk.obj) = counts(bulk.obj)/gene_length*1e3
    }
    return(list(ref.obj = ref.obj,bulk.obj = bulk.obj,common_genes = common_genes))
  }else{
    return(ref.obj)
  }

}



