


#'@param datax a list of singlecellexperiment objects.
deconference_multi_ref = function(ref.obj,bulk.obj,tau2=NULL,cell_types=NULL,sigma2=NULL,
                          est_sigma2=TRUE,meta_var='adjust',meta_mode='smooth',correction=FALSE,cellsize_est='glm',
                          marker_gene = NULL,
                          hc.type = 'hc3',w = 1){

  library(abind)
  n_ref = length(ref.obj)

  all_X_array = c()
  all_Vg_array = c()

  S_olss = c()
  S_glms = c()
  for(i in 1:n_ref){

    out_array = getXV_array(ref.obj[[i]],tau2=tau2,cell_types = cell_types,indis=NULL)
    all_X_array = abind(all_X_array,out_array$X_array)
    all_Vg_array = abind(all_Vg_array,out_array$Vg_array)
    S_olss = rbind(S_olss,out_array$S_ols)
    S_glms = rbind(S_glms,out_array$S_glm)
  }

  design.mat = getXV(all_X_array,all_Vg_array,S_olss=S_olss,S_glms=S_glms,sigma2=sigma2,
                   est_sigma2=est_sigma2,meta_var=meta_var,meta_mode=meta_mode,cell_types=cell_types)


  if(cellsize_est=='ols'){
    out = estimation_func(y=counts(bulk.obj),X=design.mat$X,Vg=design.mat$Vg,design.mat$Sigma,marker_gene=marker_gene,
                          w=w,hc.type=hc.type,correction=correction,S=design.mat$S_ols)
  }
  if(cellsize_est=='glm'){
    out = estimation_func(y=counts(bulk.obj),X=design.mat$X,Vg=design.mat$Vg,design.mat$Sigma,marker_gene=marker_gene,
                          w=w,hc.type=hc.type,correction=correction,S=design.mat$S_glm)
  }

  return(out)

}




#'@param ref.obj a list of singlecellexperiment object, each with cell_type and individual
#'@param genes external gene annotation file
multi_ref_proc = function(ref.obj,bulk.obj = NULL,genes=NULL,
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



getXV_array = function(Y,cell_type_idx=NULL,indi_idx=NULL,estimator='separate',eps=0,
                       tau2=NULL,cell_types = NULL,indis=NULL){

  #browser()

  if(is.null(cell_type_idx)&is.null(indi_idx)&class(Y) == "SingleCellExperiment"){
    cell_type_idx = Y$cell_type
    indi_idx = Y$individual
    Y = counts(Y)
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
    indi_design.mat = scRef1_proc(Yi,cell_type_i,estimator=estimator,tau2=tau2,cell_types=cell_types)
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


  S_mat = round(S_mat)
  rownames(S_mat) = indis
  colnames(S_mat) = cell_types
  S = colMeans(S_mat,na.rm = TRUE)
  S = S/S[1]

  # ## method 2: glm
  #
  S_mat_dataframe = data.frame(y = c(S_mat),
                               indi = factor(rep(indis,ncol(S_mat))),
                               type = factor(rep(cell_types,each = nrow(S_mat)),levels = cell_types))
  suppressWarnings({fit = try(MASS::glm.nb(y~.,S_mat_dataframe),silent = TRUE)})

  suppressWarnings({if(class(fit)=='try-error'){

    fit = glm(y~.,S_mat_dataframe,family = 'poisson')

  }})

  S_glm = S
  S_glm[which(!is.nan(S))] = c(1,exp(fit$coefficients[-c(1:nrow(S_mat))]))
  names(S_glm) = cell_types

  return(list(X_array = X_array,Vg_array = Vg_array,S_ols=S,S_glm=S_glm))

}


getXV = function(X_array,Vg_array,S_olss=NULL,S_glms=NULL,sigma2=NULL,
                 est_sigma2=TRUE,meta_var='adjust',meta_mode='smooth',diag_cov=FALSE,cell_types){

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

  if(!is.null(S_glms)){

    S_glm = apply(rbind(S_glms),2,mean,na.rm=T)
    S_ols = apply(rbind(S_olss),2,mean,na.rm=T)

    return(list(X=X,Vg=Vg,Sigma=Sigma,
                S_ols=S_ols,S_glm=S_glm,
                X_array=X_array,Vg_array=Vg_array))
  }else{
    return(list(X=X,Vg=Vg,Sigma=Sigma,
                X_array=X_array,Vg_array=Vg_array))
  }



}







