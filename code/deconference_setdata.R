set_data_decon = function(y=NULL,Y,
                          ref_type = 'multi_sc',
                          marker_gene = NULL,
                          cell_type_idx=NULL,indi_idx = NULL,cell_types=NULL,
                          tau2 = NULL,sigma2 = NULL,
                          w=NULL,gene_thresh=0.05,max_count_quantile=0.99){

  if(!is.null(y)){
    y = cbind(y)
    common_genes = intersect(rownames(y),rownames(Y))
    y_geneidx = match(common_genes,rownames(y))
    Y_geneidx = match(common_genes,rownames(Y))
    y = y[y_geneidx,]
    Y = Y[Y_geneidx,]
  }


  if(!is.null(marker_gene)){
    if(!is.null(y)){
      marker_gene = intersect(marker_gene,common_genes)
      gene_idx = match(marker_gene,rownames(y))
      y = y[gene_idx,]
    }
    gene_idx = match(marker_gene,rownames(Y))
    # omit.idx = which(is.na(gene_idx))
    # if(length(omit.idx)>0){
    #   gene_idx = gene_idx[-omit.idx]
    # }

    #browser()

    Y = Y[gene_idx,]

  }


  if(length(w)==1){
    w = rep(nrow(Y),w)
  }


  if(ref_type == 'bulk'){

    rm.gene = which(rowSums(Y)==0)

    if(length(rm.gene)!=0){
      Y = Y[-rm.gene,]
      if(!is.null(w)){
        w = w[-rm.gene]
      }
      if(!is.null(y)){
        y = y[-rm.gene,]
      }

    }

  }else{

    if(is.null(cell_type_idx)&is.null(indi_idx)&class(Y) == "SingleCellExperiment"){
      cell_type_idx = Y$cell_type
      indi_idx = Y$individual
      Y = counts(Y)
    }

    if(is.null(cell_types)){
      if(is.factor(cell_type_idx)){
        cell_types = levels(cell_type_idx)
      }else{
        cell_types = levels(as.factor(cell_type_idx))
      }
    }
    K = length(cell_types)

    # pick up cells that match cell types we are interested in
    temp_idx = which(cell_type_idx%in%cell_types)
    Y = Y[,temp_idx]
    cell_type_idx = cell_type_idx[temp_idx]
    indi_idx = indi_idx[temp_idx]

    if(gene_thresh<1){
      gene_thresh = round(gene_thresh*ncol(Y))
    }

    rm.gene.low = which(rowSums(Y!=0)<gene_thresh)

    # remove union of genes that expressed more than max_count_quantle each cell type

    rm.gene.high = c()
    for(k in 1:K){

      cell_k_idx = which(cell_type_idx==cell_types[k])
      if(length(cell_k_idx)!=0){
        gene_counts = rowSums(Y[,cell_k_idx])
        rm.gene.high = c(rm.gene.high,which(gene_counts>quantile(gene_counts,max_count_quantile)))
      }

    }

    rm.gene = unique(c(rm.gene.low,rm.gene.high))

    if(length(rm.gene)!=0){
      Y = Y[-rm.gene,]
      if(!is.null(sigma2)){
        sigma2 = sigma2[-rm.gene,]
      }
      if(!is.null(tau2)){
        tau2  = tau2[-rm.gene,]
      }
      if(!is.null(w)){
        w = w[-rm.gene]
      }
      if(!is.null(y)){
        y = y[-rm.gene,]
      }
    }


    # remove cells without expression
    rm.cell = which(colSums(Y)==0)
    if(length(rm.cell)!=0){
      Y = Y[,-rm.cell]
      cell_type_idx = cell_type_idx[-rm.cell]
      if(!is.null(indi_idx)){
        indi_idx = indi_idx[-rm.cell]
      }
    }
  }

  return(list(y=y,Y=Y,ref_type=ref_type,
              cell_type_idx=cell_type_idx,indi_idx=indi_idx,
              tau2=tau2,sigma2=sigma2,w=w,cell_types=cell_types,
              genes = rownames(Y)))

}


###########################
###########################

#'@title process bulk reference data
#'@return design matrix X and its variance matrix Vg
bulkRef_proc = function(Y){
  X = apply(Y,2,function(x){x/sum(x)})
  U_inv = diag(c(1/colSums(Y)))
  Vg = X%*%U_inv
  return(list(X=X,Vg=Vg,Sigma=NULL))
}

#'@title process single cell reference data of one individual
#'@param estimator aggregate or separate
#'@return design matrix X and its variance matrix Vg, and estimate cell size
scRef1_proc = function(Y,cell_type_idx,estimator='separate',tau2=NULL,cell_types = NULL){

  G = nrow(Y)
  if(is.null(cell_types)){

    if(is.factor(cell_type_idx)){
      cell_types = levels(cell_type_idx)
    }else{
      cell_types = levels(as.factor(cell_type_idx))
    }

  }


  K = length(cell_types)

  #browser()

  ## what if this individual have no cell type k?

  X = matrix(nrow=G,ncol=K)
  Vg = matrix(nrow=G,ncol=K)
  S = c()
  S_var = c()
  # need to estimate var(X)


  #browser()
  for(k in 1:K){
    cell_type = cell_types[k]
    cell_idx = which(cell_type_idx==cell_type)
    Nk = length(cell_idx)
    if(Nk==0){
      X[,k] = NA
      Vg[,k] = NA
      S[k] = NA
      S_var[k] = NA
    }else{
      Yk = Y[,cell_idx,drop=FALSE]
      S[k] = sum(Yk)/Nk
      S_var[k] = var(colSums(Yk))/Nk

      if(estimator=="aggregate"){

        X[,k] = rowSums(Yk)/sum(Yk)
        if(is.null(tau2)){
          Vg[,k] = 1/(sum(Yk))^2*rowSums((Yk-t(c(colSums(Yk))*t(X[,k,drop=FALSE]%*%t(rep(1,Nk)))))^2)
          #Vg[,k] = 1/(sum(Yk))^2*rowSums((Yk-X[,k,drop=FALSE]%*%t(rep(1,Nk))%*%diag(colSums(Yk)))^2)
        }else{
          Vg[,k] = X[,k]/sum(Yk)+tau2[,k]*sum(colSums(Yk)^2)/(sum(Yk))^2
        }

      }else if(estimator=='separate'){
        Xk = apply(Yk,2,function(z){z/sum(z)})
        # Xk = Yk%*%diag(c(1/colSums(Yk)))
        X[,k] = rowMeans(Xk)
        if(is.null(tau2)){
          Vg[,k] = 1/(Nk^2)*rowSums((Xk - X[,k,drop=FALSE]%*%t(rep(1,Nk)))^2)
        }else{
          Vg[,k] = tau2[,k]/Nk + X[,k]/Nk^2*sum(1/colSums(Yk))
        }
      }

    }
  }

  colnames(X) = cell_types
  colnames(Vg) = cell_types

  return(list(X=X,Vg=Vg,S=S,S_var=S_var,Sigma=NULL))
}

#'@title process single cell data from multiple subject.
#'@param eps a gene might have no expression in any cells of an individual so it's relative expression's standard error will be (0+eps)
#'@param est_sigma2 Indicate whether estiamte sigma^2, the variance of gene relative expresison across individuals. If yes, a PM method will be used; If not, will directly use sample variance-covariance matrix(if diag_cov, then use the diagnal of sample cov).
#'@param meta_var variance of hat{mu}, either 'plug_in' or 'adjust'
#'@param meta_mode 'universal': one sigma^2 for all X;'by_celltype': one sigma^2 for each cell type; 'local': one sigma^2 for each gene and cell type; 'smooth': fit a local polynomial to mean-var per cell type.
#'@return design matrix X and its variance matrix Vg
scRef_multi_proc = function(Y,cell_type_idx,indi_idx,estimator='separate',eps=0,
                            est_sigma2=TRUE,sigma2=NULL,tau2=NULL,meta_var='plug_in',
                            meta_mode = "universal",
                            cell_types = NULL,indis=NULL,diag_cov=FALSE){

  ##multiple individual single cell reference


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
  S_mat_var = matrix(nrow=NI,ncol=K)

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
    S_mat_var[i,] = indi_design.mat$S_var
  }


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

  # browser()

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
  fit = try(MASS::glm.nb(y~.,S_mat_dataframe),silent = TRUE)
  if(class(fit)=='try-error'){

    fit = glm(y~.,S_mat_dataframe,family = 'poisson')

  }

  S_glm = S
  S_glm[which(!is.nan(S))] = c(1,exp(fit$coefficients[-c(1:nrow(S_mat))]))
  names(S_glm) = cell_types

  return(list(X=X,Vg=Vg,Sigma=Sigma,
              S=S,S_mat=S_mat,S_glm=S_glm,
              X_array=X_array,Vg_array=Vg_array))

}


getXV_array = function(Y,cell_type_idx=NULL,indi_idx=NULL,estimator='separate',eps=0,
                       tau2=NULL,cell_types = NULL,indis=NULL){

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
  suppressWarnings(fit = try(MASS::glm.nb(y~.,S_mat_dataframe),silent = TRUE))

  suppressWarnings(if(class(fit)=='try-error'){

    fit = glm(y~.,S_mat_dataframe,family = 'poisson')

  })

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






