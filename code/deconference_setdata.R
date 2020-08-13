set_data_decon = function(y,Y,
                          ref_type = 'multi_sc',
                          #marker_gene = NULL,
                          cell_type_idx=NULL,indi_idx = NULL,
                          tau2 = NULL,sigma2 = NULL,
                          w=NULL,gene_thresh=10){

  y = cbind(y)

  # if(!is.null(marker_gene)){
  #   gene_idx = match(marker_gene,rownames(Y))
  #   omit.idx = which(is.na(gene_idx))
  #   if(length(omit.idx)>0){
  #     gene_idx = gene_idx[-omit.idx]
  #   }
  #
  #   #browser()
  #
  #   Y = Y[gene_idx,]
  #   if(!is.null(w)){
  #     w = w[gene_idx]
  #   }
  #   y = y[gene_idx,]
  #
  #   if(!is.null(sigma2)){
  #     sigma2 = sigma2[gene_idx,]
  #   }
  #   if(!is.null(tau2)){
  #     tau2  = tau2[gene_idx,]
  #   }
  # }



  if(ref_type == 'bulk'){

    rm.gene = which(rowSums(Y)==0)

    if(length(rm.gene)!=0){
      Y = Y[-rm.gene,]
      if(!is.null(w)){
        w = w[-rm.gene]
      }
      y = y[-rm.gene,]
    }

  }else{

    if(gene_thresh<1){
      gene_thresh = round(gene_thresh*ncol(Y))
    }

    rm.gene = which(rowSums(Y!=0)<gene_thresh)
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
      y = y[-rm.gene,]
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

  return(list(y=y,Y=Y,ref_type=ref_type,cell_type_idx=cell_type_idx,indi_idx=indi_idx,tau2=tau2,sigma2=sigma2,w=w))

}

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
scRef1_proc = function(Y,cell_type_idx,estimator='separate',tau2=NULL){

  G = nrow(Y)
  if(is.factor(cell_type_idx)){
    cell_types = levels(cell_type_idx)
  }else{
    cell_types = levels(as.factor(cell_type_idx))
  }

  K = length(cell_types)

  #browser()

  ## what if this individual have no cell type k?

  X = matrix(nrow=G,ncol=K)
  Vg = matrix(nrow=G,ncol=K)
  S = c()
  # need to estimate var(X)


  for(k in 1:K){
    cell_type = cell_types[k]
    cell_idx = which(cell_type_idx==cell_type)
    Nk = length(cell_idx)
    if(Nk==0){
      X[,k] = NA
      Vg[,k] = NA
      S[k] = NA
    }else{
      Yk = Y[,cell_idx,drop=FALSE]
      S[k] = sum(Yk)/Nk

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

  # if(estimator=='aggregate'){
  #
  #
  # }
  #
  # if(estimator=='separate'){
  #   for(k in 1:K){
  #     cell_idx = which(cell_type_idx==k)
  #     Nk = length(cell_idx)
  #     Yk = Y[,cell_idx]
  #     S[k] = sum(Yk)/Nk
  #     Xk = Yk%*%diag(c(1/colSums(Yk)))
  #     X[,k] = rowMeans(Xk)
  #     if(is.null(tau2)){
  #       Vg[,k] = 1/(Nk^2)*rowSums((Xk - X[,k,drop=FALSE]%*%t(rep(1,Nk)))^2)
  #     }else{
  #       Vg[,k] = tau2[,k]/Nk + X[,k]/Nk^2*sum(1/colSums(Yk))
  #     }
  #   }
  # }

  return(list(X=X,Vg=Vg,S=S,Sigma=NULL))
}

#'@title process single cell data from multiple subject.
#'@param eps a gene might have no expression in any cells of an individual so it's relative expression's standard error will be (0+eps)
#'@param est_sigma2 Indicate whether estiamte sigma^2, the variance of gene relative expresison across individuals. If yes, a PM method will be used; If not, will directly use sample variance-covariance matrix.
#'@param meta_var variance of hat{mu}, either 'plug_in' or 'adjust'
#'@param meta_mode 'universal': one sigma^2 for all X;'by_celltype': one sigma^2 for each cell type; 'local': one sigma^2 for each gene and cell type.
#'@return design matrix X and its variance matrix Vg
scRef_multi_proc = function(Y,cell_type_idx,indi_idx,estimator='separate',eps=0,
                            est_sigma2=TRUE,sigma2=NULL,tau2=NULL,meta_var='plug_in',meta_mode = "universal"){

  ##multiple individual single cell reference


  G = nrow(Y)
  K = length(unique(cell_type_idx))
  if(is.factor(indi_idx)){
    indis = levels(indi_idx)
  }else{
    indis = levels(as.factor(indi_idx))
  }

  if(is.factor(cell_type_idx)){
    cell_types = levels(cell_type_idx)
  }else{
    cell_types = levels(as.factor(cell_type_idx))
  }

  NI = length(indis)


  # first, for each individual, obtain X and Vg
  # then, for each g and k, perform PM method


  X_array = array(dim = c(G,K,NI))
  Vg_array = array(dim = c(G,K,NI))
  S_mat = matrix(nrow=NI,ncol=K)

  for(i in 1:NI){
    indi = indis[i]
    indi_cell_idx = which(indi_idx==indi)
    Yi = Y[,indi_cell_idx]
    cell_type_i = cell_type_idx[indi_cell_idx]
    indi_design.mat = scRef1_proc(Yi,cell_type_i,estimator=estimator,tau2=tau2)
    #browser()
    #print(dim(indi_design.mat$X))
    X_array[,,i] = indi_design.mat$X
    Vg_array[,,i] = indi_design.mat$Vg
    if(eps!=0){
      (Vg_array[,,i])[(Vg_array[,,i])==0] = eps
    }
    S_mat[i,] = indi_design.mat$S
  }

  # sigma2 is known
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
      Vg = t(apply(X_array,c(1),function(z){(cov(t(z),use = 'pairwise.complete.obs'))}))/NI
      Sigma=NULL
    }

  }

  #browser()

  colnames(X) = cell_types
  #colnames(Vg) = cell_types


  # estimate cell size

  ## method 1: ols

  S = colMeans(S_mat,na.rm = TRUE)
  S = S/S[1]

  ## method 2: glm

  # S_mat_dataframe = data.frame(y = c(S_mat),
  #                              indi = as.factor(rep(1:nrow(S_mat),ncol(S_mat))),
  #                              type = as.factor(rep(1:ncol(S_mat),each = nrow(S_mat))))
  # fit = glm.nb(y~.,S_mat_dataframe)

  return(list(X=X,Vg=Vg,Sigma=Sigma,S=S))

}
