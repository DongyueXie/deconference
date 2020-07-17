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
#'@return design matrix X and its variance matrix Vg
scRef1_proc = function(Y,cell_type_idx,estimator='separate',tau2=NULL){

  G = nrow(Y)
  K = length(unique(cell_type_idx))

  #browser()

  X = matrix(nrow=G,ncol=K)
  Vg = matrix(nrow=G,ncol=K)
  # need to estimate var(X)

  if(estimator=='aggregate'){
    for(k in 1:K){
      cell_idx = which(cell_type_idx==k)
      Nk = length(cell_idx)
      Yk = Y[,cell_idx]
      X[,k] = rowSums(Yk)/sum(Yk)
      if(is.null(tau2)){
        Vg[,k] = 1/(sum(Yk))^2*rowSums((Yk-X[,k,drop=FALSE]%*%t(rep(1,Nk))%*%diag(colSums(Yk)))^2)
      }else{
        Vg[,k] = X[,k]/sum(Yk)+tau2[,k]*sum(colSums(Yk)^2)/(sum(Yk))^2
      }
    }

  }

  if(estimator=='separate'){
    for(k in 1:K){
      cell_idx = which(cell_type_idx==k)
      Nk = length(cell_idx)
      Yk = Y[,cell_idx]
      Xk = Yk%*%diag(c(1/colSums(Yk)))
      X[,k] = rowMeans(Xk)
      if(is.null(tau2)){
        Vg[,k] = 1/(Nk^2)*rowSums((Xk - X[,k,drop=FALSE]%*%t(rep(1,Nk)))^2)
      }else{
        Vg[,k] = tau2[,k]/Nk + X[,k]/Nk^2*sum(1/colSums(Yk))
      }
    }
  }

  return(list(X=X,Vg=Vg,Sigma=NULL))
}

#'@title process single cell data from multiple subject.
#'@param eps a gene might have no expression in any cells of an individual so it's relative expression's standard error will be (0+eps)
#'@param est_sigma2 Indicate whether estiamte sigma^2, the variance of gene relative expresison across individuals. If yes, a PM method will be used; If not, will directly use sample variance-covariance matrix.
#'@return design matrix X and its variance matrix Vg
scRef_multi_proc = function(Y,cell_type_idx,indi_idx,estimator='separate',eps=0,
                            est_sigma2=TRUE,sigma2=NULL,tau2=NULL,meta_var='plug_in'){

  ##multiple individual single cell reference


  G = nrow(Y)
  K = length(unique(cell_type_idx))
  NI = length(unique(indi_idx))


  # first, for each individual, obtain X and Vg
  # then, for each g and k, perform PM method


  X_array = array(dim = c(G,K,NI))
  Vg_array = array(dim = c(G,K,NI))

  for(i in 1:NI){
    indi_cell_idx = which(indi_idx==i)
    Yi = Y[,indi_cell_idx]
    cell_type_i = cell_type_idx[indi_cell_idx]
    indi_design.mat = scRef1_proc(Yi,cell_type_i,estimator=estimator,tau2=tau2)
    #browser()
    #print(dim(indi_design.mat$X))
    X_array[,,i] = indi_design.mat$X
    Vg0 = indi_design.mat$Vg
    Vg0[Vg0==0] = eps
    Vg_array[,,i] = Vg0
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
        pm_out = PMmeta(x,v,sigma2=sigma2[g,k])
        X[g,k] = pm_out$mu_hat
        Vg[g,k] = pm_out$var_mu_hat
      }
    }
    Sigma = sigma2
  }else{

    if(est_sigma2){
      X = matrix(nrow = G,ncol = K)
      Vg = matrix(nrow = G,ncol = K)
      Sigma = matrix(nrow = G,ncol = K)

      for(g in 1:G){
        for(k in 1:K){
          x = X_array[g,k,]
          v = Vg_array[g,k,]
          pm_out = PMmeta(x,v,meta_var=meta_var)
          X[g,k] = pm_out$mu_hat
          Vg[g,k] = pm_out$var_mu_hat
          Sigma[g,k] = pm_out$sigma2
        }
      }

    }else{
      X = apply(X_array,c(1,2),mean)
      #browser()
      Vg = t(apply(X_array,c(1),function(z){(cov(t(z)))}))/NI
      Sigma=NULL
    }

  }

  return(list(X=X,Vg=Vg,Sigma=Sigma))

}
