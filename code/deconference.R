#'@title deconvolution inference
#'@description There are three cases: 1. bulk reference data, then input Y (gene by cell type);
#'2. single-cell reference data from the same individual, then input Y (gene by cell) and cell_type_idx indicator;
#'3. single-cell reference data from multiple-subject, then input Y (gene by cell), cell_type_idx and individual index
#'@param y a vector of bulk sample; or matrix whose columns are bulk samples
#'@param Y a reference count matrix, bulk ref: gene by cell type; single cell ref: gene by cells
#'@param ref_type bulk, sc, multi_sc
#'@param cell_type_idx a numerical vector indicates cell types of Y if using single cell reference.
#'@param indi_idx a vector indicates individuals
#'@param tau2 varaince of gene expression across single cells, gene by cell type matrix; if NULL, will be estimated.
#'@param sigma2 variance of gene expression across individuals, gene by cell type matrix.
#'@param w gene weights, length(w) = number of genes. will be adjusted to sum to G.
#'@param x_estimator separate or aggregate
#'@param est_sigma2 whether estimate sigma^2, the variance of X across individuals
#'@param meta_var 'plug_in', or 'adjust'
#'@param adjust whether adjust for uncertainty in reference matrix
#'@param alpha significance level
#'@param correction whether perform fuller's small sample correction.
#'@param a alpha in the Fuller's small sample correction
#'@param eps adjust of zero variane if a gene has no expression observed in one cell type
#'@param gene_thresh remove genes that appear in less than number of  cells
#'@return a list from estimation_func


deconference = function(data.obj,
                        marker_gene = NULL,
                        x_estimator = 'separate',
                        est_sigma2 = FALSE,
                        meta_var = 'plug_in',
                        a=0, correction=FALSE,eps=0){

  ref_type = data.obj$ref_type
  w = data.obj$w
  y = data.obj$y
  Y = data.obj$Y

  #browser()

  if(ref_type=='bulk'){

    design.mat = bulkRef_proc(Y)

  }else{

    cell_type_idx = data.obj$cell_type_idx
    tau2 = data.obj$tau2

    if(ref_type=='sc'){

    design.mat = scRef1_proc(Y,cell_type_idx,estimator=x_estimator,tau2=tau2)
   }else if(ref_type=='multi_sc'){
     indi_idx = data.obj$indi_idx
     sigma2 = data.obj$sigma2

    # multiple individual single cell reference samples, estimate sigma^2
    design.mat = scRef_multi_proc(Y,cell_type_idx,indi_idx,estimator=x_estimator,tau2=tau2,
                                  sigma2=sigma2,est_sigma2 = est_sigma2,eps=eps,meta_var=meta_var)
   }else{
     stop("unspported reference type")
   }

  }

  #browser()

  out = estimation_func(y=y,X=design.mat$X,Vg=design.mat$Vg,design.mat$Sigma,marker_gene=marker_gene,
                        w=w,a=a,correction=correction)
  return(out)
}


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


#'@title Q-statistics - (n-1) in Meta-analysis
#'@param sigma2 random effect variance
#'@param x effects of studies, a vector
#'@param v variance of x, a vector
#'@return Q value - (n-1)

Q_stats = function(sigma2,x,v){
  sum(1/(sigma2+v)*(x-(sum(x/(v+sigma2)))/(sum(1/(sigma2+v))))^2) - (length(x)-1)
}

#'@title PM meta-analysis random effecr estimator
#'@param x effects of studies, a vector
#'@param v variance of x, a vector
#'@param ub upper bound of random effect variance sigma^2
#'@return estimated random effect variance sigma^2, mu_hat and var(mu_hat)
PMmeta = function(x,v,ub=0.25,sigma2=NULL,meta_var='plug_in'){

  l = length(x)
  # do not use v=0
  rm.v = which(v==0)


  if(length(rm.v)<(l-1)){
    if(length(rm.v)!=0){
      x = x[-rm.v]
      v = v[-rm.v]
    }

    if(is.null(sigma2)){
      Q0 = Q_stats(0,x,v)
      Qub = Q_stats(ub,x,v)
      if(Q0>0){
        if(Qub<0){
          sigma2 = uniroot(Q_stats,c(0,ub),x=x,v=v)$root
        }else{
          sigma2 = ub
        }
      }else{
        sigma2 = 0
      }
    }


    w = 1/(v+sigma2)
    w = w/sum(w)
    mu_hat = sum(w*x)
    #browser()
    if(meta_var=='adjust'){
      var_mu_hat = sum(w^2/(1-w)*(x-mu_hat)^2)
    }
    if(meta_var=='plug_in'){
      var_mu_hat = 1/sum(1/(sigma2+v))
    }

    return(list(sigma2=sigma2,mu_hat=mu_hat,var_mu_hat=var_mu_hat))
  }else{

    if(is.null(sigma2)){
      return(list(sigma2=0,mu_hat=mean(x),var_mu_hat=mean(v)/l))
    }else{
      return(list(sigma2=sigma2,mu_hat=mean(x),var_mu_hat=mean(v+sigma2)/l))
    }


  }


}

#'@title unadjusted method, either using ols or sandwich for varaince estimation
unadjusted_lm = function(y,X,w=NULL){
  G = nrow(X)
  K = ncol(X)
  y = cbind(y)
  nb = ncol(y)

  if(is.null(w)){
    w  = 1/rowSums(X)
  }
  w = w/sum(w)*G

  Xw = X*sqrt(w)
  yw = cbind(y*sqrt(w))
  A = t(Xw)%*%Xw
  A_inv = solve(A)
  Hat_mat = Xw%*%A_inv%*%t(Xw)
  beta_tilde_hat = pmax(A_inv%*%t(Xw)%*%yw,0)
  beta_tilde_hat = cbind(beta_tilde_hat)
  beta_hat = apply(beta_tilde_hat,2,function(z){z/sum(z)})

  J = matrix(0,nrow=(nb*K),ncol=(nb*K))
  for(i in 1:nb){
    J[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = J_sum2one(beta_tilde_hat[,i],K)
  }

  # perform ols estimator of variance
  resid_var = t(yw)%*%(diag(G)-Hat_mat)%*%yw/(G-K)
  covb = kronecker(resid_var,A_inv)


  asyV = (J)%*%covb%*%t(J)


  beta_se = sqrt(diag(asyV))
  beta_se = matrix(beta_se,ncol=nb)

  ols.out = list(cov_beta_tilde_hat = covb,
                 beta_se = beta_se,
                 cov_beta_hat = asyV)

  # perform sandwich estimator of variance

  Sigma = matrix(0,nrow=nb*K,ncol=nb*K)
  Q_inv = matrix(0,nrow=nb*K,ncol=nb*K)
  Sigma_ii = matrix(0,nrow=nb*K,ncol=nb*K)

  for(i in 1:nb){
    Q_inv[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = A_inv
    for(j in i:nb){

      Sigma_ij = crossprod(Xw*c(Xw%*%beta_tilde_hat[,i,drop=FALSE])-Xw*c(yw[,i]),
                           Xw*c(Xw%*%beta_tilde_hat[,j,drop=FALSE])-Xw*c(yw[,j]))

      Sigma[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)] = Sigma_ij
      if(j==i){
        Sigma_ii[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = Sigma_ij
      }
    }
  }
  Sigma = Sigma+t(Sigma)-Sigma_ii
  covb = Q_inv%*%Sigma%*%Q_inv

  asyV = (J)%*%covb%*%t(J)


  beta_se = sqrt(diag(asyV))
  beta_se = matrix(beta_se,ncol=nb)

  sand.out = list(cov_beta_tilde_hat = covb,
                  beta_se = beta_se,
                  cov_beta_hat = asyV)








  return(list(beta_tilde_hat=beta_tilde_hat,
              beta_hat=beta_hat,
              sand.out=sand.out,
              ols.out=ols.out))
}



#'@param y a vector of bulk sample
#'@param X reference matrix, estimated relative expression, gene by cell
#'@param Vg variance matrix of X
#'@param Sigma variance of true X among individuals
#'@param marker_gene
#'@param w gene weights
#'@param alpha significance level
#'@param a alpha in the Fuller's small sample correction
#'@param correction whether perform fuller's small sample correction.


estimation_func = function(y,X,Vg,Sigma=NULL,marker_gene = NULL,w=NULL,a=0,correction=FALSE){

  #browser()


  if(!is.null(marker_gene)){
    gene_idx = match(marker_gene,rownames(X))
    omit.idx = which(is.na(gene_idx))
    if(length(omit.idx)>0){
      gene_idx = gene_idx[-omit.idx]
    }

    #browser()

    X = X[gene_idx,]
    Vg = Vg[gene_idx,]
    if(!is.null(w)){
      w = w[gene_idx]
    }
    y = y[gene_idx,]

    if(!is.null(Sigma)){
      Sigma = Sigma[gene_idx,]
    }

  }

  G = nrow(X)
  K = ncol(X)
  nb = ncol(y)
  if(is.null(nb)){
    nb=1
  }


  if(is.null(w)){
    if(is.null(Sigma)){
      w  = 1/rowSums(X)
    }else{
      w  = 1/(rowSums(X)/K+rowSums(Sigma)/K^2)
    }

  }
  w = w/sum(w)*G

  input = list(X=X,Vg=Vg,y=y,w=w,Sigma=Sigma)

  Xw = X*sqrt(w)
  yw = cbind(y*sqrt(w))
  A = t(Xw)%*%Xw

  if(is.matrix(Vg)){
    Vgw = Vg*w
    if(ncol(Vg)==K^2){
      V = matrix(c(colSums(Vgw)),ncol = K)
    }else if(ncol(Vg)==K){
      V = diag(c(colSums(Vgw)))
    }else{
      stop('check dimension of Vg')
    }
  }else if(is.array(Vg)){
    Vgw = Vg*rep(w,each=K*K)
    V = rowSums(Vgw,dims = 2)
  }else{
    stop("Vg input type not allowed")
  }


  ## Fuller's correction for negative matrix

  if(correction){
    M = crossprod(cbind(yw,Xw))
    B = diag(c(1/colSums(yw),1/diag(V)))

    # M = A
    # B = diag(1/diag(V))

    #print((B))
    #lambda = min(eigen(A%*%diag(1/colSums(Yuu)))$values)
    #print(eigen(A%*%U)$values)
    lambda = min(eigen(M%*%B)$values)
    #print(lambda)
    #print(paste('est',lambda))
    #print(paste('value',det(A-lambda*U_inv)))
    if(lambda>(1+1/G)){
      A = A - (1-a/G)*V
    }else{
      A = A - (lambda-1/G-a/G)*V
    }
  }else{
    A = A - V
  }
  A_inv = solve(A)

  beta_tilde_hat = pmax(A_inv%*%t(Xw)%*%yw,0)

  #browser()

  Sigma = matrix(0,nrow=nb*K,ncol=nb*K)
  Q_inv = matrix(0,nrow=nb*K,ncol=nb*K)
  Sigma_ii = matrix(0,nrow=nb*K,ncol=nb*K)

  #browser()



  for(i in 1:nb){
    Q_inv[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = A_inv
    for(j in i:nb){

      if(is.matrix(Vg)){
        if(ncol(Vg)==K^2){
          Vbi = t(apply(Vgw,1,function(z){v = matrix(z,ncol=K);v%*%beta_tilde_hat[,i]}))
          Vbj = t(apply(Vgw,1,function(z){v = matrix(z,ncol=K);v%*%beta_tilde_hat[,j]}))
        }
        if(ncol(Vg)==K){
          Vbi = (Vgw)%*%diag(c(beta_tilde_hat[,i]))
          Vbj = (Vgw)%*%diag(c(beta_tilde_hat[,j]))
        }
      }else if(is.array(Vg)){
        Vbi = t(apply(Vgw,3,function(z){z%*%beta_tilde_hat[,i]}))
        Vbj = t(apply(Vgw,3,function(z){z%*%beta_tilde_hat[,j]}))
      }


      Sigma_ij = crossprod(Xw*c(Xw%*%beta_tilde_hat[,i,drop=FALSE])-Vbi-Xw*c(yw[,i]),
                           Xw*c(Xw%*%beta_tilde_hat[,j,drop=FALSE])-Vbj-Xw*c(yw[,j]))
      Sigma[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)] = Sigma_ij
      if(j==i){
        Sigma_ii[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = Sigma_ij
      }
    }
  }
  Sigma = Sigma+t(Sigma)-Sigma_ii
  covb = Q_inv%*%Sigma%*%Q_inv



  # # only one bulk data
  # if(is.null(nb) | nb==1){
  #   #asymptotic variance
  #   Sigma = crossprod(Xw*c(Xw%*%beta_tilde_hat)-(Vgw)%*%diag(beta_tilde_hat)-Xw*c(yw))
  #   #Sigma = Sigma/G
  #   #Q_inv = G*A_inv
  #   #covb = Q_inv%*%Sigma%*%Q_inv/G
  #   covb = A_inv%*%Sigma
  # }else{
  #
  #
  # }


  beta_tilde_hat = cbind(beta_tilde_hat)


  # delta method

  # covb is a (nb*K)*(nb*K) cov matrix
  # formulate Jacobian matrix

  J = matrix(0,nrow=(nb*K),ncol=(nb*K))
  for(i in 1:nb){
    J[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = J_sum2one(beta_tilde_hat[,i],K)
  }

  asyV = (J)%*%covb%*%t(J)

  beta_hat = apply(beta_tilde_hat,2,function(z){z/sum(z)})
  beta_se = sqrt(diag(asyV))
  beta_se = matrix(beta_se,ncol=nb)

  return(list(beta_tilde_hat=beta_tilde_hat,
              beta_hat=beta_hat,
              beta_se=beta_se,
              cov_beta_hat = asyV,
              input = input))



  # if(is.null(nb) | nb==1){
  #   J = J_sum2one(beta_tilde_hat,K)
  #
  #   asyV = (J)%*%covb%*%t(J)
  #
  #   beta_hat = beta_tilde_hat/sum(beta_tilde_hat)
  #   beta_se = sqrt(diag(asyV))
  #
  #   ci = rbind(pmax(0,beta_hat - qnorm(1-alpha/2)*beta_se), pmin(1,beta_hat + qnorm(1-alpha/2)*beta_se))
  #
  #   return(list(beta_tilde_hat=beta_tilde_hat,
  #               beta_hat=beta_hat,
  #               beta_se=beta_se,
  #               ci=ci,
  #               cov_beta_hat = asyV,
  #               input = input))
  # }else{
  #
  # }


}



#'@title Jacobian matrix of sum-to-1 scale function
#'@param b beta_tilde_hat
#'@param K length of beta_tilde_hat
#'@return Jacobian matrix
J_sum2one = function(b,K){
  J = - (b)%*%t(rep(1,K))
  diag(J) = diag(J) + sum(b)
  J = J/sum(b)^2
  J
}

