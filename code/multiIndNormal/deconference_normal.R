#'@title deconvolution inference
#'@description There are three cases: 1. bulk reference data, then input Y (gene by cell type);
#'2. single-cell reference data from the same individual, then input Y (gene by cell) and cell_type indicator;
#'3. single-cell reference data from multiple-subject, then input Y (gene by cell), cell_type and individual index
#'@param y a vector of bulk sample
#'@param Y a reference count matrix, bulk ref: gene by cell type; single cell ref: gene by cells
#'@param cell_type a numerical vector indicates cell types of Y if using single cell reference.
#'@param indi_idx a vector indicate individuals
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
#'@return a list from adjusted_estimator


deconference = function(y,Y,
                        cell_type=NULL,
                        indi_idx = NULL,
                        tau2 = NULL,
                        sigma2 = NULL,
                        w=NULL,
                        x_estimator = 'separate',
                        est_sigma2 = TRUE,
                        meta_var = 'plug_in',
                        adjust=TRUE,
                        alpha=0.05,
                        a=0, correction=FALSE,eps=1e-5){

  G = length(y)
  if(is.null(cell_type)){
    # reference samples are bulk data
    design.mat = bulkRef_proc(Y)
  }else if(is.null(indi_idx)){
    # same individual single cell reference
    design.mat = scRef1_proc(Y,cell_type,estimator=x_estimator,tau2=tau2)
  }else if(est_sigma2){
    # multiple individual single cell reference samples, estimate sigma^2
    design.mat = scRef_multi_proc(Y,cell_type,indi_idx,estimator=x_estimator,
                                  sigma2=sigma2,est_sigma2 = est_sigma2,eps=eps,meta_var=meta_var)
  }else{
    # multiple individual single cell reference samples, do not estimate sigma^2
    design.mat = scRef_multi_proc(Y,cell_type,indi_idx,estimator=x_estimator,est_sigma2 = FALSE,eps=eps)
  }

  #print(colSums(Y))
  #browser()
  out = adjusted_estimator(y,design.mat$X,design.mat$Vg,w,adjust,alpha,a,correction)
  if(est_sigma2){out$sigma2 = design.mat$Sigma}
  return(out)
}

#'@title process bulk reference data
#'@return design matrix X and its variance matrix Vg
bulkRef_proc = function(Y){
  X = apply(Y,2,function(x){x/sum(x)})
  U_inv = diag(c(1/colSums(Y)))
  Vg = X%*%U_inv
  return(list(X=X,Vg=Vg))
}

#'@title process single cell reference data of one individual
#'@param estimator aggregate or separate
#'@return design matrix X and its variance matrix Vg
scRef1_proc = function(Y,cell_type,estimator='separate',tau2=NULL){

  G = nrow(Y)
  K = length(unique(cell_type))

  rm.idx = which(colSums(Y)==0)
  if(length(rm.idx)!=0){
    Y = Y[,-rm.idx]
    cell_type = cell_type[-rm.idx]
  }

  #browser()

  X = matrix(nrow=G,ncol=K)
  Vg = matrix(nrow=G,ncol=K)
  # need to estimate var(X)

  if(estimator=='aggregate'){
    for(k in 1:K){
      cell_idx = which(cell_type==k)
      Nk = length(cell_idx)
      Yk = Y[,cell_idx]
      X[,k] = rowMeans(Yk)
      if(is.null(tau2)){
        Vg[,k] = 1/(Nk)^2*rowSums((Yk-X[,k,drop=FALSE]%*%t(rep(1,Nk)))^2)
      }else{
        Vg[,k] = (tau2+sc_noise_var)/nk
      }
    }

  }

  if(estimator=='separate'){
    for(k in 1:K){
      cell_idx = which(cell_type==k)
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

  return(list(X=X,Vg=Vg))
}

#'@title process single cell data from multiple subject.
#'@param eps a gene might have no expression in any cells of an individual so it's relative expression's standard error will be (0+eps)
#'@param est_sigma2 Indicate whether estiamte sigma^2, the variance of gene relative expresison across individuals. If yes, a PM method will be used; If not, will directly use sample variance-covariance matrix.
#'@return design matrix X and its variance matrix Vg
scRef_multi_proc = function(Y,cell_type,indi_idx,estimator='separate',eps=1e-5,
                            est_sigma2=TRUE,sigma2=NULL,meta_var='plug_in'){

  ##multiple individual single cell reference


  G = nrow(Y)
  K = length(unique(cell_type))
  NI = length(unique(indi_idx))

  rm.idx = which(colSums(Y)==0)
  if(length(rm.idx)!=0){
    Y = Y[,-rm.idx]
    cell_type = cell_type[-rm.idx]
    indi_idx = indi_idx[-rm.idx]
  }

  # first, for each individual, obtain X and Vg
  # then, for each g and k, perform PM method


  X_array = array(dim = c(G,K,NI))
  Vg_array = array(dim = c(G,K,NI))

  for(i in 1:NI){
    indi_cell_idx = which(indi_idx==i)
    Yi = Y[,indi_cell_idx]
    cell_type_i = cell_type[indi_cell_idx]
    indi_design.mat = scRef1_proc(Yi,cell_type_i,estimator=estimator)
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
        pm_out = PMmeta(x,v,sigma2=sigma2[g,k])
        X[g,k] = pm_out$mu_hat
        Vg[g,k] = pm_out$var_mu_hat
      }
    }
    return(list(X=X,Vg=Vg))
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
      return(list(X=X,Vg=Vg,Sigma=Sigma))
    }else{
      X = apply(X_array,c(1,2),mean)
      Vg = t(apply(X_array,c(1),function(z){diag(cov(t(z)))}))/NI
      return(list(X=X,Vg=Vg))
    }

  }

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
PMmeta = function(x,v,ub=10,sigma2=NULL,meta_var='plug_in'){

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
  if(meta_var=='adjust'){
    var_mu_hat = sum(w^2/(1-w)*(x-mu_hat)^2)
  }
  if(meta_var=='plug_in'){
    var_mu_hat = 1/sum(1/(sigma2+v))
  }

  return(list(sigma2=sigma2,mu_hat=mu_hat,var_mu_hat=var_mu_hat))
}




#'@param y a vector of bulk sample
#'@param X reference matrix, estimated relative expression, gene by cell
#'@param Vg variance matrix of X
#'@param w gene weights
#'@param adjust whether adjust for uncertainty in refence matrix
#'@param alpha significance level
#'@param a alpha in the Fuller's small sample correction
#'@param correction whether perform fuller's small sample correction.

adjusted_estimator = function(y,X,Vg,w=NULL,
                              adjust=TRUE,
                              alpha=0.05,
                              a=0, correction=FALSE){

  G = nrow(X)
  K = ncol(X)
  if(is.null(w)){
    w  = rep(1,G)
  }else{
    w = w/sum(w)*G
  }

  Xw = X*sqrt(w)
  yw = y*sqrt(w)
  Vgw = Vg*w
  A = t(Xw)%*%Xw
  V = diag(c(colSums(Vgw)))
  if(adjust){
    ## Fuller's correction for negative matrix
    if(correction){
      M = crossprod(cbind(yw,Xw))
      B = diag(c(1/sum(w*y),1/diag(V)))
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
    #browser()
    beta_tilde_hat = pmax(c(A_inv%*%t(Xw)%*%yw),0)
    #asymptotic variance
    Sigma = crossprod(Xw*c(Xw%*%beta_tilde_hat)-(Vgw)%*%diag(beta_tilde_hat)-Xw*c(yw))
    Sigma = Sigma/G
    Q_inv = G*A_inv
    covb = Q_inv%*%Sigma%*%Q_inv/G
  }else{
    #print(y);print(X);print(w)

    #c
    lm_fit = lm(y~.-1,data.frame(y=y,X=X),weights = w)

    beta_tilde_hat = pmax(as.numeric(coefficients(lm_fit)),0)
    covb = vcov(lm_fit)
  }


  # delta method
  J = - (beta_tilde_hat)%*%t(rep(1,K))
  diag(J) = diag(J) + sum(beta_tilde_hat)
  J = J/sum(beta_tilde_hat)^2
  asyV = (J)%*%covb%*%t(J)

  beta_hat = beta_tilde_hat/sum(beta_tilde_hat)
  beta_se = sqrt(diag(asyV))

  ci = rbind(pmax(0,beta_hat - qnorm(1-alpha/2)*beta_se), pmin(1,beta_hat + qnorm(1-alpha/2)*beta_se))

  return(list(beta_tilde_hat=beta_tilde_hat,
              beta_hat=beta_hat,
              beta_se=beta_se,
              ci=ci,
              input = list(X=X,Vg=Vg)))
}




