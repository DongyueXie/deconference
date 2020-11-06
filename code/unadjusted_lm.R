#'@title unadjusted method, either using ols or sandwich for variance estimation
unadjusted_lm = function(y,X,w=NULL,groups=NULL){
  G = nrow(X)
  K = ncol(X)
  y = cbind(y)
  nb = ncol(y)

  if(!is.null(groups)){
    if(is.factor(groups)){
      group_name = levels(groups)
    }else{
      group_name = levels(as.factor(groups))
    }

    group1_idx = which(groups==group_name[1])
    a = c()
    a[group1_idx] = 1/length(group1_idx)

    group2_idx = which(groups==group_name[2])
    a[group2_idx] = -1/length(group2_idx)
  }

  if(is.null(w)){
    w  = rep(1,G)
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


  if(!is.null(groups)){

    diff_group = rowMeans(beta_hat[,group1_idx,drop=F]) - rowMeans(beta_hat[,group2_idx,drop=F])

    # N_indi = length(group1_idx) + length(group2_idx)

    V_tilde = 0

    idx = c(group1_idx,group2_idx)
    for(i in idx){
      for(j in idx){
        V_tilde = V_tilde + a[i]*a[j]*asyV[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)]
      }
    }

    z_score = diff_group/sqrt(diag(V_tilde))

    p_value = (1-pnorm(abs(z_score)))*2

    ols.out = list(cov_beta_tilde_hat = covb,
                   beta_se = beta_se,
                   cov_beta_hat = asyV,

                   diff_se = sqrt(diag(V_tilde)),
                   z_score=z_score,
                   p_value=p_value)

  }else{

    diff_group = NULL
    ols.out = list(cov_beta_tilde_hat = covb,
                   beta_se = beta_se,
                   cov_beta_hat = asyV
                   )
  }



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


  if(!is.null(groups)){

    # N_indi = length(group1_idx) + length(group2_idx)

    V_tilde = 0

    idx = c(group1_idx,group2_idx)
    for(i in idx){
      for(j in idx){
        V_tilde = V_tilde + a[i]*a[j]*asyV[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)]
      }
    }

    z_score = diff_group/sqrt(diag(V_tilde))

    p_value = (1-pnorm(abs(z_score)))*2

    sand.out = list(cov_beta_tilde_hat = covb,
                    beta_se = beta_se,
                    cov_beta_hat = asyV,
                   diff_se = sqrt(diag(V_tilde)),
                   z_score=z_score,
                   p_value=p_value)

  }else{

    sand.out = list(cov_beta_tilde_hat = covb,
                    beta_se = beta_se,
                    cov_beta_hat = asyV)

  }




  # perform hc3

  Sigma = matrix(0,nrow=nb*K,ncol=nb*K)
  #Q_inv = matrix(0,nrow=nb*K,ncol=nb*K)
  Sigma_ii = matrix(0,nrow=nb*K,ncol=nb*K)

  for(i in 1:nb){
    #Q_inv[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = A_inv

    Hi = t(t(X%*%A_inv%*%t(X))*w)
    hi = diag(Hi)
    ri = y[,i] - Hi%*%y[,i]

    for(j in i:nb){

      if(j==i){

        Sigma_ij = crossprod(c(ri)/(1-pmax(pmin(hi,1-1/G),0))*w*X)

        Sigma_ii[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = Sigma_ij
      }else{
        Hj = t(t(X%*%A_inv%*%t(X))*w)
        hj = diag(Hj)
        rj = y[,j] - Hi%*%y[,j]

        Sigma_ij = crossprod(c(ri)/(1-pmax(pmin(hi,1-1/G),0))*w*X,
                             c(rj)/(1-pmax(pmin(hj,1-1/G),0))*w*X)
      }


      Sigma[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)] = Sigma_ij

    }
  }
  Sigma = Sigma+t(Sigma)-Sigma_ii
  covb = Q_inv%*%Sigma%*%Q_inv

  asyV = (J)%*%covb%*%t(J)


  beta_se = sqrt(diag(asyV))
  beta_se = matrix(beta_se,ncol=nb)

  if(!is.null(groups)){

    # N_indi = length(group1_idx) + length(group2_idx)

    V_tilde = 0

    idx = c(group1_idx,group2_idx)
    for(i in idx){
      for(j in idx){
        V_tilde = V_tilde + a[i]*a[j]*asyV[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)]
      }
    }

    z_score = diff_group/sqrt(diag(V_tilde))

    p_value = (1-pnorm(abs(z_score)))*2

    sand.out.hc3 = list(cov_beta_tilde_hat = covb,
                        beta_se = beta_se,
                        cov_beta_hat = asyV,
                    diff_se = sqrt(diag(V_tilde)),
                    z_score=z_score,
                    p_value=p_value)

  }else{

    sand.out.hc3 = list(cov_beta_tilde_hat = covb,
                        beta_se = beta_se,
                        cov_beta_hat = asyV)

  }


  rownames(beta_hat) = colnames(X)
  rownames(beta_tilde_hat) = colnames(X)





  return(list(beta_tilde_hat=beta_tilde_hat,
              beta_hat=beta_hat,
              diff_group=diff_group,
              sand.out=sand.out,
              ols.out=ols.out,
              sand.out.hc3=sand.out.hc3))
}

