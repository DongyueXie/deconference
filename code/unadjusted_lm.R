#'@title unadjusted method, either using ols or sandwich for variance estimation
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

