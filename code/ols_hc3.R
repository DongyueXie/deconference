
ols_hc3 = function(y,X,w=NULL,groups=NULL,calc_var=TRUE){
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

  beta_hat = pmax(A_inv%*%t(Xw)%*%yw,0)
  beta_hat = cbind(beta_hat)
  p_hat = apply(beta_hat,2,function(z){z/sum(z)})
  rownames(p_hat) = colnames(X)
  rownames(beta_hat) = colnames(X)

  if(calc_var){

    res = yw - Xw%*%beta_hat
    h = rowSums((X%*%A_inv)*X)*w
    res.hc = res/(1-h)

    J = matrix(0,nrow=(nb*K),ncol=(nb*K))
    for(i in 1:nb){
      J[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = J_sum2one(beta_hat[,i],K)
    }


    score_mat = matrix(nrow=G,ncol=K*nb)
    for(i in 1:nb){
      score_i = c(res.hc[,i])*X
      score_mat[,((i-1)*K+1):(i*K)] = score_i
    }
    #
    #
    #   Sigma = matrix(0,nrow=nb*K,ncol=nb*K)
    #   Q_inv = matrix(0,nrow=nb*K,ncol=nb*K)
    #   Sigma_ii = matrix(0,nrow=nb*K,ncol=nb*K)
    #
    #   for(i in 1:nb){
    #     Q_inv[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = A_inv
    #
    #     #Hi = t(t(X%*%A_inv%*%t(X))*w)
    #     #hi = diag(Hi)
    #
    #     ri = res[,i]
    #
    #     for(j in i:nb){
    #
    #       if(j==i){
    #
    #         Sigma_ij = crossprod(c(ri)/(1-pmax(pmin(h,1-1/G),0))*w*X)
    #
    #         Sigma_ii[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = Sigma_ij
    #       }else{
    #
    #         rj = res[,j]
    #
    #         Sigma_ij = crossprod(c(ri)/(1-pmax(pmin(h,1-1/G),0))*w*X,
    #                              c(rj)/(1-pmax(pmin(h,1-1/G),0))*w*X)
    #       }
    #
    #
    #       Sigma[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)] = Sigma_ij
    #
    #     }
    #   }
    #   Sigma = Sigma+t(Sigma)-Sigma_ii
    Sigma = crossprod(score_mat)
    Sigma = (Sigma+t(Sigma))/2
    Q_inv = kronecker(diag(nb),A_inv)
    covb = Q_inv%*%Sigma%*%Q_inv

    asyV = (J)%*%covb%*%t(J)


    p_hat_se = sqrt(diag(asyV))
    p_hat_se = matrix(p_hat_se,ncol=nb)

    if(!is.null(groups)){

      V_tilde = 0

      idx = c(group1_idx,group2_idx)
      for(i in idx){
        for(j in idx){
          V_tilde = V_tilde + a[i]*a[j]*asyV[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)]
        }
      }

      diff_hat = rowMeans(p_hat[,group1_idx,drop=F]) - rowMeans(p_hat[,group2_idx,drop=F])
      z_score = diff_hat/sqrt(diag(V_tilde))

      p_value = (1-pnorm(abs(z_score)))*2

      return(list(p_hat=p_hat,
                  p_hat_se = p_hat_se,
                  diff_hat=diff_hat,
                  diff_hat_se = sqrt(diag(V_tilde)),
                  diff_z_score=z_score,
                  diff_p_value=p_value))

    }else{
      return(list(p_hat=p_hat,
                  p_hat_se = p_hat_se))
    }

  }else{
    return(list(p_hat=p_hat,
                beta_hat = beta_hat))
  }


}

