#'@param y a vector of bulk sample
#'@param X reference matrix, estimated relative expression, gene by cell
#'@param Vg variance matrix of X
#'@param X_var_pop variance of true X among individuals
#'@param marker_gene
#'@param w gene weights
#'@param hc.type hc3, hc2, or hc0
#'@param alpha significance level
#'@param a alpha in the Fuller's small sample correction
#'@param correction whether perform fuller's small sample correction.
#'@param S cell size


estimation_func = function(y,X,Vg,X_var_pop=NULL,
                           #marker_gene = NULL,
                           w=NULL,
                           hc.type='hc3',a=ncol(X)+4,
                           correction=TRUE,S=NULL,calc_cov=TRUE,
                           verbose=FALSE){

  #browser()

  #print(X)

  if(length(w)==1){
    w = rep(w,nrow(X))
  }

  # if(!is.null(marker_gene)){
  #   gene_idx = match(marker_gene,rownames(X))
  #   omit.idx = which(is.na(gene_idx))
  #   if(length(omit.idx)>0){
  #     gene_idx = gene_idx[-omit.idx]
  #   }
  #
  #   #browser()
  #
  #   X = X[gene_idx,]
  #   Vg = Vg[gene_idx,]
  #   if(!is.null(w)){
  #     w = w[gene_idx]
  #   }
  #   y = y[gene_idx,]
  #
  #   if(!is.null(X_var_pop)){
  #     X_var_pop = X_var_pop[gene_idx,]
  #   }
  #
  # }

  G = nrow(X)
  K = ncol(X)



  # nb is the number of bulk sample
  nb = ncol(y)
  if(is.null(nb)){
    nb=1
  }


  # generate weights

  if(is.null(w)){
    if(is.null(X_var_pop)){
      w  = 1/(rowSums(X)/K+rowSums(Vg)/K^2)
    }else{
      w  = 1/(rowSums(X)/K+rowSums(X_var_pop)/K^2+rowSums(Vg)/K^2)
    }
  }
  w = w/sum(w)*G



  ## take cell size into account, and scale it to O(G)
  if(is.null(S)){
    S = rep(1,K)
  }

  # X, Vg not scaled by G
  input = list(X=X,Vg=Vg,y=y,w=w,Sigma=X_var_pop,S=S)

  # scale X by G, mainly for numerical stability
  #X = X*G
  #Vg = Vg*G^2

  #input = list(X=X,Vg=Vg,y=y,w=w,Sigma=X_var_pop,S=S)

  # add library size to X
  X = X%*%diag(S)
  Xw = X*sqrt(w)
  yw = cbind(y*sqrt(w))
  A = t(Xw)%*%Xw
  d_Vg = ncol(Vg)
  Vw = Vg*w
  if(d_Vg==K^2){
    V = matrix(c(colSums(Vw*c(S%*%t(S)))),ncol = K)
  }else if(d_Vg==K){
    Vw = Vw%*%diag(S^2)
    V = diag(c(colSums(Vw)))
  }else{
    stop('check dimension of Vg')
  }



  # else if(is.array(Vg)){
  #   Vw = Vg*rep(w,each=K*K)
  #   V = rowSums(Vw,dims = 2)
  # }else{
  #   stop("Vg input type not allowed")
  # }



  #browser()
  if(verbose){
    message("estimating proportions")
  }
  beta_tilde_hat = pmax(solve(A-V)%*%t(Xw)%*%yw,0)
  ## Fuller's correction
  if(correction){

    if(verbose){
      message("performing fuller's correction")
    }

    #Z = cbind(yw,Xw)
    #m_zz = t(Z)%*%Z/G
    m_xy = t(Xw)%*%yw
    lambda = c()
    A_inv = array(dim=c(K,K,nb))

    S = matrix(0,nrow = K+1,ncol=K+1)
    S[-1,-1] = V

    #browser()
    for(ib in 1:nb){

      #print(A)

      #browser()
      #lambda[ib] = min(Re(eigen((A-m_xy[,ib]%*%t(m_xy[,ib])/sum(yw[,ib]^2))%*%solve(V))$values))

      Si = S
      Si[1,1] = 1/sum(yw[,ib])
      lambda[ib] = min(eigen(crossprod(cbind(yw[,ib],Xw))%*%solve(Si))$values)



      if(lambda[ib]>(1+1/G)){
        A_inv[,,ib] = solve(A - (1-a/G)*V)
      }else{
        A_inv[,,ib] = solve(A - (lambda[ib]-1/G-a/G)*V)
      }

      beta_tilde_hat[,ib] = pmax(A_inv[,,ib]%*%t(Xw)%*%yw[,ib],0)

    }
  }else{
    lambda = rep(1,nb)
    A_inv = array(solve(A-V),dim=c(K,K,nb))
  }

    # M = crossprod(cbind(yw,Xw))
    # B = diag(c(1/colSums(yw),1/diag(V)))
    # lambda = min(eigen(M%*%B)$values)

  #   if(lambda>(1+1/G)){
  #     A = A - (1-a/G)*V
  #   }else{
  #     A = A - (lambda-1/G-a/G)*V
  #   }
  # }else{
  #   A = A - V
  # }
  # A_inv = solve(A)
  #
  # beta_tilde_hat = pmax(A_inv%*%t(Xw)%*%yw,0)

  #browser()

  Sigma = matrix(0,nrow=nb*K,ncol=nb*K)
  Q_inv = matrix(0,nrow=nb*K,ncol=nb*K)
  Sigma_ii = matrix(0,nrow=nb*K,ncol=nb*K)

  #browser()


  if(verbose){
    message("calculating covariance matrix")
  }

  for(i in 1:nb){

    if(verbose){
      message(paste("...covaraince",i,'to',nb))
    }

    Q_inv[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = A_inv[,,i]

    if(d_Vg==K^2){
      Vbi = t(apply(Vw,1,function(z){v = matrix(z,ncol=K);v%*%beta_tilde_hat[,i]}))*pmin(lambda,1)[i]
    }
    if(d_Vg==K){
      Vbi = (Vw)%*%diag(c(beta_tilde_hat[,i]))*pmin(lambda,1)[i]
    }

    Hi = t(t(X%*%A_inv[,,i]%*%t(X))*w)
    hi = diag(Hi)
    ri = y[,i] - Hi%*%y[,i]


    for(j in i:nb){

      if(j==i){

        Vbj = Vbi

        if(hc.type == 'hc0'){
          Sigma_ij = crossprod(c(ri)*w*X+Vbi)
        }else if(hc.type == 'hc2'){
          Sigma_ij = crossprod(c(ri)/sqrt(1-pmax(pmin(hi,1-1/G),0))*w*X+Vbi)
        }else if(hc.type == 'hc3'){
          Sigma_ij = crossprod(c(ri)/(1-pmax(pmin(hi,1-1/G),0))*w*X+Vbi)
        }

        Sigma_ii[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = Sigma_ij

      }else{

        if(calc_cov){

          if(d_Vg==K^2){
            Vbj = t(apply(Vw,1,function(z){v = matrix(z,ncol=K);v%*%beta_tilde_hat[,j]}))*pmin(lambda,1)[j]
          }else if(d_Vg==K){
            Vbj = (Vw)%*%diag(c(beta_tilde_hat[,j]))*pmin(lambda,1)[j]
          }

          Hj = t(t(X%*%A_inv[,,j]%*%t(X))*w)
          hj = diag(Hj)
          rj = y[,j] - Hi%*%y[,j]
          if(hc.type == 'hc0'){
            Sigma_ij = crossprod(c(ri)*w*X+Vbi,c(rj)*w*X+Vbj)
          }else if(hc.type == 'hc2'){
            Sigma_ij = crossprod(c(ri)/sqrt(1-pmax(pmin(hi,1-1/G),0))*w*X+Vbi,c(rj)/sqrt(1-pmax(pmin(hj,1-1/G),0))*w*X+Vbj)
          }else if(hc.type == 'hc3'){
            Sigma_ij = crossprod(c(ri)/(1-pmax(pmin(hi,1-1/G),0))*w*X+Vbi,c(rj)/(1-pmax(pmin(hj,1-1/G),0))*w*X+Vbj)
          }

        }

      }

      Sigma[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)] = Sigma_ij

    }
  }
  Sigma = Sigma+t(Sigma)-Sigma_ii
  covb = Q_inv%*%Sigma%*%Q_inv

  beta_tilde_hat = cbind(beta_tilde_hat)


  # delta method

  # covb is a (nb*K)*(nb*K) cov matrix
  # formulate Jacobian matrix

  if(verbose){
    message("performing delta method")
  }

  J = matrix(0,nrow=(nb*K),ncol=(nb*K))
  for(i in 1:nb){
    J[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = J_sum2one(beta_tilde_hat[,i],K)
  }

  asyV = (J)%*%covb%*%t(J)

  beta_hat = apply(beta_tilde_hat,2,function(z){z/sum(z)})
  beta_se = sqrt(diag(asyV))
  beta_se = matrix(beta_se,ncol=nb)

  if(verbose){
    message("done")
  }

  #browser()
  # rownames(beta_hat) = colnames(X)

  return(list(beta_hat=beta_hat,
              beta_se=beta_se,
              cov_beta_hat = asyV,
              beta_tilde_hat=beta_tilde_hat,
              cov_beta_tilde_hat = covb,
              Sigma = Sigma,
              Sigma_ii = Sigma_ii,
              J=J,
              Q_inv = Q_inv,
              lambda=lambda,
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










#'@title another version, only one lambda for correction, and scale everything
#'@param y a vector of bulk sample
#'@param X reference matrix, estimated relative expression, gene by cell
#'@param Vg variance matrix of X
#'@param X_var_pop variance of true X among individuals
#'@param marker_gene
#'@param w gene weights
#'@param hc.type hc3, hc2, or hc0
#'@param alpha significance level
#'@param a alpha in the Fuller's small sample correction
#'@param correction whether perform fuller's small sample correction.
#'@param S cell size


estimation_func2 = function(y,X,Vg,X_var_pop=NULL,
                           w=NULL,
                           hc.type='hc3',a=ncol(X)+4,
                           correction=TRUE,S=NULL,calc_cov=TRUE,
                           verbose=FALSE){

  #browser()

  #print(X)

  if(length(w)==1){
    w = rep(w,nrow(X))
  }
  G = nrow(X)
  K = ncol(X)

  # nb is the number of bulk sample
  nb = ncol(y)
  if(is.null(nb)){
    nb=1
  }

  # generate weights

  if(is.null(w)){
    if(is.null(X_var_pop)){
      w  = 1/(rowSums(X)/K+rowSums(Vg)/K^2)
    }else{
      w  = 1/(rowSums(X)/K+rowSums(X_var_pop)/K^2+rowSums(Vg)/K^2)
    }
  }
  w = w/sum(w)*G

  ## take cell size into account, and scale it to O(G)
  if(is.null(S)){
    S = rep(1,K)
  }

  # X, Vg scaled by G
  input = list(X=X,Vg=Vg,y=y,w=w,Sigma=X_var_pop,S=S)

  #browser()

  #scale y
  vy_temp = apply(y,2,function(z){z*(length(z)/sum(z))^2})
  y = apply(y,2,function(z){z/sum(z)*length(z)})

  # scale X by G, mainly for numerical stability
  #X = X*G
  #Vg = Vg*G^2

  # add library size to X
  X = X%*%diag(S)
  Xw = X*sqrt(w)
  yw = cbind(y*sqrt(w))
  Vw = Vg*w
  A = t(Xw)%*%Xw
  d_Vg = ncol(Vg)
  if(d_Vg==K^2){
    V = matrix(c(colSums(Vw)),ncol = K)
  }else if(d_Vg==K){
    #Vw = Vw%*%diag(S^2)
    #Vw = Vw*(S^2)
    V = diag(c(colSums(Vw)))
  }else{
    stop('check dimension of Vg')
  }
  V = diag(S)%*%V%*%diag(S)

  #browser()
  if(verbose){
    message("estimating proportions")
  }
  Q = (A-V)/G
  Qinv = solve(Q)
  beta_tilde_hat = pmax(Qinv%*%t(Xw)%*%yw,0)/G
  ## Fuller's correction
  if(correction){

    if(verbose){
      message("performing fuller's correction")
    }

    y_temp = rowSums(yw)

    Z = cbind(y_temp,Xw)

    S = matrix(0,nrow = K+1,ncol=K+1)
    S[1,1] = sum(rowSums(vy_temp)*w)
    S[-1,-1] = V

    #browser()

    lambda0 = min(eigen(solve(S)%*%crossprod(Z))$values)

    if(lambda0>(1+1/G)){
      lambda = (1-a/G)
    }else{
      lambda = (lambda0-1/G-a/G)
    }
    Q = (A - lambda*V)/G
    Qinv = solve(Q)
    beta_tilde_hat = pmax(Qinv%*%t(Xw)%*%yw,0)/G

  }else{
    lambda = 1
    lambda0 = 1
  }

  #browser()

  Q_inv = kronecker(diag(nb),Qinv)
  H = t(t(X%*%Qinv%*%t(X))*w)/G
  # Sigma = matrix(0,nrow=nb*K,ncol=nb*K)
  # #Q_inv = matrix(0,nrow=nb*K,ncol=nb*K)
  # Sigma_ii = matrix(0,nrow=nb*K,ncol=nb*K)
  #
  # #browser()
  #
  #
  # if(verbose){
  #   message("calculating covariance matrix")
  # }
  #
  #
  # #H = t(t(X%*%Qinv%*%t(X))*w)/G
  # h = diag(H)
  # h = pmax(pmin(h,1-1/G),0)
  # for(i in 1:nb){
  #
  #   if(verbose){
  #     message(paste("...covaraince",i,'to',nb))
  #   }
  #
  #   #Q_inv[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = Qinv
  #
  #   if(d_Vg==K^2){
  #     Vbi = t(apply(Vw,1,function(z){v = matrix(z,ncol=K);v%*%beta_tilde_hat[,i]}))*lambda
  #   }
  #   if(d_Vg==K){
  #     Vbi = (Vw)%*%diag(c(beta_tilde_hat[,i]))*lambda
  #   }
  #   ri = y[,i] - H%*%y[,i]
  #
  #
  #   for(j in i:nb){
  #
  #     if(j==i){
  #
  #       if(hc.type == 'hc0'){
  #         Sigma_ij = crossprod(c(ri)*w*X+Vbi)
  #       }else if(hc.type == 'hc2'){
  #         Sigma_ij = crossprod(c(ri)/sqrt(1-h)*w*X+Vbi)
  #       }else if(hc.type == 'hc3'){
  #         Sigma_ij = crossprod(c(ri)/(1-h)*w*X+Vbi)
  #       }
  #
  #       Sigma_ii[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = Sigma_ij
  #
  #     }else{
  #
  #       if(calc_cov){
  #
  #         if(d_Vg==K^2){
  #           Vbj = t(apply(Vw,1,function(z){v = matrix(z,ncol=K);v%*%beta_tilde_hat[,j]}))*lambda
  #         }
  #         if(d_Vg==K){
  #           Vbj = (Vw)%*%diag(c(beta_tilde_hat[,j]))*lambda
  #         }
  #
  #         rj = y[,j] - H%*%y[,j]
  #         if(hc.type == 'hc0'){
  #           Sigma_ij = crossprod(c(ri)*w*X+Vbi,c(rj)*w*X+Vbj)
  #         }else if(hc.type == 'hc2'){
  #           Sigma_ij = crossprod(c(ri)/sqrt(1-h)*w*X+Vbi,c(rj)/sqrt(1-h)*w*X+Vbj)
  #         }else if(hc.type == 'hc3'){
  #           Sigma_ij = crossprod(c(ri)/(1-h)*w*X+Vbi,c(rj)/(1-h)*w*X+Vbj)
  #         }
  #
  #       }
  #
  #     }
  #
  #     Sigma[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)] = Sigma_ij
  #
  #   }
  # }
  # Sigma = (Sigma+t(Sigma)-Sigma_ii)/G/G
  Sigma = get_SIGMA(y=y,X=X,w=w,beta=beta_tilde_hat,V=Vw,H=H,nb=nb,G=G,K=K,lambda=lambda,verbose=verbose,calc_cov=calc_cov,hc.type=hc.type)
  covb = Q_inv%*%Sigma%*%Q_inv

  beta_tilde_hat = cbind(beta_tilde_hat)


  # delta method

  # covb is a (nb*K)*(nb*K) cov matrix
  # formulate Jacobian matrix

  if(verbose){
    message("performing delta method")
  }

  J = matrix(0,nrow=(nb*K),ncol=(nb*K))
  for(i in 1:nb){
    J[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = J_sum2one(beta_tilde_hat[,i],K)
  }

  asyV = (J)%*%covb%*%t(J)

  beta_hat = apply(beta_tilde_hat,2,function(z){z/sum(z)})
  beta_se = sqrt(diag(asyV))
  beta_se = matrix(beta_se,ncol=nb)

  if(verbose){
    message("done")
  }

  #browser()
  # rownames(beta_hat) = colnames(X)

  return(list(beta_hat=beta_hat,
              beta_se=beta_se,
              cov_beta_hat = asyV,
              beta_tilde_hat=beta_tilde_hat,
              cov_beta_tilde_hat = covb,
              Sigma = Sigma,
              #Sigma_ii = Sigma_ii/G^2,
              J=J,
              Q_inv = Qinv,
              Q=Q,
              V=V,
              A=A,
              #h=h,
              H=H,
              lambda=lambda0,
              input = input))

}




get_SIGMA = function(y,X,w,beta,V,H,nb,G,K,lambda,verbose,calc_cov,hc.type){



  Sigma = matrix(0,nrow=nb*K,ncol=nb*K)
  #Q_inv = matrix(0,nrow=nb*K,ncol=nb*K)
  Sigma_ii = matrix(0,nrow=nb*K,ncol=nb*K)

  #browser()


  if(verbose){
    message("calculating covariance matrix")
  }



  h = diag(H)
  h = pmax(pmin(h,1-1/G),0)
  for(i in 1:nb){

    if(verbose){
      message(paste("...covaraince",i,'to',nb))
    }

    if(ncol(V)==K^2){
      Vbi = t(apply(V,1,function(z){v = matrix(z,ncol=K);v%*%beta[,i]}))*lambda
    }
    if(ncol(V)==K){
      Vbi = (V)%*%diag(c(beta[,i]))*lambda
    }
    ri = y[,i] - H%*%y[,i]


    for(j in i:nb){

      if(j==i){

        if(hc.type == 'hc0'){
          Sigma_ij = crossprod(c(ri)*w*X+Vbi)
        }else if(hc.type == 'hc2'){
          Sigma_ij = crossprod(c(ri)/sqrt(1-h)*w*X+Vbi)
        }else if(hc.type == 'hc3'){
          Sigma_ij = crossprod(c(ri)/(1-h)*w*X+Vbi)
        }

        Sigma_ii[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = Sigma_ij

      }else{

        if(calc_cov){

          if(ncol(V)==K^2){
            Vbj = t(apply(V,1,function(z){v = matrix(z,ncol=K);v%*%beta[,j]}))*lambda
          }
          if(ncol(V)==K){
            Vbj = (V)%*%diag(c(beta[,j]))*lambda
          }

          rj = y[,j] - H%*%y[,j]
          if(hc.type == 'hc0'){
            Sigma_ij = crossprod(c(ri)*w*X+Vbi,c(rj)*w*X+Vbj)
          }else if(hc.type == 'hc2'){
            Sigma_ij = crossprod(c(ri)/sqrt(1-h)*w*X+Vbi,c(rj)/sqrt(1-h)*w*X+Vbj)
          }else if(hc.type == 'hc3'){
            Sigma_ij = crossprod(c(ri)/(1-h)*w*X+Vbi,c(rj)/(1-h)*w*X+Vbj)
          }

        }

      }

      Sigma[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)] = Sigma_ij

    }
  }
  Sigma = (Sigma+t(Sigma)-Sigma_ii)/G/G

  Sigma

}





