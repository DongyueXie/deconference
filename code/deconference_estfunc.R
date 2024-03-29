


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
#'@param X_array array of reference matrix from each ref individual
#'@param use.weight.ref whether use reference matrices to calculate weights
#'@param p.for.weight what p to use for weight calculation
#'@param p.true the true cell type proportion
#'@param R correlation matrix of genes
#'@param true.beta true beta before normalizing to 1
#'@param asyV.pos check if asymptotic Variance is positive definite
#'@param Q.pos check if Q is positive definite
#'@param nfold if using jacknife, the number of fold to be used.
#'@param folds if using jackknife, the folds to be used. folds looks like 111222333.
#'@param R01 the 0-1 correlation matrix
#'@param calc_var whether calculate the covariance matrices.
#'@param use.QP whether use quadratic programming to impose nonnegative constraint

estimation_func2 = function(y,X,Vg,X_var_pop=NULL,
                           w=NULL,
                           hc.type='hc3',
                           a=ncol(X)+4,
                           correction=FALSE,
                           S=NULL,
                           #calc_cov=TRUE,
                           verbose=FALSE,
                           #X_array = NULL,
                           #use.weight.ref = FALSE,
                           #p.for.weight = "equal",
                           #p.true = NULL,
                           #cor.idx=NULL,
                           true.beta = NULL,
                           centeringXY = FALSE,
                           Sigma.pos = FALSE,
                           Q.pos = FALSE,
                           V_tilde.pos = FALSE,
                           only.scale.pos.res=FALSE,
                           only.add.pos.res = FALSE,
                           nfold=10,
                           folds = NULL,
                           #use_all_pair_for_cov = FALSE,
                           #b=NULL,
                           #b_sd=NULL,
                           R01=NULL,
                           groups=NULL,
                           two_group_method = 'asymptotic',
                           calc_var=TRUE,
                           use.QP = FALSE
                           ){

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

    if(use.weight.ref){

      if(p.for.weight == 'true'){

      }else if(p.for.weight == 'equal'){


        Sigma_sample = c()
        for(i in 1:dim(X_array)[3]){
          Sigma_sample = cbind(Sigma_sample,(X_array[,,i] - X)%*%rep(1/K,K))
        }
        cov.Xb = Rfast::cova(t(Sigma_sample))
        w = 1/diag(cov.Xb)

      }else if(p.for.weight == 'estimated'){

      }

    }else{

      if(is.null(X_var_pop)){
        w  = 1/(rowSums(X)/K+rowSums(Vg)/K^2)
      }else{
        w  = 1/(rowSums(X)/K+rowSums(X_var_pop)/K^2+rowSums(Vg)/K^2)
      }

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


  # if(!is.null(true.beta)){
  #   true.beta = t(t(true.beta)*G/c(colSums(y)))
  # }
  # y = apply(y,2,function(z){z/sum(z)*length(z)})

  #browser()

  # scale X by G, mainly for numerical stability
  #X = X*G
  #Vg = Vg*G^2

  # add library size to X
  X = X%*%diag(S)
  Xw = X*sqrt(w)
  yw = cbind(y*sqrt(w))
  if(centeringXY){
    yw = apply(yw,2,scale,scale=FALSE)
    Xw = apply(Xw,2,scale,scale=FALSE)
  }
  Vw = Vg*w
  A = crossprod(Xw)
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
  Q = (A-V)
  if(Q.pos){
    #print(isSymmetric(Q))
    Q = (Q+t(Q))/2
    Q = make.pos.def(Q)
  }
  Qinv = solve(Q)


  if(use.QP){
    beta_tilde_hat = matrix(nrow=K,ncol=nb)
    for(b in 1:nb){
      QP.out = quadprog::solve.QP(Q,t(Xw)%*%yw[,b],diag(K),rep(0,K))
      beta_tilde_hat[,b] = QP.out$solution
    }

  }else{
    beta_tilde_hat = pmax(Qinv%*%t(Xw)%*%yw,0)
  }


  #browser()
  ## Fuller's correction
  if(correction){

    #scale y
    vy_temp = apply(y,2,function(z){z*(length(z)/sum(z))^2})

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

  if(calc_var){

    Q_inv = kronecker(diag(nb),Qinv)
    #H = t(t(X%*%Qinv%*%t(X))*w)/G
    #h = rowSums((X%*%Qinv)*X)*w
    h = rowSums((Xw%*%Qinv)*Xw)



    # Sigma = get_SIGMA2(y=yw,X=Xw,beta=if(is.null(true.beta)){beta_tilde_hat}else{true.beta},
    #                   V=Vw,h=h,nb=nb,
    #                   G=G,K=K,lambda=lambda,verbose=verbose,
    #                   calc_cov=calc_cov,hc.type=hc.type,
    #                   cor.idx=cor.idx,
    #                   only.scale.pos.res=only.scale.pos.res,
    #                   only.add.pos.res=only.add.pos.res,
    #                   nfold=nfold,
    #                   folds=folds,
    #                   use_all_pair_for_cov=use_all_pair_for_cov,
    #                   b=b,
    #                   b_sd=b_sd,
    #                   R01=R01)


    Sigma = get_SIGMA3(y=yw,X=Xw,beta=if(is.null(true.beta)){beta_tilde_hat}else{true.beta},
                       V=Vw,h=h,nb=nb,
                       G=G,K=K,verbose=verbose,
                       hc.type=hc.type,
                       only.add.pos.res=only.add.pos.res,
                       nfold=nfold,
                       folds=folds,
                       R01=R01)
    score_mat = Sigma$score_mat
    Sigma = Sigma$Sigma
    if(Sigma.pos){
      #print(isSymmetric(Sigma))
      Sigma = make.pos.def(Sigma)
    }
    covb = Q_inv%*%Sigma%*%Q_inv

    #browser()

    # if(is.null(R)){
    #   browser()
    # }

    #browser()

    # if(!is.null(R)){
    #   print(paste('adjusted for corr:',sum(diag(Sigma))))
    # }else{
    #   print(paste('not adjusted for corr:',sum(diag(Sigma))))
    # }

    #print("Looking at sum of Q_invSignaQ_inv")

    # if(!is.null(R)){
    #   print(paste('adjusted for corr:',sum(covb)*1000))
    # }else{
    #   print(paste('not adjusted for corr:',sum(covb)*1000))
    # }

    #print("Looking at diagonal of Q_invSignaQ_inv")

    # if(!is.null(R)){
    #   print(paste('adjusted for corr:',sum(diag(covb))))
    # }else{
    #   print(paste('not adjusted for corr:',sum(diag(covb))))
    # }

    #browser()




    # delta method

    # covb is a (nb*K)*(nb*K) cov matrix
    # formulate Jacobian matrix

    if(verbose){
      message("performing delta method")
    }

    beta_tilde_hat = cbind(beta_tilde_hat)
    rownames(beta_tilde_hat) = colnames(X)
    J = matrix(0,nrow=(nb*K),ncol=(nb*K))

    if(is.null(true.beta)){

      for(i in 1:nb){
        J[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = J_sum2one(beta_tilde_hat[,i],K)
      }

    }else{

      for(i in 1:nb){
        J[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = J_sum2one(true.beta[,i],K)
      }

    }



    asyV = (J)%*%covb%*%t(J)
    # if(asyV.pos){
    #   diag(asyV) = pmax(diag(asyV),0)
    # }

    # if(!is.null(R)){
    #   print(paste('adjusted for corr:',sum(asyV)))
    # }else{
    #   print(paste('not adjusted for corr:',sum(asyV)))
    # }

    # if(!is.null(R)){
    #   print(paste('adjusted for corr:',sum(diag(asyV))))
    # }else{
    #   print(paste('not adjusted for corr:',sum(diag(asyV))))
    # }

    #browser()

    p_hat = apply(beta_tilde_hat,2,function(z){z/sum(z)})
    p_hat_se = sqrt(diag(asyV))
    p_hat_se = matrix(p_hat_se,ncol=nb)
    rownames(p_hat) = colnames(X)
    rownames(p_hat_se) = colnames(X)

    if(!is.null(groups)){
      two_group_res = get_two_group_res(score_mat,J,Q_inv,groups,K,p_hat,asyV,two_group_method,R01,V_tilde.pos)
    }else{
      two_group_res = NULL
    }

    if(verbose){
      message("done")
    }

    #browser()




    return(list(p_hat=p_hat,
                p_hat_se=p_hat_se,
                p_hat_cov = asyV,
                beta_hat=beta_tilde_hat,
                beta_hat_se = matrix(sqrt(diag(covb)),ncol=nb),
                beta_hat_cov = covb,
                Sigma = Sigma,
                two_group_res=two_group_res
                #Sigma_ii = Sigma_ii/G^2,
                #J=J,
                #Q_inv = Q_inv,
                #Q=Q,
                #V=V,
                #A=A,
                #h=h,
                #H=H,
                #lambda=lambda0,
                #input = input
    ))

  }else{

    p_hat = apply(beta_tilde_hat,2,function(z){z/sum(z)})
    return(list(p_hat=p_hat,
                beta_hat=beta_tilde_hat
                #Sigma_ii = Sigma_ii/G^2,
                #J=J,
                #Q_inv = Q_inv,
                #Q=Q,
                #V=V,
                #A=A,
                #h=h,
                #H=H,
                #lambda=lambda0,
                #input = input
    ))

  }



}

#'@title two group test function
#'@param socre_mat score matri, N by K
#'@param J Jacobian matrix
#'@param Q_inv inverse of Q
#'@param groups group labels
#'@param K number of cell types
#'@param p_hat estimated proportions
#'@param p_hat_cov covaraince matrix of p_hat
#'@param R01 0-1 correlation matrix
get_two_group_res = function(score_mat,J,Q_inv,groups,K,p_hat,p_hat_cov,method,R01,V_tilde.pos){
  nb = length(groups)
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

  diff_group = rowMeans(p_hat[,group1_idx,drop=F]) - rowMeans(p_hat[,group2_idx,drop=F])
  if(method=='asymptotic'){
    temp = t(J%*%Q_inv%*%t(score_mat))
    temp2 = 0
    for(i in 1:nb){
      temp2 = temp2 + a[i]*temp[,((i-1)*K+1):(i*K)]
    }
    if(!is.null(R01)){
      V_tilde = crossprod(temp2,R01)%*%temp2
    }else{
      V_tilde = crossprod(temp2)
    }

  }
  if(method=='var(delta_hat)'){
    V_tilde = 0
    idx = c(group1_idx,group2_idx)
    for(i in idx){
      for(j in idx){
        V_tilde = V_tilde + a[i]*a[j]*p_hat_cov[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)]
      }
    }
  }
  if(V_tilde.pos){
    #print(isSymmetric(V_tilde))
    V_tilde = (V_tilde+t(V_tilde))/2
    V_tilde = make.pos.def(V_tilde)
  }
  z_score = diff_group/sqrt(diag(V_tilde))
  p_value = (1-pnorm(abs(z_score)))*2
  return(list(p_value = p_value,z_score=z_score,V_tilde=V_tilde,
              diff_hat=diff_group,diff_hat_se = sqrt(diag(V_tilde)),a=a))
}


# new function, faster by exploring matrix operation
get_SIGMA3 = function(y,X,beta,V,h,nb,G,K,verbose,hc.type,
                      only.add.pos.res,nfold,folds,
                      R01){

  res = y - X%*%beta
  h = pmax(pmin(abs(h),1-1/G),0)
  if(hc.type == 'hc0'){
    res.hc = res
  }else if(hc.type == 'hc2'){
    res.hc = res/sqrt(1-h)
  }else if(hc.type == 'hc3'){
    res.hc = res/(1-h)
  }else if(hc.type=='jackknife'){
    res.hc = get_jack_res(y,X,V,nfold=nfold,folds=folds)
  }else if(hc.type=='jackknife_indep'){
    res.hc = get_jack_res_indep(y,X,V,nfold=nfold,folds=folds,R01=R01)
  }


  d_V = ncol(V)
  # formulate score matrix, G by nk
  score_mat = matrix(nrow=G,ncol=K*nb)
  for(i in 1:nb){

    if(d_V==K^2){
      #Vbi = t(apply(V,1,function(z){v = matrix(z,ncol=K);v%*%beta[,i]}))
      Vbi = V%*%matrix(c(rep(c(beta[,i],rep(0,K^2)),K-1),beta[,i]),ncol=K)
    }
    if(d_V==K){
      Vbi = (V)%*%diag(c(beta[,i]))
    }
    score_i = (c(res.hc[,i])*X+Vbi)
    score_mat[,((i-1)*K+1):(i*K)] = score_i

  }

  if(verbose){
    message("calculating covariance matrix")
  }

  if(!is.null(R01)){
    if(!only.add.pos.res){
      # R01 01 cor mat, G by G
      # score_mat G by nK
      Sigma = crossprod(score_mat,R01)%*%score_mat

    }else{
      #browser()
      R01_idx = summary(R01)
      # obtain average matrix A
      #i_idx = R01_idx$i
      #dp = diff(R01@p)
      #j_idx = R01_idx$j
      if(verbose){
        message("averaging A over individuals")
      }
      res_prod = 0
      for(i in 1:nb){
        res_prod = res_prod + (res.hc[R01_idx$i,i]*res.hc[R01_idx$j,i])
      }
      res_prod = res_prod>0
      print(sum(res_prod))
      # res_prod = 1
      # for(i in 1:nb){
      #   res_prod = res_prod * (res.hc[R01_idx$i,i]*res.hc[R01_idx$j,i])>0
      # }

      # res.hc G by n residual matrix
      # length(R01_idx$i)
      #idx = rowSums((res.hc[R01_idx$i,]*res.hc[R01_idx$j,]))>0
      # mean.res = rowMeans(res.hc)
      # res_prod = (mean.res[R01_idx$i]*mean.res[R01_idx$j])>0
      A = sparseMatrix(i=R01_idx$i[res_prod],j = R01_idx$j[res_prod],dims = c(G,G))
      Sigma = crossprod(score_mat,A)%*%score_mat
    }
  }else{
    Sigma = crossprod(score_mat)
  }


  Sigma = (Sigma+t(Sigma))/2
  return(list(Sigma=Sigma,score_mat=score_mat))


}




#'@description allow correlations among samples, and use yw and Xw, not X and y
#'#'@param only.scale.pos.res only apply hc adjustment to positive empirical covariance, when accounting for correlation?
get_SIGMA2 = function(y,X,beta,V,h,nb,G,K,lambda,verbose,calc_cov,hc.type,cor.idx,
                      only.scale.pos.res,only.add.pos.res,nfold,folds,use_all_pair_for_cov,
                      b,b_sd,R01){



  Sigma = matrix(0,nrow=nb*K,ncol=nb*K)
  Sigma_ii = matrix(0,nrow=nb*K,ncol=nb*K)

  #browser()

  res = y - X%*%beta


  h = pmax(pmin(abs(h),1-1/G),0)
  if(hc.type == 'hc0'){
    res.hc = res
  }else if(hc.type == 'hc2'){
    res.hc = res/sqrt(1-h)
  }else if(hc.type == 'hc3'){
    res.hc = res/(1-h)
  }else if(hc.type=='jackknife'){
    res.hc = get_jack_res(y,X,V,nfold=nfold,folds=folds)
  }else if(hc.type=='jackknife_indep'){
    res.hc = get_jack_res_indep(y,X,V,nfold=nfold,folds=folds,R01=R01)
  }





  if(verbose){
    message("calculating covariance matrix")
  }



  for(i in 1:nb){

    if(verbose){
      message(paste("...covariance",i,'to',nb))
    }

    if(ncol(V)==K^2){
      Vbi = t(apply(V,1,function(z){v = matrix(z,ncol=K);v%*%beta[,i]}))*lambda
    }
    if(ncol(V)==K){
      Vbi = (V)%*%diag(c(beta[,i]))*lambda
    }
    ri = res.hc[,i]
    score.temp = (c(ri)*X+Vbi)

    for(j in i:nb){
      Sigma_ij = matrix(0,nrow=K,ncol=K)
      if(j==i){
        if(use_all_pair_for_cov){
          Sigma_ij = tcrossprod(colSums(score.temp))
        }else{
          if(!is.null(R01)){
            #browser()
            if(only.add.pos.res){


              ## Check if the result is correct

              ri.sign = sign(ri)
              # tcrossprod(ri.sign) * R01
              #idx = R@x>0
              #tic()
              i_idx = R01@i+1
              dp = diff(R01@p)
              j_idx = rep(seq_along(dp),dp)
              idx = (ri.sign[i_idx]*ri.sign[j_idx])>0
              # to calculate over all individuals:
              # Res: residual matrix, G by n
              # idx = rowSums((Res[i_idx,]*Res[j_idx,]))>0
              A = sparseMatrix(i=i_idx[idx],j = j_idx[idx],dims = c(G,G))
              Sigma_ij = Sigma_ij = t(score.temp)%*%A%*%score.temp
              #toc()




              #########################
              # temp = t(R01*ri.sign)*ri.sign
              # idx = temp@x>0
              # i_idx = temp@i[idx]+1
              # dp = diff(temp@p)
              # j_idx = rep(seq_along(dp),dp)
              # j_idx = j_idx[idx]
              # A = sparseMatrix(i=i_idx,j = j_idx,dims = c(G,G))
              # Sigma_ij = Sigma_ij = t(score.temp)%*%A%*%score.temp
              ###########################
              #browser()
              #Sigma_ij = t(score.temp)%*%((ri>0)*t(R01*(ri>0)))%*%score.temp + t(score.temp)%*%((ri<0)*t(R01*(ri<0)))%*%score.temp
            }else{
              Sigma_ij = t(score.temp)%*%R01%*%score.temp
            }
          }else{
            Sigma_ij = crossprod(score.temp)
          }
        }
        #browser()
        Sigma_ii[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = as.matrix(Sigma_ij)



          # if(!is.null(cor.idx)){
          #
          #   # l1.temp = (c(ri)*X+Vbi)[cor.idx[,1],,drop=FALSE]
          #   # l2.temp = (c(ri)*X+Vbi)[cor.idx[,2],,drop=FALSE]
          #   # Sigma_ij = Sigma_ij + crossprod(l1.temp,l2.temp) + crossprod(l2.temp,l1.temp)
          #
          #   if(only.scale.pos.res){
          #     # find pos res prod
          #
          #     score.temp.no.hc = (c(res[,i])*X+Vbi)
          #     # ecov = 0
          #     # for(cc in 1:dim(cor.idx)[1]){
          #     #   res.prod = ri[cor.idx[cc,][1]] * ri[cor.idx[cc,][2]]
          #     #   if(res.prod>0){
          #     #     ecov = ecov + tcrossprod(score.temp[cor.idx[cc,][1],],score.temp[cor.idx[cc,][2],])
          #     #   }else{
          #     #     ecov = ecov + tcrossprod(score.temp.no.hc[cor.idx[cc,][1],],
          #     #                              score.temp.no.hc[cor.idx[cc,][2],])
          #     #   }
          #     # }
          #
          #     res.prod = ri[cor.idx[,1]] * ri[cor.idx[,2]]
          #     s.idx = which(res.prod>0)
          #
          #     Sigma_ij = Sigma_ij+
          #       crossprod((score.temp)[cor.idx[s.idx,1],,drop=FALSE],(score.temp)[cor.idx[s.idx,2],,drop=FALSE])+
          #       crossprod((score.temp.no.hc)[cor.idx[-s.idx,1],,drop=FALSE],(score.temp.no.hc)[cor.idx[-s.idx,2],,drop=FALSE])
          #
          #
          #
          #
          #   }else if(only.add.pos.res){
          #
          #     #browser()
          #
          #     s.idx = which((ri[cor.idx[,1]] * ri[cor.idx[,2]])>0)
          #     #
          #     # idx.temp = split(s.idx,1:100)
          #     # cc.temp = lapply(idx.temp,function(idx){crossprod((score.temp)[cor.idx[idx,1],,drop=FALSE],
          #     #                                                   (score.temp)[cor.idx[idx,2],,drop=FALSE])})
          #     # cc.temp = Reduce('+',cc.temp)
          #     # Sigma_ij = Sigma_ij + cc.temp
          #
          #     Sigma_ij = Sigma_ij+
          #       crossprod((score.temp)[cor.idx[s.idx,1],,drop=FALSE],(score.temp)[cor.idx[s.idx,2],,drop=FALSE])
          #   }else{
          #     # idx.temp = split(1:dim(cor.idx)[1],1:100)
          #     # cc.temp = lapply(idx.temp,function(idx){crossprod((score.temp)[cor.idx[idx,1],,drop=FALSE],
          #     #                                                   (score.temp)[cor.idx[idx,2],,drop=FALSE])})
          #     # cc.temp = Reduce('+',cc.temp)
          #     # Sigma_ij = Sigma_ij + cc.temp
          #     # browser()
          #     Sigma_ij = Sigma_ij + crossprod((score.temp)[cor.idx[,1],,drop=FALSE],
          #                                     (score.temp)[cor.idx[,2],,drop=FALSE])
          #   }

        # }
        #
        #
        #
        # }

        ###########
        # if(!is.null(R)){
        #   print(paste('adjusted for corr:',sum(Sigma_ij)))
        # }else{
        #   print(paste('not adjusted for corr:',sum(Sigma_ij)))
        # }

        ############



      }else{

        # need to revise and add different options to calc variance
        if(calc_cov){

          if(ncol(V)==K^2){
            Vbj = t(apply(V,1,function(z){v = matrix(z,ncol=K);v%*%beta[,j]}))*lambda
          }
          if(ncol(V)==K){
            Vbj = (V)%*%diag(c(beta[,j]))*lambda
          }

          rj = res.hc[,j]
          score.temp.j = (c(rj)*X+Vbj)

          if(!is.null(R01)){
            if(only.add.pos.res){
              Sigma_ij = t(score.temp)%*%(t(R01*(ri>0))*(rj>0))%*%(score.temp.j)
              + t(score.temp)%*%(t(R01*(ri<0))*(rj<0))%*%(score.temp.j)
            }else{
              Sigma_ij = t(score.temp)%*%R01%*%(score.temp.j)
            }

          }else{
            Sigma_ij = crossprod(score.temp,score.temp.j)
          }

          # if(!is.null(cor.idx)){
          #   # Sigma_ij = Sigma_ij
          #   # + crossprod((c(ri)*X+Vbi)[cor.idx[,1],,drop=FALSE],(c(rj)*X+Vbj)[cor.idx[,2],,drop=FALSE])
          #   # + crossprod((c(ri)*X+Vbi)[cor.idx[,2],,drop=FALSE],(c(rj)*X+Vbj)[cor.idx[,1],,drop=FALSE])
          #
          #   Sigma_ij = Sigma_ij
          #   + crossprod((c(ri)*X+Vbi)[cor.idx[,1],,drop=FALSE],(c(rj)*X+Vbj)[cor.idx[,2],,drop=FALSE])
          # }
          Sigma[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)] = as.matrix(Sigma_ij)
        }

      # }
    }
    }
  }

  Sigma = (Sigma+t(Sigma)+Sigma_ii)

  Sigma

}


get_jack_res = function(y,X,V,nfold = 10,folds=NULL){
  n = nrow(y)
  n_bulk = ncol(y)
  K = ncol(X)

  d_Vg = ncol(V)

  #browser()


  res = matrix(nrow=n,ncol=n_bulk)
  if(is.null(folds)){
    folds = cut(1:n,breaks = nfold,labels = F)
  }else{
    nfold = length(table(folds))
  }


  # for the genes in each fold, estimate the beta hat use the genes in the rest folds
  for(f in 1:nfold){
    idx = which(folds==f)
    X.temp = X[-idx,]
    y.temp = y[-idx,]
    if(d_Vg==K^2){
      V.temp = matrix(c(colSums(V[-idx,])),ncol = K)
    }else if(d_Vg==K){
      V.temp = diag(c(colSums(V[-idx,])))
    }else{
      stop('check dimension of V')
    }

    bhat = pmax(solve(crossprod(X.temp)-V.temp)%*%t(X.temp)%*%y.temp,0)
    res[idx,] = y[idx,] - X[idx,]%*%bhat
  }

  res

}

#
# get_jack_res2 = function(y,X,V,nfold = 10,folds=NULL,R01){
#   n = nrow(y)
#   n_bulk = ncol(y)
#   K = ncol(X)
#
#   d_Vg = ncol(V)
#
#   #browser()
#
#
#   res = matrix(nrow=n,ncol=n_bulk)
#   if(is.null(folds)){
#     folds = cut(1:n,breaks = nfold,labels = F)
#   }else{
#     nfold = length(table(folds))
#   }
#
#
#   for(f in 1:nfold){
#     idx = which(folds==f)
#
#     temp = rowSums(R01[idx, ])
#     index = idx[which(temp == 0)]
#     if(length(index)<(length(idx)/2)){
#       index = idx[(order(temp,decreasing = F))[1:(length(idx)/2)]]
#     }
#
#     X.temp = X[index,]
#     y.temp = y[index,]
#     if(d_Vg==K^2){
#       V.temp = matrix(c(colSums(V[index,])),ncol = K)
#     }else if(d_Vg==K){
#       V.temp = diag(c(colSums(V[index,])))
#     }else{
#       stop('check dimension of V')
#     }
#
#     bhat = pmax(solve(crossprod(X.temp)-V.temp)%*%t(X.temp)%*%y.temp,0)
#     res[idx,] = y[idx,] - X[idx,]%*%bhat
#   }
#
#   res
#
# }
#
#
# get_jack_res3 = function(y,X,V,nfold = 10,folds=NULL,R01){
#   n = nrow(y)
#   n_bulk = ncol(y)
#   K = ncol(X)
#
#   d_Vg = ncol(V)
#
#   #browser()
#
#
#   res = matrix(nrow=n,ncol=n_bulk)
#   if(is.null(folds)){
#     folds = cut(1:n,breaks = nfold,labels = F)
#   }else{
#     nfold = length(table(folds))
#   }
#
#
#   for(f in 1:nfold){
#     idx = which(folds==f)
#
#     temp = rowSums(R01[-idx, ])
#     index = ((1:n)[-idx])[which(temp == 0)]
#     if(length(index)<(length(temp)/2)){
#       index = ((1:n)[-idx])[(order(temp,decreasing = F))[1:(length(temp)/2)]]
#     }
#
#     X.temp = X[index,]
#     y.temp = y[index,]
#     if(d_Vg==K^2){
#       V.temp = matrix(c(colSums(V[index,])),ncol = K)
#     }else if(d_Vg==K){
#       V.temp = diag(c(colSums(V[index,])))
#     }else{
#       stop('check dimension of V')
#     }
#
#     bhat = pmax(solve(crossprod(X.temp)-V.temp)%*%t(X.temp)%*%y.temp,0)
#     res[idx,] = y[idx,] - X[idx,]%*%bhat
#   }
#
#   res
#
# }


get_jack_res_indep = function(y,X,V,nfold = 10,folds=NULL,R01){
  n = nrow(y)
  n_bulk = ncol(y)
  K = ncol(X)

  d_Vg = ncol(V)

  #browser()


  res = matrix(nrow=n,ncol=n_bulk)
  if(is.null(folds)){
    folds = cut(1:n,breaks = nfold,labels = F)
  }else{
    nfold = length(table(folds))
  }


  for(f in 1:nfold){
    idx = which(folds==f)

    temp = colSums(R01[idx, ])
    #index = ((1:n)[-idx])[which(temp == 0)]
    index = which(temp==0)
    #print(length(index))
    if(length(index)<(n/2)){
      index = ((1:n)[-idx])[(order(temp[-idx],decreasing = F))[1:((n-length(idx))/2)]]
    }

    X.temp = X[index,]
    y.temp = y[index,]
    if(d_Vg==K^2){
      V.temp = matrix(c(colSums(V[index,])),ncol = K)
    }else if(d_Vg==K){
      V.temp = diag(c(colSums(V[index,])))
    }else{
      stop('check dimension of V')
    }

    bhat = pmax(solve(crossprod(X.temp)-V.temp)%*%t(X.temp)%*%y.temp,0)
    res[idx,] = y[idx,] - X[idx,]%*%bhat
  }

  res

}

#'@title make a symmetric matrix pos.def
#'@description by setting the negative evs to be a small number
#'@param
make.pos.def = function(X,s=100){
  attempt = try(chol(X),silent = T)
  if(class(attempt)[1]=="try-error"){
    X.evd = eigen(X)
    D = X.evd$values
    if(D[1]<0){
      stop('matrix is negative definite')
    }
    D[D<=0] = D[1]/s
    X = X.evd$vectors%*%diag(D)%*%t(X.evd$vectors)
  }
  X
}


# get_jack_res_b = function(y,X,V,nfold = 10,folds=NULL,b,b_sd){
#   n = nrow(y)
#   n_bulk = ncol(y)
#   K = ncol(X)
#
#   d_Vg = ncol(V)
#
#
#
#   res = matrix(nrow=n,ncol=n_bulk)
#   if(is.null(folds)){
#     folds = cut(1:n,breaks = nfold,labels = F)
#   }
#
#
#   for(f in 1:nfold){
#     idx = which(folds==f)
#     X.temp = X[-idx,]
#     y.temp = y[-idx,]
#     if(d_Vg==K^2){
#       V.temp = matrix(c(colSums(V[-idx,])),ncol = K)
#     }else if(d_Vg==K){
#       V.temp = diag(c(colSums(V[-idx,])))
#     }else{
#       stop('check dimension of V')
#     }
#
#     if(!is.null(b)){
#       bhat = (b+rnorm(length(b),0,b_sd))%*%t(rep(1,n_bulk))
#     }else{
#       bhat = pmax(solve(crossprod(X.temp)-V.temp)%*%t(X.temp)%*%y.temp,0)
#     }
#
#
#     res[idx,] = y[idx,] - X[idx,]%*%bhat
#   }
#
#   res
#
# }






####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################







#'
#'
#' #'@param y a vector of bulk sample
#' #'@param X reference matrix, estimated relative expression, gene by cell
#' #'@param Vg variance matrix of X
#' #'@param X_var_pop variance of true X among individuals
#' #'@param marker_gene
#' #'@param w gene weights
#' #'@param hc.type hc3, hc2, or hc0
#' #'@param alpha significance level
#' #'@param a alpha in the Fuller's small sample correction
#' #'@param correction whether perform fuller's small sample correction.
#' #'@param S cell size
#'
#'
#' estimation_func = function(y,X,Vg,X_var_pop=NULL,
#'                            #marker_gene = NULL,
#'                            w=NULL,
#'                            hc.type='hc3',a=ncol(X)+4,
#'                            correction=TRUE,S=NULL,calc_cov=TRUE,
#'                            verbose=FALSE){
#'
#'   #browser()
#'
#'   #print(X)
#'
#'   if(length(w)==1){
#'     w = rep(w,nrow(X))
#'   }
#'
#'   # if(!is.null(marker_gene)){
#'   #   gene_idx = match(marker_gene,rownames(X))
#'   #   omit.idx = which(is.na(gene_idx))
#'   #   if(length(omit.idx)>0){
#'   #     gene_idx = gene_idx[-omit.idx]
#'   #   }
#'   #
#'   #   #browser()
#'   #
#'   #   X = X[gene_idx,]
#'   #   Vg = Vg[gene_idx,]
#'   #   if(!is.null(w)){
#'   #     w = w[gene_idx]
#'   #   }
#'   #   y = y[gene_idx,]
#'   #
#'   #   if(!is.null(X_var_pop)){
#'   #     X_var_pop = X_var_pop[gene_idx,]
#'   #   }
#'   #
#'   # }
#'
#'   G = nrow(X)
#'   K = ncol(X)
#'
#'
#'
#'   # nb is the number of bulk sample
#'   nb = ncol(y)
#'   if(is.null(nb)){
#'     nb=1
#'   }
#'
#'
#'   # generate weights
#'
#'   if(is.null(w)){
#'     if(is.null(X_var_pop)){
#'       w  = 1/(rowSums(X)/K+rowSums(Vg)/K^2)
#'     }else{
#'       w  = 1/(rowSums(X)/K+rowSums(X_var_pop)/K^2+rowSums(Vg)/K^2)
#'     }
#'   }
#'   w = w/sum(w)*G
#'
#'
#'
#'   ## take cell size into account, and scale it to O(G)
#'   if(is.null(S)){
#'     S = rep(1,K)
#'   }
#'
#'   # X, Vg not scaled by G
#'   input = list(X=X,Vg=Vg,y=y,w=w,Sigma=X_var_pop,S=S)
#'
#'   # scale X by G, mainly for numerical stability
#'   #X = X*G
#'   #Vg = Vg*G^2
#'
#'   #input = list(X=X,Vg=Vg,y=y,w=w,Sigma=X_var_pop,S=S)
#'
#'   # add library size to X
#'   X = X%*%diag(S)
#'   Xw = X*sqrt(w)
#'   yw = cbind(y*sqrt(w))
#'   A = t(Xw)%*%Xw
#'   d_Vg = ncol(Vg)
#'   Vw = Vg*w
#'   if(d_Vg==K^2){
#'     V = matrix(c(colSums(Vw*c(S%*%t(S)))),ncol = K)
#'   }else if(d_Vg==K){
#'     Vw = Vw%*%diag(S^2)
#'     V = diag(c(colSums(Vw)))
#'   }else{
#'     stop('check dimension of Vg')
#'   }
#'
#'
#'
#'   # else if(is.array(Vg)){
#'   #   Vw = Vg*rep(w,each=K*K)
#'   #   V = rowSums(Vw,dims = 2)
#'   # }else{
#'   #   stop("Vg input type not allowed")
#'   # }
#'
#'
#'
#'   #browser()
#'   if(verbose){
#'     message("estimating proportions")
#'   }
#'   beta_tilde_hat = pmax(solve(A-V)%*%t(Xw)%*%yw,0)
#'   ## Fuller's correction
#'   if(correction){
#'
#'     if(verbose){
#'       message("performing fuller's correction")
#'     }
#'
#'     #Z = cbind(yw,Xw)
#'     #m_zz = t(Z)%*%Z/G
#'     m_xy = t(Xw)%*%yw
#'     lambda = c()
#'     A_inv = array(dim=c(K,K,nb))
#'
#'     S = matrix(0,nrow = K+1,ncol=K+1)
#'     S[-1,-1] = V
#'
#'     #browser()
#'     for(ib in 1:nb){
#'
#'       #print(A)
#'
#'       #browser()
#'       #lambda[ib] = min(Re(eigen((A-m_xy[,ib]%*%t(m_xy[,ib])/sum(yw[,ib]^2))%*%solve(V))$values))
#'
#'       Si = S
#'       Si[1,1] = 1/sum(yw[,ib])
#'       lambda[ib] = min(eigen(crossprod(cbind(yw[,ib],Xw))%*%solve(Si))$values)
#'
#'
#'
#'       if(lambda[ib]>(1+1/G)){
#'         A_inv[,,ib] = solve(A - (1-a/G)*V)
#'       }else{
#'         A_inv[,,ib] = solve(A - (lambda[ib]-1/G-a/G)*V)
#'       }
#'
#'       beta_tilde_hat[,ib] = pmax(A_inv[,,ib]%*%t(Xw)%*%yw[,ib],0)
#'
#'     }
#'   }else{
#'     lambda = rep(1,nb)
#'     A_inv = array(solve(A-V),dim=c(K,K,nb))
#'   }
#'
#'   # M = crossprod(cbind(yw,Xw))
#'   # B = diag(c(1/colSums(yw),1/diag(V)))
#'   # lambda = min(eigen(M%*%B)$values)
#'
#'   #   if(lambda>(1+1/G)){
#'   #     A = A - (1-a/G)*V
#'   #   }else{
#'   #     A = A - (lambda-1/G-a/G)*V
#'   #   }
#'   # }else{
#'   #   A = A - V
#'   # }
#'   # A_inv = solve(A)
#'   #
#'   # beta_tilde_hat = pmax(A_inv%*%t(Xw)%*%yw,0)
#'
#'   #browser()
#'
#'   Sigma = matrix(0,nrow=nb*K,ncol=nb*K)
#'   Q_inv = matrix(0,nrow=nb*K,ncol=nb*K)
#'   Sigma_ii = matrix(0,nrow=nb*K,ncol=nb*K)
#'
#'   #browser()
#'
#'
#'   if(verbose){
#'     message("calculating covariance matrix")
#'   }
#'
#'   for(i in 1:nb){
#'
#'     if(verbose){
#'       message(paste("...covaraince",i,'to',nb))
#'     }
#'
#'     Q_inv[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = A_inv[,,i]
#'
#'     if(d_Vg==K^2){
#'       Vbi = t(apply(Vw,1,function(z){v = matrix(z,ncol=K);v%*%beta_tilde_hat[,i]}))*pmin(lambda,1)[i]
#'     }
#'     if(d_Vg==K){
#'       Vbi = (Vw)%*%diag(c(beta_tilde_hat[,i]))*pmin(lambda,1)[i]
#'     }
#'
#'     Hi = t(t(X%*%A_inv[,,i]%*%t(X))*w)
#'     hi = diag(Hi)
#'     ri = y[,i] - Hi%*%y[,i]
#'
#'
#'     for(j in i:nb){
#'
#'       if(j==i){
#'
#'         Vbj = Vbi
#'
#'         if(hc.type == 'hc0'){
#'           Sigma_ij = crossprod(c(ri)*w*X+Vbi)
#'         }else if(hc.type == 'hc2'){
#'           Sigma_ij = crossprod(c(ri)/sqrt(1-pmax(pmin(hi,1-1/G),0))*w*X+Vbi)
#'         }else if(hc.type == 'hc3'){
#'           Sigma_ij = crossprod(c(ri)/(1-pmax(pmin(hi,1-1/G),0))*w*X+Vbi)
#'         }
#'
#'         Sigma_ii[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = Sigma_ij
#'
#'       }else{
#'
#'         if(calc_cov){
#'
#'           if(d_Vg==K^2){
#'             Vbj = t(apply(Vw,1,function(z){v = matrix(z,ncol=K);v%*%beta_tilde_hat[,j]}))*pmin(lambda,1)[j]
#'           }else if(d_Vg==K){
#'             Vbj = (Vw)%*%diag(c(beta_tilde_hat[,j]))*pmin(lambda,1)[j]
#'           }
#'
#'           Hj = t(t(X%*%A_inv[,,j]%*%t(X))*w)
#'           hj = diag(Hj)
#'           rj = y[,j] - Hi%*%y[,j]
#'           if(hc.type == 'hc0'){
#'             Sigma_ij = crossprod(c(ri)*w*X+Vbi,c(rj)*w*X+Vbj)
#'           }else if(hc.type == 'hc2'){
#'             Sigma_ij = crossprod(c(ri)/sqrt(1-pmax(pmin(hi,1-1/G),0))*w*X+Vbi,c(rj)/sqrt(1-pmax(pmin(hj,1-1/G),0))*w*X+Vbj)
#'           }else if(hc.type == 'hc3'){
#'             Sigma_ij = crossprod(c(ri)/(1-pmax(pmin(hi,1-1/G),0))*w*X+Vbi,c(rj)/(1-pmax(pmin(hj,1-1/G),0))*w*X+Vbj)
#'           }
#'
#'         }
#'
#'       }
#'
#'       Sigma[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)] = Sigma_ij
#'
#'     }
#'   }
#'   Sigma = Sigma+t(Sigma)-Sigma_ii
#'   covb = Q_inv%*%Sigma%*%Q_inv
#'
#'   beta_tilde_hat = cbind(beta_tilde_hat)
#'
#'
#'   # delta method
#'
#'   # covb is a (nb*K)*(nb*K) cov matrix
#'   # formulate Jacobian matrix
#'
#'   if(verbose){
#'     message("performing delta method")
#'   }
#'
#'   J = matrix(0,nrow=(nb*K),ncol=(nb*K))
#'   for(i in 1:nb){
#'     J[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = J_sum2one(beta_tilde_hat[,i],K)
#'   }
#'
#'   asyV = (J)%*%covb%*%t(J)
#'
#'   beta_hat = apply(beta_tilde_hat,2,function(z){z/sum(z)})
#'   beta_se = sqrt(diag(asyV))
#'   beta_se = matrix(beta_se,ncol=nb)
#'
#'   if(verbose){
#'     message("done")
#'   }
#'
#'   #browser()
#'   # rownames(beta_hat) = colnames(X)
#'
#'   return(list(beta_hat=beta_hat,
#'               beta_se=beta_se,
#'               cov_beta_hat = asyV,
#'               beta_tilde_hat=beta_tilde_hat,
#'               cov_beta_tilde_hat = covb,
#'               Sigma = Sigma,
#'               Sigma_ii = Sigma_ii,
#'               J=J,
#'               Q_inv = Q_inv,
#'               lambda=lambda,
#'               input = input))
#'
#'
#'
#'   # if(is.null(nb) | nb==1){
#'   #   J = J_sum2one(beta_tilde_hat,K)
#'   #
#'   #   asyV = (J)%*%covb%*%t(J)
#'   #
#'   #   beta_hat = beta_tilde_hat/sum(beta_tilde_hat)
#'   #   beta_se = sqrt(diag(asyV))
#'   #
#'   #   ci = rbind(pmax(0,beta_hat - qnorm(1-alpha/2)*beta_se), pmin(1,beta_hat + qnorm(1-alpha/2)*beta_se))
#'   #
#'   #   return(list(beta_tilde_hat=beta_tilde_hat,
#'   #               beta_hat=beta_hat,
#'   #               beta_se=beta_se,
#'   #               ci=ci,
#'   #               cov_beta_hat = asyV,
#'   #               input = input))
#'   # }else{
#'   #
#'   # }
#'
#'
#' }
#'
#'
#'
#'
#' #'@description allow correlations among samples.
#' get_SIGMA = function(y,X,w,beta,V,h,nb,G,K,lambda,verbose,calc_cov,hc.type,R){
#'
#'
#'
#'   Sigma = matrix(0,nrow=nb*K,ncol=nb*K)
#'   #Q_inv = matrix(0,nrow=nb*K,ncol=nb*K)
#'   Sigma_ii = matrix(0,nrow=nb*K,ncol=nb*K)
#'
#'   res = y - X%*%beta
#'
#'   #browser()
#'
#'   if(!is.null(R)){
#'     cor.idx = which(R!=0,arr.ind = T)
#'     cor.idx = t(apply(cor.idx,1,sort))
#'     cor.idx = cor.idx[!duplicated(cor.idx),]
#'     cor.idx = cor.idx[(cor.idx[,1]!=cor.idx[,2]),]
#'   }
#'
#'
#'   if(verbose){
#'     message("calculating covariance matrix")
#'   }
#'
#'
#'   h = pmax(pmin(h,1-1/G),0)
#'   for(i in 1:nb){
#'
#'     if(verbose){
#'       message(paste("...covaraince",i,'to',nb))
#'     }
#'
#'     if(ncol(V)==K^2){
#'       Vbi = t(apply(V,1,function(z){v = matrix(z,ncol=K);v%*%beta[,i]}))*lambda
#'     }
#'     if(ncol(V)==K){
#'       Vbi = (V)%*%diag(c(beta[,i]))*lambda
#'     }
#'     ri = res[,i]
#'
#'
#'     for(j in i:nb){
#'
#'       if(j==i){
#'
#'         if(hc.type == 'hc0'){
#'           Sigma_ij = crossprod(c(ri)*w*X+Vbi)
#'         }else if(hc.type == 'hc2'){
#'           Sigma_ij = crossprod(c(ri)/sqrt(1-h)*w*X+Vbi)
#'         }else if(hc.type == 'hc3'){
#'           Sigma_ij = crossprod(c(ri)/(1-h)*w*X+Vbi)
#'         }
#'
#'         if(!is.null(R)){
#'
#'           l1.temp = (c(ri)*w*X+Vbi)[cor.idx[,1],,drop=FALSE]
#'           l2.temp = (c(ri)*w*X+Vbi)[cor.idx[,2],,drop=FALSE]
#'           Sigma_ij = Sigma_ij + crossprod(l1.temp,l2.temp) + crossprod(l2.temp,l1.temp)
#'
#'         }
#'
#'         ###########
#'         # if(!is.null(R)){
#'         #   print(paste('adjusted for corr:',sum(Sigma_ij)))
#'         # }else{
#'         #   print(paste('not adjusted for corr:',sum(Sigma_ij)))
#'         # }
#'
#'         ############
#'
#'         Sigma_ii[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = Sigma_ij
#'
#'       }else{
#'
#'         if(calc_cov){
#'
#'           if(ncol(V)==K^2){
#'             Vbj = t(apply(V,1,function(z){v = matrix(z,ncol=K);v%*%beta[,j]}))*lambda
#'           }
#'           if(ncol(V)==K){
#'             Vbj = (V)%*%diag(c(beta[,j]))*lambda
#'           }
#'
#'           rj = res[,j]
#'           if(hc.type == 'hc0'){
#'             Sigma_ij = crossprod(c(ri)*w*X+Vbi,c(rj)*w*X+Vbj)
#'           }else if(hc.type == 'hc2'){
#'             Sigma_ij = crossprod(c(ri)/sqrt(1-h)*w*X+Vbi,c(rj)/sqrt(1-h)*w*X+Vbj)
#'           }else if(hc.type == 'hc3'){
#'             Sigma_ij = crossprod(c(ri)/(1-h)*w*X+Vbi,c(rj)/(1-h)*w*X+Vbj)
#'           }
#'
#'           if(!is.null(R)){
#'             Sigma_ij = Sigma_ij
#'             + crossprod((c(ri)*w*X+Vbi)[cor.idx[,1],,drop=FALSE],(c(rj)*w*X+Vbj)[cor.idx[,2],,drop=FALSE])
#'             + crossprod((c(ri)*w*X+Vbi)[cor.idx[,2],,drop=FALSE],(c(rj)*w*X+Vbj)[cor.idx[,1],,drop=FALSE])
#'           }
#'
#'         }
#'
#'       }
#'
#'       Sigma[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)] = Sigma_ij
#'
#'     }
#'   }
#'   Sigma = (Sigma+t(Sigma)-Sigma_ii)
#'
#'   Sigma
#'
#' }
#'
#'
#'
