

#'@title generalized least square estimator
#'@param y a vector of bulk sample
#'@param Y A matrix, gene by cell type
#'@param w weights for gene
#'@param adjust whether adjust for uncertainty in Yr
#'@param alpha significance level
#'@param a alpha in the Fuller's small sample correction
#'@param correction whether perform fuller's small sample correction.
#'@return p_tilde_hat: before delta method; p_hat: estimate of proportions; p_se: standard error of estimator; ci: confidence intervals

poi_gls = function(y,Y,w=NULL,
                   adjust=TRUE,
                   alpha=0.05,a=0,
                   correction=FALSE){

  G = nrow(Y)
  K = ncol(Y)

  U = diag(c(colSums(Y)))
  U_inv = diag(c(1/colSums(Y)))
  Yu = apply(Y,2,function(x){x/sum(x)})
  Yuu = Yu%*%U_inv
  A = t(Yu)%*%Yu
  if(adjust){

    ## diag(colSums(Yuu)) = U_inv

    ## Fuller's correction for negative matrix

    if(correction){
      M = t(cbind(y,Yu))%*%cbind(y,Yu)
      B = diag(c(1/sum(y),colSums(Y)))
      #print((B))
      #lambda = min(eigen(A%*%diag(1/colSums(Yuu)))$values)
      #print(eigen(A%*%U)$values)
      lambda = min(eigen(M%*%B)$values)
      #print(lambda)
      #print(paste('est',lambda))
      #print(paste('value',det(A-lambda*U_inv)))
      if(lambda>(1+1/G)){
        A = A - (1-a/G)*diag(c(colSums(Yuu)))
      }else{
        A = A - (lambda-1/G-a/G)*diag(c(colSums(Yuu)))
      }
    }else{
      A = A - diag(c(colSums(Yuu)))
    }
    p_tilde_hat = solve(A)%*%t(Yu)%*%y
  }else{
    p_tilde_hat = solve(A)%*%t(Yu)%*%y
  }

  #print(solve(A))

  ## truncate p_tilde_hat
  p_tilde_hat = pmax(p_tilde_hat,0)

  #Q = 0
  Sigma=0

  if(adjust){
    for(i in 1:G){
      ag = tcrossprod(Yu[i,])-diag(Yuu[i,])
      #Q = Q + ag
      nabla = (ag%*%p_tilde_hat-y[i]*Yu[i,])
      Sigma = Sigma + tcrossprod(nabla)
    }
  }else{
    for(i in 1:G){
      ag = tcrossprod(Yu[i,])
      #Q = Q + ag
      nabla = (ag%*%p_tilde_hat-y[i]*Yu[i,])
      Sigma = Sigma + tcrossprod(nabla)
    }
  }

  Q = A/G
  Sigma = Sigma/G

  J = - (p_tilde_hat)%*%t(rep(1,K))
  diag(J) = diag(J) + sum(p_tilde_hat)
  J = J/sum(p_tilde_hat)^2

  asyV = ((J)%*%solve(Q)%*%Sigma%*%solve(Q)%*%t(J))/G

  p_hat = p_tilde_hat/sum(p_tilde_hat)

  p_se = sqrt(diag(asyV))

  ci_left = pmax(0,p_hat - qnorm(1-alpha/2)*p_se)
  ci_right = pmin(1,p_hat + qnorm(1-alpha/2)*p_se)

  return(list(p_tilde_hat=p_tilde_hat,
              p_hat=p_hat,
              p_se=p_se,
              ci_left=ci_left,
              ci_right=ci_right))

}





#'@title This simulation averages over parameter b.

#'@param ref ref data matrix, Gene by cell type
#'@param s cell size
#'@param Ng number of genes to use in simulation, gene will be resampled each round
#'@param b cell type proportions
#'@param bulk_lib_size bulk sample library size
#'@param ref_lib_size reference sample library size, here reference data are also bulk data
#'@param nreps number of repetition.
#'@param alpha significancel level
#'@return a list with estimates and standard errors.
#'@return 1. true proportion; 2. mean of estimates; 3. standard error of estimates; 4. coverage.
#'@return input parameters


simu_study = function(ref,s,Ng,b,
                      bulk_lib_size = 50,
                      ref_lib_size = 30,
                      nreps=100,alpha=0.05,a=0,correction=FALSE){

  G = nrow(ref)
  K = ncol(ref)


  if(missing(s)){
    s = rep(1,K)
  }

  b = b/sum(b)
  p = (b*s)/sum(b*s)


  est_adj = matrix(nrow=nreps,ncol=K)
  se_adj = matrix(nrow=nreps,ncol=K)

  est_unadj = matrix(nrow=nreps,ncol=K)
  se_unadj = matrix(nrow=nreps,ncol=K)



  for(rep in 1:nreps){

    #Obtain Theta
    ref_rep = ref[sample(1:G,Ng),]
    Theta = apply(ref_rep,2,function(z){z/sum(z)})

    # bulk data gene relative expression.
    Xb = Theta%*%diag(s*(Ng/G))%*%b
    thetab = Xb/sum(Xb)

    #reference data
    #Cr = rpois(K,ref_lib_size)
    #if(prod(Cr)==0){Cr[which(Cr==0)]=ref_lib_size}
    Cr = rep(ref_lib_size,K)
    U = diag(Cr*Ng)
    Y = matrix(rpois(Ng*K,Theta%*%U),ncol=K)
    #bulk data
    y = rpois(Ng,rpois(1,bulk_lib_size)*Ng*thetab)


    fit_adj = poi_gls(y,Y,adjust = TRUE,alpha=alpha,a=a,correction=correction)
    est_adj[rep,] = fit_adj$p_hat
    se_adj[rep,] = fit_adj$p_se

    fit_unadj = poi_gls(y,Y,adjust = FALSE,alpha=alpha)
    est_unadj[rep,] = fit_unadj$p_hat
    se_unadj[rep,] = fit_unadj$p_se

  }


  ########
  #1. mean squared error, and standard error
  mean_est_adj = apply(est_adj,2,mean)
  mse_adj = apply((est_adj - rep(1,nreps)%*%t(p))^2,2,mean)
  se_est_adj = apply(est_adj,2,sd)
  mean_est_unadj = apply(est_unadj,2,mean)
  mse_unadj = apply((est_unadj - rep(1,nreps)%*%t(p))^2,2,mean)
  se_est_unadj = apply(est_unadj,2,sd)
  #mean_se_adj = apply(se_adj,2,mean)
  #mean_se_unadj = apply(se_unadj,2,mean)

  ###################
  #2. coverage
  ci_l = est_adj - qnorm(1-alpha/2)*se_adj
  ci_r = est_adj + qnorm(1-alpha/2)*se_adj
  covergae_adj = ((rep(1,nreps)%*%t(p))>=ci_l) & ((rep(1,nreps)%*%t(p))<=ci_r)
  covergae_adj=apply(covergae_adj,2,mean)

  ci_l = est_unadj - qnorm(1-alpha/2)*se_unadj
  ci_r = est_unadj + qnorm(1-alpha/2)*se_unadj
  covergae_unadj = ((rep(1,nreps)%*%t(p))>=ci_l) & ((rep(1,nreps)%*%t(p))<=ci_r)
  covergae_unadj=apply(covergae_unadj,2,mean)

  return(list(p=p,
              mean_est_adj=mean_est_adj,
              mse_adj=mse_adj,
              se_est_adj=se_est_adj,
              mean_est_unadj=mean_est_unadj,
              mse_unadj=mse_unadj,
              se_est_unadj=se_est_unadj,
              covergae_adj=covergae_adj,
              covergae_unadj=covergae_unadj,
              se_adj = se_adj,
              se_unadj = se_unadj,
              #mean_se_adj=mean_se_adj,
              #mean_se_unadj=mean_se_unadj,
              est_adj = est_adj,
              est_unadj = est_unadj,
              simu_param = list(ref=ref,s=s,Ng=Ng,b=b,
                                ref_lib_size = ref_lib_size,
                                bulk_lib_size = bulk_lib_size,
                                nreps=nreps,alpha=alpha)))

}

####################################
############## test ################
####################################

set.seed(12345)
G = 300
K = 4
b = 1:K
b = b/sum(b)
library(gtools)
X = t(rdirichlet(K,rep(1,G)))

s=rep(1,K)
results = simu_study(X,s,G,b,ref_lib_size = 0.5,correction = F)

