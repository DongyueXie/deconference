
#'@description generate multivariate log-normal random variable.
#'@param n number of samples
#'@param mu mean of the log normal distributed variable
#'@param s vector, standard deviations of the log normal distributed variable
#'@param R correlation matrix
sim_MLN = function(n,mu,s,R){
  n.mu = log(mu^2/sqrt(mu^2+s^2))
  n.s = sqrt(log(1+s^2/mu^2))
  n.Sigma = t(R*n.s)*n.s
  exp(mvnfast::rmvn(n,mu = n.mu,sigma = n.Sigma))
}


is.identity = function(X){
  (sum(X)==nrow(X))
}


#'@description For simplicity, we only simulate data at individual level. - only generate X from U

#'@param R correlation matrix, sparse form
simu_corr_simple = function(ref,
                            b,
                            R,
                            sigma2,
                            nreps = 100,
                            bulk_lib_size = 500,
                            n_indi = 8,
                            printevery=10,
                            verbose=FALSE,
                            alpha=0.05,
                            groups = c(rep(1,ncol(b)/2),rep(2,ncol(b)/2)),
                            centeringXY = FALSE,
                            true.beta.for.Sigma=FALSE,
                            calc_cov = T){

  G = nrow(ref)
  K = ncol(ref)
  n_bulk = ncol(b)
  b = apply(b,2,function(z){z/sum(z)})
  is.indep = is.identity(R)

  # est_adj = matrix(nrow=nreps,ncol=n_bulk*K)
  # se_adj = matrix(nrow=nreps,ncol=n_bulk*K)
  #
  # est_unadj = matrix(nrow=nreps,ncol=n_bulk*K)
  # se_unadj_hc0 = matrix(nrow=nreps,ncol=n_bulk*K)
  # #se_unadj_cv = matrix(nrow=nreps,ncol=n_bulk*K)
  # se_unadj_hc3 = matrix(nrow=nreps,ncol=n_bulk*K)

  data_sparsity = c()

  est_adj = matrix(nrow=nreps,ncol=n_bulk*K)
  se_adj_hc0 = matrix(nrow=nreps,ncol=n_bulk*K)
  se_adj_hc3 = matrix(nrow=nreps,ncol=n_bulk*K)

  beta_hat = matrix(nrow=nreps,ncol=n_bulk*K)
  beta_se_adj_hc0 = matrix(nrow=nreps,ncol=n_bulk*K)
  beta_se_adj_hc3 = matrix(nrow=nreps,ncol=n_bulk*K)

  diff_adj = matrix(nrow=nreps,ncol=K)
  diff_adj_se_hc0 = matrix(nrow=nreps,ncol=K)
  diff_adj_p_hc0 = matrix(nrow = nreps,ncol=K)

  diff_adj_se_hc3 = matrix(nrow=nreps,ncol=K)
  diff_adj_p_hc3 = matrix(nrow = nreps,ncol=K)

  est_unadj = matrix(nrow=nreps,ncol=n_bulk*K)
  se_unadj_hc0 = matrix(nrow=nreps,ncol=n_bulk*K)
  diff_unadj = matrix(nrow=nreps,ncol=K)
  diff_unadj_se_hc0 = matrix(nrow=nreps,ncol=K)
  diff_unadj_p_hc0 = matrix(nrow = nreps,ncol=K)


  se_unadj_cv = matrix(nrow=nreps,ncol=n_bulk*K)
  diff_unadj_se_cv = matrix(nrow=nreps,ncol=K)
  diff_unadj_p_cv = matrix(nrow = nreps,ncol=K)

  se_unadj_hc3 = matrix(nrow=nreps,ncol=n_bulk*K)
  diff_unadj_se_hc3 = matrix(nrow=nreps,ncol=K)
  diff_unadj_p_hc3 = matrix(nrow = nreps,ncol=K)

  gene_names = rownames(ref)
  true_betas = matrix(nrow=nreps,ncol=n_bulk*K)


  ## pre calculate MLN and generate independent normals



  n.ref = matrix(nrow=G,ncol=K)
  n.Sigma.chol = list()
  n.Sigma = matrix(nrow=G,ncol=K)
  for(k in 1:K){
    n.ref[,k] = log(ref[,k]^2/sqrt(ref[,k]^2+sigma2[,k]))
    n.s = sqrt(log(1+sigma2[,k]^2/ref[,k]^2))
    if(!is.indep){
      n.Sigma.chol[[k]] = t(n.s*t(chol(R)))
    }else{
      n.Sigma[,k] = n.s^2
    }
  }


  if(!is.indep){
    cor.idx = which(R!=0,arr.ind = T)
    #cor.idx = t(apply(cor.idx,1,sort))
    #cor.idx = cor.idx[!duplicated(cor.idx),]
    cor.idx = cor.idx[(cor.idx[,1]!=cor.idx[,2]),]
  }else{
    cor.idx = NULL
  }


  for(reps in 1:nreps){

    if(reps%%printevery==0){print(sprintf("running %d (out of %d)",reps,nreps))}



    # generate individual reference matrices


    X_array = array(dim=c(G,K,n_indi+n_bulk))

    for(k in 1:K){
      if(is.indep){
        X_array[,k,] = exp(matrix(rnorm(G*(n_indi+n_bulk),n.ref[,k],sqrt(n.Sigma[,k])),ncol=n_indi+n_bulk))
      }else{
        X_array[,k,] = t(exp(mvnfast::rmvn(n_indi+n_bulk,mu = n.ref[,k],sigma = n.Sigma.chol[[k]],isChol = TRUE)))
      }
    }

    X_array_bulk = X_array[,,1:n_bulk]
    X_array = X_array[,,-(1:n_bulk)]



    # generate bulk data

    # mb = 0
    # for(k in 1:K){
    #   if(!is.identity(R)){
    #     mb = mb + b[k,]*exp(mvnfast::rmvn(n_bulk,mu = n.ref[,k],sigma = n.Sigma.chol[[k]],isChol = TRUE))
    #   }else{
    #     mb = mb + b[k,]*exp(mvnfast::rmvn(n_bulk,mu = n.ref[,k],sigma = n.Sigma.chol[[k]],isChol = TRUE))
    #   }
    #
    #   #mb = mb + b[k,]*sim_MLN(n_bulk,ref[,k],sqrt(sigma2[,k]),R)
    # }

    mb = lapply(1:n_bulk,function(i){X_array_bulk[,,i]%*%b[,i]})

    mb = do.call(cbind,mb)
    true.beta = t(t(b)*c(apply(mb,2,function(z){bulk_lib_size*G/sum(z)})))
    true_betas[reps,] = c(true.beta)
    thetab = apply(mb,2,function(z){z/sum(z)})

    # mb = matrix(nrow=G,ncol=n_bulk)
    # for(i in 1:n_bulk){
    #
    #   X_b = matrix(nrow=G,ncol=K)
    #   for(k in 1:K){
    #     X_b[,k] = sim_MLN(1,ref[,k],sqrt(sigma2[,k]),R)
    #   }
    #
    #   mb[,i] = X_b%*%b[,i]
    #
    #
    # }
    # thetab = apply(mb,2,function(z){z/sum(z)})


    #browser()
    y = matrix(rpois(G*n_bulk,bulk_lib_size*G*thetab),nrow=G)
    rownames(y) = gene_names

    #bulks = SingleCellExperiment(assays = list(counts = y),
    #                             colData = DataFrame(individual = 1:n_bulk))



    ## generate reference matrices X

    # X_array = array(dim=c(G,K,n_indi))
    #
    # for(k in 1:K){
    #   x_k = sim_MLN(n_indi,ref[,k],sqrt(sigma2[,k]),R)
    #   X_array[,k,] = t(x_k)
    # }



    # fit model

    #browser()

    X = apply(X_array,c(1,2),mean,na.rm=TRUE)
    V = t(apply(X_array,c(1),function(z){(cov(t(z),use = 'pairwise.complete.obs'))}))/n_indi


    fit.adj.hc0 = estimation_func2(y=y,X=X,Vg=V,
                                   w=1,hc.type='hc0',correction=FALSE,
                                   calc_cov=calc_cov,verbose=verbose,
                                   cor.idx=cor.idx,
                                   centeringXY=centeringXY,
                                   true.beta = if(true.beta.for.Sigma){true.beta}else{NULL})

    # fit.adj.hc3 = estimation_func2(y=y,X=X,Vg=V,
    #                                w=1,hc.type='hc3',correction=FALSE,
    #                                calc_cov=calc_cov,verbose=verbose,
    #                                cor.idx=cor.idx,
    #                                centeringXY=centeringXY,
    #                                true.beta = if(true.beta.for.Sigma){true.beta}else{NULL})

    fit.unadj.hc0 = estimation_func2(y=y,X=X,Vg=V,
                                     w=1,hc.type='hc0',correction=FALSE,
                                     calc_cov=calc_cov,verbose=verbose,
                                     cor.idx=NULL,
                                     centeringXY=centeringXY,
                                     true.beta = if(true.beta.for.Sigma){true.beta}else{NULL})



    est_adj[reps,] = c(fit.adj.hc0$beta_hat)
    se_adj_hc0[reps,] = c(fit.adj.hc0$beta_se)
    #se_adj_hc3[reps,] = c(fit.adj.hc3$beta_se)
    diff_out_hc0 = two_group_test(fit.adj.hc0,groups)
    #diff_out_hc3 = two_group_test(fit.adj.hc3,groups)

    beta_hat[reps,] = c(fit.adj.hc0$beta_tilde_hat)
    beta_se_adj_hc0[reps,] = c(fit.adj.hc0$beta_tilde_se)
    #beta_se_adj_hc3[reps,] = c(fit.adj.hc3$beta_tilde_se)

    #browser()

    diff_adj[reps,] = c(diff_out_hc0$diff_group)
    diff_adj_se_hc0[reps,] = c(diff_out_hc0$diff_se)
    diff_adj_p_hc0[reps,] = c(diff_out_hc0$p_value)

    #diff_adj_se_hc3[reps,] = c(diff_out_hc3$diff_se)
    #diff_adj_p_hc3[reps,] = c(diff_out_hc3$p_value)


    #est_unadj[reps,] = c(fit_unadj$beta_hat)
    se_unadj_hc0[reps,] = c(fit.unadj.hc0$beta_se)
    #diff_unadj[reps,] = c(fit_unadj$diff_group)
    diff_out_hc0_unadj = two_group_test(fit.unadj.hc0,groups)
    diff_unadj_se_hc0[reps,] = c(diff_out_hc0_unadj$diff_se)
    diff_unadj_p_hc0[reps,] = c(diff_out_hc0_unadj$p_value)

    #se_unadj_cv[reps,] = c(fit_unadj$ols.out$beta_se)
    #diff_unadj_se_cv[reps,] = c(fit_unadj$ols.out$diff_se)
    #diff_unadj_p_cv[reps,] = fit_unadj$ols.out$p_value

    #se_unadj_hc3[reps,] = c(fit_unadj$sand.out.hc3$beta_se)
    #diff_unadj_se_hc3[reps,] = c(fit_unadj$sand.out.hc3$diff_se)
    #diff_unadj_p_hc3[reps,] = fit_unadj$sand.out.hc3$p_value



    #beta_hat_cov_adj[[reps]] = fit_adj$cov_beta_hat
    #beta_hat_cov_unadj[[reps]] = fit_unadj$sand.out$cov_beta_hat
    #beta_hat_cov_unadj_cv[[reps]] = fit_unadj$ols.out$cov_beta_hat
    #beta_hat_cov_unadj_hc3[[reps]] = fit_unadj$sand.out.hc3$cov_beta_hat

  }


  ########
  #1. mean squared error, and standard error
  mean_est_adj = apply(est_adj,2,mean)
  mse_adj = apply((est_adj - rep(1,nreps)%*%t(c(b)))^2,2,mean)
  se_est_adj = apply(est_adj,2,sd)
  #mean_est_unadj = apply(est_unadj,2,mean)
  #mse_unadj = apply((est_unadj - rep(1,nreps)%*%t(c(b)))^2,2,mean)
  #se_est_unadj = apply(est_unadj,2,sd)
  #mean_se_adj = apply(se_adj,2,mean)
  #mean_se_unadj = apply(se_unadj,2,mean)

  ###################
  #2. coverage

  ## adj

  ## hc0
  ci_l = est_adj - qnorm(1-alpha/2)*se_adj_hc0
  ci_r = est_adj + qnorm(1-alpha/2)*se_adj_hc0
  coverage_adj_hc0 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_adj_hc0=apply(coverage_adj_hc0,2,mean,na.rm=TRUE)

  # ## hc3
  # ci_l = est_adj - qnorm(1-alpha/2)*se_adj_hc3
  # ci_r = est_adj + qnorm(1-alpha/2)*se_adj_hc3
  # coverage_adj_hc3 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  # coverage_adj_hc3=apply(coverage_adj_hc3,2,mean,na.rm=TRUE)


  # unadj
  ## sand
  ci_l = est_adj - qnorm(1-alpha/2)*se_unadj_hc0
  ci_r = est_adj + qnorm(1-alpha/2)*se_unadj_hc0
  coverage_unadj_hc0 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_unadj_hc0=apply(coverage_unadj_hc0,2,mean,na.rm=TRUE)


  ## cv
  # ci_l = est_unadj - qnorm(1-alpha/2)*se_unadj_cv
  # ci_r = est_unadj + qnorm(1-alpha/2)*se_unadj_cv
  # coverage_unadj_cv = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  # coverage_unadj_cv=apply(coverage_unadj_cv,2,mean)

  ## hc3
  # ci_l = est_unadj - qnorm(1-alpha/2)*se_unadj_hc3
  # ci_r = est_unadj + qnorm(1-alpha/2)*se_unadj_hc3
  # coverage_unadj_hc3 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  # coverage_unadj_hc3=apply(coverage_unadj_hc3,2,mean)

  ## diff

  ci_l = diff_adj - qnorm(1-alpha/2)*diff_adj_se_hc0
  ci_r = diff_adj + qnorm(1-alpha/2)*diff_adj_se_hc0
  coverage_diff_adj_hc0 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
  coverage_diff_adj_hc0=apply(coverage_diff_adj_hc0,2,mean,na.rm=TRUE)

  # ci_l = diff_adj - qnorm(1-alpha/2)*diff_adj_se_hc3
  # ci_r = diff_adj + qnorm(1-alpha/2)*diff_adj_se_hc3
  # coverage_diff_adj_hc3 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
  # coverage_diff_adj_hc3=apply(coverage_diff_adj_hc3,2,mean,na.rm=TRUE)

  ci_l = diff_adj - qnorm(1-alpha/2)*diff_unadj_se_hc0
  ci_r = diff_adj + qnorm(1-alpha/2)*diff_unadj_se_hc0
  coverage_diff_unadj_hc0 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
  coverage_diff_unadj_hc0=apply(coverage_diff_unadj_hc0,2,mean,na.rm=TRUE)

  # ci_l = diff_unadj - qnorm(1-alpha/2)*diff_unadj_se_cv
  # ci_r = diff_unadj + qnorm(1-alpha/2)*diff_unadj_se_cv
  # coverage_diff_unadj_cv = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
  # coverage_diff_unadj_cv = apply(coverage_diff_unadj_cv,2,mean)
  #
  # ci_l = diff_unadj - qnorm(1-alpha/2)*diff_unadj_se_hc3
  # ci_r = diff_unadj + qnorm(1-alpha/2)*diff_unadj_se_hc3
  # coverage_diff_unadj_hc3 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
  # coverage_diff_unadj_hc3 = apply(coverage_diff_unadj_hc3,2,mean)

  power_adj_hc0 = apply(diff_adj_p_hc0<=alpha,2,mean,na.rm=TRUE)
  # power_adj_hc3 = apply(diff_adj_p_hc3<=alpha,2,mean,na.rm=TRUE)
  power_unadj_hc0 = apply(diff_unadj_p_hc0<=alpha,2,mean,na.rm=TRUE)
  #power_unadj_cv = apply(diff_unadj_p_cv<=alpha,2,mean)
  #power_unadj_hc3 = apply(diff_unadj_p_hc3<=alpha,2,mean)


  return(list(b=b,
              coverage_diff_adj_hc0=coverage_diff_adj_hc0,
              #coverage_diff_adj_hc3=coverage_diff_adj_hc3,
              coverage_diff_unadj_hc0 = coverage_diff_unadj_hc0,
              #coverage_diff_unadj_cv = coverage_diff_unadj_cv,
              #coverage_diff_unadj_hc3 = coverage_diff_unadj_hc3,
              #beta_hat_cov_adj=beta_hat_cov_adj,
              #beta_hat_cov_unadj = beta_hat_cov_unadj,
              #beta_hat_cov_unadj_cv = beta_hat_cov_unadj_cv,
              #beta_hat_cov_unadj_hc3 = beta_hat_cov_unadj_hc3,
              power_adj_hc0=power_adj_hc0,
              #power_adj_hc3=power_adj_hc3,
              power_unadj_hc0=power_unadj_hc0,
              #power_unadj_cv=power_unadj_cv,
              #power_unadj_hc3=power_unadj_hc3,

              diff_adj_p_hc0 = diff_adj_p_hc0,
              #diff_adj_p_hc3 = diff_adj_p_hc3,

              diff_unadj_p_hc0=diff_unadj_p_hc0,
              #diff_unadj_p_cv=diff_unadj_p_cv,
              #diff_unadj_p_hc3=diff_unadj_p_hc3,

              diff_adj=diff_adj,
              #diff_unadj=diff_unadj,

              diff_adj_se_hc0=diff_adj_se_hc0,
              #diff_adj_se_hc3=diff_adj_se_hc3,

              diff_unadj_se_hc0=diff_unadj_se_hc0,
              #diff_unadj_se_cv=diff_unadj_se_cv,
              #diff_unadj_se_hc3=diff_unadj_se_hc3,



              mean_est_adj=mean_est_adj,
              mse_adj=mse_adj,
              se_est_adj=se_est_adj,

              # mean_est_unadj=mean_est_unadj,
              # mse_unadj=mse_unadj,
              # se_est_unadj=se_est_unadj,

              coverage_adj_hc0=coverage_adj_hc0,
              #coverage_adj_hc3=coverage_adj_hc3,
              coverage_unadj_hc0=coverage_unadj_hc0,
              #coverage_unadj_cv=coverage_unadj_cv,
              #coverage_unadj_hc3=coverage_unadj_hc3,

              se_adj_hc0 = se_adj_hc0,
              #se_adj_hc3 = se_adj_hc3,
              se_unadj_hc0 = se_unadj_hc0,
              #se_unadj_cv = se_unadj_cv,
              est_adj = est_adj,

              beta_hat = beta_hat,
              beta_se_adj_hc0=beta_se_adj_hc0,
              #beta_se_adj_hc3=beta_se_adj_hc3,
              true_betas = true_betas,
              #est_unadj = est_unadj,
              simu_param = list(ref=ref,b=b,

                                bulk_lib_size = bulk_lib_size,
                                nreps=nreps,alpha=alpha)))



}


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

















#'@param ref reference matrix
#'@param b K by #bulk matrix
#'

simu_correlation = function(ref,b,
                            R,
                          nreps = 100,
                          bulk_lib_size = 500,
                          sc_lib_size = 0.1,
                          nk = 100,
                          tau2,
                          sigma2,
                          expr_dist = 'log-normal',
                          x_estimator = 'separate',
                          n_indi = 8,
                          est_pop_var = FALSE,
                          meta_var='adjust',
                          meta_mode = 'local',
                          filter.gene = FALSE,
                          known_var=FALSE,
                          printevery=10,
                          verbose=FALSE,
                          alpha=0.05,
                          groups = c(rep(1,ncol(b)/2),rep(2,ncol(b)/2))){

  G = nrow(ref)
  K = ncol(ref)
  n_bulk = ncol(b)
  b = apply(b,2,function(z){z/sum(z)})


  # est_adj = matrix(nrow=nreps,ncol=n_bulk*K)
  # se_adj = matrix(nrow=nreps,ncol=n_bulk*K)
  #
  # est_unadj = matrix(nrow=nreps,ncol=n_bulk*K)
  # se_unadj_hc0 = matrix(nrow=nreps,ncol=n_bulk*K)
  # #se_unadj_cv = matrix(nrow=nreps,ncol=n_bulk*K)
  # se_unadj_hc3 = matrix(nrow=nreps,ncol=n_bulk*K)

  data_sparsity = c()

  est_adj = matrix(nrow=nreps,ncol=n_bulk*K)
  se_adj_hc0 = matrix(nrow=nreps,ncol=n_bulk*K)
  se_adj_hc3 = matrix(nrow=nreps,ncol=n_bulk*K)

  diff_adj = matrix(nrow=nreps,ncol=K)
  diff_adj_se_hc0 = matrix(nrow=nreps,ncol=K)
  diff_adj_p_hc0 = matrix(nrow = nreps,ncol=K)

  diff_adj_se_hc3 = matrix(nrow=nreps,ncol=K)
  diff_adj_p_hc3 = matrix(nrow = nreps,ncol=K)

  est_unadj = matrix(nrow=nreps,ncol=n_bulk*K)
  se_unadj_hc0 = matrix(nrow=nreps,ncol=n_bulk*K)
  diff_unadj = matrix(nrow=nreps,ncol=K)
  diff_unadj_se_hc0 = matrix(nrow=nreps,ncol=K)
  diff_unadj_p_hc0 = matrix(nrow = nreps,ncol=K)


  se_unadj_cv = matrix(nrow=nreps,ncol=n_bulk*K)
  diff_unadj_se_cv = matrix(nrow=nreps,ncol=K)
  diff_unadj_p_cv = matrix(nrow = nreps,ncol=K)

  se_unadj_hc3 = matrix(nrow=nreps,ncol=n_bulk*K)
  diff_unadj_se_hc3 = matrix(nrow=nreps,ncol=K)
  diff_unadj_p_hc3 = matrix(nrow = nreps,ncol=K)

  gene_names = rownames(ref)


  if(expr_dist=='log-normal'){
    mu = log(ref^2/sqrt(ref^2+sigma2))
    v = log(sigma2/ref^2+1)
  }

  if(expr_dist=='gamma'){
    shape = ref^2/sigma2
    rate=ref/sigma2
  }

  for(reps in 1:nreps){

    if(reps%%printevery==0){print(sprintf("running %d (out of %d)",reps,nreps))}


    # generate bulk data


    mb = 0
    for(k in 1:K){
      mb = mb + b[k,]*sim_MLN(n_bulk,ref[,k],sqrt(sigma2[,k]),R)
    }
    mb = t(mb)



    thetab = apply(mb,2,function(z){z/sum(z)})



    # mb = matrix(nrow=G,ncol=n_bulk)
    # for(i in 1:n_bulk){
    #
    #  X_b = matrix(nrow=G,ncol=K)
    #  for(k in 1:K){
    #    X_b[,k] = sim_MLN(1,ref[,k],sqrt(sigma2[,k]),R)
    #  }
    #
    #   mb[,i] = X_b%*%b[,i]
    #
    #
    # }
    # thetab = apply(mb,2,function(z){z/sum(z)})


    y = matrix(rpois(G*n_bulk,bulk_lib_size*G*thetab),nrow=G)
    rownames(y) = gene_names



    #bulks = SingleCellExperiment(assays = list(counts = y),
    #                             colData = DataFrame(individual = 1:n_bulk))



    ## generate single cell data

    Y = matrix(nrow=G,ncol=nk*n_indi*K)
    indi_idx = rep(1:n_indi,each = K*nk)
    cell_type = c()
    Cr = c()
    for(i in 1:n_indi){
      #print(i)
      indi_cell_idx = which(indi_idx==i)
      indi_celltype_idx = rep(1:K,each=nk)

      X_i = matrix(nrow=G,ncol=K)
      for(k in 1:K){
        X_i[,k] = sim_MLN(1,ref[,k],sqrt(sigma2[,k]),R)
      }

      Y_i = matrix(nrow=G,ncol=K*nk)
      Cr_i = c()
      for(k in 1:K){
        cc = which(indi_celltype_idx==k)
        ag = X_i[,k]%*%t(rep(1,nk))
        #Cr_i[cc] = rnbinom(nk,sc_lib_size*G,0.5)+1
        Cr_i[cc] = sc_lib_size*G


        if(expr_dist=='gamma'){
          atau2 = tau2[,k]%*%t(rep(1,nk))
          Y_i[,cc] = matrix(rgamma(G*nk,ag^2/atau2,rate=ag/atau2),ncol=nk)
        }

        if(expr_dist=='log-normal'){
          atau2 = tau2[,k]%*%t(rep(1,nk))
          Y_i[,cc] = exp(matrix(rnorm(G*nk,log(ag^2/sqrt(ag^2+atau2)),sqrt(log(atau2/ag^2+1))),ncol=nk))
        }


      }
      cell_type = c(cell_type,indi_celltype_idx)
      Y[,indi_cell_idx] = Y_i
      Cr = c(Cr,Cr_i)
    }

    Y = matrix(rpois(G*nk*n_indi*K,t(t(Y)*Cr/median(colSums(Y)))),ncol=nk*n_indi*K)

    rownames(Y) = gene_names


    # fit model

    data.obj = set_data_decon(y,Y,
                              ref_type = 'multi_sc',
                              indi_idx=indi_idx,
                              cell_type_idx=cell_type,
                              filter.gene=filter.gene)

    #browser()
    if(known_var){
      dmat = scRef_multi_proc(data.obj$Y,data.obj$cell_type_idx,data.obj$indi_idx,
                              estimator=x_estimator,est_sigma2 = est_pop_var,
                              meta_var=meta_var,
                              meta_mode=meta_mode,
                              sigma2 = sigma2,
                              tau2=tau2/G^2)
      #browser()
    }else{
      dmat = scRef_multi_proc(data.obj$Y,data.obj$cell_type_idx,data.obj$indi_idx,
                              estimator=x_estimator,est_sigma2 = est_pop_var,
                              meta_var=meta_var,
                              meta_mode=meta_mode)
    }


    fit.adj.hc0 = estimation_func2(y=data.obj$y,X=dmat$X,Vg=dmat$Vg,dmat$Sigma,
                                   w=1,hc.type='hc0',correction=FALSE,
                                   S=dmat$S,calc_cov=TRUE,verbose=verbose,
                                   ref_weights = FALSE,R=R)

    fit.adj.hc3 = estimation_func2(y=data.obj$y,X=dmat$X,Vg=dmat$Vg,dmat$Sigma,
                                   w=1,hc.type='hc3',correction=FALSE,
                                   S=dmat$S,calc_cov=TRUE,verbose=verbose,
                                   ref_weights = FALSE,R=R)

    fit.unadj.hc0 = estimation_func2(y=data.obj$y,X=dmat$X,Vg=dmat$Vg,dmat$Sigma,
                                 w=1,hc.type='hc0',correction=FALSE,
                                 S=dmat$S,calc_cov=TRUE,verbose=verbose,
                                 ref_weights = FALSE,R=NULL)


    data_sparsity[reps] = sum(Y==0)/prod(dim(Y))


    est_adj[reps,] = c(fit.adj.hc0$beta_hat)
    se_adj_hc0[reps,] = c(fit.adj.hc0$beta_se)
    se_adj_hc3[reps,] = c(fit.adj.hc3$beta_se)
    diff_out_hc0 = two_group_test(fit.adj.hc0,groups)
    diff_out_hc3 = two_group_test(fit.adj.hc3,groups)

    #browser()

    diff_adj[reps,] = c(diff_out_hc0$diff_group)
    diff_adj_se_hc0[reps,] = c(diff_out_hc0$diff_se)
    diff_adj_p_hc0[reps,] = c(diff_out_hc0$p_value)

    diff_adj_se_hc3[reps,] = c(diff_out_hc3$diff_se)
    diff_adj_p_hc3[reps,] = c(diff_out_hc3$p_value)


    #est_unadj[reps,] = c(fit_unadj$beta_hat)
    se_unadj_hc0[reps,] = c(fit.unadj.hc0$beta_se)
    #diff_unadj[reps,] = c(fit_unadj$diff_group)
    diff_out_hc0_unadj = two_group_test(fit.unadj.hc0,groups)
    diff_unadj_se_hc0[reps,] = c(diff_out_hc0_unadj$diff_se)
    diff_unadj_p_hc0[reps,] = c(diff_out_hc0_unadj$p_value)

    #se_unadj_cv[reps,] = c(fit_unadj$ols.out$beta_se)
    #diff_unadj_se_cv[reps,] = c(fit_unadj$ols.out$diff_se)
    #diff_unadj_p_cv[reps,] = fit_unadj$ols.out$p_value

    #se_unadj_hc3[reps,] = c(fit_unadj$sand.out.hc3$beta_se)
    #diff_unadj_se_hc3[reps,] = c(fit_unadj$sand.out.hc3$diff_se)
    #diff_unadj_p_hc3[reps,] = fit_unadj$sand.out.hc3$p_value



    #beta_hat_cov_adj[[reps]] = fit_adj$cov_beta_hat
    #beta_hat_cov_unadj[[reps]] = fit_unadj$sand.out$cov_beta_hat
    #beta_hat_cov_unadj_cv[[reps]] = fit_unadj$ols.out$cov_beta_hat
    #beta_hat_cov_unadj_hc3[[reps]] = fit_unadj$sand.out.hc3$cov_beta_hat

  }


  ########
  #1. mean squared error, and standard error
  mean_est_adj = apply(est_adj,2,mean)
  mse_adj = apply((est_adj - rep(1,nreps)%*%t(c(b)))^2,2,mean)
  se_est_adj = apply(est_adj,2,sd)
  #mean_est_unadj = apply(est_unadj,2,mean)
  #mse_unadj = apply((est_unadj - rep(1,nreps)%*%t(c(b)))^2,2,mean)
  #se_est_unadj = apply(est_unadj,2,sd)
  #mean_se_adj = apply(se_adj,2,mean)
  #mean_se_unadj = apply(se_unadj,2,mean)

  ###################
  #2. coverage

  ## adj

  ## hc0
  ci_l = est_adj - qnorm(1-alpha/2)*se_adj_hc0
  ci_r = est_adj + qnorm(1-alpha/2)*se_adj_hc0
  coverage_adj_hc0 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_adj_hc0=apply(coverage_adj_hc0,2,mean)

  ## hc3
  ci_l = est_adj - qnorm(1-alpha/2)*se_adj_hc3
  ci_r = est_adj + qnorm(1-alpha/2)*se_adj_hc3
  coverage_adj_hc3 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_adj_hc3=apply(coverage_adj_hc3,2,mean)


  # unadj
  ## sand
  ci_l = est_adj - qnorm(1-alpha/2)*se_unadj_hc0
  ci_r = est_adj + qnorm(1-alpha/2)*se_unadj_hc0
  coverage_unadj_hc0 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_unadj_hc0=apply(coverage_unadj_hc0,2,mean)


  ## cv
  # ci_l = est_unadj - qnorm(1-alpha/2)*se_unadj_cv
  # ci_r = est_unadj + qnorm(1-alpha/2)*se_unadj_cv
  # coverage_unadj_cv = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  # coverage_unadj_cv=apply(coverage_unadj_cv,2,mean)

  ## hc3
  # ci_l = est_unadj - qnorm(1-alpha/2)*se_unadj_hc3
  # ci_r = est_unadj + qnorm(1-alpha/2)*se_unadj_hc3
  # coverage_unadj_hc3 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  # coverage_unadj_hc3=apply(coverage_unadj_hc3,2,mean)

  ## diff

  ci_l = diff_adj - qnorm(1-alpha/2)*diff_adj_se_hc0
  ci_r = diff_adj + qnorm(1-alpha/2)*diff_adj_se_hc0
  coverage_diff_adj_hc0 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
  coverage_diff_adj_hc0=apply(coverage_diff_adj_hc0,2,mean)

  ci_l = diff_adj - qnorm(1-alpha/2)*diff_adj_se_hc3
  ci_r = diff_adj + qnorm(1-alpha/2)*diff_adj_se_hc3
  coverage_diff_adj_hc3 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
  coverage_diff_adj_hc3=apply(coverage_diff_adj_hc3,2,mean)

  ci_l = diff_adj - qnorm(1-alpha/2)*diff_unadj_se_hc0
  ci_r = diff_adj + qnorm(1-alpha/2)*diff_unadj_se_hc0
  coverage_diff_unadj_hc0 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
  coverage_diff_unadj_hc0=apply(coverage_diff_unadj_hc0,2,mean)

  # ci_l = diff_unadj - qnorm(1-alpha/2)*diff_unadj_se_cv
  # ci_r = diff_unadj + qnorm(1-alpha/2)*diff_unadj_se_cv
  # coverage_diff_unadj_cv = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
  # coverage_diff_unadj_cv = apply(coverage_diff_unadj_cv,2,mean)
  #
  # ci_l = diff_unadj - qnorm(1-alpha/2)*diff_unadj_se_hc3
  # ci_r = diff_unadj + qnorm(1-alpha/2)*diff_unadj_se_hc3
  # coverage_diff_unadj_hc3 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
  # coverage_diff_unadj_hc3 = apply(coverage_diff_unadj_hc3,2,mean)

  power_adj_hc0 = apply(diff_adj_p_hc0<=alpha,2,mean)
  power_adj_hc3 = apply(diff_adj_p_hc3<=alpha,2,mean)
  power_unadj_hc0 = apply(diff_unadj_p_hc0<=alpha,2,mean)
  #power_unadj_cv = apply(diff_unadj_p_cv<=alpha,2,mean)
  #power_unadj_hc3 = apply(diff_unadj_p_hc3<=alpha,2,mean)


  return(list(b=b,
              coverage_diff_adj_hc0=coverage_diff_adj_hc0,
              coverage_diff_adj_hc3=coverage_diff_adj_hc3,
              coverage_diff_unadj_hc0 = coverage_diff_unadj_hc0,
              #coverage_diff_unadj_cv = coverage_diff_unadj_cv,
              #coverage_diff_unadj_hc3 = coverage_diff_unadj_hc3,
              #beta_hat_cov_adj=beta_hat_cov_adj,
              #beta_hat_cov_unadj = beta_hat_cov_unadj,
              #beta_hat_cov_unadj_cv = beta_hat_cov_unadj_cv,
              #beta_hat_cov_unadj_hc3 = beta_hat_cov_unadj_hc3,
              power_adj_hc0=power_adj_hc0,
              power_adj_hc3=power_adj_hc3,
              power_unadj_hc0=power_unadj_hc0,
              #power_unadj_cv=power_unadj_cv,
              #power_unadj_hc3=power_unadj_hc3,

              diff_adj_p_hc0 = diff_adj_p_hc0,
              diff_adj_p_hc3 = diff_adj_p_hc3,

              diff_unadj_p_hc0=diff_unadj_p_hc0,
              #diff_unadj_p_cv=diff_unadj_p_cv,
              #diff_unadj_p_hc3=diff_unadj_p_hc3,

              diff_adj=diff_adj,
              #diff_unadj=diff_unadj,

              diff_adj_se_hc0=diff_adj_se_hc0,
              diff_adj_se_hc3=diff_adj_se_hc3,

              diff_unadj_se_hc0=diff_unadj_se_hc0,
              #diff_unadj_se_cv=diff_unadj_se_cv,
              #diff_unadj_se_hc3=diff_unadj_se_hc3,

              data_sparsity = data_sparsity,

              mean_est_adj=mean_est_adj,
              mse_adj=mse_adj,
              se_est_adj=se_est_adj,

              # mean_est_unadj=mean_est_unadj,
              # mse_unadj=mse_unadj,
              # se_est_unadj=se_est_unadj,

              coverage_adj_hc0=coverage_adj_hc0,
              coverage_adj_hc3=coverage_adj_hc3,
              coverage_unadj_hc0=coverage_unadj_hc0,
              #coverage_unadj_cv=coverage_unadj_cv,
              #coverage_unadj_hc3=coverage_unadj_hc3,

              se_adj_hc0 = se_adj_hc0,
              se_adj_hc3 = se_adj_hc3,
              se_unadj_hc0 = se_unadj_hc0,
              #se_unadj_cv = se_unadj_cv,
              #se_unadj_hc3 = se_unadj_hc3,
              est_adj = est_adj,
              #est_unadj = est_unadj,
              simu_param = list(ref=ref,b=b,
                                sc_lib_size = sc_lib_size,
                                bulk_lib_size = bulk_lib_size,nk=nk,
                                nreps=nreps,alpha=alpha)))



}


