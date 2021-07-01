
source('code/simulation/get_cor_pairs.R')

#'@description For simplicity, we only simulate data at individual level. - only generate X from U;

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
                            alpha.cor = 0.1,
                            groups = c(rep(1,ncol(b)/2),rep(2,ncol(b)/2)),
                            centeringXY = FALSE,
                            true.beta.for.Sigma=FALSE,
                            calc_cov = T,
                            est_cor = TRUE,
                            only.scale.pos.res = FALSE,
                            n_bulk_for_cor = 100,
                            cor_method = 'testing',
                            nfold = 10){

  is.identity = function(X){
    (sum(X)==nrow(X))
  }

  genp = function(K){
    p = runif(K)
    p/sum(p)
  }

  G = nrow(ref)
  K = ncol(ref)
  n_bulk = ncol(b)
  b = apply(b,2,function(z){z/sum(z)})
  is.indep = (is.identity(R))|(is.null(R))

  # est_adj = matrix(nrow=nreps,ncol=n_bulk*K)
  # se_adj = matrix(nrow=nreps,ncol=n_bulk*K)
  #
  # est_unadj = matrix(nrow=nreps,ncol=n_bulk*K)
  # se_unadj_hc0 = matrix(nrow=nreps,ncol=n_bulk*K)
  # #se_unadj_cv = matrix(nrow=nreps,ncol=n_bulk*K)
  # se_unadj_hc3 = matrix(nrow=nreps,ncol=n_bulk*K)


  ## 2. adjust for variance of U_hat, do not adjust for correlation, then hc0, hc2, hc3.
  ## 3. adjust for varaince of U_hat, adjust for correlation, then hc0, hc2, hc3.

  data_sparsity = c()

  p_hat = matrix(nrow=nreps,ncol=n_bulk*K)
  p_hat_se_adj_hc0 = matrix(nrow=nreps,ncol=n_bulk*K)
  p_hat_se_adj_hc2 = matrix(nrow=nreps,ncol=n_bulk*K)
  p_hat_se_adj_hc3 = matrix(nrow=nreps,ncol=n_bulk*K)
  p_hat_se_adj_jack = matrix(nrow=nreps,ncol=n_bulk*K)

  beta_hat = matrix(nrow=nreps,ncol=n_bulk*K)
  beta_se_adj_hc0 = matrix(nrow=nreps,ncol=n_bulk*K)
  beta_se_adj_hc2 = matrix(nrow=nreps,ncol=n_bulk*K)
  beta_se_adj_hc3 = matrix(nrow=nreps,ncol=n_bulk*K)
  beta_se_adj_jack = matrix(nrow=nreps,ncol=n_bulk*K)

  diff_hat = matrix(nrow=nreps,ncol=K)
  diff_adj_se_hc0 = matrix(nrow=nreps,ncol=K)
  diff_adj_se_hc2 = matrix(nrow=nreps,ncol=K)
  diff_adj_se_hc3 = matrix(nrow=nreps,ncol=K)
  diff_adj_se_jack = matrix(nrow=nreps,ncol=K)
  diff_adj_p_hc0 = matrix(nrow = nreps,ncol=K)
  diff_adj_p_hc2 = matrix(nrow = nreps,ncol=K)
  diff_adj_p_hc3 = matrix(nrow = nreps,ncol=K)
  diff_adj_p_jack = matrix(nrow = nreps,ncol=K)

  p_hat_se_unadj_hc0 = matrix(nrow=nreps,ncol=n_bulk*K)
  p_hat_se_unadj_hc2 = matrix(nrow=nreps,ncol=n_bulk*K)
  p_hat_se_unadj_hc3 = matrix(nrow=nreps,ncol=n_bulk*K)

  diff_unadj_se_hc0 = matrix(nrow=nreps,ncol=K)
  diff_unadj_se_hc2 = matrix(nrow=nreps,ncol=K)
  diff_unadj_se_hc3 = matrix(nrow=nreps,ncol=K)
  diff_unadj_p_hc0 = matrix(nrow = nreps,ncol=K)
  diff_unadj_p_hc2 = matrix(nrow = nreps,ncol=K)
  diff_unadj_p_hc3 = matrix(nrow = nreps,ncol=K)



  gene_names = rownames(ref)
  true_betas = matrix(nrow=nreps,ncol=n_bulk*K)


  ## pre calculate MLN and generate independent normal

  n.ref = matrix(nrow=G,ncol=K)
  n.Sigma.chol = list()
  n.Sigma = matrix(nrow=G,ncol=K)
  chol.R = chol(R)
  for(k in 1:K){
    n.ref[,k] = log(ref[,k]^2/sqrt(ref[,k]^2+sigma2[,k]))
    n.s = sqrt(log(1+sigma2[,k]/ref[,k]^2))
    if(!is.indep){
      n.Sigma.chol[[k]] = t(n.s*t(chol.R))
    }else{
      n.Sigma[,k] = n.s^2
    }
  }


  if(!est_cor){
    if(!is.indep){
      cor.idx = which(R!=0,arr.ind = T)
      cor.idx = cor.idx[(cor.idx[,1]!=cor.idx[,2]),]
    }else{
      cor.idx = NULL
    }
    n_bulk_for_cor = 0
  }

  cor.idx.all = list()




  for(reps in 1:nreps){

    if(reps%%printevery==0){print(sprintf("running %d (out of %d)",reps,nreps))}



    # generate individual reference matrices


    n.temp = n_indi+n_bulk+n_bulk_for_cor
    X_array = array(dim=c(G,K,n.temp))

    for(k in 1:K){
      if(is.indep){
        X_array[,k,] = exp(matrix(rnorm(G*n.temp,n.ref[,k],sqrt(n.Sigma[,k])),ncol=n.temp))
      }else{
        X_array[,k,] = t(exp(mvnfast::rmvn(n.temp,mu = n.ref[,k],sigma = n.Sigma.chol[[k]],isChol = TRUE)))
      }
    }

    #browser()

    X_array_bulk = X_array[,,1:n_bulk]

    X_array_ref = X_array[,,(n_bulk+1):(n_bulk+n_indi)]

    X_array_bulk_for_cor = X_array[,,-(1:(n_bulk+n_indi))]



    mb = lapply(1:n_bulk,function(i){X_array_bulk[,,i]%*%b[,i]})
    mb = do.call(cbind,mb)
    true.beta = t(t(b)*c(apply(mb,2,function(z){bulk_lib_size*G/sum(z)})))
    true_betas[reps,] = c(true.beta)
    thetab = apply(mb,2,function(z){z/sum(z)})


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


    ## get correlated pairs

    if(est_cor){

      # generate bulk data for estimating correlation

      mb = apply(X_array_bulk_for_cor,3,function(z){z%*%genp(K)})
      thetab = apply(mb,2,function(z){z/sum(z)})
      bulk_for_cor = matrix(rpois(G*n_bulk_for_cor,bulk_lib_size*G*thetab),nrow=G)
      rownames(bulk_for_cor) = gene_names

      cor.idx = get_cor_pairs2(bulk_for_cor,alpha=alpha.cor,method=cor_method)

      cor.idx.all[[reps]] = cor.idx

    }

    # fit model

    #browser()

    X = apply(X_array_ref,c(1,2),mean,na.rm=TRUE)
    V = t(apply(X_array_ref,c(1),function(z){(cov(t(z),use = 'complete.obs'))}))/n_indi


    fit.adj.hc0 = estimation_func2(y=y,X=X,Vg=V,
                                   w=1,hc.type='hc0',correction=FALSE,
                                   calc_cov=calc_cov,verbose=verbose,
                                   cor.idx=cor.idx,
                                   centeringXY=centeringXY,
                                   true.beta = if(true.beta.for.Sigma){true.beta}else{NULL},
                                   only.scale.pos.res=only.scale.pos.res)

    fit.adj.hc2 = estimation_func2(y=y,X=X,Vg=V,
                                   w=1,hc.type='hc2',correction=FALSE,
                                   calc_cov=calc_cov,verbose=verbose,
                                   cor.idx=cor.idx,
                                   centeringXY=centeringXY,
                                   true.beta = if(true.beta.for.Sigma){true.beta}else{NULL},
                                   only.scale.pos.res=only.scale.pos.res)

    fit.adj.hc3 = estimation_func2(y=y,X=X,Vg=V,
                                   w=1,hc.type='hc3',correction=FALSE,
                                   calc_cov=calc_cov,verbose=verbose,
                                   cor.idx=cor.idx,
                                   centeringXY=centeringXY,
                                   true.beta = if(true.beta.for.Sigma){true.beta}else{NULL},
                                   only.scale.pos.res=only.scale.pos.res)

    fit.adj.jack = estimation_func2(y=y,X=X,Vg=V,
                                   w=1,hc.type='jackknife',correction=FALSE,
                                   calc_cov=calc_cov,verbose=verbose,
                                   cor.idx=cor.idx,
                                   centeringXY=centeringXY,
                                   true.beta = if(true.beta.for.Sigma){true.beta}else{NULL},
                                   only.scale.pos.res=only.scale.pos.res,
                                   nfold = nfold)

    fit.unadj.hc0 = estimation_func2(y=y,X=X,Vg=V,
                                     w=1,hc.type='hc0',correction=FALSE,
                                     calc_cov=calc_cov,verbose=verbose,
                                     cor.idx=NULL,
                                     centeringXY=centeringXY,
                                     true.beta = if(true.beta.for.Sigma){true.beta}else{NULL},
                                     only.scale.pos.res=only.scale.pos.res)
    fit.unadj.hc2 = estimation_func2(y=y,X=X,Vg=V,
                                     w=1,hc.type='hc2',correction=FALSE,
                                     calc_cov=calc_cov,verbose=verbose,
                                     cor.idx=NULL,
                                     centeringXY=centeringXY,
                                     true.beta = if(true.beta.for.Sigma){true.beta}else{NULL},
                                     only.scale.pos.res=only.scale.pos.res)
    fit.unadj.hc3 = estimation_func2(y=y,X=X,Vg=V,
                                     w=1,hc.type='hc3',correction=FALSE,
                                     calc_cov=calc_cov,verbose=verbose,
                                     cor.idx=NULL,
                                     centeringXY=centeringXY,
                                     true.beta = if(true.beta.for.Sigma){true.beta}else{NULL},
                                     only.scale.pos.res=only.scale.pos.res)



    p_hat[reps,] = c(fit.adj.hc0$beta_hat)
    p_hat_se_adj_hc0[reps,] = c(fit.adj.hc0$beta_se)
    p_hat_se_adj_hc2[reps,] = c(fit.adj.hc2$beta_se)
    p_hat_se_adj_hc3[reps,] = c(fit.adj.hc3$beta_se)
    p_hat_se_adj_jack[reps,] = c(fit.adj.jack$beta_se)



    beta_hat[reps,] = c(fit.adj.hc0$beta_tilde_hat)
    beta_se_adj_hc0[reps,] = c(fit.adj.hc0$beta_tilde_se)
    beta_se_adj_hc2[reps,] = c(fit.adj.hc2$beta_tilde_se)
    beta_se_adj_hc3[reps,] = c(fit.adj.hc3$beta_tilde_se)
    beta_se_adj_jack[reps,] = c(fit.adj.jack$beta_tilde_se)


    #browser()

    diff_out_hc0 = two_group_test(fit.adj.hc0,groups)
    diff_hat[reps,] = c(diff_out_hc0$diff_group)
    diff_adj_se_hc0[reps,] = c(diff_out_hc0$diff_se)
    diff_adj_p_hc0[reps,] = c(diff_out_hc0$p_value)

    diff_out_hc2 = two_group_test(fit.adj.hc2,groups)
    diff_adj_se_hc2[reps,] = c(diff_out_hc2$diff_se)
    diff_adj_p_hc2[reps,] = c(diff_out_hc2$p_value)

    diff_out_hc3 = two_group_test(fit.adj.hc3,groups)
    diff_adj_se_hc3[reps,] = c(diff_out_hc3$diff_se)
    diff_adj_p_hc3[reps,] = c(diff_out_hc3$p_value)

    diff_out_jack = two_group_test(fit.adj.jack,groups)
    diff_adj_se_jack[reps,] = c(diff_out_jack$diff_se)
    diff_adj_p_jack[reps,] = c(diff_out_jack$p_value)


    #est_unadj[reps,] = c(fit_unadj$beta_hat)
    p_hat_se_unadj_hc0[reps,] = c(fit.unadj.hc0$beta_se)
    p_hat_se_unadj_hc2[reps,] = c(fit.unadj.hc2$beta_se)
    p_hat_se_unadj_hc3[reps,] = c(fit.unadj.hc3$beta_se)
    #diff_unadj[reps,] = c(fit_unadj$diff_group)

    diff_out_hc0_unadj = two_group_test(fit.unadj.hc0,groups)
    diff_unadj_se_hc0[reps,] = c(diff_out_hc0_unadj$diff_se)
    diff_unadj_p_hc0[reps,] = c(diff_out_hc0_unadj$p_value)

    diff_out_hc2_unadj = two_group_test(fit.unadj.hc2,groups)
    diff_unadj_se_hc2[reps,] = c(diff_out_hc2_unadj$diff_se)
    diff_unadj_p_hc2[reps,] = c(diff_out_hc2_unadj$p_value)

    diff_out_hc3_unadj = two_group_test(fit.unadj.hc3,groups)
    diff_unadj_se_hc3[reps,] = c(diff_out_hc3_unadj$diff_se)
    diff_unadj_p_hc3[reps,] = c(diff_out_hc3_unadj$p_value)

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
  mean_p_hat = apply(p_hat,2,mean)
  mse_adj = apply((p_hat - rep(1,nreps)%*%t(c(b)))^2,2,mean)
  se_p_hat = apply(p_hat,2,sd)


  ###################
  #2. coverage

  ## adj

  ## hc0
  ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_adj_hc0
  ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_adj_hc0
  coverage_adj_hc0 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_adj_hc0=apply(coverage_adj_hc0,2,mean,na.rm=TRUE)

  ## hc2
  ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_adj_hc2
  ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_adj_hc2
  coverage_adj_hc2 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_adj_hc2=apply(coverage_adj_hc2,2,mean,na.rm=TRUE)

  ## hc3
  ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_adj_hc3
  ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_adj_hc3
  coverage_adj_hc3 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_adj_hc3=apply(coverage_adj_hc3,2,mean,na.rm=TRUE)

  ## jack
  ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_adj_jack
  ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_adj_jack
  coverage_adj_jack = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_adj_jack=apply(coverage_adj_jack,2,mean,na.rm=TRUE)


  # unadj

  ## HC0
  ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_unadj_hc0
  ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_unadj_hc0
  coverage_unadj_hc0 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_unadj_hc0=apply(coverage_unadj_hc0,2,mean,na.rm=TRUE)

  ## HC2
  ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_unadj_hc2
  ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_unadj_hc2
  coverage_unadj_hc2 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_unadj_hc2=apply(coverage_unadj_hc2,2,mean,na.rm=TRUE)

  ## HC3
  ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_unadj_hc3
  ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_unadj_hc3
  coverage_unadj_hc3 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_unadj_hc3=apply(coverage_unadj_hc3,2,mean,na.rm=TRUE)


  # diff

  ## adjust

  ### hc0
  ci_l = diff_hat - qnorm(1-alpha/2)*diff_adj_se_hc0
  ci_r = diff_hat + qnorm(1-alpha/2)*diff_adj_se_hc0
  coverage_diff_adj_hc0 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
  coverage_diff_adj_hc0=apply(coverage_diff_adj_hc0,2,mean,na.rm=TRUE)

  ### hc2
  ci_l = diff_hat - qnorm(1-alpha/2)*diff_adj_se_hc2
  ci_r = diff_hat + qnorm(1-alpha/2)*diff_adj_se_hc2
  coverage_diff_adj_hc2 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc2$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc2$a)))<=ci_r)
  coverage_diff_adj_hc2=apply(coverage_diff_adj_hc2,2,mean,na.rm=TRUE)

  ### hc3
  ci_l = diff_hat - qnorm(1-alpha/2)*diff_adj_se_hc3
  ci_r = diff_hat + qnorm(1-alpha/2)*diff_adj_se_hc3
  coverage_diff_adj_hc3 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc3$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc3$a)))<=ci_r)
  coverage_diff_adj_hc3=apply(coverage_diff_adj_hc3,2,mean,na.rm=TRUE)

  ### jack
  ci_l = diff_hat - qnorm(1-alpha/2)*diff_adj_se_jack
  ci_r = diff_hat + qnorm(1-alpha/2)*diff_adj_se_jack
  coverage_diff_adj_jack = ((rep(1,nreps)%*%t(c(b%*%diff_out_jack$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_jack$a)))<=ci_r)
  coverage_diff_adj_jack=apply(coverage_diff_adj_jack,2,mean,na.rm=TRUE)



  ## unadjust

  ### hc0
  ci_l = diff_hat - qnorm(1-alpha/2)*diff_unadj_se_hc0
  ci_r = diff_hat + qnorm(1-alpha/2)*diff_unadj_se_hc0
  coverage_diff_unadj_hc0 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
  coverage_diff_unadj_hc0=apply(coverage_diff_unadj_hc0,2,mean,na.rm=TRUE)

  ### hc2
  ci_l = diff_hat - qnorm(1-alpha/2)*diff_unadj_se_hc2
  ci_r = diff_hat + qnorm(1-alpha/2)*diff_unadj_se_hc2
  coverage_diff_unadj_hc2 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc2$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc2$a)))<=ci_r)
  coverage_diff_unadj_hc2=apply(coverage_diff_unadj_hc2,2,mean,na.rm=TRUE)

  ### hc3
  ci_l = diff_hat - qnorm(1-alpha/2)*diff_unadj_se_hc3
  ci_r = diff_hat + qnorm(1-alpha/2)*diff_unadj_se_hc3
  coverage_diff_unadj_hc3 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc3$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc3$a)))<=ci_r)
  coverage_diff_unadj_hc3=apply(coverage_diff_unadj_hc3,2,mean,na.rm=TRUE)



  power_adj_hc0 = apply(diff_adj_p_hc0<=alpha,2,mean,na.rm=TRUE)
  power_adj_hc2 = apply(diff_adj_p_hc2<=alpha,2,mean,na.rm=TRUE)
  power_adj_hc3 = apply(diff_adj_p_hc3<=alpha,2,mean,na.rm=TRUE)
  power_adj_jack = apply(diff_adj_p_jack<=alpha,2,mean,na.rm=TRUE)
  power_unadj_hc0 = apply(diff_unadj_p_hc0<=alpha,2,mean,na.rm=TRUE)
  power_unadj_hc2 = apply(diff_unadj_p_hc2<=alpha,2,mean,na.rm=TRUE)
  power_unadj_hc3 = apply(diff_unadj_p_hc3<=alpha,2,mean,na.rm=TRUE)


  return(list(b=b,
              coverage_diff_adj_hc0=coverage_diff_adj_hc0,
              coverage_diff_adj_hc2=coverage_diff_adj_hc2,
              coverage_diff_adj_hc3=coverage_diff_adj_hc3,
              coverage_diff_adj_jack=coverage_diff_adj_jack,
              coverage_diff_unadj_hc0 = coverage_diff_unadj_hc0,
              coverage_diff_unadj_hc2 = coverage_diff_unadj_hc2,
              coverage_diff_unadj_hc3 = coverage_diff_unadj_hc3,

              power_adj_hc0=power_adj_hc0,
              power_adj_hc2=power_adj_hc2,
              power_adj_hc3=power_adj_hc3,
              power_adj_jack=power_adj_jack,
              power_unadj_hc0=power_unadj_hc0,
              power_unadj_hc2=power_unadj_hc2,
              power_unadj_hc3=power_unadj_hc3,

              diff_adj_p_hc0 = diff_adj_p_hc0,
              diff_adj_p_hc2 = diff_adj_p_hc2,
              diff_adj_p_hc3 = diff_adj_p_hc3,
              diff_adj_p_jack = diff_adj_p_jack,

              diff_unadj_p_hc0=diff_unadj_p_hc0,
              diff_unadj_p_hc2=diff_unadj_p_hc2,
              diff_unadj_p_hc3=diff_unadj_p_hc3,

              diff_hat=diff_hat,
              #diff_unadj=diff_unadj,

              diff_adj_se_hc0=diff_adj_se_hc0,
              diff_adj_se_hc2=diff_adj_se_hc2,
              diff_adj_se_hc3=diff_adj_se_hc3,
              diff_adj_se_jack=diff_adj_se_jack,

              diff_unadj_se_hc0=diff_unadj_se_hc0,
              diff_unadj_se_hc2=diff_unadj_se_hc2,
              diff_unadj_se_hc3=diff_unadj_se_hc3,



              mean_p_hat=mean_p_hat,
              mse_adj=mse_adj,
              se_p_hat=se_p_hat,

              # mean_est_unadj=mean_est_unadj,
              # mse_unadj=mse_unadj,
              # se_est_unadj=se_est_unadj,

              coverage_adj_hc0=coverage_adj_hc0,
              coverage_adj_hc2=coverage_adj_hc2,
              coverage_adj_hc3=coverage_adj_hc3,
              coverage_adj_jack=coverage_adj_jack,
              coverage_unadj_hc0=coverage_unadj_hc0,
              coverage_unadj_hc2=coverage_unadj_hc2,
              coverage_unadj_hc3=coverage_unadj_hc3,

              p_hat_se_adj_hc0 = p_hat_se_adj_hc0,
              p_hat_se_adj_hc2 = p_hat_se_adj_hc2,
              p_hat_se_adj_hc3 = p_hat_se_adj_hc3,
              p_hat_se_adj_jack = p_hat_se_adj_jack,

              p_hat_se_unadj_hc0 = p_hat_se_unadj_hc0,
              p_hat_se_unadj_hc2 = p_hat_se_unadj_hc2,
              p_hat_se_unadj_hc3 = p_hat_se_unadj_hc3,

              p_hat = p_hat,

              beta_hat = beta_hat,
              beta_se_adj_hc0=beta_se_adj_hc0,
              beta_se_adj_hc2=beta_se_adj_hc2,
              beta_se_adj_hc3=beta_se_adj_hc3,
              beta_se_adj_jack=beta_se_adj_jack,
              cor.idx.all=cor.idx.all,

              true_betas = true_betas,
              simu_param = list(b=b,
                                bulk_lib_size = bulk_lib_size,
                                nreps=nreps,alpha=alpha,nfold=nfold)))



}


