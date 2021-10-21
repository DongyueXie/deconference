
############## Run Simulation for manuscript ###############
source('code/deconference_main.R')
source('code/simulation/get_cor_pairs.R')
devtools::load_all('D://githubs/MuSiC')

#'@description Generate data: based on real Xin data, generate uncorrelated or correlated individual reference matrices.
#'Compare methods: ols, adjust for measurement error, adjust for both measurement error and correlation, adjust for both measurement error and correlation and add weights.
#'For simplicity, we only simulate data at individual level. - only generate X from U;

#'@param ref population gene expression matrix U, of dimension G by K
#'@param p1 first group mean
#'@param p2 second group mean
#'@param n_bulk number of bulk data
#'@param dirichlet whether generate random p
#'@param R correlation matrix, sparse form
#'@param sigma2 biological variance matrix across individuals, of dimension G by K
#'@param alpha significance level of confidence interval
#'@param alpha.cor fdr level of detecting correlations
#'@param groups indicates two groups of bulk individual
#'@param n_bulk_for_cor number of bulk data generated for inferring correlations
simu_study     =   function(ref,
                            sigma2,
                            p1,
                            p2,
                            n_bulk=100,
                            dirichlet = TRUE,
                            dirichlet.scale = 10,
                            R=NULL,
                            nreps = 100,
                            bulk_lib_size = 500,
                            n_ref = 10,
                            printevery=10,
                            verbose=FALSE,
                            alpha=0.05,
                            alpha.cor = 0.5,
                            groups = c(rep(1,n_bulk/2),rep(2,n_bulk/2)),
                            centeringXY = FALSE,
                            true.beta.for.Sigma=FALSE,
                            est_cor = TRUE,
                            n_bulk_for_cor = 100,
                            cor_method = 'testing',
                            add.music=TRUE,
                            make.cov.pos = FALSE){

  is.identity = function(X){
    if(!is.null(X)){
      (sum(X)==nrow(X))
    }else{
      FALSE
    }

  }

  genp = function(K){
    p = runif(K)
    p/sum(p)
  }

  G = nrow(ref)
  K = ncol(ref)
  is.indep = (is.identity(R))|(is.null(R))




  p_hat_ols = array(dim = c(K,n_bulk,nreps))
  p_hat_ols_se = array(dim = c(K,n_bulk,nreps))

  p_hat = array(dim = c(K,n_bulk,nreps))
  p_hat_se = array(dim = c(K,n_bulk,nreps))
  p_hat_se_cor = array(dim = c(K,n_bulk,nreps))
  #p_hat_se_cor_cv = array(dim = c(K,n_bulk,nreps))

  p_hat_weight = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cor = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cor_cv = array(dim = c(K,n_bulk,nreps))

  diff_hat_ols = matrix(nrow=nreps,ncol=K)
  diff_hat_ols_se = matrix(nrow=nreps,ncol=K)

  diff_hat = matrix(nrow=nreps,ncol=K)
  diff_hat_se = matrix(nrow=nreps,ncol=K)
  diff_hat_se_cor = matrix(nrow=nreps,ncol=K)
  #diff_hat_se_cor_cv = matrix(nrow=nreps,ncol=K)

  diff_hat_weight = matrix(nrow=nreps,ncol=K)
  diff_hat_weight_se = matrix(nrow=nreps,ncol=K)
  diff_hat_weight_se_cor = matrix(nrow=nreps,ncol=K)
  diff_hat_weight_se_cor_cv = matrix(nrow=nreps,ncol=K)


  gene_names = rownames(ref)
  true_betas = matrix(nrow=nreps,ncol=n_bulk*K)

  p_hat_music = array(dim = c(K,n_bulk,nreps))

  ## pre calculate MLN and generate independent normal

  norm.ref = matrix(nrow=G,ncol=K)
  norm.Sigma.chol = list()
  norm.Sigma = matrix(nrow=G,ncol=K)
  if(!is.indep){
    chol.R = chol(R)
  }
  for(k in 1:K){
    norm.ref[,k] = log(ref[,k]^2/sqrt(ref[,k]^2+sigma2[,k]))
    norm.s = sqrt(log(1+sigma2[,k]/ref[,k]^2))
    if(!is.indep){
      norm.Sigma.chol[[k]] = t(norm.s*t(chol.R))
    }else{
      norm.Sigma[,k] = norm.s^2
    }
  }


  if(!est_cor){
    if(!is.indep){
      cor.idx = which(R!=0,arr.ind = T)
      R01 = sparseMatrix(i = cor.idx[,1],j = cor.idx[,2],dims=c(G,G))
      cor.idx = cor.idx[(cor.idx[,1]!=cor.idx[,2]),]

    }else{
      cor.idx = NULL
      R01 = NULL
    }
  }else{

    X_array = array(dim=c(G,K,n_bulk_for_cor))
    for(k in 1:K){
      if(is.indep){
        X_array[,k,] = exp(matrix(rnorm(G*n_bulk_for_cor,norm.ref[,k],sqrt(norm.Sigma[,k])),ncol=n_bulk_for_cor))
      }else{
        X_array[,k,] = t(exp(mvnfast::rmvn(n_bulk_for_cor,mu = norm.ref[,k],sigma = norm.Sigma.chol[[k]],isChol = TRUE)))
      }
    }
    mb = apply(X_array,3,function(z){z%*%genp(K)})
    thetab = apply(mb,2,function(z){z/sum(z)})
    bulk_for_cor = matrix(rpois(G*n_bulk_for_cor,bulk_lib_size*G*thetab),nrow=G)
    bulk_for_cor = apply(bulk_for_cor,2,function(z){z/sum(z)*1e6})
    rownames(bulk_for_cor) = gene_names
    cor.idx = get_cor_pairs2(bulk_for_cor,alpha=alpha.cor,method=cor_method)
    R01 = sparseMatrix(i=cor.idx[,1],j=cor.idx[,2],dims=c(G,G))
    diag(R01) = 1

  }

  # all_fit = list()
  true_p = array(dim = c(K,n_bulk,nreps))

  for(reps in 1:nreps){

    if(reps%%printevery==0){print(sprintf("running %d (out of %d)",reps,nreps))}

    ## generate group p's
    nb1 = table(groups)[1]
    nb2 = n_bulk - nb1
    if(dirichlet){
      b1 = t(gtools::rdirichlet(nb1,p1*dirichlet.scale))
      b2 = t(gtools::rdirichlet(nb2,p2*dirichlet.scale))
      b = cbind(b1,b2)
    }else{
      b = cbind(matrix(p1,ncol=nb1,nrow=K),matrix(p2,ncol=nb2,nrow=K))
    }
    true_p[,,reps] = b

    ## generate individual reference matrices

    n.temp = n_ref+n_bulk
    X_array = array(dim=c(G,K,n.temp))

    for(k in 1:K){
      if(is.indep){
        X_array[,k,] = exp(matrix(rnorm(G*n.temp,norm.ref[,k],sqrt(norm.Sigma[,k])),ncol=n.temp))
      }else{
        X_array[,k,] = t(exp(mvnfast::rmvn(n.temp,mu = norm.ref[,k],sigma = norm.Sigma.chol[[k]],isChol = TRUE)))
      }
    }

    #browser()

    X_array_bulk = X_array[,,1:n_bulk]

    X_array_ref = X_array[,,(n_bulk+1):(n_bulk+n_ref)]




    mb = lapply(1:n_bulk,function(i){X_array_bulk[,,i]%*%b[,i]})
    mb = do.call(cbind,mb)
    true.beta = t(t(b)*c(apply(mb,2,function(z){bulk_lib_size*G/sum(z)})))
    true_betas[reps,] = c(true.beta)
    thetab = apply(mb,2,function(z){z/sum(z)})


    #browser()
    y = matrix(rpois(G*n_bulk,bulk_lib_size*G*thetab),nrow=G)
    rownames(y) = gene_names


    # fit model

    #browser()

    X = apply(X_array_ref,c(1,2),mean,na.rm=TRUE)
    V = t(apply(X_array_ref,c(1),function(z){(cov(t(z),use = 'complete.obs'))}))/n_ref



    # ols fit
    fit.ols = ols_hc3(y,X,groups=groups)


    # adjust for measurement error, do not adjust for correlation
    fit.err = estimation_func2(y=y,
                               X=X,
                               Vg=V,
                               w=1,
                               hc.type='hc3',
                               correction=FALSE,
                               verbose=verbose,
                               R01=NULL,
                               true.beta = if(true.beta.for.Sigma){true.beta}else{NULL},
                               groups = groups,
                               Q.pos = make.cov.pos,
                               Sigma.pos = make.cov.pos,
                               V_tilde.pos = make.cov.pos)

    fit.err.cor = estimation_func2(y=y,
                               X=X,
                               Vg=V,
                               w=1,
                               hc.type='hc3',
                               correction=FALSE,
                               verbose=verbose,
                               R01=R01,
                               true.beta = if(true.beta.for.Sigma){true.beta}else{NULL},
                               groups = groups,
                               Q.pos = make.cov.pos,
                               Sigma.pos = make.cov.pos,
                               V_tilde.pos = make.cov.pos)



    # fit.err.cor.cv = estimation_func2(y=y,
    #                                X=X,
    #                                Vg=V,
    #                                w=1,
    #                                hc.type='jackknife_indep',
    #                                correction=FALSE,
    #                                verbose=verbose,
    #                                R01=R01,
    #                                true.beta = if(true.beta.for.Sigma){true.beta}else{NULL})

    fit.vash = vashr::vash(sqrt(rowSums(V)),df=n_ref-1)
    w = 1/(fit.vash$sd.post)^2

    fit.err.weight = estimation_func2(y=y,
                                      X=X,
                                      Vg=V,
                                      w=w,
                                      hc.type='hc3',
                                      correction=FALSE,
                                      verbose=verbose,
                                      R01=NULL,
                                      true.beta = if(true.beta.for.Sigma){true.beta}else{NULL},
                                      groups = groups,
                                      Q.pos = make.cov.pos,
                                      Sigma.pos = make.cov.pos,
                                      V_tilde.pos = make.cov.pos)

    fit.err.cor.weight = estimation_func2(y=y,
                                   X=X,
                                   Vg=V,
                                   w=w,
                                   hc.type='hc3',
                                   correction=FALSE,
                                   verbose=verbose,
                                   R01=R01,
                                   true.beta = if(true.beta.for.Sigma){true.beta}else{NULL},
                                   groups = groups,
                                   Q.pos = make.cov.pos,
                                   Sigma.pos = make.cov.pos,
                                   V_tilde.pos = make.cov.pos)

    fit.err.cor.weight.cv = estimation_func2(y=y,
                                          X=X,
                                          Vg=V,
                                          w=w,
                                          hc.type='jackknife_indep',
                                          correction=FALSE,
                                          verbose=verbose,
                                          R01=R01,
                                          true.beta = if(true.beta.for.Sigma){true.beta}else{NULL},
                                          groups = groups,
                                          Q.pos = make.cov.pos,
                                          Sigma.pos = make.cov.pos,
                                          V_tilde.pos = make.cov.pos)

    if(add.music){
      Sigma.music = t(apply(X_array_ref,c(1),function(z){diag(cov(t(z),use = 'complete.obs'))}))
      Est.prop.allgene = NULL
      for(bb in 1:n_bulk){
        fit.music = music.basic(y[,bb],X,S=1,Sigma.music,iter.max=1000,nu=1e-4,eps=0.01)
        Est.prop.allgene = cbind(Est.prop.allgene, fit.music$p.weight)
      }
      p_hat_music[,,reps] = Est.prop.allgene
    }
#
#     temp = list()
#     temp[[1]] = fit.ols
#     temp[[2]] = fit.err
#     temp[[3]] = fit.err.cor
#     temp[[4]] = fit.err.cor.weight
#     all_fit[[reps]] = temp

    # temp = list()
    # temp$fit.ols =fit.ols
    # temp$fit.err =fit.err
    # temp$fit.err.cor = fit.err.cor
    # #temp$fit.err.cor.cv = fit.err.cor.cv
    # temp$fit.err.cor.weight = fit.err.cor.weight
    # temp$fit.err.cor.weight.cv = fit.err.cor.weight.cv
    # all_fit[[reps]] = temp

    p_hat_ols[,,reps] = fit.ols$p_hat
    p_hat_ols_se[,,reps] = fit.ols$p_hat_se

    p_hat[,,reps] = fit.err$p_hat
    p_hat_se[,,reps] = fit.err$p_hat_se
    p_hat_se_cor[,,reps] = fit.err.cor$p_hat_se
    #p_hat_se_cor_cv[,,reps] = fit.err.cor.cv$p_hat_se

    p_hat_weight[,,reps] = fit.err.cor.weight$p_hat
    p_hat_weight_se[,,reps] = fit.err.weight$p_hat_se
    p_hat_weight_se_cor[,,reps] = fit.err.cor.weight$p_hat_se
    p_hat_weight_se_cor_cv[,,reps] = fit.err.cor.weight.cv$p_hat_se

    diff_hat_ols[reps,] = fit.ols$diff_hat
    diff_hat_ols_se[reps,] = fit.ols$diff_hat_se

    #diff_test = two_group_test(fit.err,groups)
    diff_hat[reps,] = fit.err$two_group_res$diff_hat
    diff_hat_se[reps,] = fit.err$two_group_res$diff_hat_se

    #diff_test_cor = two_group_test(fit.err.cor,groups)
    diff_hat_se_cor[reps,] = fit.err.cor$two_group_res$diff_hat_se

    #diff_test_cor_cv = two_group_test(fit.err.cor.cv,groups)
    #diff_hat_se_cor_cv[reps,] = c(diff_test_cor_cv$diff_se)

    #diff_test_cor_weight = two_group_test(fit.err.cor.weight,groups)
    diff_hat_weight[reps,] = fit.err.cor.weight$two_group_res$diff_hat
    diff_hat_weight_se_cor[reps,] = fit.err.cor.weight$two_group_res$diff_hat_se


    diff_hat_weight_se[reps,] = fit.err.weight$two_group_res$diff_hat_se

    #diff_test_cor_weight_cv = two_group_test(fit.err.cor.weight.cv,groups)
    diff_hat_weight_se_cor_cv[reps,] = fit.err.cor.weight.cv$two_group_res$diff_hat_se

  }

  return(list(p_hat_ols = p_hat_ols,
              p_hat_ols_se = p_hat_ols_se,

              p_hat = p_hat,
              p_hat_se = p_hat_se,
              p_hat_se_cor = p_hat_se_cor,
              #p_hat_se_cor_cv = p_hat_se_cor_cv,

              p_hat_weight = p_hat_weight,
              p_hat_weight_se = p_hat_weight_se,
              p_hat_weight_se_cor = p_hat_weight_se_cor,
              p_hat_weight_se_cor_cv = p_hat_weight_se_cor_cv,

              diff_hat_ols = diff_hat_ols,
              diff_hat_ols_se = diff_hat_ols_se,

              diff_hat = diff_hat,
              diff_hat_se = diff_hat_se,
              diff_hat_se_cor = diff_hat_se_cor,
              #diff_hat_se_cor_cv = diff_hat_se_cor_cv,

              diff_hat_weight = diff_hat_weight,
              diff_hat_weight_se = diff_hat_weight_se,
              diff_hat_weight_se_cor = diff_hat_weight_se_cor,
              diff_hat_weight_se_cor_cv = diff_hat_weight_se_cor_cv,

              p_hat_music = p_hat_music,
              true_p = true_p,

              input = list(p1=p1,p2=p2,groups=groups)))





}

calc_ci = function(b,V,alpha=0.05,eps=1e-7,logit=TRUE){
  if(logit){
    b = b+eps
    b.temp = log(b/(1-b))
    J = diag(1/(b*(1-b)))
    V.temp = J%*%V%*%J
    ci.l = b.temp-qnorm(1-alpha/2)*sqrt(diag(V.temp))
    ci.r = b.temp+qnorm(1-alpha/2)*sqrt(diag(V.temp))
    return(cbind(1/(1+exp(-ci.l)),1/(1+exp(-ci.r))))
  }else{
    ci.l = b-qnorm(1-alpha/2)*sqrt(diag(V))
    ci.r = b+qnorm(1-alpha/2)*sqrt(diag(V))
    return(cbind(ci.l,ci.r))
  }

}

#
#
# ########
# #1. mean squared error, and standard error
# mean_p_hat = apply(p_hat,2,mean)
# mse_adj = apply((p_hat - rep(1,nreps)%*%t(c(b)))^2,2,mean)
# se_p_hat = apply(p_hat,2,sd)
#
#
# ###################
# #2. coverage
#
# ## adj
#
# ## hc0
# ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_adj_hc0
# ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_adj_hc0
# coverage_adj_hc0 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
# coverage_adj_hc0=apply(coverage_adj_hc0,2,mean,na.rm=TRUE)
#
# ## hc2
# ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_adj_hc2
# ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_adj_hc2
# coverage_adj_hc2 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
# coverage_adj_hc2=apply(coverage_adj_hc2,2,mean,na.rm=TRUE)
#
# ## hc3
# ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_adj_hc3
# ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_adj_hc3
# coverage_adj_hc3 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
# coverage_adj_hc3=apply(coverage_adj_hc3,2,mean,na.rm=TRUE)
#
# ## jack
# ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_adj_jack
# ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_adj_jack
# coverage_adj_jack = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
# coverage_adj_jack=apply(coverage_adj_jack,2,mean,na.rm=TRUE)
#
#
# # unadj
#
# ## HC0
# ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_unadj_hc0
# ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_unadj_hc0
# coverage_unadj_hc0 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
# coverage_unadj_hc0=apply(coverage_unadj_hc0,2,mean,na.rm=TRUE)
#
# ## HC2
# ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_unadj_hc2
# ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_unadj_hc2
# coverage_unadj_hc2 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
# coverage_unadj_hc2=apply(coverage_unadj_hc2,2,mean,na.rm=TRUE)
#
# ## HC3
# ci_l = p_hat - qnorm(1-alpha/2)*p_hat_se_unadj_hc3
# ci_r = p_hat + qnorm(1-alpha/2)*p_hat_se_unadj_hc3
# coverage_unadj_hc3 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
# coverage_unadj_hc3=apply(coverage_unadj_hc3,2,mean,na.rm=TRUE)
#
#
# # diff
#
# ## adjust
#
# ### hc0
# ci_l = diff_hat - qnorm(1-alpha/2)*diff_adj_se_hc0
# ci_r = diff_hat + qnorm(1-alpha/2)*diff_adj_se_hc0
# coverage_diff_adj_hc0 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
# coverage_diff_adj_hc0=apply(coverage_diff_adj_hc0,2,mean,na.rm=TRUE)
#
# ### hc2
# ci_l = diff_hat - qnorm(1-alpha/2)*diff_adj_se_hc2
# ci_r = diff_hat + qnorm(1-alpha/2)*diff_adj_se_hc2
# coverage_diff_adj_hc2 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc2$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc2$a)))<=ci_r)
# coverage_diff_adj_hc2=apply(coverage_diff_adj_hc2,2,mean,na.rm=TRUE)
#
# ### hc3
# ci_l = diff_hat - qnorm(1-alpha/2)*diff_adj_se_hc3
# ci_r = diff_hat + qnorm(1-alpha/2)*diff_adj_se_hc3
# coverage_diff_adj_hc3 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc3$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc3$a)))<=ci_r)
# coverage_diff_adj_hc3=apply(coverage_diff_adj_hc3,2,mean,na.rm=TRUE)
#
# ### jack
# ci_l = diff_hat - qnorm(1-alpha/2)*diff_adj_se_jack
# ci_r = diff_hat + qnorm(1-alpha/2)*diff_adj_se_jack
# coverage_diff_adj_jack = ((rep(1,nreps)%*%t(c(b%*%diff_out_jack$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_jack$a)))<=ci_r)
# coverage_diff_adj_jack=apply(coverage_diff_adj_jack,2,mean,na.rm=TRUE)
#
#
#
# ## unadjust
#
# ### hc0
# ci_l = diff_hat - qnorm(1-alpha/2)*diff_unadj_se_hc0
# ci_r = diff_hat + qnorm(1-alpha/2)*diff_unadj_se_hc0
# coverage_diff_unadj_hc0 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc0$a)))<=ci_r)
# coverage_diff_unadj_hc0=apply(coverage_diff_unadj_hc0,2,mean,na.rm=TRUE)
#
# ### hc2
# ci_l = diff_hat - qnorm(1-alpha/2)*diff_unadj_se_hc2
# ci_r = diff_hat + qnorm(1-alpha/2)*diff_unadj_se_hc2
# coverage_diff_unadj_hc2 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc2$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc2$a)))<=ci_r)
# coverage_diff_unadj_hc2=apply(coverage_diff_unadj_hc2,2,mean,na.rm=TRUE)
#
# ### hc3
# ci_l = diff_hat - qnorm(1-alpha/2)*diff_unadj_se_hc3
# ci_r = diff_hat + qnorm(1-alpha/2)*diff_unadj_se_hc3
# coverage_diff_unadj_hc3 = ((rep(1,nreps)%*%t(c(b%*%diff_out_hc3$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out_hc3$a)))<=ci_r)
# coverage_diff_unadj_hc3=apply(coverage_diff_unadj_hc3,2,mean,na.rm=TRUE)
#
#
#
# power_adj_hc0 = apply(diff_adj_p_hc0<=alpha,2,mean,na.rm=TRUE)
# power_adj_hc2 = apply(diff_adj_p_hc2<=alpha,2,mean,na.rm=TRUE)
# power_adj_hc3 = apply(diff_adj_p_hc3<=alpha,2,mean,na.rm=TRUE)
# power_adj_jack = apply(diff_adj_p_jack<=alpha,2,mean,na.rm=TRUE)
# power_unadj_hc0 = apply(diff_unadj_p_hc0<=alpha,2,mean,na.rm=TRUE)
# power_unadj_hc2 = apply(diff_unadj_p_hc2<=alpha,2,mean,na.rm=TRUE)
# power_unadj_hc3 = apply(diff_unadj_p_hc3<=alpha,2,mean,na.rm=TRUE)
#
#
# return(list(b=b,
#             coverage_diff_adj_hc0=coverage_diff_adj_hc0,
#             coverage_diff_adj_hc2=coverage_diff_adj_hc2,
#             coverage_diff_adj_hc3=coverage_diff_adj_hc3,
#             coverage_diff_adj_jack=coverage_diff_adj_jack,
#             coverage_diff_unadj_hc0 = coverage_diff_unadj_hc0,
#             coverage_diff_unadj_hc2 = coverage_diff_unadj_hc2,
#             coverage_diff_unadj_hc3 = coverage_diff_unadj_hc3,
#
#             power_adj_hc0=power_adj_hc0,
#             power_adj_hc2=power_adj_hc2,
#             power_adj_hc3=power_adj_hc3,
#             power_adj_jack=power_adj_jack,
#             power_unadj_hc0=power_unadj_hc0,
#             power_unadj_hc2=power_unadj_hc2,
#             power_unadj_hc3=power_unadj_hc3,
#
#             diff_adj_p_hc0 = diff_adj_p_hc0,
#             diff_adj_p_hc2 = diff_adj_p_hc2,
#             diff_adj_p_hc3 = diff_adj_p_hc3,
#             diff_adj_p_jack = diff_adj_p_jack,
#
#             diff_unadj_p_hc0=diff_unadj_p_hc0,
#             diff_unadj_p_hc2=diff_unadj_p_hc2,
#             diff_unadj_p_hc3=diff_unadj_p_hc3,
#
#             diff_hat=diff_hat,
#             #diff_unadj=diff_unadj,
#
#             diff_adj_se_hc0=diff_adj_se_hc0,
#             diff_adj_se_hc2=diff_adj_se_hc2,
#             diff_adj_se_hc3=diff_adj_se_hc3,
#             diff_adj_se_jack=diff_adj_se_jack,
#
#             diff_unadj_se_hc0=diff_unadj_se_hc0,
#             diff_unadj_se_hc2=diff_unadj_se_hc2,
#             diff_unadj_se_hc3=diff_unadj_se_hc3,
#
#
#
#             mean_p_hat=mean_p_hat,
#             mse_adj=mse_adj,
#             se_p_hat=se_p_hat,
#
#             # mean_est_unadj=mean_est_unadj,
#             # mse_unadj=mse_unadj,
#             # se_est_unadj=se_est_unadj,
#
#             coverage_adj_hc0=coverage_adj_hc0,
#             coverage_adj_hc2=coverage_adj_hc2,
#             coverage_adj_hc3=coverage_adj_hc3,
#             coverage_adj_jack=coverage_adj_jack,
#             coverage_unadj_hc0=coverage_unadj_hc0,
#             coverage_unadj_hc2=coverage_unadj_hc2,
#             coverage_unadj_hc3=coverage_unadj_hc3,
#
#             p_hat_se_adj_hc0 = p_hat_se_adj_hc0,
#             p_hat_se_adj_hc2 = p_hat_se_adj_hc2,
#             p_hat_se_adj_hc3 = p_hat_se_adj_hc3,
#             p_hat_se_adj_jack = p_hat_se_adj_jack,
#
#             p_hat_se_unadj_hc0 = p_hat_se_unadj_hc0,
#             p_hat_se_unadj_hc2 = p_hat_se_unadj_hc2,
#             p_hat_se_unadj_hc3 = p_hat_se_unadj_hc3,
#
#             p_hat = p_hat,
#
#             beta_hat = beta_hat,
#             beta_se_adj_hc0=beta_se_adj_hc0,
#             beta_se_adj_hc2=beta_se_adj_hc2,
#             beta_se_adj_hc3=beta_se_adj_hc3,
#             beta_se_adj_jack=beta_se_adj_jack,
#             cor.idx.all=cor.idx.all,
#
#             true_betas = true_betas,
#             simu_param = list(b=b,
#                               bulk_lib_size = bulk_lib_size,
#                               nreps=nreps,alpha=alpha,nfold=nfold)))
#





ols_hc3 = function(y,X,w=NULL,groups=NULL){
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
  res = yw - Xw%*%beta_hat

  J = matrix(0,nrow=(nb*K),ncol=(nb*K))
  for(i in 1:nb){
    J[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = J_sum2one(beta_hat[,i],K)
  }




  Sigma = matrix(0,nrow=nb*K,ncol=nb*K)
  Q_inv = matrix(0,nrow=nb*K,ncol=nb*K)
  Sigma_ii = matrix(0,nrow=nb*K,ncol=nb*K)

  for(i in 1:nb){
    Q_inv[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = A_inv

    #Hi = t(t(X%*%A_inv%*%t(X))*w)
    #hi = diag(Hi)


    h = rowSums((X%*%A_inv)*X)*w
    ri = res[,i]

    for(j in i:nb){

      if(j==i){

        Sigma_ij = crossprod(c(ri)/(1-pmax(pmin(h,1-1/G),0))*w*X)

        Sigma_ii[((i-1)*K+1):(i*K),((i-1)*K+1):(i*K)] = Sigma_ij
      }else{

        rj = res[,j]

        Sigma_ij = crossprod(c(ri)/(1-pmax(pmin(h,1-1/G),0))*w*X,
                             c(rj)/(1-pmax(pmin(h,1-1/G),0))*w*X)
      }


      Sigma[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)] = Sigma_ij

    }
  }
  Sigma = Sigma+t(Sigma)-Sigma_ii
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


  }


  return(list(p_hat=p_hat,
              p_hat_se = p_hat_se,
              diff_hat=diff_hat,
              diff_hat_se = sqrt(diag(V_tilde)),
              diff_z_score=z_score,
              diff_p_value=p_value))




}


