


#'@param ref reference matrix
#'@param b K by #bulk matrix
#'

simu_studyX = function(ref,b,
                       bulk_lib_size = 500,
                       w.mode = 'equal',
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
                       known_var=FALSE){

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


  for(reps in 1:1){

    #if(reps%%printevery==0){print(sprintf("running %d (out of %d)",reps,nreps))}
    gene_names = rownames(ref)


    if(expr_dist=='log-normal'){
      mu = log(ref^2/sqrt(ref^2+sigma2))
      v = log(sigma2/ref^2+1)
    }

    if(expr_dist=='gamma'){
      shape = ref^2/sigma2
      rate=ref/sigma2
    }




    # generate bulk data

    mb = matrix(nrow=G,ncol=n_bulk)
    for(i in 1:n_bulk){

      if(expr_dist=='gamma'){
        X_b = matrix(rgamma(G*K,shape,rate),ncol=K)
      }

      if(expr_dist=='log-normal'){
        X_b = exp(matrix(rnorm(G*K,mu,sqrt(v)),ncol=K))
      }

      mb[,i] = X_b%*%b[,i]


    }
    thetab = apply(mb,2,function(z){z/sum(z)})


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

      if(expr_dist=='gamma'){
        X_i = matrix(rgamma(G*K,shape,rate),ncol=K)
      }
      if(expr_dist=='log-normal'){
        X_i = exp(matrix(rnorm(G*K,mu,sqrt(v)),ncol=K))
      }

      Y_i = matrix(nrow=G,ncol=K*nk)
      Cr_i = c()
      for(k in 1:K){
        cc = which(indi_celltype_idx==k)
        ag = X_i[,k]%*%t(rep(1,nk))
        Cr_i[cc] = rnbinom(nk,sc_lib_size*G,0.5)+1


        if(expr_dist=='gamma'){
          atau2 = tau2[,k]%*%t(rep(1,nk))
          Y_i[,cc] = matrix(rgamma(G*nk,ag^2/atau2,rate=ag/atau2),ncol=nk)
        }

        if(expr_dist=='log-normal'){
          atau2 = tau2[,k]%*%t(rep(1,nk))
          # print(k)
          # print(dim(log(ag^2/sqrt(ag^2+atau2))))
          # print(dim(log(atau2/ag^2+1)))
          # if(k==2){
          #   browser()
          # }
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




    fit.adj.hc0 = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v=TRUE,hc.type='hc0')
    fit.adj.hc3 = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v=TRUE,hc.type='hc3')
    fit.unadj.hc0 = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v  =FALSE,hc.type='hc0')
    fit.unadj.hc3 = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v=FALSE,hc.type='hc3')
    fit.unadj.const = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v=FALSE,hc.type='const')



  }

  return(list(se_adj_hc0 = fit.adj.hc0$beta_hat_se,
              se_adj_hc3 = fit.adj.hc3$beta_hat_se,
              se_unadj_hc0 = fit.unadj.hc0$beta_hat_se,
              se_unadj_hc3 = fit.unadj.hc3$beta_hat_se,
              se_unadj_const = fit.unadj.const$beta_hat_se,
              est_adj = fit.adj.hc0$beta_hat,
              est_unadj = fit.unadj.hc0$beta_hat))

}




#'@param ref reference matrix
#'@param b K by 1 vector
#'

simu_studyX2 = function(ref,b,
                        nreps = 100,
                       bulk_lib_size = 500,
                       w.mode = 'equal',
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
                       known_var=FALSE){

  G = nrow(ref)
  K = ncol(ref)
  #n_bulk = ncol(b)
  #b = apply(b,2,function(z){z/sum(z)})
  b = b/sum(b)


  est_adj = matrix(nrow=nreps,ncol=K)
  se_adj_hc0 = matrix(nrow=nreps,ncol=K)
  se_adj_hc3 = matrix(nrow=nreps,ncol=K)

  est_unadj = matrix(nrow=nreps,ncol=K)
  se_unadj_hc0 = matrix(nrow=nreps,ncol=K)
  se_unadj_const = matrix(nrow=nreps,ncol=K)
  se_unadj_hc3 = matrix(nrow=nreps,ncol=K)

  data_sparsity = c()


  if(expr_dist=='log-normal'){
    mu = log(ref^2/sqrt(ref^2+sigma2))
    v = log(sigma2/ref^2+1)
  }

  if(expr_dist=='gamma'){
    shape = ref^2/sigma2
    rate=ref/sigma2
  }



  n_bulk = 1
  for(reps in 1:nreps){

    #if(reps%%printevery==0){print(sprintf("running %d (out of %d)",reps,nreps))}
    gene_names = rownames(ref)

    # generate bulk data

    if(expr_dist=='gamma'){
      X_b = matrix(rgamma(G*K,shape,rate),ncol=K)
    }

    if(expr_dist=='log-normal'){
      X_b = exp(matrix(rnorm(G*K,mu,sqrt(v)),ncol=K))
    }

    mb = X_b%*%b

    thetab = apply(mb,2,function(z){z/sum(z)})

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

      if(expr_dist=='gamma'){
        X_i = matrix(rgamma(G*K,shape,rate),ncol=K)
      }
      if(expr_dist=='log-normal'){
        X_i = exp(matrix(rnorm(G*K,mu,sqrt(v)),ncol=K))
      }

      Y_i = matrix(nrow=G,ncol=K*nk)
      Cr_i = c()
      for(k in 1:K){
        cc = which(indi_celltype_idx==k)
        ag = X_i[,k]%*%t(rep(1,nk))
        Cr_i[cc] = rnbinom(nk,sc_lib_size*G,0.5)+1


        if(expr_dist=='gamma'){
          atau2 = tau2[,k]%*%t(rep(1,nk))
          Y_i[,cc] = matrix(rgamma(G*nk,ag^2/atau2,rate=ag/atau2),ncol=nk)
        }

        if(expr_dist=='log-normal'){
          atau2 = tau2[,k]%*%t(rep(1,nk))
          # print(k)
          # print(dim(log(ag^2/sqrt(ag^2+atau2))))
          # print(dim(log(atau2/ag^2+1)))
          # if(k==2){
          #   browser()
          # }
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




    fit.adj.hc0 = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v=TRUE,hc.type='hc0')
    fit.adj.hc3 = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v=TRUE,hc.type='hc3')
    fit.unadj.hc0 = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v  =FALSE,hc.type='hc0')
    fit.unadj.hc3 = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v=FALSE,hc.type='hc3')
    fit.unadj.const = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v=FALSE,hc.type='const')

    est_adj[reps,] = fit.adj.hc0$beta_hat
    est_unadj[reps,] = fit.unadj.hc0$beta_hat

    se_adj_hc0[reps,] =fit.adj.hc0$beta_hat_se
    se_adj_hc3[reps,] =fit.adj.hc3$beta_hat_se

    se_unadj_hc0[reps,] = fit.unadj.hc0$beta_hat_se
    se_unadj_const[reps,] = fit.unadj.const$beta_hat_se
    se_unadj_hc3[reps,] = fit.unadj.hc3$beta_hat_se



  }

  return(list(se_adj_hc0 = se_adj_hc0,
              se_adj_hc3 = se_adj_hc3,
              se_unadj_hc0 = se_unadj_hc0,
              se_unadj_hc3 = se_unadj_hc3,
              se_unadj_const = se_unadj_const,
              est_adj = est_adj,
              est_unadj = est_unadj))

}





calc_ci = function(b,V,alpha=0.05,eps=1e-7,trans=TRUE){
  if(trans){
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

#'@param bhat nbulk*K matrix
calc_coverage = function(bhat,b_se,alpha,b,trans = TRUE){
  nb = nrow(bhat)
  K = ncol(bhat)
  cover = matrix(nrow=nb,ncol=K)
  for(i in 1:nb){
    ci.temp = calc_ci(bhat[i,],diag(b_se[i,]^2),alpha=alpha,trans=trans)
    cover[i,] = (b[i,]>=ci.temp[,1])&(b[i,]<=ci.temp[,2])
  }
  mean(cover)
}

get_summary = function(simu.res,alpha=0.05,b,trans=TRUE){

  error = c(rmse(simu.res$est_adj,t(b)),rmse(simu.res$est_unadj,t(b)))
  error = rbind(error)
  colnames(error) = c("adj",'ols')
  rownames(error) = 'rmse'
  print(round(error,4))

  cov_adj = c(calc_coverage(simu.res$est_adj,simu.res$se_adj_hc0,alpha,t(b),trans=trans),
              calc_coverage(simu.res$est_adj,simu.res$se_adj_hc3,alpha,t(b),trans=trans),
              0)
  cov_unadj = c(calc_coverage(simu.res$est_unadj,simu.res$se_unadj_hc0,alpha,t(b),trans=trans),
                calc_coverage(simu.res$est_unadj,simu.res$se_unadj_hc3,alpha,t(b),trans=trans),
                calc_coverage(simu.res$est_unadj,simu.res$se_unadj_const,alpha,t(b),trans=trans))
  covs = cbind(cov_adj,cov_unadj)
  colnames(covs) = c("adj",'ols')
  rownames(covs) = c("hc0","hc3",'const.var')
  print(covs)
}


#'@description b is a vector now
#'
get_summary2 = function(simu.res,alpha=0.05,b,trans=TRUE){

  nreps = nrow(simu.res$est_adj)
  b.mat = rep(1,nreps)%*%t(b)
  error = c(rmse(simu.res$est_adj,b.mat),rmse(simu.res$est_unadj,b.mat))
  error = rbind(error)
  colnames(error) = c("adj",'ols')
  rownames(error) = 'rmse'
  print(round(error,4))

  cov_adj = c(calc_coverage(simu.res$est_adj,simu.res$se_adj_hc0,alpha,b.mat,trans=trans),
              calc_coverage(simu.res$est_adj,simu.res$se_adj_hc3,alpha,b.mat,trans=trans),
              0)
  cov_unadj = c(calc_coverage(simu.res$est_unadj,simu.res$se_unadj_hc0,alpha,b.mat,trans=trans),
                calc_coverage(simu.res$est_unadj,simu.res$se_unadj_hc3,alpha,b.mat,trans=trans),
                calc_coverage(simu.res$est_unadj,simu.res$se_unadj_const,alpha,b.mat,trans=trans))
  covs = cbind(cov_adj,cov_unadj)
  colnames(covs) = c("adj",'ols')
  rownames(covs) = c("hc0","hc3",'const.var')
  print(covs)
}





