########################################
############ simulation function for 1. bulk reference 2. same individual single cell referece 3. multi indi sc
########################################


#'@title Simulation study function of 1. bulk as reference and 2. one individual single cell reference. This simulation averages over parameter b.
#'@param ref ref data matrix, Gene by cell type, relative expression
#'@param Ng number of genes to use in simulation, gene will be resampled each round
#'@param b cell type proportions
#'@param ref_type 'sc' for single cell ref; 'bulk' for bulk ref, multi_sc
#'@param hc.type type of sandwich estimator
#'@param tau2 a G by K matrix, entry(g,k) is the variance of gene g in cell type k, across cells
#'@param tau2known whether \tau^2, the variance of gene expr across cell, is known
#'@param sigma2 a G by K matrix, entry(g,k) is the variance of gene g in cell type k,across individuals
#'@param weight default, optimal, equal, external
#'@param marker_gene marker gene index
#'@param bulk_lib_size bulk sample library size
#'@param sc_lib_size single cell ref library size
#'@param ref_lib_size bulk ref data library size
#'@param nk number of single cells of each cell type
#'@param x_estimator estimate the reference matrix from single cell data, separate or aggregate
#'@param nreps number of repetition.
#'@param alpha significance level
#'@param a fuller's correction parameter a
#'@param correction whether perform fullers correction.
#'@param s cell size
#'@param add_y2 whether add another bulk sample and perform hypothesis testing.
#'@param b2 the second bulk sample's cell type proportion
#'@return a list with estimates and standard errors.
#'@return 1. true proportion; 2. mean of estimates; 3. standard error of estimates; 4. coverage.
#'@return input parameters


simu_study = function(ref,Ng=nrow(ref),b,ref_type='sc',
                      hc.type = 'hc3',
                      marker_gene = NULL,
                      weight = 'default',
                      bulk_lib_size = 50,
                      sc_lib_size = 0.1,
                      ref_lib_size = 50,
                      nk = 100,
                      ref_scale = Ng,
                      same_indi = FALSE,
                      mean_to_var_sigma = 10,
                      mean_to_var_tau = 10,
                      tau2 = matrix(runif(nrow(ref)*length(b),0.5/(nrow(ref))^2,1/(nrow(ref))^2),nrow=nrow(ref)),
                      tau2known=FALSE,
                      x_estimator = 'aggregate',
                      nreps=100,alpha=0.05,a=0,
                      correction=TRUE,s,printevery=10,
                      add_y2 = TRUE, b2=NULL,gene_thresh=10,
                      n_indi = 6,sigma2known=FALSE,
                      sigma2 = matrix(runif(nrow(ref)*length(b),0.5/(nrow(ref))^3,1/(nrow(ref))^3),nrow=nrow(ref)),
                      est_pop_var = FALSE,meta_var='adjust',meta_mode = 'one',eps=0){

  G = nrow(ref)
  K = ncol(ref)

  ref = apply(ref,2,function(z){z/sum(z)})


  if(missing(s)){
    s = rep(1,K)
  }

  b = b/sum(b)
  p = (b*s)/sum(b*s)

  if(length(tau2)==1){
    tau2 = matrix(tau2,nrow=Ng,ncol=K)
  }

  if(length(sigma2)==1){
    sigma2 = matrix(sigma2,nrow=Ng,ncol=K)
  }

  if(add_y2 & !is.null(b2)){
    b2 = b2/sum(b2)
    p2 = (b2*s)/sum(b2*s)
  }else{
    p2 = NULL
  }


  if(add_y2){
    est_adj = matrix(nrow=nreps,ncol=2*K)
    se_adj = matrix(nrow=nreps,ncol=2*K)
    diff_adj = matrix(nrow=nreps,ncol=K)
    diff_adj_se = matrix(nrow=nreps,ncol=K)
    diff_adj_p = matrix(nrow = nreps,ncol=K)

    est_unadj = matrix(nrow=nreps,ncol=2*K)
    se_unadj = matrix(nrow=nreps,ncol=2*K)
    diff_unadj = matrix(nrow=nreps,ncol=K)
    diff_unadj_se = matrix(nrow=nreps,ncol=K)
    diff_unadj_p = matrix(nrow = nreps,ncol=K)


    se_unadj_cv = matrix(nrow=nreps,ncol=2*K)

    diff_unadj_se_cv = matrix(nrow=nreps,ncol=K)
    diff_unadj_p_cv = matrix(nrow = nreps,ncol=K)
  }else{
    est_adj = matrix(nrow=nreps,ncol=K)
    se_adj = matrix(nrow=nreps,ncol=K)

    est_unadj = matrix(nrow=nreps,ncol=K)
    se_unadj = matrix(nrow=nreps,ncol=K)


    se_unadj_cv = matrix(nrow=nreps,ncol=K)
  }




  data_sparsity = c()
  beta_hat_cov_adj = list()
  beta_hat_cov_unadj = list()
  beta_hat_cov_unadj_cv = list()

  sigma2_hat = list()# unly use when ref is multi-sc and est_sigma2 IS TRURE

  for(reps in 1:nreps){


    if(reps%%printevery==0){print(sprintf("done %d (out of %d)",reps,nreps))}

    #Obtain X
    if(Ng==G){
      ref_rep = ref
    }else{
      ref_rep = ref[sample(1:G,Ng),]
    }

    gene_names = rownames(ref_rep)

    Theta = apply(ref_rep,2,function(z){z/sum(z)})*ref_scale

    if(same_indi){

      # bulk data gene relative expression.

      mb = Theta%*%diag(s*(Ng/G))%*%b
      thetab = mb/sum(mb)

      if(add_y2){
        mb2 = Theta%*%diag(s*(Ng/G))%*%b2
        thetab2 = mb2/sum(mb2)
      }


    }else{

      # bulk data gene relative expression.
      ## generate theta_i for bulk individual
      if(!is.null(mean_to_var_sigma)){
        X_b = matrix(rgamma(Ng*K,Theta*mean_to_var_sigma,rate=mean_to_var_sigma),ncol=K)
      }else{
        X_b = matrix(rgamma(Ng*K,Theta^2/sigma2,rate=Theta/sigma2),ncol=K)
      }

      mb = X_b%*%diag(s*(Ng/G))%*%b
      thetab = mb/sum(mb)

      #mb = Theta%*%diag(s*(Ng/G))%*%b
      #thetab = mb/sum(mb)

      if(add_y2){

        if(!is.null(mean_to_var_sigma)){
          X_b2 = matrix(rgamma(Ng*K,Theta*mean_to_var_sigma,rate=mean_to_var_sigma),ncol=K)
        }else{
          X_b2 = matrix(rgamma(Ng*K,Theta^2/sigma2,rate=Theta/sigma2),ncol=K)
        }
        mb2 = X_b2%*%diag(s*(Ng/G))%*%b2
        thetab2 = mb2/sum(mb2)

        #mb2 = Theta%*%diag(s*(Ng/G))%*%b2
        #thetab2 = mb2/sum(mb2)
      }

    }

    y = rpois(Ng,bulk_lib_size*Ng*thetab)

    if(add_y2){
      y2 = rpois(Ng,bulk_lib_size*Ng*thetab2)
      y = cbind(y,y2)
    }


    #browser()

    if(weight == 'equal'){
      w = rep(1,Ng)
    }else if(weight == 'default'){
      w = NULL
    }else if(weight == 'external'){
      w = rpois(Ng,bulk_lib_size*Ng*thetab)
    }


    if(ref_type=='multi_sc'){

      Cr = rnbinom(nk*n_indi*K,sc_lib_size*Ng,0.5)+1

      Y = matrix(nrow=Ng,ncol=nk*n_indi*K)
      indi_idx = rep(1:n_indi,each = K*nk)
      cell_type = c()
      for(i in 1:n_indi){
        #print(i)
        indi_cell_idx = which(indi_idx==i)
        indi_celltype_idx = rep(1:K,each=nk)

        if(!is.null(mean_to_var_sigma)){
          X_i = matrix(rgamma(Ng*K,Theta*mean_to_var_sigma,rate=mean_to_var_sigma),ncol=K)
        }else{
          X_i = matrix(rgamma(Ng*K,Theta^2/sigma2,rate=Theta/sigma2),ncol=K)
        }

        Y_i = matrix(nrow=Ng,ncol=K*nk)
        for(k in 1:K){
          cc = which(indi_celltype_idx==k)
          ag = X_i[,k]%*%t(rep(1,nk))


          if(!is.null(mean_to_var_tau)){
            Y_i[,cc] = matrix(rgamma(Ng*nk,ag*mean_to_var_tau,rate=mean_to_var_tau),ncol=nk)
          }else{
            atau2 = tau2[,k]%*%t(rep(1,nk))
            Y_i[,cc] = matrix(rgamma(Ng*nk,ag^2/atau2,rate=ag/atau2),ncol=nk)
          }


        }
        cell_type = c(cell_type,indi_celltype_idx)
        Y[,indi_cell_idx] = Y_i
      }

      #Y = matrix(rnorm(Ng*n_indi*K*nk,Y,sqrt(sc_noise_var)),ncol=n_indi*K*nk)
      #y = rnorm(Ng,Xb,sqrt(bulk_noise_var))

      Y = matrix(rpois(Ng*nk*n_indi*K,t(t(Y)*Cr/median(colSums(Y)))),ncol=nk*n_indi*K)

      rownames(Y) = gene_names


      if(sigma2known){
        if(tau2known){
          data.obj = set_data_decon(y,Y,
                                    ref_type = 'multi_sc',
                                    #marker_gene=marker_gene,
                                    indi_idx=indi_idx,
                                    cell_type_idx=cell_type,
                                    tau2 = tau2,sigma2 = sigma2,
                                    w=w,gene_thresh=gene_thresh)

        }else{

          data.obj = set_data_decon(y,Y,
                                    ref_type = 'multi_sc',
                                    #marker_gene=marker_gene,
                                    indi_idx=indi_idx,
                                    cell_type_idx=cell_type,
                                    tau2 = NULL,sigma2 = sigma2,
                                    w=w,gene_thresh=gene_thresh)

        }
      }else{
        data.obj = set_data_decon(y,Y,
                                  ref_type = 'multi_sc',
                                  #marker_gene=marker_gene,
                                  indi_idx=indi_idx,
                                  cell_type_idx=cell_type,
                                  tau2 = NULL,sigma2 = NULL,
                                  w=w,gene_thresh=gene_thresh)

      }


      fit_adj = deconference(data.obj,marker_gene=marker_gene,hc.type=hc.type,
                             est_pop_var=est_pop_var,x_estimator=x_estimator,
                             meta_var=meta_var,meta_mode=meta_mode,
                             correction=correction,eps=eps)
      sigma2_hat[[reps]] = fit_adj$input$Sigma



      fit_unadj = unadjusted_lm(fit_adj$input$y,fit_adj$input$X,w=fit_adj$input$w)


      #browser()

    }



    #reference data
    if(ref_type=='sc'){

      #t1=Sys.time()

      if(!same_indi){

        if(!is.null(mean_to_var_sigma)){
          X_i = matrix(rgamma(Ng*K,Theta*mean_to_var_sigma,rate=mean_to_var_sigma),ncol=K)
        }else{
          X_i = matrix(rgamma(Ng*K,Theta^2/sigma2,rate=Theta/sigma2),ncol=K)
        }

      }else{
          X_i = Theta
        }


      cell_type = rep(1:K,each=nk)
      Cr = rnbinom(nk*K,sc_lib_size*Ng,0.5)+1
      Y = matrix(nrow=Ng,ncol=nk*K)
      #tau2 = matrix(nrow = Ng,ncol = K)
      for(k in 1:K){
        cell_idx = which(cell_type==k)
        #Y[,cell_idx] = t(rdirichlet(nk,Theta[,k]*Ng))
        # Y[,cell_idx] = matrix(rgamma(Ng*nk,Theta[,k]%*%t(rep(1,nk))*snr,rate=1),ncol=nk)/snr
        ag = X_i[,k]%*%t(rep(1,nk))
        atau2 = tau2[,k]%*%t(rep(1,nk))

        if(!is.null(mean_to_var_tau)){
          Y[,cell_idx] = matrix(rgamma(Ng*nk,ag*mean_to_var_tau,rate=mean_to_var_tau),ncol=nk)
        }else{
          atau2 = tau2[,k]%*%t(rep(1,nk))
          Y[,cell_idx] = matrix(rgamma(Ng*nk,ag^2/atau2,
                                       rate=ag/atau2),ncol=nk)
        }


        # tau2[,k] = Theta[,k]*(1-Theta[,k]) / (Ng+1)
      }
      #Y = matrix(rpois(Ng*nk*K,Y%*%diag(Cr)),ncol=nk*K)
      Y = matrix(rpois(Ng*nk*K,t(t(Y)*Cr/median(colSums(Y)))),ncol=nk*K)
      rownames(Y) = gene_names

      #t2=Sys.time()

      #print(paste('generate data takes',t2-t1))



      if(tau2known){
        data.obj = set_data_decon(y,Y,
                                  ref_type = ref_type,
                                  #marker_gene=marker_gene,
                                  cell_type_idx=cell_type,
                                  tau2 = tau2,sigma2 = sigma2,
                                  w=w,gene_thresh=gene_thresh)

      }else{
        data.obj = set_data_decon(y,Y,
                                  ref_type = ref_type,
                                  #marker_gene=marker_gene,
                                  cell_type_idx=cell_type,
                                  tau2 = NULL,sigma2 = sigma2,
                                  w=w,gene_thresh=gene_thresh)

      }

      fit_adj = deconference(data.obj,marker_gene=marker_gene,hc.type=hc.type,
                             x_estimator=x_estimator,a=a,correction=correction,eps=eps)
      fit_unadj = unadjusted_lm(fit_adj$input$y,fit_adj$input$X,w=fit_adj$input$w)


      #t3 = Sys.time()
      #print(paste('fit model takes',t3-t2))
    }
    if(ref_type=='bulk'){

      #t1=Sys.time()


      if(!same_indi){

        if(!is.null(mean_to_var_sigma)){
          X_i = matrix(rgamma(Ng*K,Theta*mean_to_var_sigma,rate=mean_to_var_sigma),ncol=K)
        }else{
          X_i = matrix(rgamma(Ng*K,Theta^2/sigma2,rate=Theta/sigma2),ncol=K)
        }

      }else{
        X_i = Theta
      }


      Cr = rpois(K,ref_lib_size*Ng)+1
      #Cr = rep(ref_lib_size,K)
      U = diag(Cr)
      #Theta = apply(ref_rep,2,function(z){z/sum(z)})
      Y = matrix(rpois(Ng*K,X_i%*%U/median(colSums(X_i))),ncol=K)
      rownames(Y) = gene_names

      #t2=Sys.time()

      #print(paste('generate data takes',t2-t1))


      data.obj = set_data_decon(y,Y,ref_type = ref_type,
                                #marker_gene=marker_gene,
                                w=w)

      fit_adj = deconference(data.obj,marker_gene=marker_gene,hc.type=hc.type,
                             correction=correction,eps=eps)
      fit_unadj = unadjusted_lm(fit_adj$input$y,fit_adj$input$X,w=fit_adj$input$w)

      #t3 = Sys.time()
      #print(paste('fit model takes',t3-t2))
    }


    data_sparsity[reps] = sum(Y==0)/prod(dim(Y))

    if(add_y2){
      est_adj[reps,] = c(fit_adj$beta_hat)
      se_adj[reps,] = c(fit_adj$beta_se)
      diff_adj[reps,] = fit_adj$beta_hat[,1] - fit_adj$beta_hat[,2]
      diff_adj_se[reps,] = sqrt(fit_adj$beta_se[,1]^2 + fit_adj$beta_se[,2]^2 - 2*diag(fit_adj$cov_beta_hat[1:K,(K+1):(2*K)]))
      diff_adj_p[reps,] = (1-pnorm(abs(diff_adj[reps,])/diff_adj_se[reps,]))*2


      est_unadj[reps,] = c(fit_unadj$beta_hat)
      se_unadj[reps,] = c(fit_unadj$sand.out$beta_se)
      diff_unadj[reps,] = fit_unadj$beta_hat[,1] - fit_unadj$beta_hat[,2]
      diff_unadj_se[reps,] = sqrt(fit_unadj$sand.out$beta_se[,1]^2 + fit_unadj$sand.out$beta_se[,2]^2 - 2*diag(fit_unadj$sand.out$cov_beta_hat[1:K,(K+1):(2*K)]))
      diff_unadj_p[reps,] = (1-pnorm(abs(diff_unadj[reps,])/diff_unadj_se[reps,]))*2

      se_unadj_cv[reps,] = c(fit_unadj$ols.out$beta_se)
      diff_unadj_se_cv[reps,] = sqrt(fit_unadj$ols.out$beta_se[,1]^2 + fit_unadj$ols.out$beta_se[,2]^2 - 2*diag(fit_unadj$ols.out$cov_beta_hat[1:K,(K+1):(2*K)]))
      diff_unadj_p_cv[reps,] = (1-pnorm(abs(diff_unadj[reps,])/diff_unadj_se_cv[reps,]))*2



      beta_hat_cov_adj[[reps]] = fit_adj$cov_beta_hat
      beta_hat_cov_unadj[[reps]] = fit_unadj$sand.out$cov_beta_hat
      beta_hat_cov_unadj_cv[[reps]] = fit_unadj$ols.out$cov_beta_hat
    }else{
      est_adj[reps,] = fit_adj$beta_hat
      se_adj[reps,] = fit_adj$beta_se

      est_unadj[reps,] = fit_unadj$beta_hat
      se_unadj[reps,] = fit_unadj$sand.out$beta_se
      se_unadj_cv[reps,] = fit_unadj$ols.out$beta_se
    }



    #browser()

  }


  ########
  #1. mean squared error, and standard error
  mean_est_adj = apply(est_adj,2,mean)
  mse_adj = apply((est_adj - rep(1,nreps)%*%t(c(p,p2)))^2,2,mean)
  se_est_adj = apply(est_adj,2,sd)
  mean_est_unadj = apply(est_unadj,2,mean)
  mse_unadj = apply((est_unadj - rep(1,nreps)%*%t(c(p,p2)))^2,2,mean)
  se_est_unadj = apply(est_unadj,2,sd)
  #mean_se_adj = apply(se_adj,2,mean)
  #mean_se_unadj = apply(se_unadj,2,mean)

  ###################
  #2. coverage

  ci_l = est_adj - qnorm(1-alpha/2)*se_adj
  ci_r = est_adj + qnorm(1-alpha/2)*se_adj
  covergae_adj = ((rep(1,nreps)%*%t(c(p,p2)))>=ci_l) & ((rep(1,nreps)%*%t(c(p,p2)))<=ci_r)
  covergae_adj=apply(covergae_adj,2,mean)


  ## sand
  ci_l = est_unadj - qnorm(1-alpha/2)*se_unadj
  ci_r = est_unadj + qnorm(1-alpha/2)*se_unadj
  covergae_unadj = ((rep(1,nreps)%*%t(c(p,p2)))>=ci_l) & ((rep(1,nreps)%*%t(c(p,p2)))<=ci_r)
  covergae_unadj=apply(covergae_unadj,2,mean)


  ## cv
  ci_l = est_unadj - qnorm(1-alpha/2)*se_unadj_cv
  ci_r = est_unadj + qnorm(1-alpha/2)*se_unadj_cv
  covergae_unadj_cv = ((rep(1,nreps)%*%t(c(p,p2)))>=ci_l) & ((rep(1,nreps)%*%t(c(p,p2)))<=ci_r)
  covergae_unadj_cv=apply(covergae_unadj_cv,2,mean)




  ###################
  #3. compare two individual bulk data
  if(add_y2){

    ci_l = diff_adj - qnorm(1-alpha/2)*diff_adj_se
    ci_r = diff_adj + qnorm(1-alpha/2)*diff_adj_se
    covergae_diff_adj = ((rep(1,nreps)%*%t(c(p-p2)))>=ci_l) & ((rep(1,nreps)%*%t(c(p-p2)))<=ci_r)
    covergae_diff_adj=apply(covergae_diff_adj,2,mean)

    ci_l = diff_unadj - qnorm(1-alpha/2)*diff_unadj_se
    ci_r = diff_unadj + qnorm(1-alpha/2)*diff_unadj_se
    covergae_diff_unadj = ((rep(1,nreps)%*%t(c(p-p2)))>=ci_l) & ((rep(1,nreps)%*%t(c(p-p2)))<=ci_r)
    covergae_diff_unadj=apply(covergae_diff_unadj,2,mean)

    ci_l = diff_unadj - qnorm(1-alpha/2)*diff_unadj_se_cv
    ci_r = diff_unadj + qnorm(1-alpha/2)*diff_unadj_se_cv
    covergae_diff_unadj_cv = ((rep(1,nreps)%*%t(c(p-p2)))>=ci_l) & ((rep(1,nreps)%*%t(c(p-p2)))<=ci_r)
    covergae_diff_unadj_cv = apply(covergae_diff_unadj_cv,2,mean)

    power_adj = apply(diff_adj_p<=alpha,2,mean)
    power_unadj = apply(diff_unadj_p<=alpha,2,mean)
    power_unadj_cv = apply(diff_unadj_p_cv<=alpha,2,mean)


    return(list(p=p,p2=p2,
                covergae_diff_adj=covergae_diff_adj,
                covergae_diff_unadj = covergae_diff_unadj,
                covergae_diff_unadj_cv = covergae_diff_unadj_cv,
                beta_hat_cov_adj=beta_hat_cov_adj,
                beta_hat_cov_unadj = beta_hat_cov_unadj,
                beta_hat_cov_unadj_cv = beta_hat_cov_unadj_cv,
                power_adj=power_adj,
                power_unadj=power_unadj,
                power_unadj_cv=power_unadj_cv,
                diff_adj_p = diff_adj_p,
                diff_unadj_p=diff_unadj_p,
                diff_unadj_p_cv=diff_unadj_p_cv,
                diff_adj=diff_adj,
                diff_unadj=diff_unadj,
                diff_adj_se=diff_adj_se,
                diff_unadj_se=diff_unadj_se,
                diff_unadj_se_cv=diff_unadj_se_cv,
                data_sparsity = data_sparsity,

                mean_est_adj=mean_est_adj,
                mse_adj=mse_adj,
                se_est_adj=se_est_adj,

                mean_est_unadj=mean_est_unadj,
                mse_unadj=mse_unadj,
                se_est_unadj=se_est_unadj,

                covergae_adj=covergae_adj,
                covergae_unadj=covergae_unadj,
                covergae_unadj_cv=covergae_unadj_cv,

                se_adj = se_adj,
                se_unadj = se_unadj,
                se_unadj_cv = se_unadj_cv,
                est_adj = est_adj,
                est_unadj = est_unadj,
                sigma2_hat=sigma2_hat,
                simu_param = list(ref=ref,s=s,Ng=Ng,b=b,ref_type=ref_type,tau2known=tau2known,
                                  sc_lib_size = sc_lib_size,ref_lib_size=ref_lib_size,
                                  bulk_lib_size = bulk_lib_size,nk=nk,
                                  nreps=nreps,alpha=alpha,a=a,correction=correction)))

  }else{
    return(list(p=p,
                mean_est_adj=mean_est_adj,
                mse_adj=mse_adj,
                se_est_adj=se_est_adj,
                mean_est_unadj=mean_est_unadj,
                mse_unadj=mse_unadj,
                se_est_unadj=se_est_unadj,
                covergae_adj=covergae_adj,
                covergae_unadj=covergae_unadj,
                covergae_unadj_cv=covergae_unadj_cv,
                se_adj = se_adj,
                se_unadj = se_unadj,
                se_unadj_cv = se_unadj_cv,
                est_adj = est_adj,
                est_unadj = est_unadj,
                sigma2_hat=sigma2_hat,
                simu_param = list(ref=ref,s=s,Ng=Ng,b=b,ref_type=ref_type,tau2known=tau2known,
                                  sc_lib_size = sc_lib_size,ref_lib_size=ref_lib_size,
                                  bulk_lib_size = bulk_lib_size,nk=nk,
                                  nreps=nreps,alpha=alpha,a=a,correction=correction)))
  }





}


############## test ################
# set.seed(12345)
#
# G = 1000
# K = 4
# b = 1:K
# b = b/sum(b)
# b2 = b
# library(gtools)
# ref = t(rdirichlet(K,rep(1,G)))
# test = simu_study(ref,G,b,ref_type='bulk',nreps = 100,b2=b2)
#
# test = simu_study(ref,G,b,ref_type='sc',nreps = 100,b2=b2)
#
# test = simu_study(ref,G,b,ref_type='bulk',nreps = 100,b2=b2)

