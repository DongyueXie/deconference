
########################################
############ simulation function for multiple individuals
########################################






### multiple individual #####
#'@description sigma2 is the gene relative expresison variance across individuals.
#'@param n_indi number of individuals
#'@param sigma2known assume sigma2 known or not
#'@param sigma2 value of sigma2.
#'@param est_sigma2 whether estimate sigma2 or directly use sample var/cov.
simu_study_multiInd = function(ref,Ng,b,ref_type='sc',
                               n_indi = 6,
                               bulk_lib_size = 50,
                               sc_lib_size = 0.2,
                               #ref_lib_size = 5,
                               tau2known = FALSE,
                               tau2 = matrix(runif(Ng*length(b),1/(Ng*10),1/(Ng*10)),nrow=Ng),
                               #snr=3,
                               x_estimator='aggregate',
                               nk = 100,
                               sigma2known=FALSE,
                               sigma2 = matrix(runif(Ng*length(b),1/(Ng*100),1/(Ng*100)),nrow=Ng),
                               est_sigma2 = FALSE,
                               meta_var='adjust',
                               addw=TRUE,
                               eps=0,
                               nreps=100,alpha=0.05,a=0,
                               correction=FALSE,s,printevery=10,
                               add_y2 = TRUE, b2=NULL){

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
  }else{
    est_adj = matrix(nrow=nreps,ncol=K)
    se_adj = matrix(nrow=nreps,ncol=K)

    est_unadj = matrix(nrow=nreps,ncol=K)
    se_unadj = matrix(nrow=nreps,ncol=K)
  }



  sigma2_hat = list()

  data_sparsity = c()
  beta_hat_cov = list()


  for(reps in 1:nreps){

    if(reps%%printevery==0){print(sprintf("done %d (out of %d)",reps,nreps))}

    #Obtain X
    ref_rep = ref[sample(1:G,Ng),]
    Theta = apply(ref_rep,2,function(z){z/sum(z)})

    ##### bulk data gene relative expression.

    ## generate theta_i for bulk individual

    X_i = matrix(rgamma(Ng*K,Theta^2/sigma2,rate=Theta/sigma2),ncol=K)



    # X_i = matrix(rgamma(Ng*K,Theta*snr,rate=1),ncol=K)/snr

    # X_i = matrix(nrow=Ng,ncol=K)
    # for(k in 1:K){
    #   X_i[,k] = c(rdirichlet(1,Theta[,k]*Ng))
    # }

    #mb = Theta%*%diag(s*(Ng/G))%*%b
    mb = X_i%*%diag(s*(Ng/G))%*%b
    Xb = mb/sum(mb)

    if(add_y2){
      X_i2 = matrix(rgamma(Ng*K,Theta^2/sigma2,rate=Theta/sigma2),ncol=K)
      mb2 = X_i2%*%diag(s*(Ng/G))%*%b2
      Xb2 = mb2/sum(mb2)
    }

    if(addw){
      w = runif(Ng)
    }else{
      w=NULL
    }

    #reference data
    if(ref_type=='sc'){

      #### data generation steps #####
      # 1. for each cell type, draw gene relative expression for each individual
      # 2. for each cell type of each individual, draw single cells gene relative expression and cell library size.


      ###########################################
      #
      ## cell library size
      Cr = rnbinom(nk*n_indi*K,sc_lib_size*Ng,0.5)+1

      Y = matrix(nrow=Ng,ncol=nk*n_indi*K)
      indi_idx = rep(1:n_indi,each = K*nk)
      cell_type = c()
      for(i in 1:n_indi){
        #print(i)
        indi_cell_idx = which(indi_idx==i)
        indi_celltype_idx = rep(1:K,each=nk)
        X_i = matrix(rgamma(Ng*K,Theta^2/sigma2,rate=Theta/sigma2),ncol=K)
        Y_i = matrix(nrow=Ng,ncol=K*nk)
        for(k in 1:K){
          cc = which(indi_celltype_idx==k)
          ag = X_i[,k]%*%t(rep(1,nk))
          atau2 = tau2[,k]%*%t(rep(1,nk))
          Y_i[,cc] = matrix(rgamma(Ng*nk,ag^2/atau2,rate=ag/atau2),ncol=nk)
        }
        cell_type = c(cell_type,indi_celltype_idx)
        Y[,indi_cell_idx] = Y_i
      }

      #Y = matrix(rnorm(Ng*n_indi*K*nk,Y,sqrt(sc_noise_var)),ncol=n_indi*K*nk)
      #y = rnorm(Ng,Xb,sqrt(bulk_noise_var))

      Y = matrix(rpois(Ng*nk*n_indi*K,t(t(Y)*Cr)),ncol=nk*n_indi*K)


      y = rpois(Ng,rpois(1,bulk_lib_size)*Ng*Xb)

      if(add_y2){
        y2 = rpois(Ng,rpois(1,bulk_lib_size)*Ng*Xb2)
        y = cbind(y,y2)
      }

      data_sparsity[reps] = sum(Y==0)/prod(dim(Y))
      ############################################

      # cell_type = rep(1:K,each=(nk*n_indi))
      #
      # ## cell library size
      # Cr = rnbinom(nk*n_indi*K,sc_lib_size*Ng,0.5)+1
      #
      # Y = matrix(nrow=Ng,ncol=nk*n_indi*K)
      # #tau2 = matrix(nrow = Ng,ncol = K)
      # for(k in 1:K){
      #   cell_idx = which(cell_type==k)
      #   Yik = c()
      #   for(i in 1:n_indi){
      #     Theta_ik = rgamma(Ng,Theta[,k]^2/sigma2,rate=Theta[,k]/sigma2)
      #     Yik = matrix(rgamma(Ng*nk,(Theta_ik%*%t(rep(1,nk)))^2/tau2,
      #                         rate=(Theta_ik%*%t(rep(1,nk)))/tau2),ncol=nk)
      #     # Theta_ik = rgamma(Ng,Theta[,k]*snr,rate=1)/snr
      #     # Yik = matrix(rgamma(Ng*nk,Theta_ik%*%t(rep(1,nk))*snr,rate=1),ncol=nk)/snr
      #     #Theta_ik = c(rdirichlet(1,Theta[,k]*Ng))
      #     #Yik = cbind(Yik,t(rdirichlet(nk,Theta_ik*Ng)))
      #   }
      #   Y[,cell_idx] = Yik
      #   #tau2[,k] = Theta[,k]*(1-Theta[,k]) / (Ng+1)
      # }
      # Y = matrix(rpois(Ng*nk*n_indi*K,t(t(Y)*Cr)),ncol=nk*n_indi*K)
      # y = rpois(Ng,rpois(1,bulk_lib_size)*Ng*Xb)
      # indi_idx = rep(rep(1:n_indi,each=nk),K)

      #browser()
      if(sigma2known){
        if(tau2known){
          fit_adj = deconference(y,Y,ref_type='multi_sc',cell_type=cell_type,indi_idx=indi_idx,w=w,est_sigma2=est_sigma2,
                                 sigma2=matrix(sigma2,nrow=Ng,ncol=K),x_estimator=x_estimator,
                                 tau2 = matrix(tau2,nrow=Ng,ncol=K),meta_var=meta_var,
                                 adjust = TRUE,alpha=alpha,a=a,correction=correction,eps=eps)
        }else{
          fit_adj = deconference(y,Y,ref_type='multi_sc',cell_type=cell_type,indi_idx=indi_idx,w=w,est_sigma2=est_sigma2,
                                 sigma2=matrix(sigma2,nrow=Ng,ncol=K),meta_var=meta_var,
                                 tau2=NULL,
                                 x_estimator=x_estimator,
                                 adjust = TRUE,alpha=alpha,a=a,correction=correction,eps=eps)
        }
      }else{
        fit_adj = deconference(y,Y,ref_type='multi_sc',cell_type=cell_type,indi_idx=indi_idx,w=w,est_sigma2=est_sigma2,
                               sigma2=NULL,tau2=NULL,meta_var=meta_var,x_estimator=x_estimator,
                               adjust = TRUE,alpha=alpha,a=a,correction=correction,eps=eps)
        sigma2_hat[[reps]] = fit_adj$sigma2
      }

      fit_unadj = estimation_func(fit_adj$input$y,fit_adj$input$X,fit_adj$input$Vg,w=fit_adj$input$w,
                                                adjust=FALSE,
                                                alpha=alpha)

      #fit_unadj = deconference(y,Y,cell_type=cell_type,indi_idx=indi_idx,w=w,adjust = FALSE,alpha=alpha)

    }
    # if(ref_type=='bulk'){
    #   Cr = rpois(K,ref_lib_size*Ng)
    #   if(prod(Cr)==0){Cr[which(Cr==0)]=ref_lib_size*Ng}
    #   #Cr = rep(ref_lib_size,K)
    #   U = diag(Cr)
    #   Y = matrix(rpois(Ng*K,Theta%*%U),ncol=K)
    #   #bulk data
    #   y = rpois(Ng,rpois(1,bulk_lib_size)*Ng*Xb)
    #   fit_adj = deconference(y,Y,w=w,adjust = TRUE,alpha=alpha,a=a,correction=correction)
    #   fit_unadj = deconference(y,Y,w=w,adjust = FALSE,alpha=alpha)
    #
    # }

    if(add_y2){
      est_adj[reps,] = c(fit_adj$beta_hat)
      se_adj[reps,] = c(fit_adj$beta_se)

      diff_adj[reps,] = fit_adj$beta_hat[,1] - fit_adj$beta_hat[,2]
      diff_adj_se[reps,] = sqrt(fit_adj$beta_se[,1]^2 + fit_adj$beta_se[,2]^2 - 2*diag(fit_adj$cov_beta_hat[1:K,(K+1):(2*K)]))
      #diff_adj_p[reps,] = fit_adj$p_value_diff


      est_unadj[reps,] = c(fit_unadj$beta_hat)
      se_unadj[reps,] = c(fit_unadj$beta_se)

      diff_unadj[reps,] = fit_unadj$beta_hat[,1] - fit_unadj$beta_hat[,2]
      diff_unadj_se[reps,] = sqrt(fit_unadj$beta_se[,1]^2 + fit_unadj$beta_se[,2]^2 - 2*diag(fit_unadj$cov_beta_hat[1:K,(K+1):(2*K)]))
      #diff_unadj_p[reps,] = fit_unadj$p_value_diff

      beta_hat_cov[[reps]] = fit_adj$cov_beta_hat
    }else{
      est_adj[reps,] = fit_adj$beta_hat
      se_adj[reps,] = fit_adj$beta_se

      est_unadj[reps,] = fit_unadj$beta_hat
      se_unadj[reps,] = fit_unadj$beta_se
    }
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

  ci_l = est_unadj - qnorm(1-alpha/2)*se_unadj
  ci_r = est_unadj + qnorm(1-alpha/2)*se_unadj
  covergae_unadj = ((rep(1,nreps)%*%t(c(p,p2)))>=ci_l) & ((rep(1,nreps)%*%t(c(p,p2)))<=ci_r)
  covergae_unadj=apply(covergae_unadj,2,mean)

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

    #power_adj = apply(diff_adj_p<=alpha,2,mean)
    #power_unadj = apply(diff_unadj_p<=alpha,2,mean)


    return(list(p=p,p2=p2,
                sigma2_hat = sigma2_hat,
                covergae_diff_adj=covergae_diff_adj,
                covergae_diff_unadj = covergae_diff_unadj,
                beta_hat_cov=beta_hat_cov,
                #power_adj=power_adj,
                #power_unadj=power_unadj,
                diff_adj_p = diff_adj_p,
                diff_unadj_p=diff_unadj_p,
                diff_adj=diff_adj,
                diff_unadj=diff_unadj,
                diff_adj_se=diff_adj_se,
                diff_unadj_se=diff_unadj_se,

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
                est_adj = est_adj,
                est_unadj = est_unadj,
                data_sparsity = data_sparsity,
                simu_param = list(ref=ref,s=s,Ng=Ng,b=b,ref_type=ref_type,addw=addw,
                                  sc_lib_size = sc_lib_size,
                                  bulk_lib_size = bulk_lib_size,
                                  nk=nk,
                                  tau2=tau2,sigma2=sigma2,sigma2known=sigma2known,
                                  nreps=nreps,alpha=alpha,a=a,correction=correction)))

  }else{
    return(list(p=p,
                sigma2_hat = sigma2_hat,
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
                est_adj = est_adj,
                est_unadj = est_unadj,
                data_sparsity = data_sparsity,
                simu_param = list(ref=ref,s=s,Ng=Ng,b=b,ref_type=ref_type,addw=addw,
                                  sc_lib_size = sc_lib_size,
                                  bulk_lib_size = bulk_lib_size,
                                  nk=nk,
                                  tau2=tau2,sigma2=sigma2,sigma2known=sigma2known,
                                  nreps=nreps,alpha=alpha,a=a,correction=correction)))
  }

}
