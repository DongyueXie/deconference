
########################################
############ simulation function for multiple individuals:
############ simulate normal data and check if codes are correct.
########################################






### multiple individual #####
#'@description sigma2 is the gene relative expresison variance across individuals.
#'@param n_indi number of individuals
#'@param sigma2known assume sigma2 known or not
#'@param sigma2 value of sigma2. If specified, assume all the same cross genes and cell type.
#'@param est_sigma2 whether estimate sigma2 or directly use sample var/cov.
simu_study_multiInd_normal = function(ref,Ng,b,ref_type='sc',
                                      n_indi = 6,
                                      bulk_noise_var = 1,
                                      sc_noise_var = 1,
                                      #ref_lib_size = 5,
                                      tau2known = FALSE,
                                      tau2 = 1/Ng^2,
                                      #snr=3,
                                      x_estimator = 'aggregate',
                                      nk = 100,
                                      sigma2known=TRUE,
                                      sigma2 = 1/Ng^2,
                                      est_sigma2 = TRUE,
                                      meta_var='adjust',
                                      addw=TRUE,
                                      nreps=100,alpha=0.05,a=0,
                                      correction=FALSE,s,printevery=10){

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




  sigma2_hat = list()


  for(reps in 1:nreps){

    if(reps%%printevery==0){print(sprintf("done %d (out of %d)",reps,nreps))}

    #Obtain X
    #ref_rep = ref[sample(1:G,Ng),]
    #Theta = ref_rep
    Theta = ref

    ##### bulk data gene relative expression.

    ## generate theta_i for bulk individual

    X_ib = matrix(rnorm(Ng*K,Theta,sqrt(sigma2)),ncol=K)
    # X_i = matrix(rgamma(G*K,Theta^2/sigma2,rate=Theta/sigma2),ncol=K)

    # X_i = matrix(rgamma(Ng*K,Theta*snr,rate=1),ncol=K)/snr

    # X_i = matrix(nrow=Ng,ncol=K)
    # for(k in 1:K){
    #   X_i[,k] = c(rdirichlet(1,Theta[,k]*Ng))
    # }

    #mb = Theta%*%diag(s*(Ng/G))%*%b
    mb = X_ib%*%diag(s*(Ng/G))%*%b
    Xb = mb

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
      Y = matrix(nrow=Ng,ncol=nk*n_indi*K)
      indi_idx = rep(1:n_indi,each = K*nk)
      cell_type = c()
      for(i in 1:n_indi){
        #print(i)
        indi_cell_idx = which(indi_idx==i)
        indi_celltype_idx = rep(1:K,each=nk)
        X_i = matrix(rnorm(G*K,Theta,sqrt(sigma2)),ncol=K)
        Y_i = matrix(nrow=Ng,ncol=K*nk)
        for(k in 1:K){
          cc = which(indi_celltype_idx==k)
          Y_i[,cc] = matrix(rnorm(Ng*nk,X_i[,k]%*%t(rep(1,nk)),sqrt(tau2)),ncol=nk)
          #Y_i[,cc] = t(MASS::mvrnorm(nk,X_i[,k],diag(Ng)*sqrt(tau2)))
        }
        cell_type = c(cell_type,indi_celltype_idx)
        Y[,indi_cell_idx] = Y_i
      }

      Y = matrix(rnorm(Ng*n_indi*K*nk,Y,sqrt(sc_noise_var)),ncol=n_indi*K*nk)
      y = rnorm(Ng,Xb,sqrt(bulk_noise_var))
      ############################################



      # ############################################
      # cell_type = rep(1:K,each=(nk*n_indi))
      #
      # ## cell library size
      #
      # # Cr = rnbinom(nk*n_indi*K,sc_lib_size*Ng,0.5)+1
      #
      # Y = matrix(nrow=Ng,ncol=nk*n_indi*K)
      # #tau2 = matrix(nrow = Ng,ncol = K)
      # for(k in 1:K){
      #   cell_idx = which(cell_type==k)
      #   Yik = c()
      #   for(i in 1:n_indi){
      #     #Theta_ik = rgamma(Ng,Theta[,k]^2/sigma2,rate=Theta[,k]/sigma2)
      #     Theta_ik = rnorm(Ng,Theta[,k],sqrt(sigma2))
      #     Yik = matrix(rnorm(Ng*nk,Theta_ik%*%t(rep(1,nk)),sqrt(tau2)),ncol=nk)
      #     # Theta_ik = rgamma(Ng,Theta[,k]*snr,rate=1)/snr
      #     # Yik = matrix(rgamma(Ng*nk,Theta_ik%*%t(rep(1,nk))*snr,rate=1),ncol=nk)/snr
      #     #Theta_ik = c(rdirichlet(1,Theta[,k]*Ng))
      #     #Yik = cbind(Yik,t(rdirichlet(nk,Theta_ik*Ng)))
      #   }
      #   Y[,cell_idx] = Yik
      #   #tau2[,k] = Theta[,k]*(1-Theta[,k]) / (Ng+1)
      # }
      # Y = matrix(rnorm(Ng*nk*n_indi*K,Y,sqrt(sc_noise_var)),ncol=nk*n_indi*K)
      #
      # #Xb = Theta%*%diag(s*(Ng/G))%*%b
      # y = rnorm(Ng,Xb,sqrt(bulk_noise_var))
      #
      # indi_idx = rep(rep(1:n_indi,each=nk),K)
      # #############################################################################
      #browser()
      if(sigma2known){
        if(tau2known){
          #browser()
          fit_adj = deconference(y,Y,cell_type,indi_idx,w=w,est_sigma2=est_sigma2,
                                 sigma2=matrix(sigma2,nrow=Ng,ncol=K),x_estimator=x_estimator,
                                 tau2 = matrix(tau2,nrow=Ng,ncol=K),meta_var=meta_var,
                                 adjust = TRUE,alpha=alpha,a=a,correction=correction)
        }else{
          fit_adj = deconference(y,Y,cell_type,indi_idx,w=w,est_sigma2=est_sigma2,x_estimator=x_estimator,
                                 sigma2=matrix(sigma2,nrow=Ng,ncol=K),meta_var=meta_var,
                                 adjust = TRUE,alpha=alpha,a=a,correction=correction)
        }
      }else{
        fit_adj = deconference(y,Y,cell_type,indi_idx,w=w,est_sigma2=est_sigma2,
                               sigma2=NULL,meta_var=meta_var,x_estimator=x_estimator,
                               adjust = TRUE,alpha=alpha,a=a,correction=correction)
        sigma2_hat[[reps]] = fit_adj$sigma2
      }

      fit_unadj = deconference(y,Y,cell_type,indi_idx,w=w,adjust = FALSE,alpha=alpha,x_estimator=x_estimator,)

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

    est_adj[reps,] = fit_adj$beta_hat
    se_adj[reps,] = fit_adj$beta_se


    est_unadj[reps,] = fit_unadj$beta_hat
    se_unadj[reps,] = fit_unadj$beta_se
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
              simu_param = list(ref=ref,s=s,Ng=Ng,b=b,ref_type=ref_type,addw=addw,

                                nk=nk,
                                tau2=tau2,sigma2=sigma2,sigma2known=sigma2known,
                                nreps=nreps,alpha=alpha,a=a,correction=correction)))

}





###############################

