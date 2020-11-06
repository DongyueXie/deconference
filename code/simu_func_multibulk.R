########################################
############ simulation function for 1. bulk reference 2. same individual single cell referece 3. multi indi sc
########################################


#'@title Simulation study function of 1. bulk as reference and 2. one individual single cell reference. 3. multi indi sc as ref; multiple bulk sample
#'@param ref ref data matrix, Gene by cell type, relative expression
#'@param Ng number of genes to use in simulation, gene will be resampled each round
#'@param b cell type proportions, a K by n_sample matrix.
#'@param ref_type 'sc' for single cell ref; 'bulk' for bulk ref; multi_sc
#'@param hc.type type of sandwich estimator
#'@param tau2 a G by K matrix, entry(g,k) is the variance of gene g in cell type k, across cells
#'@param tau2known whether \tau^2, the variance of gene expr across cell, is known
#'@param sigma2 a G by K matrix, entry(g,k) is the variance of gene g in cell type k,across individuals
#'@param weight default, optimal, equal, external
#'@param marker_gene marker gene index
#'@param bulk_lib_size bulk sample library size
#'@param sc_lib_size single cell ref library size
#'@param ref_bulklib_size bulk ref data library size
#'@param nk number of single cells of each cell type
#'@param x_estimator estimate the reference matrix from single cell data, separate or aggregate
#'@param nreps number of repetition.
#'@param alpha significance level
#'@param a fuller's correction parameter a
#'@param correction whether perform fullers correction.
#'@param s cell size
#'@param s_type known, equal or estimate
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
                      ref_bulklib_size = 50,
                      nk = 100,
                      ref_scale = Ng,
                      same_indi = FALSE,
                      mean_to_var_sigma = 10,
                      mean_to_var_tau = 10,
                      tau2 = matrix(runif(nrow(ref)*length(b),0.5/(nrow(ref))^2,1/(nrow(ref))^2),nrow=nrow(ref)),
                      tau2known=FALSE,
                      x_estimator = 'aggregate',
                      nreps=100,alpha=0.05,a=0,
                      correction=TRUE,
                      s=c(1.00,1.13,1.13,0.94,1.37,0.78),
                      s_type = 'known',
                      printevery=10,
                      gene_thresh=10,
                      n_indi = 10,
                      sigma2 = matrix(runif(nrow(ref)*length(b),0.5/(nrow(ref))^3,1/(nrow(ref))^3),nrow=nrow(ref)),
                      sigma2known=FALSE,
                      est_pop_var = FALSE,meta_var='adjust',meta_mode = 'universal',eps=0,
                      groups = c(rep(1,ncol(b)/2),rep(2,ncol(b)/2))){

  G = nrow(ref)
  K = ncol(ref)
  n_bulk = ncol(b)

  ref = apply(ref,2,function(z){z/sum(z)})


  if(missing(s)){
    s = rep(1,K)
  }

  b = apply(b,2,function(z){z/sum(z)})
  p = apply(b,2,function(z){(z*s)/sum(z*s)})

  if(length(tau2)==1){
    tau2 = matrix(tau2,nrow=Ng,ncol=K)
  }

  if(length(sigma2)==1){
    sigma2 = matrix(sigma2,nrow=Ng,ncol=K)
  }

  est_adj = matrix(nrow=nreps,ncol=n_bulk*K)
  se_adj = matrix(nrow=nreps,ncol=n_bulk*K)
  diff_adj = matrix(nrow=nreps,ncol=K)
  diff_adj_se = matrix(nrow=nreps,ncol=K)
  diff_adj_p = matrix(nrow = nreps,ncol=K)

  est_unadj = matrix(nrow=nreps,ncol=n_bulk*K)
  se_unadj = matrix(nrow=nreps,ncol=n_bulk*K)
  diff_unadj = matrix(nrow=nreps,ncol=K)
  diff_unadj_se = matrix(nrow=nreps,ncol=K)
  diff_unadj_p = matrix(nrow = nreps,ncol=K)


  se_unadj_cv = matrix(nrow=nreps,ncol=n_bulk*K)
  diff_unadj_se_cv = matrix(nrow=nreps,ncol=K)
  diff_unadj_p_cv = matrix(nrow = nreps,ncol=K)

  se_unadj_hc3 = matrix(nrow=nreps,ncol=n_bulk*K)
  diff_unadj_se_hc3 = matrix(nrow=nreps,ncol=K)
  diff_unadj_p_hc3 = matrix(nrow = nreps,ncol=K)









  data_sparsity = c()
  beta_hat_cov_adj = list()
  beta_hat_cov_unadj = list()
  beta_hat_cov_unadj_cv = list()
  beta_hat_cov_unadj_hc3 = list()

  sigma2_hat = list()# unly use when ref is multi-sc and est_sigma2 IS TRURE

  for(reps in 1:nreps){


    if(reps%%printevery==0){print(sprintf("running %d (out of %d)",reps,nreps))}

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

      mb = Theta%*%diag(s)%*%b



    }else{

      # bulk data gene relative expression.
      ## generate theta_i for bulk individual
      mb = matrix(nrow=Ng,ncol=n_bulk)
      for(i in 1:n_bulk){

        if(!is.null(mean_to_var_sigma)){
          X_b = matrix(rgamma(Ng*K,Theta*mean_to_var_sigma,rate=mean_to_var_sigma),ncol=K)
        }else{
          X_b = matrix(rgamma(Ng*K,Theta^2/sigma2,rate=Theta/sigma2),ncol=K)
        }

        mb[,i] = X_b%*%diag(s)%*%b[,i]


      }
      thetab = apply(mb,2,function(z){z/sum(z)})


    }

    y = matrix(rpois(Ng*n_bulk,bulk_lib_size*Ng*thetab),nrow=Ng)




    #browser()

    if(weight == 'equal'){
      w = rep(1,Ng)
    }else if(weight == 'default'){
      w = NULL
    }else if(weight == 'external'){
      w = rpois(Ng,bulk_lib_size*Ng*rowMeans(thetab))
    }


    if(ref_type=='multi_sc'){



      #Cr = rnbinom(nk*n_indi*K,sc_lib_size*Ng,0.5)+1

      Y = matrix(nrow=Ng,ncol=nk*n_indi*K)
      indi_idx = rep(1:n_indi,each = K*nk)
      cell_type = c()
      Cr = c()
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
        Cr_i = c()
        for(k in 1:K){
          cc = which(indi_celltype_idx==k)
          ag = X_i[,k]%*%t(rep(1,nk))
          Cr_i[cc] = rnbinom(nk,sc_lib_size*Ng*s[k],0.5)+1


          if(!is.null(mean_to_var_tau)){
            Y_i[,cc] = matrix(rgamma(Ng*nk,ag*mean_to_var_tau,rate=mean_to_var_tau),ncol=nk)
          }else{
            atau2 = tau2[,k]%*%t(rep(1,nk))
            Y_i[,cc] = matrix(rgamma(Ng*nk,ag^2/atau2,rate=ag/atau2),ncol=nk)
          }


        }
        cell_type = c(cell_type,indi_celltype_idx)
        Y[,indi_cell_idx] = Y_i
        Cr = c(Cr,Cr_i)
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

      browser()

      fit_adj = deconference(data.obj,marker_gene=marker_gene,hc.type=hc.type,
                             est_pop_var=est_pop_var,x_estimator=x_estimator,
                             meta_var=meta_var,meta_mode=meta_mode,
                             correction=correction,eps=eps)
      sigma2_hat[[reps]] = fit_adj$input$Sigma



      fit_unadj = unadjusted_lm(fit_adj$input$y,fit_adj$input$X,w=fit_adj$input$w,groups)


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
                             x_estimator=x_estimator,correction=correction,eps=eps)
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


      Cr = rpois(K,ref_bulklib_size*Ng)+1
      #Cr = rep(ref_bulklib_size,K)
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


    est_adj[reps,] = c(fit_adj$beta_hat)
    se_adj[reps,] = c(fit_adj$beta_se)
    diff_out = two_group_test(fit_adj,groups)

    #browser()

    diff_adj[reps,] = c(diff_out$diff_group)
    diff_adj_se[reps,] = c(diff_out$diff_se)
    diff_adj_p[reps,] = c(diff_out$p_value)


    est_unadj[reps,] = c(fit_unadj$beta_hat)
    se_unadj[reps,] = c(fit_unadj$sand.out$beta_se)
    diff_unadj[reps,] = c(fit_unadj$diff_group)
    diff_unadj_se[reps,] = c(fit_unadj$sand.out$diff_se)
    diff_unadj_p[reps,] = fit_unadj$sand.out$p_value

    se_unadj_cv[reps,] = c(fit_unadj$ols.out$beta_se)
    diff_unadj_se_cv[reps,] = c(fit_unadj$ols.out$diff_se)
    diff_unadj_p_cv[reps,] = fit_unadj$ols.out$p_value

    se_unadj_hc3[reps,] = c(fit_unadj$sand.out.hc3$beta_se)
    diff_unadj_se_hc3[reps,] = c(fit_unadj$sand.out.hc3$diff_se)
    diff_unadj_p_hc3[reps,] = fit_unadj$sand.out.hc3$p_value



    beta_hat_cov_adj[[reps]] = fit_adj$cov_beta_hat
    beta_hat_cov_unadj[[reps]] = fit_unadj$sand.out$cov_beta_hat
    beta_hat_cov_unadj_cv[[reps]] = fit_unadj$ols.out$cov_beta_hat
    beta_hat_cov_unadj_hc3[[reps]] = fit_unadj$sand.out.hc3$cov_beta_hat



  }


  ########
  #1. mean squared error, and standard error
  mean_est_adj = apply(est_adj,2,mean)
  mse_adj = apply((est_adj - rep(1,nreps)%*%t(c(b)))^2,2,mean)
  se_est_adj = apply(est_adj,2,sd)
  mean_est_unadj = apply(est_unadj,2,mean)
  mse_unadj = apply((est_unadj - rep(1,nreps)%*%t(c(b)))^2,2,mean)
  se_est_unadj = apply(est_unadj,2,sd)
  #mean_se_adj = apply(se_adj,2,mean)
  #mean_se_unadj = apply(se_unadj,2,mean)

  ###################
  #2. coverage

  ci_l = est_adj - qnorm(1-alpha/2)*se_adj
  ci_r = est_adj + qnorm(1-alpha/2)*se_adj
  coverage_adj = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_adj=apply(coverage_adj,2,mean)


  ## sand
  ci_l = est_unadj - qnorm(1-alpha/2)*se_unadj
  ci_r = est_unadj + qnorm(1-alpha/2)*se_unadj
  coverage_unadj = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_unadj=apply(coverage_unadj,2,mean)


  ## cv
  ci_l = est_unadj - qnorm(1-alpha/2)*se_unadj_cv
  ci_r = est_unadj + qnorm(1-alpha/2)*se_unadj_cv
  coverage_unadj_cv = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_unadj_cv=apply(coverage_unadj_cv,2,mean)

  ## hc3
  ci_l = est_unadj - qnorm(1-alpha/2)*se_unadj_hc3
  ci_r = est_unadj + qnorm(1-alpha/2)*se_unadj_hc3
  coverage_unadj_hc3 = ((rep(1,nreps)%*%t(c(b)))>=ci_l) & ((rep(1,nreps)%*%t(c(b)))<=ci_r)
  coverage_unadj_hc3=apply(coverage_unadj_hc3,2,mean)

  ## diff

  ci_l = diff_adj - qnorm(1-alpha/2)*diff_adj_se
  ci_r = diff_adj + qnorm(1-alpha/2)*diff_adj_se
  coverage_diff_adj = ((rep(1,nreps)%*%t(c(b%*%diff_out$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out$a)))<=ci_r)
  coverage_diff_adj=apply(coverage_diff_adj,2,mean)

  ci_l = diff_unadj - qnorm(1-alpha/2)*diff_unadj_se
  ci_r = diff_unadj + qnorm(1-alpha/2)*diff_unadj_se
  coverage_diff_unadj = ((rep(1,nreps)%*%t(c(b%*%diff_out$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out$a)))<=ci_r)
  coverage_diff_unadj=apply(coverage_diff_unadj,2,mean)

  ci_l = diff_unadj - qnorm(1-alpha/2)*diff_unadj_se_cv
  ci_r = diff_unadj + qnorm(1-alpha/2)*diff_unadj_se_cv
  coverage_diff_unadj_cv = ((rep(1,nreps)%*%t(c(b%*%diff_out$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out$a)))<=ci_r)
  coverage_diff_unadj_cv = apply(coverage_diff_unadj_cv,2,mean)

  ci_l = diff_unadj - qnorm(1-alpha/2)*diff_unadj_se_hc3
  ci_r = diff_unadj + qnorm(1-alpha/2)*diff_unadj_se_hc3
  coverage_diff_unadj_hc3 = ((rep(1,nreps)%*%t(c(b%*%diff_out$a)))>=ci_l) & ((rep(1,nreps)%*%t(c(b%*%diff_out$a)))<=ci_r)
  coverage_diff_unadj_hc3 = apply(coverage_diff_unadj_hc3,2,mean)

  power_adj = apply(diff_adj_p<=alpha,2,mean)
  power_unadj = apply(diff_unadj_p<=alpha,2,mean)
  power_unadj_cv = apply(diff_unadj_p_cv<=alpha,2,mean)
  power_unadj_hc3 = apply(diff_unadj_p_hc3<=alpha,2,mean)


  return(list(b,
              coverage_diff_adj=coverage_diff_adj,
              coverage_diff_unadj = coverage_diff_unadj,
              coverage_diff_unadj_cv = coverage_diff_unadj_cv,
              coverage_diff_unadj_hc3 = coverage_diff_unadj_hc3,
              beta_hat_cov_adj=beta_hat_cov_adj,
              beta_hat_cov_unadj = beta_hat_cov_unadj,
              beta_hat_cov_unadj_cv = beta_hat_cov_unadj_cv,
              beta_hat_cov_unadj_hc3 = beta_hat_cov_unadj_hc3,
              power_adj=power_adj,
              power_unadj=power_unadj,
              power_unadj_cv=power_unadj_cv,
              power_unadj_hc3=power_unadj_hc3,
              diff_adj_p = diff_adj_p,
              diff_unadj_p=diff_unadj_p,
              diff_unadj_p_cv=diff_unadj_p_cv,
              diff_unadj_p_hc3=diff_unadj_p_hc3,
              diff_adj=diff_adj,
              diff_unadj=diff_unadj,
              diff_adj_se=diff_adj_se,
              diff_unadj_se=diff_unadj_se,
              diff_unadj_se_cv=diff_unadj_se_cv,
              diff_unadj_se_hc3=diff_unadj_se_hc3,
              data_sparsity = data_sparsity,

              mean_est_adj=mean_est_adj,
              mse_adj=mse_adj,
              se_est_adj=se_est_adj,

              mean_est_unadj=mean_est_unadj,
              mse_unadj=mse_unadj,
              se_est_unadj=se_est_unadj,

              coverage_adj=coverage_adj,
              coverage_unadj=coverage_unadj,
              coverage_unadj_cv=coverage_unadj_cv,
              coverage_unadj_hc3=coverage_unadj_hc3,

              se_adj = se_adj,
              se_unadj = se_unadj,
              se_unadj_cv = se_unadj_cv,
              se_unadj_hc3 = se_unadj_hc3,
              est_adj = est_adj,
              est_unadj = est_unadj,
              sigma2_hat=sigma2_hat,
              simu_param = list(ref=ref,s=s,Ng=Ng,b=b,ref_type=ref_type,tau2known=tau2known,
                                sc_lib_size = sc_lib_size,ref_bulklib_size=ref_bulklib_size,
                                bulk_lib_size = bulk_lib_size,nk=nk,
                                nreps=nreps,alpha=alpha,correction=correction)))






}


############## test ################
set.seed(12345)
G = 1000
K = 6
b1 = c(0.05,0.05,0.05,0.05,0.1,0.7)
b2 = c(0.05,0.06,0.08,0.1,0.2,0.51)
n_bulk = 10
b = cbind(t(rdirichlet(n_bulk/2,b1*K)),t(rdirichlet(n_bulk/2,b2*K)))
library(gtools)
ref = t(rdirichlet(K,rep(1,G)))
test = simu_study(ref,G,b,
                  ref_type='multi_sc',
                  nreps = 1,
                  sc_lib_size = 0.1,
                  printevery = 1,same_indi = F,mean_to_var_sigma = 1/3,mean_to_var_tau = 1/3,tau2 = NULL,sigma2 = NULL,
                  tau2known=F,sigma2known = F,est_pop_var = T,correction = T,
                  weight = 'equal',hc.type = 'hc3',n_indi = 10,meta_mode = 'smooth')
