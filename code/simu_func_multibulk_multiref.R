########################################
############ simulation function for multi indi sc ref and multi bulk deconv
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
#'@param sc_lib_size_UMI single cell ref library size for UMI, sc_lib_size*G is sequencing depth; a vector gives different ls for different UMI datasets
#'@param sc_lib_size_nonUMI single cell ref library size for nonUMI; a vector gives different ls for different nonUMI datasets
#'@param n_ref_UMI number of reference dataset from UMI
#'@param n_ref_nonUMI number of reference dataset from UMI
#'@param expr_dist the distribution of gene relative expression, to draw from, either gamma or log-normal
#'@param nk number of single cells of each cell type
#'@param ref_scale a scaling scalar times gene relative expression
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


simu_study = function(ref,b,
                      hc.type = 'hc3',
                      marker_gene = NULL,
                      weight = 'equal',
                      bulk_lib_size = 100,
                      sc_lib_size_UMI = 0.1,
                      sc_lib_size_nonUMI = 20,
                      n_ref_UMI = 1,
                      n_ref_nonUMI = 1,
                      gene_length = NULL,
                      nk = 50,
                      #ref_scale = 1,
                      same_indi = FALSE,
                      mean_to_var_sigma = 10,
                      mean_to_var_tau_UMI = 10,
                      mean_to_var_tau_nonUMI = 1,
                      zero_inflated_nonUMI = TRUE,
                      tau2_UMI = matrix(runif(nrow(ref)*length(b),0.5/(nrow(ref))^2,0.8/(nrow(ref))^2),nrow=nrow(ref)),
                      tau2_nonUMI = matrix(runif(nrow(ref)*length(b),0.8/(nrow(ref))^2,1/(nrow(ref))^2),nrow=nrow(ref)),
                      tau2known=FALSE,
                      sigma2 = matrix(runif(nrow(ref)*length(b),0.5/(nrow(ref))^3,1/(nrow(ref))^3),nrow=nrow(ref)),
                      sigma2known=FALSE,
                      full_sigma2 = FALSE,
                      expr_dist = 'log-normal',
                      x_estimator = 'separate',
                      cellsize_est = 'glm',
                      nreps=100,alpha=0.05,
                      #a=0,
                      correction=FALSE,
                      s=NULL,
                      #s_type = 'unknown',
                      printevery=10,
                      n_indi = 8,
                      est_pop_var = TRUE,
                      meta_var='adjust',
                      meta_mode = 'local',
                      eps=0,
                      Ng=nrow(ref),
                      groups = c(rep(1,ncol(b)/2),rep(2,ncol(b)/2))){

  G = nrow(ref)
  K = ncol(ref)
  n_bulk = ncol(b)
  if(is.null(gene_length)){
    gene_length = rep(1,G)
  }else{
    gene_length = gene_length/sum(gene_length)*G
  }


  ref = apply(ref,2,function(z){z/sum(z)})*G
  if(is.null(rownames(ref))){
    rownames(ref) = 1:nrow(ref)
  }


  if(is.null(s)){
    s = rep(1,K)
  }

  if(length(nk)==1){
    nk = rep(nk,K)
  }

  b = apply(b,2,function(z){z/sum(z)})
  p = apply(b,2,function(z){(z*s)/sum(z*s)})

  if(length(tau2_UMI)==1){
    tau2_UMI = matrix(tau2_UMI,nrow=Ng,ncol=K)
  }
  if(length(tau2_nonUMI)==1){
    tau2_nonUMI = matrix(tau2_nonUMI,nrow=Ng,ncol=K)
  }

  if(length(sigma2)==1){
    sigma2 = matrix(sigma2,nrow=Ng,ncol=K)
  }

  # if(!full_sigma2){
  #
  #   # if using log normal distribution, pre-calculate mu and var needed, for generating individual X
  #   if(expr_dist=='log-normal'){
  #     if(!is.null(mean_to_var_sigma)){
  #       mu = log(ref^2/sqrt(ref^2+ref/mean_to_var_sigma))
  #       v = log(ref/mean_to_var_sigma/ref^2+1)
  #
  #     }else{
  #       mu = log(ref^2/sqrt(ref^2+sigma2))
  #       v = log(sigma2/ref^2+1)
  #
  #     }
  #   }
  #
  #   if(expr_dist=='gamma'){
  #     if(!is.null(mean_to_var_sigma)){
  #       shape = ref*mean_to_var_sigma
  #       rate=mean_to_var_sigma
  #
  #     }else{
  #       shape = ref^2/sigma2
  #       rate=ref/sigma2
  #     }
  #   }
  #
  # }

  if(weight == 'equal'){
    w = rep(1,Ng)
  }else{
    w = NULL
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


  data_sparsity = matrix(nrow=nreps,ncol=(n_ref_UMI+n_ref_nonUMI))
  X_hat = list()
  S_hat = matrix(nrow=nreps,ncol=K)
  y_all = list()
  beta_hat_cov_adj = list()
  beta_hat_cov_unadj = list()
  beta_hat_cov_unadj_cv = list()
  beta_hat_cov_unadj_hc3 = list()

  cell_size_est = matrix(nrow=nreps,ncol=K)

  sigma2_hat = list()# unly use when ref is multi-sc and est_sigma2 IS TRUE

  for(reps in 1:nreps){


    if(reps%%printevery==0){print(sprintf("running %d (out of %d)",reps,nreps))}

    ## generate simulated dataset
    #print(s)
    simu_data = generate_simu_data(ref=ref,b=b,
                                              bulk_lib_size = bulk_lib_size,
                                              sc_lib_size_UMI = sc_lib_size_UMI,
                                              sc_lib_size_nonUMI = sc_lib_size_nonUMI,
                                              n_ref_UMI = n_ref_UMI,
                                              n_ref_nonUMI = n_ref_nonUMI,
                                              gene_length = gene_length,
                                              nk = nk,
                                              same_indi = same_indi,
                                              mean_to_var_sigma = mean_to_var_sigma,
                                              mean_to_var_tau_UMI = mean_to_var_tau_UMI,
                                              mean_to_var_tau_nonUMI = mean_to_var_tau_nonUMI,
                                              zero_inflated_nonUMI = zero_inflated_nonUMI,
                                              tau2_UMI = tau2_UMI,
                                              tau2_nonUMI = tau2_nonUMI,
                                              sigma2 = sigma2,
                                              full_sigma2 = full_sigma2,
                                              expr_dist = expr_dist,
                                              s=s,
                                              n_indi = n_indi,
                                              Ng=Ng,
                                              groups = groups)




    #print('running regression')
    fit_adj = deconference_multi_ref(simu_data$refs,simu_data$bulks,cell_types=1:K,sigma2=NULL,tau2=NULL,
                                                est_sigma2=est_pop_var,meta_var=meta_var,meta_mode=meta_mode,
                                     gene_length_adjust = TRUE,
                                     gene_length=gene_length,protocol = c(rep('UMI',n_ref_UMI),rep('nonUMI',n_ref_nonUMI)),
                                                correction=correction,cellsize_est=cellsize_est,
                                                #marker_gene = NULL,
                                                hc.type = hc.type,w = w)

    #browser()
    cell_size_est[reps,] = fit_adj$input$S
    sigma2_hat[[reps]] = fit_adj$input$Sigma



    #print('running unadj regression')
    #browser()
    fit_unadj = unadjusted_lm(fit_adj$input$y,fit_adj$input$X%*%diag(fit_adj$input$S),w=fit_adj$input$w,groups)

    # if(reps==27){
    #   browser()
    # }


    #print('record sparsity')

    data_sparsity[reps,] = unlist(lapply(simu_data$refs, function(Y){sum(counts(Y)==0)/prod(dim(counts(Y)))}))
    X_hat[[reps]] = fit_adj$input$X
    S_hat[reps,] = fit_adj$input$S
    y_all[[reps]] = fit_adj$input$y


    #print('recording regression estimations')

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


  return(list(b=b,
              X_hat=X_hat,
              S_hat=S_hat,
              y_all=y_all,
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
              cell_size_est=cell_size_est,

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
              simu_param = list(ref=ref,s=s,Ng=Ng,b=b,bulk_lib_size = bulk_lib_size,nk=nk,
                                nreps=nreps,alpha=alpha,correction=correction,gene_length=gene_length)))

}


# m is a matrix, each row is a mean
# V is a covariance matrix
draw_mvlnorm = function(M,V){
  n = nrow(M)
  p = ncol(M)
  X = matrix(nrow=n,ncol=p)

  for(i in 1:n){

    mu = log(M[i,]^2/sqrt(diag(V)+M[i,]^2))
    ccov = log(1+V/(M[i,]%*%t(M[i,])))
    temp = try(MASS::mvrnorm(1,mu,ccov),silent = TRUE)
    if(class(temp)=='try-error'){
      ccov = ccov + (-min(eigen(ccov)$values)+1e-8)*diag(p)
      temp = MASS::mvrnorm(1,mu,ccov)
    }
    X[i,] = temp
  }
  exp(X)
}



generate_simu_data = function(ref,b,
                              bulk_lib_size = 100,
                              sc_lib_size_UMI = 0.1,
                              sc_lib_size_nonUMI = 20,
                              n_ref_UMI = 1,
                              n_ref_nonUMI = 1,
                              gene_length = NULL,
                              nk = 50,
                              same_indi = FALSE,
                              mean_to_var_sigma = 10,
                              mean_to_var_tau_UMI = 10,
                              mean_to_var_tau_nonUMI = 1,
                              zero_inflated_nonUMI = TRUE,
                              tau2_UMI = matrix(runif(nrow(ref)*length(b),0.5/(nrow(ref))^2,0.8/(nrow(ref))^2),nrow=nrow(ref)),
                              tau2_nonUMI = matrix(runif(nrow(ref)*length(b),0.8/(nrow(ref))^2,1/(nrow(ref))^2),nrow=nrow(ref)),
                              sigma2 = matrix(runif(nrow(ref)*length(b),0.5/(nrow(ref))^3,1/(nrow(ref))^3),nrow=nrow(ref)),
                              full_sigma2 = FALSE,
                              expr_dist = 'log-normal',
                              s=NULL,
                              n_indi = 8,
                              Ng=nrow(ref),
                              groups = c(rep(1,ncol(b)/2),rep(2,ncol(b)/2))){

  G = nrow(ref)
  K = ncol(ref)
  n_bulk = ncol(b)
  if(is.null(gene_length)){
    gene_length = rep(1,G)
  }else{
    gene_length = gene_length/sum(gene_length)*G
  }


  ref = apply(ref,2,function(z){z/sum(z)})*G
  if(is.null(rownames(ref))){
    rownames(ref) = 1:nrow(ref)
  }


  if(is.null(s)){
    s = rep(1,K)
  }

  if(length(nk)==1){
    nk = rep(nk,K)
  }

  b = apply(b,2,function(z){z/sum(z)})
  p = apply(b,2,function(z){(z*s)/sum(z*s)})

  if(length(tau2_UMI)==1){
    tau2_UMI = matrix(tau2_UMI,nrow=Ng,ncol=K)
  }
  if(length(tau2_nonUMI)==1){
    tau2_nonUMI = matrix(tau2_nonUMI,nrow=Ng,ncol=K)
  }

  if(length(sigma2)==1){
    sigma2 = matrix(sigma2,nrow=Ng,ncol=K)
  }

  if(!full_sigma2){

    # if using log normal distribution, pre-calculate mu and var needed, for generating individual X
    if(expr_dist=='log-normal'){
      if(!is.null(mean_to_var_sigma)){
        mu = log(ref^2/sqrt(ref^2+ref/mean_to_var_sigma))
        v = log(ref/mean_to_var_sigma/ref^2+1)

      }else{
        mu = log(ref^2/sqrt(ref^2+sigma2))
        v = log(sigma2/ref^2+1)

      }
    }

    if(expr_dist=='gamma'){
      if(!is.null(mean_to_var_sigma)){
        shape = ref*mean_to_var_sigma
        rate=mean_to_var_sigma

      }else{
        shape = ref^2/sigma2
        rate=ref/sigma2
      }
    }

  }

  gene_names = rownames(ref)
  names(gene_length) = gene_names

  #Theta = apply(ref_rep,2,function(z){z/sum(z)})*ref_scale

  if(same_indi){

    # bulk data gene relative expression.


    mb = (ref*gene_length)%*%diag(s)%*%b

  }else{

    # bulk data gene relative expression.
    ## generate theta_i for bulk individual
    mb = matrix(nrow=Ng,ncol=n_bulk)
    for(i in 1:n_bulk){

      #browser()

      if(!full_sigma2){
        if(expr_dist=='gamma'){
          X_b = matrix(rgamma(Ng*K,shape,rate),ncol=K)
        }

        if(expr_dist=='log-normal'){
          X_b = exp(matrix(rnorm(Ng*K,mu,sqrt(v)),ncol=K))
        }
      }else{
        X_b = draw_mvlnorm(ref,sigma2)
      }

      #browser()
      mb[,i] = (X_b*gene_length)%*%diag(s)%*%b[,i]

    }
  }

  thetab = apply(mb,2,function(z){z/sum(z)})

  y = matrix(rpois(Ng*n_bulk,bulk_lib_size*Ng*thetab),nrow=Ng)
  rownames(y) = gene_names

  bulks = SingleCellExperiment(assays = list(counts = y),
                               colData = DataFrame(individual = 1:n_bulk,groups=groups))




  #browser()


  ## generate multi ref data


  ### generate UMI data



  refs = list()
  for(u in 1:(n_ref_UMI+n_ref_nonUMI)){

    #print(paste('UMI:',u))

    Y = matrix(nrow=Ng,ncol=sum(nk)*n_indi)
    #browser()
    indi_idx = rep(1:n_indi,each = sum(nk))
    cell_type = c()
    Cr = c()
    if(u<=n_ref_UMI){
      sc_lib_size = sc_lib_size_UMI[u]
      mean_to_var_tau = mean_to_var_tau_UMI
      tau2 = tau2_UMI

    }else{
      sc_lib_size = sc_lib_size_nonUMI[u-n_ref_UMI]
      mean_to_var_tau = mean_to_var_tau_nonUMI
      tau2 = tau2_nonUMI
    }
    for(i in 1:n_indi){
      #print(paste('indi',i))
      indi_cell_idx = which(indi_idx==i)
      indi_celltype_idx = rep(1:K,times = nk)

      if(!full_sigma2){
        if(expr_dist=='gamma'){
          X_i = matrix(rgamma(Ng*K,shape,rate),ncol=K)
        }

        if(expr_dist=='log-normal'){
          X_i = exp(matrix(rnorm(Ng*K,mu,sqrt(v)),ncol=K))
        }
      }else{
        X_i = draw_mvlnorm(ref,sigma2)
      }

      Y_i = matrix(nrow=Ng,ncol=sum(nk))
      Cr_i = c()
      for(k in 1:K){
        cc = which(indi_celltype_idx==k)
        ag = X_i[,k]%*%t(rep(1,nk[k]))

        Cr_i[cc] = (rnbinom(nk[k],sc_lib_size*Ng*s[k],0.5)+1)/Ng
        #Cr_i[cc]  = rep(sc_lib_size*Ng,nk[k])
        #browser()

        ## if UMI or nonUMI+not zero inflated -- same
        ## if nonUMI and zero inflated -- zi


        ### z_mat indicates 0 or not
        z_mat = matrix(1,nrow=Ng,ncol=nk[k])
        if(u>n_ref_UMI&zero_inflated_nonUMI){
          tempx = c(X_i[,k])
          tempx = (tempx-mean(tempx))/sd(tempx)
          pi0 = exp(tempx)/(1+exp(tempx))
          ag = ag/pi0
          z_mat = apply(z_mat,2,function(z){rbinom(Ng,1,pi0)})
        }

        if(expr_dist=='gamma'){
          if(!is.null(mean_to_var_tau)){
            Y_i[,cc] = matrix(rgamma(Ng*nk[k],ag*mean_to_var_tau,rate=mean_to_var_tau),ncol=nk[k])
          }else{
            atau2 = tau2[,k]%*%t(rep(1,nk[k]))
            Y_i[,cc] = matrix(rgamma(Ng*nk[k],ag^2/atau2,rate=ag/atau2),ncol=nk[k])
          }
        }

        if(expr_dist=='log-normal'){
          if(!is.null(mean_to_var_tau)){
            Y_i[,cc] = exp(matrix(rnorm(Ng*nk[k],log(ag^2/sqrt(ag^2+ag/mean_to_var_tau)),sqrt(log(ag/mean_to_var_tau/ag^2+1))),ncol=nk[k]))
          }else{
            atau2 = tau2[,k]%*%t(rep(1,nk))
            Y_i[,cc] = exp(matrix(rnorm(Ng*nk[k],log(ag^2/sqrt(ag^2+atau2)),sqrt(log(atau2/ag^2+1))),ncol=nk[k]))
          }
        }

        #browser()

        Y_i[,cc] = Y_i[,cc]*z_mat




      }
      cell_type = c(cell_type,indi_celltype_idx)
      Y[,indi_cell_idx] = Y_i
      Cr = c(Cr,Cr_i)
    }

    #Y = matrix(rnorm(Ng*n_indi*K*nk,Y,sqrt(sc_noise_var)),ncol=n_indi*K*nk)
    #y = rnorm(Ng,Xb,sqrt(bulk_noise_var))


    #print('generating Y')

    #browser()

    if(u>n_ref_UMI){
      #s_all = s[match(cell_type,1:K)]
      #Y = t(t(Y*gene_length)*s_all)
      Y = Y*gene_length
    }

    #browser()

    #Y = apply(Y,2,function(z){z/sum(z)})
    Y = matrix(rpois(Ng*sum(nk)*n_indi,t(t(Y)*Cr)),ncol=sum(nk)*n_indi)

    rownames(Y) = gene_names

    refs[[u]] = SingleCellExperiment(assays = list(counts = Y),
                                     colData = DataFrame(cell_type = cell_type,individual = indi_idx))

    #print('done creating  se obj')

  }


  return(list(refs=refs,bulks=bulks,simu_param = list(ref=ref,s=s,Ng=Ng,b=b,
                                                      bulk_lib_size = bulk_lib_size,nk=nk,gene_length=gene_length)))



}


############# test ################
# set.seed(1234)
# library(gtools)
# G = 300
# #K = 6
# #b1 = c(0.05,0.05,0.05,0.05,0.1,0.7)
# #b2 = c(0.05,0.06,0.08,0.1,0.2,0.51)
# #n_bulk = 6
# K=4
# b1 = 1:4/10
# b2=b1
# b = cbind(b1,b2)
# ref = t(rdirichlet(K,rep(1,G)))
# gene_length=readRDS('data/gene_length.rds')
# s = c(1,1.5,1.3,0.7)
# source('code/deconference_main.R')
# mean_to_var_sigma = 100
# sigma2 = matrix(runif(K^2,1/mean_to_var_sigma/5,1/mean_to_var_sigma/2),ncol=K)
# sigma2 = (sigma2+t(sigma2))/2
# diag(sigma2) = 1/mean_to_var_sigma
# test = simu_study(ref,b,nreps = 10,printevery = 1,gene_length = sample(gene_length,G),bulk_lib_size = 300,s=s,
#                   expr_dist = 'log-normal',sc_lib_size_UMI = c(0.1,0.2),sc_lib_size_nonUMI = c(10,20),zero_inflated_nonUMI = T,full_sigma2 = TRUE,sigma2=sigma2,
#                   nk=100,mean_to_var_sigma = 300,mean_to_var_tau_UMI = 100,mean_to_var_tau_nonUMI = 100,n_ref_UMI = 1,n_ref_nonUMI = 0)
#
# print(test$coverage_adj)
# print(test$coverage_unadj_hc3)
#
# print(test$mean_est_adj)
# print(test$mean_est_unadj)
#
# print(test$se_est_adj)
# print(test$se_est_unadj)
#
# print(mean(test$mse_adj))
# print(mean(test$mse_unadj))
