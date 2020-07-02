bootstrap_deconv = function(y,Y,cell_type,estimator='separate',nboot=1000,adjust=TRUE,ref_type = "bulk"){
  y = cbind(y)
  if(ref_type=='bulk'){
    # reference samples are bulk data
    out = bootstrap_bulk_proc(y,Y,nboot,adjust)
  }

  if(ref_type=='sc'){
    # reference samples are bulk data
    out = bootstrap_sc_proc(y,Y,cell_type,estimator='separate',nboot,adjust)
  }
  return(out)
}

bootstrap_sc_proc = function(y,Y,cell_type,estimator='separate',nboot,adjust){

  beta_hat = c(deconference(y,Y,cell_type=cell_type,ref_type = "sc",x_estimator=estimator,adjust=adjust)$beta_hat)

  beta_hats = matrix(nrow=nboot,ncol=K*ncol(y))

  for(i in 1:nboot){
    g_idx = sample(1:G,G,replace=TRUE)
    y_i = y[g_idx,]
    Y_i = Y[g_idx,]
    out = deconference(y_i,Y_i,cell_type=cell_type,ref_type = "sc",x_estimator=estimator,adjust=adjust)
    beta_hats[i,] = c(out$beta_hat)
  }

  return(list(beta_hat=beta_hat,beta_hats=beta_hats))

}


bootstrap_bulk_proc = function(y,Y,nboot=1000,adjust=TRUE){

  G = nrow(Y)
  K = ncol(Y)

  beta_hat = c(deconference(y,Y,ref_type = "bulk",adjust=adjust)$beta_hat)

  beta_hats = matrix(nrow=nboot,ncol=K*ncol(y))


  for(i in 1:nboot){
    g_idx = sample(1:G,G,replace=TRUE)
    y_i = y[g_idx,]
    Y_i = Y[g_idx,]
    out = deconference(y_i,Y_i,ref_type = "bulk",adjust=adjust)
    beta_hats[i,] = c(out$beta_hat)
  }

  #browser()

  return(list(beta_hat=beta_hat,beta_hats=beta_hats))

}

get_ci = function(out,alpha=0.05){
  ci = apply(out$beta_hats,2,function(z){c(quantile(z,alpha/2),quantile(z,1-alpha/2))})
  ci
}

get_cover = function(ci,p){
  1*((p > ci[1,]) & (p < ci[2,]))
}

get_ci_diff = function(out,K,alpha=0.05){
  beta_diff = out$beta_hats[,1:K] - out$beta_hats[,(K+1):(2*K)]
  apply(beta_diff,2,function(z){c(quantile(z,alpha/2),quantile(z,1-alpha/2))})
}


simu_study_boot = function(ref,Ng,b,ref_type='bulk',bulk_lib_size = 50,
                           ref_lib_size = 30,
                           nreps=100,alpha=0.05,a=0,
                           correction=FALSE,s,printevery=10,
                           tau2 = matrix(runif(Ng*length(b),1/(Ng*10),1/(Ng*10)),nrow=Ng),
                           b2=NULL,nboot=1000){

  G = nrow(ref)
  K = ncol(ref)

  ref = apply(ref,2,function(z){z/sum(z)})


  if(missing(s)){
    s = rep(1,K)
  }

  b = b/sum(b)
  p = (b*s)/sum(b*s)



    b2 = b2/sum(b2)
    p2 = (b2*s)/sum(b2*s)



    est_adj = matrix(nrow=nreps,ncol=2*K)
    coverage_adj = matrix(nrow=nreps,ncol=2*K)
    coverage_diff_adj = matrix(nrow=nreps,ncol=K)
    #se_adj = matrix(nrow=nreps,ncol=2*K)
    #diff_adj = matrix(nrow=nreps,ncol=K)
    #diff_adj_se = matrix(nrow=nreps,ncol=K)
    #diff_adj_p = matrix(nrow = nreps,ncol=K)

    est_unadj = matrix(nrow=nreps,ncol=2*K)
    coverage_unadj = matrix(nrow=nreps,ncol=2*K)
    coverage_diff_unadj = matrix(nrow=nreps,ncol=K)
    #se_unadj = matrix(nrow=nreps,ncol=2*K)
    #diff_unadj = matrix(nrow=nreps,ncol=K)
    #diff_unadj_se = matrix(nrow=nreps,ncol=K)
    #diff_unadj_p = matrix(nrow = nreps,ncol=K)


  for(reps in 1:nreps){


    if(reps%%printevery==0){print(sprintf("done %d (out of %d)",reps,nreps))}




    #Obtain X
    ref_rep = ref[sample(1:G,Ng),]
    Theta = apply(ref_rep,2,function(z){z/sum(z)})

    # bulk data gene relative expression.
    mb = Theta%*%diag(s*(Ng/G))%*%b
    thetab = mb/sum(mb)


      mb2 = Theta%*%diag(s*(Ng/G))%*%b2
      thetab2 = mb2/sum(mb2)


      if(ref_type=='bulk'){


        Cr = rpois(K,ref_lib_size*Ng)+1
        #Cr = rep(ref_lib_size,K)
        U = diag(Cr)
        Y = matrix(rpois(Ng*K,Theta%*%U),ncol=K)
        #bulk data
        y = rpois(Ng,rpois(1,bulk_lib_size)*Ng*thetab)

          y2 = rpois(Ng,rpois(1,bulk_lib_size)*Ng*thetab2)
          y = cbind(y,y2)


        #t2=Sys.time()

        #print(paste('generate data takes',t2-t1))

        fit_adj = bootstrap_deconv(y,Y,nboot=nboot,adjust=TRUE,ref_type = ref_type)
        fit_unadj = bootstrap_deconv(y,Y,nboot=nboot,adjust=FALSE,ref_type = ref_type)

        #browser()

        est_adj[reps,] = fit_adj$beta_hat
        coverage_adj[reps,] = get_cover(get_ci(fit_adj,alpha=alpha),c(p,p2))
        coverage_diff_adj[reps,] = get_cover(get_ci_diff(fit_adj,K=K,alpha=alpha),(p-p2))

        est_unadj[reps,] = fit_unadj$beta_hat
        coverage_unadj[reps,] = get_cover(get_ci(fit_unadj,alpha=alpha),c(p,p2))
        coverage_diff_unadj[reps,] = get_cover(get_ci_diff(fit_unadj,K=K,alpha=alpha),(p-p2))

        #t3 = Sys.time()
        #print(paste('fit model takes',t3-t2))
      }



      if(ref_type=='sc'){

        #t1=Sys.time()


        cell_type = rep(1:K,each=nk)
        Cr = rnbinom(nk*K,sc_lib_size*Ng,0.5)+1
        Y = matrix(nrow=Ng,ncol=nk*K)
        #tau2 = matrix(nrow = Ng,ncol = K)
        for(k in 1:K){
          cell_idx = which(cell_type==k)
          #Y[,cell_idx] = t(rdirichlet(nk,Theta[,k]*Ng))
          # Y[,cell_idx] = matrix(rgamma(Ng*nk,Theta[,k]%*%t(rep(1,nk))*snr,rate=1),ncol=nk)/snr
          ag = Theta[,k]%*%t(rep(1,nk))
          atau2 = tau2[,k]%*%t(rep(1,nk))
          Y[,cell_idx] = matrix(rgamma(Ng*nk,ag^2/atau2,
                                       rate=ag/atau2),ncol=nk)
          # tau2[,k] = Theta[,k]*(1-Theta[,k]) / (Ng+1)
        }
        #Y = matrix(rpois(Ng*nk*K,Y%*%diag(Cr)),ncol=nk*K)
        Y = matrix(rpois(Ng*nk*K,t(t(Y)*Cr)),ncol=nk*K)

        y = rpois(Ng,rpois(1,bulk_lib_size)*Ng*thetab)

          y2 = rpois(Ng,rpois(1,bulk_lib_size)*Ng*thetab2)
          y = cbind(y,y2)

          fit_adj = bootstrap_deconv(y,Y,nboot=nboot,adjust=TRUE,ref_type = ref_type)
          fit_unadj = bootstrap_deconv(y,Y,nboot=nboot,adjust=FALSE,ref_type = ref_type)



      }


  }

    return(list(est_adj = est_adj,
           est_unadj = est_unadj,
           coverage_adj = apply(coverage_adj,2,mean),
           coverage_unadj = apply(coverage_unadj,2,mean),
           coverage_diff_adj = apply(coverage_diff_adj,2,mean),
           coverage_diff_unadj = apply(coverage_diff_unadj,2,mean)))




}



