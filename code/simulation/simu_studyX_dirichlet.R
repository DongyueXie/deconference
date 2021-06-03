
#'@description generate X usign dirichlet distribution

simu_studyX_dirichlet = function(ref,b,
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
                       filter.gene = FALSE){

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
    dmat = scRef_multi_proc(data.obj$Y,data.obj$cell_type_idx,data.obj$indi_idx,
                            estimator=x_estimator,est_sigma2 = est_pop_var,
                            meta_var=meta_var,
                            meta_mode=meta_mode)



    fit.adj.hc0 = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v=TRUE,hc.type='hc0')
    fit.adj.hc3 = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v=TRUE,hc.type='hc3')
    fit.unadj.hc0 = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v  =FALSE,hc.type='hc0')
    fit.unadj.hc3 = wols(dmat$X,data.obj$y,w.mode=w.mode,X_array = dmat$X_array,Vg = dmat$Vg,adj.v=FALSE,hc.type='hc3')



  }

  return(list(se_adj_hc0 = fit.adj.hc0$beta_hat_se,
              se_adj_hc3 = fit.adj.hc3$beta_hat_se,
              se_unadj_hc0 = fit.unadj.hc0$beta_hat_se,
              se_unadj_hc3 = fit.unadj.hc3$beta_hat_se,
              est_adj = fit.adj.hc0$beta_hat,
              est_unadj = fit.unadj.hc0$beta_hat))

}


