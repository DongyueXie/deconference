


simu_neuron = function(indis_ref,
                       ref.idx,
                       b,
                       cor.idx,
                       calc_cov = T,
                       bulk_lib_size = 500,
                       groups = c(rep(1,ncol(b)/2),rep(2,ncol(b)/2)),
                       centeringXY=FALSE,
                       only.scale.pos.res = FALSE,
                       only.add.pos.res = FALSE,
                       verbose = F,
                       weighted = FALSE,
                       w = NULL,
                       eb.V = FALSE){

  G = dim(indis_ref)[1]
  K = dim(indis_ref)[2]
  n_sub = dim(indis_ref)[3]
  n_ref = length(ref.idx)
  n_bulk = n_sub-n_ref
  gene_names = dimnames(indis_ref)[[1]]

  X_array_ref = indis_ref[,,ref.idx]
  X_array_bulk = indis_ref[,,-ref.idx]

  X = apply(X_array_ref,c(1,2),mean,na.rm=TRUE)
  #browser()
  if(eb.V){
    #browser()
    V.diag = t(apply(X_array_ref,c(1),function(z){diag(cov(t(z),use = 'complete.obs'))}))
    V.diag.sd = apply(V.diag,2,function(z){vashr::vash(sqrt(z),df=n_ref-1)$sd.post})
    V = matrix(nrow=G,ncol=K^2)
    for (g in 1:G){
      V[g,] = c(t(cor(t(X_array_ref[g,,]),use = 'complete.obs')*V.diag.sd[g,])*V.diag.sd[g,]/n_ref)
    }
    V[is.na(V)] = 0
  }else{
    V = t(apply(X_array_ref,c(1),function(z){(cov(t(z),use = 'complete.obs'))}))/n_ref
  }


  if(weighted){
    if(is.null(w)){
      # calc weights for each cell type, then average,
      if(eb.V){
        w = 1/(rowSums(V))
      }else{
        V.temp = t(apply(X_array_ref,c(1),function(z){(cov(t(z),use = 'complete.obs'))}))
        fit.vash = vashr::vash(sqrt(rowSums(V.temp)),df=n_ref-1)
        w = 1/(fit.vash$sd.post)^2
      }
      w = w/sum(w)*G
    }
  }else{
    w = 1
  }

  # create bulk data

  mb = lapply(1:n_bulk,function(i){X_array_bulk[,,i]%*%b[,i]})
  mb = do.call(cbind,mb)
  thetab = apply(mb,2,function(z){z/sum(z)})

  y = matrix(rpois(G*n_bulk,bulk_lib_size*G*thetab),nrow=G)
  rownames(y) = gene_names

  # ols

  fit.ols = unadjusted_lm(y,X,w=w,groups = groups)

  # adjust for measurement error, not for correlation

  fit.err.hc0 = estimation_func2(y=y,X=X,Vg=V,
                                 w=w,hc.type='hc0',correction=FALSE,
                                 calc_cov=calc_cov,verbose=verbose,
                                 cor.idx=NULL,
                                 centeringXY=centeringXY,
                                 true.beta = NULL,
                                 only.scale.pos.res=only.scale.pos.res,
                                 only.add.pos.res=only.add.pos.res)

  fit.err.hc3 = estimation_func2(y=y,X=X,Vg=V,
                                 w=w,hc.type='hc3',correction=FALSE,
                                 calc_cov=calc_cov,verbose=verbose,
                                 cor.idx=NULL,
                                 centeringXY=centeringXY,
                                 true.beta = NULL,
                                 only.scale.pos.res=only.scale.pos.res,
                                 only.add.pos.res=only.add.pos.res)

  # adjust both measurement error and correlation

  fit.err.cor.hc0 = estimation_func2(y=y,X=X,Vg=V,
                                 w=w,hc.type='hc0',correction=FALSE,
                                 calc_cov=calc_cov,verbose=verbose,
                                 cor.idx=cor.idx,
                                 centeringXY=centeringXY,
                                 true.beta = NULL,
                                 only.scale.pos.res=only.scale.pos.res,
                                 only.add.pos.res=only.add.pos.res)

  fit.err.cor.hc3 = estimation_func2(y=y,X=X,Vg=V,
                                 w=w,hc.type='hc3',correction=FALSE,
                                 calc_cov=calc_cov,verbose=verbose,
                                 cor.idx=cor.idx,
                                 centeringXY=centeringXY,
                                 true.beta = NULL,
                                 only.scale.pos.res=only.scale.pos.res,
                                 only.add.pos.res=only.add.pos.res)

  out = list(fit.ols=fit.ols,
             fit.err.hc0=fit.err.hc0,
             fit.err.hc3=fit.err.hc3,
             fit.err.cor.hc0=fit.err.cor.hc0,
             fit.err.cor.hc3=fit.err.cor.hc3,
             w=w,
             input = list(ref.idx=ref.idx,b=b))
  return(out)

}



simu_neuron_music = function(indis_ref,
                             ref.idx,
                             b,
                             bulk_lib_size = 500){

  devtools::load_all('D://githubs/MuSiC')

  G = dim(indis_ref)[1]
  K = dim(indis_ref)[2]
  n_sub = dim(indis_ref)[3]
  n_ref = length(ref.idx)
  n_bulk = n_sub-n_ref
  gene_names = dimnames(indis_ref)[[1]]

  X_array_ref = indis_ref[,,ref.idx]
  X_array_bulk = indis_ref[,,-ref.idx]

  X = apply(X_array_ref,c(1,2),mean,na.rm=TRUE)
  Sigma = t(apply(X_array_ref,c(1),function(z){diag(cov(t(z),use = 'complete.obs'))}))

  # create bulk data

  mb = lapply(1:n_bulk,function(i){X_array_bulk[,,i]%*%b[,i]})
  mb = do.call(cbind,mb)
  thetab = apply(mb,2,function(z){z/sum(z)})

  y = matrix(rpois(G*n_bulk,bulk_lib_size*G*thetab),nrow=G)
  rownames(y) = gene_names


  Est.prop.allgene = NULL
  for(nb in 1:n_bulk){
    fit.music = music.basic(y[,nb],X,S=1,Sigma,iter.max=1000,nu=1e-4,eps=0.01)
    Est.prop.allgene = cbind(Est.prop.allgene, fit.music$p.nnls)
  }

  Est.prop.allgene

}
