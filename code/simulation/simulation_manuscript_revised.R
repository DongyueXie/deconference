
############## Run Simulation for manuscript ###############
source('code/deconference_main.R')
source('code/ols_hc3.R')
source('code/CIBERSORT.R')
devtools::load_all('D://githubs/MuSiC')
library(MEAD)
library(gtools)
library(Matrix)
#'@description Generate data: based on real Xin data, generate uncorrelated or correlated individual reference matrices.
#'Compare methods: ols, adjust for measurement error, adjust for both measurement error and correlation, adjust for both measurement error and correlation and add weights.
#'For simplicity, we only simulate data at individual level. - only generate X from U;

#'@param ref,sigma2 population gene expression matrix U, of dimension G by K; and variance matrix
#'@param p cell type proportion
#'@param n_bulk number of bulk data
#'@param dirichlet whether generate random p
#'@param n_bulk_for_cor number of bulk data generated for inferring correlations
simu_study = function(ref,
                      sigma2,
                      p,
                      n_bulk=100,
                      dirichlet = TRUE,
                      dirichlet.scale = 10,
                      d=500,
                      nreps = 100,
                      bulk_lib_size = 500,
                      n_ref = 10,
                      printevery=10,
                      n_bulk_for_cor = 100,
                      gamma_shape = 5,
                      bulk_bias_sd = 0.1,
                      nfold = 10
){


  genp = function(K){
    p = runif(K)
    p/sum(p)
  }

  G = nrow(ref)
  K = ncol(ref)

  is.indep = (d==0)
  if(d!=0){
    R = matrix(0,nrow=G,ncol=G)
    for(i in 1:G){
      for(j in i:min(i+d,G)){
        R[i,j] = max(1-abs(i-j)/d,0)
      }
    }
    R = R+t(R) - diag(G)
    R = Matrix(R,sparse = TRUE)
  }









  ########
  # rmse #
  ########

  # rmse:
  ## 1. ols
  ## 2. MEAD, NO WEIGHT
  ## 3. MEAD Marker
  ## 4. MEAD
  ## 5. MuSiC
  ## 6. CIBERSORT

  p_hat_ols = array(dim = c(K,n_bulk,nreps))
  p_hat = array(dim = c(K,n_bulk,nreps))
  p_hat_marker = array(dim = c(K,n_bulk,nreps))
  p_hat_weight = array(dim = c(K,n_bulk,nreps))
  p_hat_music = array(dim = c(K,n_bulk,nreps))
  p_hat_ciber = array(dim = c(K,n_bulk,nreps))

  ############
  # Variance #
  ############

  # variance:
  ## 1. OLS
  ## 2. MEAD, no weight, no cor
  ## 3. MEAD, no corr
  ## 4. MEAD, known cor
  ## 5. MEAD, est cor 0.1
  ## 6. MEAD, est cor 0.3
  ## 7. MEAD, est cor 0.5
  ## 8. MEAD+cv, known cor
  ## 9. MEAD+cv, est cor 0.1
  ## 10. MEAD+cv, est cor 0.3
  ## 11. MEAD+cv, est cor 0.5

  p_hat_ols_se = array(dim = c(K,n_bulk,nreps))
  p_hat_se = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_knowncor = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cor01 = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cor03 = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cor05 = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cv_knowncor = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cv_cor01 = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cv_cor03 = array(dim = c(K,n_bulk,nreps))
  p_hat_weight_se_cv_cor05 = array(dim = c(K,n_bulk,nreps))




  gene_names = rownames(ref)
  celltypes = colnames(ref)
  #true_betas = matrix(nrow=nreps,ncol=n_bulk*K)
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

  # true R01
  cor.idx = which(R!=0,arr.ind = T)
  R01_true = sparseMatrix(i = cor.idx[,1],j = cor.idx[,2],dims=c(G,G))

  # estimated R01
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
  #cor.idx = get_cor_pairs2(bulk_for_cor,alpha=alpha.cor,method=cor_method)
  R01_01 = get_R01(bulk_for_cor,0.1)
  R01_03 = get_R01(bulk_for_cor,0.3)
  R01_05 = get_R01(bulk_for_cor,0.5)

  # all_fit = list()
  true_p = array(dim = c(K,n_bulk,nreps))

  for(reps in 1:nreps){

    if(reps%%printevery==0){print(sprintf("running %d (out of %d)",reps,nreps))}

    ## generate group p's

    if(dirichlet){
      b = t(gtools::rdirichlet(n_bulk,p*dirichlet.scale))
    }else{
      b = matrix(p,ncol=nb1,nrow=K)
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


    # add variance to mimic the average of cell counts
    for(r in 1:n_ref){
      for(k in 1:K){
        X_array_ref[,k,r] = rgamma(G,shape=gamma_shape,scale = X_array_ref[,k,r]/gamma_shape)
      }
    }

    # select marker gene
    ref_samples = c()
    for(k in 1:K){
      temp = X_array_ref[,k,]
      colnames(temp) = rep(colnames(ref)[k],n_ref)
      ref_samples = cbind(ref_samples,temp)
    }
    rownames(ref_samples) = gene_names
    sig_mat = build_signature_matrix_CIBERSORT(ref_samples)
    sig_mat = sig_mat[,match(celltypes,colnames(sig_mat))]
    sig_genes = rownames(sig_mat)
    marker.idx = match(sig_genes,gene_names)


    mb = lapply(1:n_bulk,function(i){X_array_bulk[,,i]%*%b[,i]})
    mb = do.call(cbind,mb)
    mb = mb*rnorm(G,1,bulk_bias_sd)
    #true.beta = t(t(b)*c(apply(mb,2,function(z){bulk_lib_size*G/sum(z)})))
    #true_betas[reps,] = c(true.beta)
    thetab = apply(mb,2,function(z){z/sum(z)})


    #browser()
    y = matrix(rpois(G*n_bulk,bulk_lib_size*G*thetab),nrow=G)
    rownames(y) = gene_names


    # fit model

    #browser()

    X = apply(X_array_ref,c(1,2),mean,na.rm=TRUE)
    rownames(X) = gene_names
    colnames(X) = 1:K
    V = t(apply(X_array_ref,c(1),function(z){(cov(t(z),use = 'complete.obs'))}))/n_ref



    # ols fit
    fit.ols = ols_hc3(y,X)
    p_hat_ols[,,reps] = fit.ols$p_hat
    p_hat_ols_se[,,reps] = fit.ols$p_hat_se


    # adjust for measurement error, do not adjust for correlation
    fit.err = MEAD_est(y=y,
                       X=X,
                       Vg=V,
                       w=1,
                       hc.type='hc3',
                       R01=NULL)
    p_hat[,,reps] = fit.err$p_hat
    p_hat_se[,,reps] = fit.err$p_hat_se


    # use marker, no weight
    fit.marker = MEAD_est(y=y[marker.idx,],
                          X=X[marker.idx,],
                          Vg = V[marker.idx,],
                          w=1,
                          calc_var = FALSE)
    p_hat_marker[,,reps] = fit.marker$p_hat


    # add weight, no cor
    fit.vash = vashr::vash(sqrt(rowSums(V)),df=n_ref-1)
    w = 1/(fit.vash$sd.post)^2

    fit.mead.nocor = MEAD_est(y=y,
                              X=X,
                              Vg=V,
                              w=w,
                              hc.type='hc3',
                              R01=NULL)
    p_hat_weight[,,reps] = fit.mead.nocor$p_hat
    p_hat_weight_se[,,reps] = fit.mead.nocor$p_hat_se


    fit.mead.knowncor = MEAD_est(y=y,
                                 X=X,
                                 Vg=V,
                                 w=w,
                                 hc.type='hc3',
                                 R01=R01_true)
    p_hat_weight_se_knowncor[,,reps] = fit.mead.knowncor$p_hat_se

    fit.mead.cor01 = MEAD_est(y=y,
                                 X=X,
                                 Vg=V,
                                 w=w,
                                 hc.type='hc3',
                                 R01=R01_01)
    p_hat_weight_se_cor01[,,reps] = fit.mead.cor01$p_hat_se


    fit.mead.cor03 = MEAD_est(y=y,
                              X=X,
                              Vg=V,
                              w=w,
                              hc.type='hc3',
                              R01=R01_03)

    p_hat_weight_se_cor03[,,reps] = fit.mead.cor03$p_hat_se

    fit.mead.cor05 = MEAD_est(y=y,
                              X=X,
                              Vg=V,
                              w=w,
                              hc.type='hc3',
                              R01=R01_05)

    p_hat_weight_se_cor05[,,reps] = fit.mead.cor05$p_hat_se


    folds = cut(1:G,breaks = nfold,labels = F)

    fit.mead.cv.knowncor = MEAD_est(y=y,
                              X=X,
                              Vg=V,
                              w=w,
                              hc.type='cv_indep',
                              R01=R01_true,
                              folds=folds)

    p_hat_weight_se_cv_knowncor[,,reps] = fit.mead.cv.knowncor$p_hat_se

    fit.mead.cv.cor01 = MEAD_est(y=y,
                                    X=X,
                                    Vg=V,
                                    w=w,
                                    hc.type='cv_indep',
                                    R01=R01_01,
                                 folds=folds)

    p_hat_weight_se_cv_cor01[,,reps] = fit.mead.cv.cor01$p_hat_se



    fit.mead.cv.cor03 = MEAD_est(y=y,
                                 X=X,
                                 Vg=V,
                                 w=w,
                                 hc.type='cv_indep',
                                 R01=R01_03,
                                 folds=folds)

    p_hat_weight_se_cv_cor03[,,reps] = fit.mead.cv.cor03$p_hat_se


    fit.mead.cv.cor05 = MEAD_est(y=y,
                                 X=X,
                                 Vg=V,
                                 w=w,
                                 hc.type='cv_indep',
                                 R01=R01_05,
                                 folds=folds)

    p_hat_weight_se_cv_cor05[,,reps] = fit.mead.cv.cor05$p_hat_se


    # run music
    Sigma.music = t(apply(X_array_ref,c(1),function(z){diag(cov(t(z),use = 'complete.obs'))}))
    Est.prop.allgene = NULL
    for(bb in 1:n_bulk){
      fit.music = music.basic(y[,bb],X,S=1,Sigma.music,iter.max=1000,nu=1e-4,eps=0.01)
      Est.prop.allgene = cbind(Est.prop.allgene, fit.music$p.weight)
    }
    p_hat_music[,,reps] = Est.prop.allgene

    # run cibersort
    fit_cibersort = CIBERSORT(sig_mat,data.frame(GeneSymbol = rownames(y),y))
    c.order = match(celltypes,colnames(fit_cibersort[,1:K]))
    p_hat_ciber[,,reps] = t(fit_cibersort[,1:K])[c.order,]


  }

  return(list(p_hat_ols = p_hat_ols,
                p_hat = p_hat,
                p_hat_marker = p_hat_marker,
                p_hat_weight = p_hat_weight,
                p_hat_music = p_hat_music,
                p_hat_ciber = p_hat_ciber,

              p_hat_ols_se = p_hat_ols_se,
              p_hat_se = p_hat_se,
              p_hat_weight_se = p_hat_weight_se,
              p_hat_weight_se_knowncor = p_hat_weight_se_knowncor,
              p_hat_weight_se_cor01 = p_hat_weight_se_cor01,
              p_hat_weight_se_cor03 = p_hat_weight_se_cor03,
              p_hat_weight_se_cor05 = p_hat_weight_se_cor05,
              p_hat_weight_se_cv_knowncor = p_hat_weight_se_cv_knowncor,
              p_hat_weight_se_cv_cor01 = p_hat_weight_se_cv_cor01,
              p_hat_weight_se_cv_cor03 = p_hat_weight_se_cv_cor03,
              p_hat_weight_se_cv_cor05 = p_hat_weight_se_cv_cor05,


                true_p = true_p,
                p = p))





}



##################################
#################################
###### test function #############
##################################
###################################
# print('Running test...')
# K = 3
# G = 100
# ref = abs(matrix(rnorm(G*K,3,3),nrow=G,ncol=K))
# rownames(ref) = 1:G
# colnames(ref) = c('a','beta','c')
# sigma2 = ref/3
# p = 1:K/(sum(1:K))
# simu_study(ref,
#                       sigma2,
#                       p,
#                       n_bulk=2,
#                       dirichlet = TRUE,
#                       dirichlet.scale = 10,
#                       d=5,
#                       nreps = 2,
#                       bulk_lib_size = 500,
#                       n_ref = 10,
#                       printevery=1,
#                       n_bulk_for_cor = 100,
#                       gamma_shape = 5,
#                       bulk_bias_sd = 0.1,
#                       nfold = 10
# )
#
# n_bulk=2
# dirichlet = TRUE
# dirichlet.scale = 10
# d=5
# nreps = 2
# bulk_lib_size = 500
# n_ref = 10
# printevery=1
# n_bulk_for_cor = 100
# gamma_shape = 5
# bulk_bias_sd = 0.1
# nfold = 10


xin = readRDS('data/pancreas/xin_ref_sigma9496.rds')
ref = xin$ref
sigma2 = xin$sigma2
set.seed(12345)
simu_out = simu_study(ref,sigma2,p=c(0.5,0.3,0.1,0.1),n_bulk = 50,dirichlet.scale=1000,printevery = 1,gamma_shape = 5)
saveRDS(simu_out,file = paste('output/manuscript/simulation/final/50bulk_500genecor_dirichlet1000_gammashape5.rds',sep=''))

#simu_out = readRDS('output/manuscript/simulation/final/50bulk_500genecor_dirichlet1000_gammashape5.rds')

set.seed(12345)
simu_out = simu_study(ref,sigma2,p=c(0.5,0.3,0.1,0.1),n_bulk = 50,dirichlet.scale=100,printevery = 1,gamma_shape = 5)
saveRDS(simu_out,file = paste('output/manuscript/simulation/final/50bulk_500genecor_dirichlet100_gammashape5.rds',sep=''))

#simu_out = readRDS('output/manuscript/simulation/final/50bulk_500genecor_dirichlet100_gammashape5.rds')

for(m in 1:6){
  print(paste(names(simu_out)[m],':',get_rmse_array(simu_out[[m]],simu_out$true_p)))
}

print(paste(names(simu_out)[7],':',mean(get_coverage_p(simu_out$p_hat_ols,simu_out[[7]],simu_out$true_p))))
print(paste(names(simu_out)[8],':',mean(get_coverage_p(simu_out$p_hat,simu_out[[8]],simu_out$true_p))))
for(m in 9:17){
  print(paste(names(simu_out)[m],':',mean(get_coverage_p(simu_out$p_hat_weight,simu_out[[m]],simu_out$true_p))))
}

################################
######## plot for CI length#####

# get_mean_coverage_per_individual = function(p_hat,p_hat_se,true_p,alpha = 0.05){
#   lower_array = p_hat - qnorm(1-alpha/2)*p_hat_se
#   lower_array = ifelse(lower_array>0,lower_array,0)
#
#   upper_array = p_hat + qnorm(1-alpha/2)*p_hat_se
#   upper_array = ifelse(upper_array>1,1,upper_array)
#
#   lower_mean = apply(lower_array,c(1,2),mean,na.rm=TRUE)
#   upper_mean = apply(upper_array,c(1,2),mean,na.rm=TRUE)
#   return(list(lower_mean=lower_mean,upper_mean = upper_mean))
# }
#
# lu = get_mean_coverage_per_individual(simu_out$p_hat_weight,simu_out$p_hat_weight_se_cv_cor05,simu_out$true_p)

simu_out = readRDS('output/manuscript/simulation/final/50bulk_500genecor_dirichlet10_gammashape5.rds')

get_coverage_for_one_rep = function(p_hat,p_hat_se,true_p,alpha = 0.05){
  lower_array = p_hat - qnorm(1-alpha/2)*p_hat_se
  lower_array = ifelse(lower_array>0,lower_array,0)

  upper_array = p_hat + qnorm(1-alpha/2)*p_hat_se
  upper_array = ifelse(upper_array>1,1,upper_array)

  return(list(lower=lower_array,upper = upper_array,true_p = true_p))
}
lu = get_coverage_for_one_rep(simu_out$p_hat_weight,simu_out$p_hat_weight_se_cv_cor01,simu_out$true_p)

library(ggplot2)
library(gridExtra)


plot_list = list()
iter=1
celltypes = c('alpha', 'beta', 'delta', 'gamma')
for(k in 1:4){

  datax = data.frame(lower = lu$lower[k,,iter],
                     upper = lu$upper[k,,iter],
                     true_p = lu$true_p[k,,iter])
  datax = datax[order(datax$true_p),]
  datax$indi = 1:50

  plot_list[[k]] = ggplot(datax,aes(x = indi,y=true_p))+
    geom_line()+
    geom_point(size = 1,color = 'red')+
    #scale_linetype_manual(values=c('twodash'))+
    #geom_hline(yintercept = 0.5,linetype='dashed', col = 'red')+
    #geom_hline(yintercept = 0.1,linetype='dashed', col = 'red')+
    #scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
    xlab("Sorted bulk samples") +
    ylab("p")+
    geom_ribbon(data=datax,aes(ymin=lower,ymax=upper),alpha=0.3)+
    ggtitle(paste('Cell type:', celltypes[k]))+
     ylim(0,1)
}

grid.arrange(grobs=plot_list,ncol=2)
# for(iter in 1:6){
#   datax = data.frame(individuals = rep(1:50,4),
#                      ci = c(lu$lower[1,,iter],lu$upper[1,,iter],lu$lower[3,,iter],lu$upper[3,,iter]),
#                      true_p = c(lu$true_p[1,,iter],lu$true_p[3,,iter]),
#                      p = rep(c(0.5,0.1),each=100),
#                      bound = rep(rep(c('lower','upper'),each=50),2))
#
#   datax$p = factor(datax$p)
#   datax$bound = factor(datax$bound)
#
#   # ggplot(datax,aes(x=individuals,y=ci,shape = p))+
#   #   geom_point(color='black')+
#   #   scale_shape_manual(values=c(3, 16))+
#   #   #scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))+
#   #   geom_hline(yintercept = 0.5,linetype='dashed', col = 'red')+
#   #   geom_hline(yintercept = 0.1,linetype='dashed', col = 'red')+
#   #   scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
#   #   xlab("Bulk samples") +
#   #   ylab("")
#   #   #labs(color='CI bound',shape = 'p')
#
#   plot_list[[iter]] = ggplot(datax,aes(x=individuals,y=ci,linetype = p,color = bound))+
#     geom_line(size = 0.8)+
#     scale_linetype_manual(values=c('twodash', 'solid'))+
#     scale_color_manual(values=c('black','black'))+
#     #geom_hline(yintercept = 0.5,linetype='dashed', col = 'red')+
#     #geom_hline(yintercept = 0.1,linetype='dashed', col = 'red')+
#     scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
#     xlab("Bulk samples") +
#     ylab("")+ guides(color = 'none')
# }





#labs(color='CI bound',shape = 'p')
###################


get_rmse = function(p_hat,b){
  K = dim(p_hat)[1]
  nb = dim(p_hat)[2]
  nreps = dim(p_hat)[3]
  rmses = c()
  for(i in 1:nb){
    err = c()
    for(j in 1:nreps){
      err[j] = sum((p_hat[,i,j]-b[,i])^2)
    }
    rmses[i] = sqrt(mean(err))
  }
  names(rmses) = paste('bulk',1:nb)
  round(rmses,3)

}

get_rmse_array = function(x,y){
  nn = dim(x)[3]
  ses = c()
  for(i in 1:nn){
    ses[i] = mean((x[,,i]-y[,,i])^2)
  }
  sqrt(mean(ses))
}

rmse = function(x,y){
  sqrt(mean((x-y)^2))
}

get_coverage_p = function(p_hat,p_hat_se,b_array){

  K = dim(p_hat)[1]
  nb = dim(p_hat)[2]
  z = array(dim = dim(p_hat))
  for(i in 1:dim(z)[3]){
    z[,,i] = (p_hat[,,i]-b_array[,,i])/p_hat_se[,,i]
  }
  crg = apply(z,c(1,2),function(z){round(mean(abs(z)<1.96,na.rm=T),3)})
  rownames(crg) = paste('cell',1:K)
  colnames(crg) = paste('bulk',1:nb)
  crg
}
