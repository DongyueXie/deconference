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



########################
###### get coverage ####
########################

K = 4
alpha.cors = c(0,0.1,0.5)
cases = c("null")
nbs = c(100)
dirichlet.scale = c(5)

coverages = c()
coverages.sd = c()

set.seed(12345)
for(case in cases){

  if(case=='null'){
    p1 = c(0.5,0.3,0.1,0.1)
    p2 = c(0.5,0.3,0.1,0.1)
  }else if(case=='all_diff'){
    p1 = c(0.15,0.2,0.45,0.2)
    p2 = c(0.1,0.1,0.3,0.5)
  }
  dif = p1-p2

  for(nb in nbs){
    for(aa in dirichlet.scale){
      for(alpha.cor in alpha.cors){

        if(alpha.cor==0){
          est_cor = FALSE
          cor.status = 'trueR'
        }else{
          est_cor=TRUE
          cor.status = paste('cor0',alpha.cor*10,sep = '')
        }

        print(paste('Running:',case,'nb=',nb,'cor:',cor.status,"dirichlet.scale=",aa))

        #simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1,alpha.cor = alpha.cor,est_cor=est_cor)
        simu_out = readRDS(file = paste('output/manuscript/simulation/simulation_',nb,'bulk_500genecor_',cor.status,'_',case,'_dirichlet',aa,'no_pd.rds',sep=''))

        #print("coverage of p")
        cc = c(mean(get_coverage_p(simu_out$p_hat_ols,simu_out$p_hat_ols_se,simu_out$true_p)),
               mean(get_coverage_p(simu_out$p_hat,simu_out$p_hat_se,simu_out$true_p)),
               mean(get_coverage_p(simu_out$p_hat_weight,simu_out$p_hat_weight_se,simu_out$true_p)),
               mean(get_coverage_p(simu_out$p_hat,simu_out$p_hat_se_cor,simu_out$true_p)),
               mean(get_coverage_p(simu_out$p_hat_weight,simu_out$p_hat_weight_se_cor,simu_out$true_p)),
               mean(get_coverage_p(simu_out$p_hat_weight,simu_out$p_hat_weight_se_cor_cv,simu_out$true_p)))
        coverages = cbind(coverages,cc)
        cc.sd = c(sd(get_coverage_p(simu_out$p_hat_ols,simu_out$p_hat_ols_se,simu_out$true_p)),
                  sd(get_coverage_p(simu_out$p_hat,simu_out$p_hat_se,simu_out$true_p)),
                  sd(get_coverage_p(simu_out$p_hat_weight,simu_out$p_hat_weight_se,simu_out$true_p)),
                  sd(get_coverage_p(simu_out$p_hat,simu_out$p_hat_se_cor,simu_out$true_p)),
                  sd(get_coverage_p(simu_out$p_hat_weight,simu_out$p_hat_weight_se_cor,simu_out$true_p)),
                  sd(get_coverage_p(simu_out$p_hat_weight,simu_out$p_hat_weight_se_cor_cv,simu_out$true_p)))
        coverages.sd = cbind(coverages.sd,cc.sd)

      }
    }
  }
}
rownames(coverages) = c('ols+hc3','mea.err+hc3','mea.err+hc3+weight','mea.err+hc3+cor','mea.err+hc3+cor+weight','mea.err+cv+cor+weight')
colnames(coverages) = c('known.cor','fdr0.1','fdr0.5')
rownames(coverages.sd) = c('ols+hc3','mea.err+hc3','mea.err+hc3+weight','mea.err+hc3+cor','mea.err+hc3+cor+weight','mea.err+cv+cor+weight')
colnames(coverages.sd) = c('known.cor','fdr0.1','fdr0.5')
round(coverages,3)
round(coverages.sd,3)

########################
#### get group diff ####
########################



########################
###### get rmse ########
########################
K = 4
nreps = 100
alpha.cors = c(0)
cases = c("all_diff")
nbs = c(100)
dirichlet.scale = c(10)
rmses = c()

set.seed(12345)
for(case in cases){

  if(case=='null'){
    p1 = c(0.5,0.3,0.1,0.1)
    p2 = c(0.5,0.3,0.1,0.1)
  }else if(case=='all_diff'){
    p1 = c(0.15,0.2,0.45,0.2)
    p2 = c(0.1,0.1,0.3,0.5)
  }
  dif = p1-p2

  for(nb in nbs){
    for(aa in dirichlet.scale){
      for(alpha.cor in alpha.cors){

        if(alpha.cor==0){
          est_cor = FALSE
          cor.status = 'trueR'
        }else{
          est_cor=TRUE
          cor.status = paste('cor0',alpha.cor*10,sep = '')
        }

        print(paste('Running:',case,'nb=',nb,'cor:',cor.status,"dirichlet.scale=",aa))

        #simu_out = simu_study(ref,b,R=A,sigma2,printevery = 1,alpha.cor = alpha.cor,est_cor=est_cor)
        simu_out = readRDS(file = paste('output/manuscript/simulation/simulation_',nb,'bulk_500genecor_',cor.status,'_',case,'_dirichlet',aa,'no_pd.rds',sep=''))
        simu_out_marker = readRDS(file = paste('output/manuscript/simulation/simulation_',nb,'bulk_500genecor_marker_meaerr_phat_',case,'_dirichlet',aa,'no_pd.rds',sep=''))

        rr = c(get_rmse_array(simu_out$p_hat_ols,simu_out$true_p),
               get_rmse_array(simu_out$p_hat,simu_out$true_p),
               get_rmse_array(simu_out$p_hat_weight,simu_out$true_p),
               get_rmse_array(simu_out$p_hat_music,simu_out$true_p),
               get_rmse_array(simu_out_marker$p_hat_weight_marker,simu_out_marker$true_p))
        rmses = cbind(rmses,rr)

      }
    }
  }
}
rownames(rmses) = c('ols','mea.err','mea.err+weight','music','mea.err+marker')
colnames(rmses) = cases
round(rmses,3)

