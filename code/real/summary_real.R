get_rmse = function(p_hat,b){
  K = dim(p_hat)[1]
  nb = dim(p_hat)[2]
  n_rep = dim(p_hat)[3]
  rmses = c()
  for(i in 1:nb){
    err = c()
    for(j in 1:n_rep){
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

get_rmse_array_cell_type = function(x,y){
  nn = dim(x)[3]
  K = dim(x)[1]
  ses = matrix(nrow=nn,ncol=K)
  for(i in 1:nn){
    ses[i,] = c(rowMeans((x[,,i]-y[,,i])^2))
  }
  sqrt(colMeans(ses))
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

library(reticulate)
np <- import("numpy")
celltypes = c("DA" ,"Epen1","Sert","FPP","P_FPP","U_Neur")
celltypes_sieve = c('DA', 'Epen1', 'FPP', 'P_FPP', 'Sert', 'U_Neur')

n_ref=11
K = 6
n_rep = 100
n_bulk = 86
celltypes = c("DA","Epen1","Sert","FPP","P_FPP","U_Neur")







dirichlet.scales = c(5)
cases = c("all_diff")
#cases = c("null")
rmses = c()

#method_list = c('my','my_unweighted','music','cibersort','rnasieve','ols')
method_list = c('my','music','cibersort','rnasieve','ols')

for(aa in dirichlet.scales){
  for(case in cases){

    if(case=='null'){
      p1 = c(0.3,0.2,0.15,0.15,0.1,0.1)
      p2 = c(0.3,0.2,0.15,0.15,0.1,0.1)
    }else if(case=='all_diff'){
      p1 = c(0.15,0.15,0.1,0.1,0.2,0.3)
      p2 = c(0.1,0.1,0.2,0.3,0.15,0.15)
    }

    print(paste('Running:',case,'dirichlet.scale=',aa,"n_bulk=",n_bulk,'n_ref=',n_ref))




    meth.order=c()

    # my fit
    if('my'%in%method_list){
      out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_",case,'.rds',sep=''))
      p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_my,p_true_array))
      meth.order = c(meth.order,'my')
    }


    # my_unweighted fit
    if('my_unweighted'%in%method_list){
      out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_unweighted_",case,'.rds',sep=''))
      p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_my,p_true_array))
      meth.order = c(meth.order,'my_unweighted')
    }


    # music
    if('music'%in%method_list){
      out = readRDS(paste("output/manuscript/real/music/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
      p_hat_array_music = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        c.order = match(celltypes,colnames(out$music_fit[[r]]$Est.prop.weighted))
        p_hat_array_music[,,r] = t(out$music_fit[[r]]$Est.prop.weighted)[c.order,]
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_music,p_true_array))
      meth.order = c(meth.order,'music')
    }



    # cibersort
    if('cibersort'%in%method_list){
      out = readRDS(paste("output/manuscript/real/cibersort/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
      p_hat_array_cibersort = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        c.order = match(celltypes,colnames(out$cibersort_fit[[r]][,1:K]))
        p_hat_array_cibersort[,,r] = t(out$cibersort_fit[[r]][,1:K])[c.order,]
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_cibersort,p_true_array))
      meth.order = c(meth.order,'cibersort')
    }



    # rna-sieve
    if('rnasieve'%in%method_list){
      p_hat_sieve = np$load(paste("output/manuscript/real/rnasieve/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'_p_hat.npy',sep=''))
      for(r in 1:n_rep){
        p_hat_sieve[,,r] = p_hat_sieve[match(celltypes,celltypes_sieve),,r]
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_sieve,p_true_array))
      meth.order = c(meth.order,'rnasieve')
    }


    # ols
    if('ols'%in%method_list){
      out = readRDS(paste("output/manuscript/real/ols/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
      p_hat_array_ols = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_ols[,,r] = out$ols_fit[[r]]$p_hat
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_ols,p_true_array))
      meth.order = c(meth.order,'ols')


    }


  }
}

rownames(rmses) = celltypes
colnames(rmses) = meth.order
round(rmses,3)

plot(rmses[,1],ylab='RMSE',ylim=range(rmses),xaxt = 'n',xlab='cell types')
axis(1, at=1:K, labels=celltypes)
lines(rmses[,2],type='p',pch=2)
lines(rmses[,3],type='p',pch=3)
lines(rmses[,4],type='p',pch = 4)
lines(rmses[,5],type='p',pch = 5)
legend('topright',c('mea.err','music','cibersort','rnasieve','ols'),pch=c(1,2,3,4,5))


### coverage of p

dirichlet.scales = c(5)
cases = c("null")
method_list = c('my','my_unweighted','rnasieve','ols')

res=c()
res2 = c()
for(aa in dirichlet.scales){
  for(case in cases){

    if(case=='null'){
      p1 = c(0.3,0.2,0.15,0.15,0.1,0.1)
      p2 = c(0.3,0.2,0.15,0.15,0.1,0.1)
    }else if(case=='all_diff'){
      p1 = c(0.15,0.15,0.1,0.1,0.2,0.3)
      p2 = c(0.1,0.1,0.2,0.3,0.15,0.15)
    }

    print(paste('Running:',case,'dirichlet.scale=',aa))

    meth.order=c()
    if('my'%in%method_list){
      out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_",case,'.rds',sep=''))
      p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
      p_hat_array_se = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
        p_hat_array_se[,,r] = out$my_fit[[r]]$p_hat_se
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      #get_rmse_array(p_hat_array_my,p_true_array)

      res = rbind(res, rowMeans(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array)))
      res2 = rbind(res2,apply(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array),1,sd))
      meth.order = c(meth.order,'my')
    }


    if('my_unweighted'%in%method_list){
      out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_unweighted_",case,'.rds',sep=''))
      p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
      p_hat_array_se = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
        p_hat_array_se[,,r] = out$my_fit[[r]]$p_hat_se
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      #get_rmse_array(p_hat_array_my,p_true_array)

      res = rbind(res, rowMeans(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array)))
      res2 = rbind(res2,apply(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array),1,sd))
      meth.order = c(meth.order,'my_unweighted')
    }

    if('ols'%in%method_list){
      out = readRDS(paste("output/manuscript/real/ols/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
      p_hat_array_ols = array(dim = c(K,n_bulk,n_rep))
      p_hat_array_se = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_ols[,,r] = out$ols_fit[[r]]$p_hat
        p_hat_array_se[,,r] = out$ols_fit[[r]]$p_hat_se
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      #get_rmse_array(p_hat_array_my,p_true_array)

      res = rbind(res, rowMeans(get_coverage_p(p_hat_array_ols,p_hat_array_se,p_true_array)))
      res2 = rbind(res2,apply(get_coverage_p(p_hat_array_ols,p_hat_array_se,p_true_array),1,sd))
      meth.order = c(meth.order,'ols')
    }

    if('rnasieve'%in%method_list){
      p_hat_sieve = np$load(paste("output/manuscript/real/rnasieve/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'_p_hat.npy',sep=''))
      p_hat_l_sieve = np$load(paste("output/manuscript/real/rnasieve/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'_p_hat_ci_l.npy',sep=''))
      p_hat_r_sieve = np$load(paste("output/manuscript/real/rnasieve/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'_p_hat_ci_r.npy',sep=''))

      for(r in 1:n_rep){
        p_hat_sieve[,,r] = p_hat_sieve[match(celltypes,celltypes_sieve),,r]
        p_hat_l_sieve[,,r] = p_hat_l_sieve[match(celltypes,celltypes_sieve),,r]
        p_hat_r_sieve[,,r] = p_hat_r_sieve[match(celltypes,celltypes_sieve),,r]
      }

      res = rbind(res, rowMeans(get_coverage_p(p_hat_sieve,(p_hat_r_sieve-p_hat_l_sieve)/1.96/2,p_true_array)))
      res2 = rbind(res2,apply(get_coverage_p(p_hat_sieve,(p_hat_r_sieve-p_hat_l_sieve)/1.96/2,p_true_array),1,sd))
      meth.order = c(meth.order,'rnasieve')
    }


  }
}

colnames(res) = celltypes
colnames(res2) = celltypes
rownames(res) = meth.order
rownames(res2) = meth.order
res
res2


### two sample t test

n_bulk = 86
n_rep = 100
dirichlet.scales = c(5)
cases = c("all_diff")

add_rnasieve = T

res=c()
res_mat = c()
res2 = c()
for(aa in dirichlet.scales){
  for(case in cases){

    if(case=='null'){
      p1 = c(0.3,0.2,0.15,0.15,0.1,0.1)
      p2 = c(0.3,0.2,0.15,0.15,0.1,0.1)
    }else if(case=='all_diff'){
      p1 = c(0.15,0.15,0.1,0.1,0.2,0.3)
      p2 = c(0.1,0.1,0.2,0.3,0.15,0.15)
    }

    print(paste('Running:',case,'dirichlet.scale=',aa))





    # my fit
    out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_",case,'.rds',sep=''))
    p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
    p_true_array = array(dim = c(K,n_bulk,n_rep))
    for(r in 1:n_rep){
      p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
      p_true_array[,,r] = out$rep_info[[r]]$b
    }


    # music
    out = readRDS(paste("output/manuscript/real/music/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
    p_hat_array_music = array(dim = c(K,n_bulk,n_rep))
    for(r in 1:n_rep){
      c.order = match(celltypes,colnames(out$music_fit[[r]]$Est.prop.weighted))
      p_hat_array_music[,,r] = t(out$music_fit[[r]]$Est.prop.weighted)[c.order,]
    }


    # cibersort
    out = readRDS(paste("output/manuscript/real/cibersort/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
    p_hat_array_cibersort = array(dim = c(K,n_bulk,n_rep))
    for(r in 1:n_rep){
      c.order = match(celltypes,colnames(out$cibersort_fit[[r]][,1:K]))
      p_hat_array_cibersort[,,r] = t(out$cibersort_fit[[r]][,1:K])[c.order,]
    }

    # rna-sieve
    if(add_rnasieve){
      p_hat_sieve = np$load(paste("output/manuscript/real/rnasieve/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'_p_hat.npy',sep=''))
      for(r in 1:n_rep){
        p_hat_sieve[,,r] = p_hat_sieve[match(celltypes,celltypes_sieve),,r]
      }
    }



    dif = p1-p2
    cover.naive0 = matrix(nrow = n_rep,ncol=K)
    sd.naive0 = matrix(nrow = n_rep,ncol=K)
    for(i in 1:n_rep){
      for(k in 1:K){
        temp = t.test(p_true_array[k,out$rep_info[[i]]$groups==1,i],p_true_array[k,out$rep_info[[i]]$groups==2,i])
        sd.naive0[i,k] = temp$stderr
        cover.naive0[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }
    mean(cover.naive0)

    cover.naive = matrix(nrow = n_rep,ncol=K)
    sd.naive = matrix(nrow = n_rep,ncol=K)
    for(i in 1:n_rep){
      for(k in 1:K){
        temp = t.test(p_hat_array_my[k,out$rep_info[[i]]$groups==1,i],p_hat_array_my[k,out$rep_info[[i]]$groups==2,i])
        sd.naive[i,k] = temp$stderr
        cover.naive[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }

    cover.naive.music = matrix(nrow = n_rep,ncol=K)
    sd.naive.music = matrix(nrow = n_rep,ncol=K)
    for(i in 1:n_rep){
      for(k in 1:K){
        temp = t.test(p_hat_array_music[k,out$rep_info[[i]]$groups==1,i],p_hat_array_music[k,out$rep_info[[i]]$groups==2,i])
        sd.naive.music[i,k] = temp$stderr
        cover.naive.music[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }

    cover.naive.cibersort = matrix(nrow = n_rep,ncol=K)
    sd.naive.cibersort = matrix(nrow = n_rep,ncol=K)
    for(i in 1:n_rep){
      for(k in 1:K){
        temp = t.test(p_hat_array_cibersort[k,out$rep_info[[i]]$groups==1,i],p_hat_array_cibersort[k,out$rep_info[[i]]$groups==2,i])
        sd.naive.cibersort[i,k] = temp$stderr
        cover.naive.cibersort[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }

    if(add_rnasieve){
      cover.naive.rnasieve = matrix(nrow = n_rep,ncol=K)
      sd.naive.rnasieve = matrix(nrow = n_rep,ncol=K)
      for(i in 1:n_rep){
        for(k in 1:K){
          temp = t.test(p_hat_sieve[k,out$rep_info[[i]]$groups==1,i],p_hat_sieve[k,out$rep_info[[i]]$groups==2,i])
          sd.naive.rnasieve[i,k] = temp$stderr
          cover.naive.rnasieve[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
        }
      }
    }else{
      cover.naive.rnasieve = NULL
    }


    # diff_hat_array_my = matrix(nrow=n_rep,ncol=K)
    # diff_hat_array_se = matrix(nrow=n_rep,ncol=K)
    # for(r in 1:n_rep){
    #   diff_hat_array_my[r,] = out$my_fit[[r]]$two_group_res$diff_hat
    #   diff_hat_array_se[r,] = out$my_fit[[r]]$two_group_res$diff_hat_se
    # }
    # temp = abs(diff_hat_array_my - rep(1,n_rep)%*%t(dif))/diff_hat_array_se
    # cover.asy.cv = temp<1.96


    cc = c(mean(cover.naive0,na.rm=TRUE),
           mean(cover.naive.music,na.rm=TRUE),
           mean(cover.naive.cibersort,na.rm=TRUE),
           mean(cover.naive.rnasieve,na.rm=TRUE),
           mean(cover.naive,na.rm=TRUE))

    res_mat = cbind(colMeans(cover.naive0),
                    colMeans(cover.naive.music),
                    colMeans(cover.naive.cibersort),
                    colMeans(cover.naive.rnasieve),
                    colMeans(cover.naive))

    res = cbind(res,cc)


  }
}



rownames(res)=(c('t+truep','t+music','t+cibersort','t+rnasieve','t+mea.err+weight'))
#colnames(res) = cases
res

rownames(res_mat) = celltypes
colnames(res_mat) = c('truep','music','cibersort','rnasieve','mea.err')

plot(res_mat[,5],ylab='Coverage',ylim=range(res_mat),xaxt = 'n',xlab='cell types')
abline(h=0.95,lty=2)
axis(1, at=1:K, labels=celltypes)
lines(res_mat[,2],type='p',pch=2)
lines(res_mat[,3],type='p',pch=3)
lines(res_mat[,4],type='p',pch = 4)
lines(res_mat[,1],type='p',pch = 5)
legend('bottomleft',c('mea.err','music','cibersort','rnasieve','truep'),pch=c(1,2,3,4,5))


diff_hat_sieve = matrix(nrow=n_rep,ncol=K)
for(i in 1:n_rep){
  diff_hat_sieve[i,] = rowMeans(p_hat_sieve[,out$rep_info[[i]]$groups==1,i]) - rowMeans(p_hat_sieve[,out$rep_info[[i]]$groups==2,i])
}
diff_true = matrix(dif,nrow=n_rep,ncol=K,byrow = T)

sqrt(colMeans(diff_hat_sieve - diff_true)^2)


#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################
#################################################



# look at rmse

get_rmse = function(p_hat,b){
  K = dim(p_hat)[1]
  nb = dim(p_hat)[2]
  n_rep = dim(p_hat)[3]
  rmses = c()
  for(i in 1:nb){
    err = c()
    for(j in 1:n_rep){
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


# rmse

out = readRDS('output/manuscript/real/neuron_ref11_rep100_dirichlet10_corfdr005.rds')

# 3 methods
K = 6
n_rep = 100
n_bulk = 86

# my fit
p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
for(r in 1:n_rep){
  p_hat_array[,,r] = out$my_fit[[r]]$p_hat
  p_hat_array_my[,,r] = out$rep_info[[r]]$b
}
get_rmse_array(p_hat_array,p_hat_array_my)

# music
p_hat_array_music = array(dim = c(K,n_bulk,n_rep))
for(r in 1:n_rep){
  p_hat_array_music[,,r] = t(out$music_fit[[r]]$Est.prop.weighted)
}
get_rmse_array(p_hat_array_music,p_true_array)

# cibersort
p_hat_array_cibersort = array(dim = c(K,n_bulk,n_rep))
for(r in 1:n_rep){
  p_hat_array_cibersort[,,r] = t(out$cibersort_fit[[r]][,1:K])
}
get_rmse_array(p_hat_array_cibersort,p_true_array)

# coverage of p

p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
p_hat_array_se = array(dim = c(K,n_bulk,n_rep))
p_true_array = array(dim = c(K,n_bulk,n_rep))
for(r in 1:n_rep){
  p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
  p_hat_array_se[,,r] = out$my_fit[[r]]$p_hat_se
  p_true_array[,,r] = out$rep_info[[r]]$b
}
get_rmse_array(p_hat_array_my,p_true_array)

mean(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array))
sd(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array))

# coverage of two sample test
dif = out$p1-out$p2
cover.naive0 = matrix(nrow = n_rep,ncol=K)
sd.naive0 = matrix(nrow = n_rep,ncol=K)
for(i in 1:n_rep){
  for(k in 1:K){
    temp = t.test(p_true_array[k,out$rep_info[[i]]$groups==1,i],p_true_array[k,out$rep_info[[i]]$groups==2,i])
    sd.naive0[i,k] = temp$stderr
    cover.naive0[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
  }
}
mean(cover.naive0)

cover.naive = matrix(nrow = n_rep,ncol=K)
sd.naive = matrix(nrow = n_rep,ncol=K)
for(i in 1:n_rep){
  for(k in 1:K){
    temp = t.test(p_hat_array_my[k,out$rep_info[[i]]$groups==1,i],p_hat_array_my[k,out$rep_info[[i]]$groups==2,i])
    sd.naive[i,k] = temp$stderr
    cover.naive[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
  }
}

cover.naive.music = matrix(nrow = n_rep,ncol=K)
sd.naive.music = matrix(nrow = n_rep,ncol=K)
for(i in 1:n_rep){
  for(k in 1:K){
    temp = t.test(p_hat_array_music[k,out$rep_info[[i]]$groups==1,i],p_hat_array_music[k,out$rep_info[[i]]$groups==2,i])
    sd.naive.music[i,k] = temp$stderr
    cover.naive.music[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
  }
}

cover.naive.cibersort = matrix(nrow = n_rep,ncol=K)
sd.naive.cibersort = matrix(nrow = n_rep,ncol=K)
for(i in 1:n_rep){
  for(k in 1:K){
    temp = t.test(p_hat_array_cibersort[k,out$rep_info[[i]]$groups==1,i],p_hat_array_cibersort[k,out$rep_info[[i]]$groups==2,i])
    sd.naive.cibersort[i,k] = temp$stderr
    cover.naive.cibersort[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
  }
}

diff_hat_array_my = matrix(nrow=n_rep,ncol=K)
diff_hat_array_se = matrix(nrow=n_rep,ncol=K)
for(r in 1:n_rep){
  diff_hat_array_my[r,] = out$my_fit[[r]]$two_group_res$diff_hat
  diff_hat_array_se[r,] = out$my_fit[[r]]$two_group_res$diff_hat_se
}
temp = abs(diff_hat_array_my - rep(1,n_rep)%*%t(dif))/diff_hat_array_se
cover.asy.cv = temp<1.96


cc = c(mean(cover.naive0),
       mean(cover.naive.music),
       mean(cover.naive),
       mean(cover.asy.cv,na.rm=TRUE))



names(cc)=(c('t+truep','t+music','t+mea.err+weight','asy.weight.hc3.cv'))
cc




