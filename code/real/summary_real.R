# look at rmse

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

# rmse

out = readRDS('output/manuscript/real/neuron_ref11_rep100_dirichlet10_corfdr005.rds')

# 3 methods
K = 6
nreps = 100
n_bulk = 86

# my fit
p_hat_array_my = array(dim = c(K,n_bulk,nreps))
p_hat_array_my = array(dim = c(K,n_bulk,nreps))
for(r in 1:nreps){
  p_hat_array[,,r] = out$my_fit[[r]]$p_hat
  p_hat_array_my[,,r] = out$rep_info[[r]]$b
}
get_rmse_array(p_hat_array,p_hat_array_my)

# music
p_hat_array_music = array(dim = c(K,n_bulk,nreps))
for(r in 1:nreps){
  p_hat_array_music[,,r] = t(out$music_fit[[r]]$Est.prop.weighted)
}
get_rmse_array(p_hat_array_music,p_true_array)

# cibersort
p_hat_array_cibersort = array(dim = c(K,n_bulk,nreps))
for(r in 1:nreps){
  p_hat_array_cibersort[,,r] = t(out$cibersort_fit[[r]][,1:K])
}
get_rmse_array(p_hat_array_cibersort,p_true_array)

# coverage of p

p_hat_array_my = array(dim = c(K,n_bulk,nreps))
p_hat_array_se = array(dim = c(K,n_bulk,nreps))
p_true_array = array(dim = c(K,n_bulk,nreps))
for(r in 1:nreps){
  p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
  p_hat_array_se[,,r] = out$my_fit[[r]]$p_hat_se
  p_true_array[,,r] = out$rep_info[[r]]$b
}
get_rmse_array(p_hat_array_my,p_true_array)

mean(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array))
sd(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array))

# coverage of two sample test
dif = out$p1-out$p2
cover.naive0 = matrix(nrow = nreps,ncol=K)
sd.naive0 = matrix(nrow = nreps,ncol=K)
for(i in 1:nreps){
  for(k in 1:K){
    temp = t.test(p_true_array[k,out$rep_info[[i]]$groups==1,i],p_true_array[k,out$rep_info[[i]]$groups==2,i])
    sd.naive0[i,k] = temp$stderr
    cover.naive0[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
  }
}
mean(cover.naive0)

cover.naive = matrix(nrow = nreps,ncol=K)
sd.naive = matrix(nrow = nreps,ncol=K)
for(i in 1:nreps){
  for(k in 1:K){
    temp = t.test(p_hat_array_my[k,out$rep_info[[i]]$groups==1,i],p_hat_array_my[k,out$rep_info[[i]]$groups==2,i])
    sd.naive[i,k] = temp$stderr
    cover.naive[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
  }
}

cover.naive.music = matrix(nrow = nreps,ncol=K)
sd.naive.music = matrix(nrow = nreps,ncol=K)
for(i in 1:nreps){
  for(k in 1:K){
    temp = t.test(p_hat_array_music[k,out$rep_info[[i]]$groups==1,i],p_hat_array_music[k,out$rep_info[[i]]$groups==2,i])
    sd.naive.music[i,k] = temp$stderr
    cover.naive.music[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
  }
}

cover.naive.cibersort = matrix(nrow = nreps,ncol=K)
sd.naive.cibersort = matrix(nrow = nreps,ncol=K)
for(i in 1:nreps){
  for(k in 1:K){
    temp = t.test(p_hat_array_cibersort[k,out$rep_info[[i]]$groups==1,i],p_hat_array_cibersort[k,out$rep_info[[i]]$groups==2,i])
    sd.naive.cibersort[i,k] = temp$stderr
    cover.naive.cibersort[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
  }
}

diff_hat_array_my = matrix(nrow=nreps,ncol=K)
diff_hat_array_se = matrix(nrow=nreps,ncol=K)
for(r in 1:nreps){
  diff_hat_array_my[r,] = out$my_fit[[r]]$two_group_res$diff_hat
  diff_hat_array_se[r,] = out$my_fit[[r]]$two_group_res$diff_hat_se
}
temp = abs(diff_hat_array_my - rep(1,nreps)%*%t(dif))/diff_hat_array_se
cover.asy.cv = temp<1.96


cc = c(mean(cover.naive0),
       mean(cover.naive.music),
       mean(cover.naive),
       mean(cover.asy.cv,na.rm=TRUE))



names(cc)=(c('t+truep','t+music','t+mea.err+weight','asy.weight.hc3.cv'))
cc









#################################################
#################################################
#################################################
#################################################
#################################################

n_ref=11
K = 6
nreps = 100
n_bulk = 86
celltypes = c("DA","Epen1","Sert","FPP","P_FPP","U_Neur")
dirichlet.scales = c(10)
cases = c("null","all_diff")

rmses = c()

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

    out = readRDS(paste("output/manuscript/real/neuron_ref",n_ref,"_rep",n_rep,"_dirichlet",aa,"_corfdr005_",case,'.rds',sep=''))



    # my fit
    p_hat_array_my = array(dim = c(K,n_bulk,nreps))
    p_true_array = array(dim = c(K,n_bulk,nreps))
    for(r in 1:nreps){
      p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
      p_true_array[,,r] = out$rep_info[[r]]$b
    }


    # music
    p_hat_array_music = array(dim = c(K,n_bulk,nreps))
    for(r in 1:nreps){
      c.order = match(celltypes,colnames(out$music_fit[[r]]$Est.prop.weighted))
      p_hat_array_music[,,r] = t(out$music_fit[[r]]$Est.prop.weighted)[c.order,]
    }


    # cibersort
    p_hat_array_cibersort = array(dim = c(K,n_bulk,nreps))
    for(r in 1:nreps){
      c.order = match(celltypes,colnames(out$cibersort_fit[[r]][,1:K]))
      p_hat_array_cibersort[,,r] = t(out$cibersort_fit[[r]][,1:K])[c.order,]
    }


    rmses = cbind(rmses,c(get_rmse_array(p_hat_array_my,p_true_array),
                          get_rmse_array(p_hat_array_music,p_true_array),
                          get_rmse_array(p_hat_array_cibersort,p_true_array)))

  }
}

rownames(rmses) = c('mea.err','music','cibersort')
colnames(rmses) = cases
round(rmses,3)



### coverage of p

dirichlet.scales = c(5,10)
cases = c("null","all_diff")

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

    out = readRDS(paste("output/manuscript/real/neuron_ref",n_ref,"_rep",n_rep,"_dirichlet",aa,"_corfdr005_",case,'.rds',sep=''))



    p_hat_array_my = array(dim = c(K,n_bulk,nreps))
    p_hat_array_se = array(dim = c(K,n_bulk,nreps))
    p_true_array = array(dim = c(K,n_bulk,nreps))
    for(r in 1:nreps){
      p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
      p_hat_array_se[,,r] = out$my_fit[[r]]$p_hat_se
      p_true_array[,,r] = out$rep_info[[r]]$b
    }
    get_rmse_array(p_hat_array_my,p_true_array)




    res = c(res,mean(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array)))
    res2 = c(res2, sd(get_coverage_p(p_hat_array_my,p_hat_array_se,p_true_array)))

  }
}


res
res2


### two sample t test

library(reticulate)
np <- import("numpy")
celltypes = c("DA" ,"Epen1","Sert","FPP","P_FPP","U_Neur")
celltypes_sieve = c('DA', 'Epen1', 'FPP', 'P_FPP', 'Sert', 'U_Neur')
dirichlet.scales = c(5)
cases = c("null","all_diff")

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

    out = readRDS(paste("output/manuscript/real/neuron_ref",n_ref,"_rep",n_rep,"_dirichlet",aa,"_corfdr005_",case,'.rds',sep=''))



    # my fit
    p_hat_array_my = array(dim = c(K,n_bulk,nreps))
    p_true_array = array(dim = c(K,n_bulk,nreps))
    for(r in 1:nreps){
      p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
      p_true_array[,,r] = out$rep_info[[r]]$b
    }


    # music
    p_hat_array_music = array(dim = c(K,n_bulk,nreps))
    for(r in 1:nreps){
      c.order = match(celltypes,colnames(out$music_fit[[r]]$Est.prop.weighted))
      p_hat_array_music[,,r] = t(out$music_fit[[r]]$Est.prop.weighted)[c.order,]
    }


    # cibersort
    p_hat_array_cibersort = array(dim = c(K,n_bulk,nreps))
    for(r in 1:nreps){
      c.order = match(celltypes,colnames(out$cibersort_fit[[r]][,1:K]))
      p_hat_array_cibersort[,,r] = t(out$cibersort_fit[[r]][,1:K])[c.order,]
    }

    # rna-sieve
    p_hat_sieve = np$load(paste("output/manuscript/real/rnasieve/neuron_ref",n_ref,"_rep",n_rep,"_dirichlet",aa,"_",case,'_p_hat.npy',sep=''))
    for(r in 1:nreps){
      p_hat_sieve[,,r] = p_hat_sieve[match(celltypes,celltypes_sieve),,r]
    }


    dif = p1-p2
    cover.naive0 = matrix(nrow = nreps,ncol=K)
    sd.naive0 = matrix(nrow = nreps,ncol=K)
    for(i in 1:nreps){
      for(k in 1:K){
        temp = t.test(p_true_array[k,out$rep_info[[i]]$groups==1,i],p_true_array[k,out$rep_info[[i]]$groups==2,i])
        sd.naive0[i,k] = temp$stderr
        cover.naive0[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }
    mean(cover.naive0)

    cover.naive = matrix(nrow = nreps,ncol=K)
    sd.naive = matrix(nrow = nreps,ncol=K)
    for(i in 1:nreps){
      for(k in 1:K){
        temp = t.test(p_hat_array_my[k,out$rep_info[[i]]$groups==1,i],p_hat_array_my[k,out$rep_info[[i]]$groups==2,i])
        sd.naive[i,k] = temp$stderr
        cover.naive[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }

    cover.naive.music = matrix(nrow = nreps,ncol=K)
    sd.naive.music = matrix(nrow = nreps,ncol=K)
    for(i in 1:nreps){
      for(k in 1:K){
        temp = t.test(p_hat_array_music[k,out$rep_info[[i]]$groups==1,i],p_hat_array_music[k,out$rep_info[[i]]$groups==2,i])
        sd.naive.music[i,k] = temp$stderr
        cover.naive.music[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }

    cover.naive.cibersort = matrix(nrow = nreps,ncol=K)
    sd.naive.cibersort = matrix(nrow = nreps,ncol=K)
    for(i in 1:nreps){
      for(k in 1:K){
        temp = t.test(p_hat_array_cibersort[k,out$rep_info[[i]]$groups==1,i],p_hat_array_cibersort[k,out$rep_info[[i]]$groups==2,i])
        sd.naive.cibersort[i,k] = temp$stderr
        cover.naive.cibersort[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }

    cover.naive.rnasieve = matrix(nrow = nreps,ncol=K)
    sd.naive.rnasieve = matrix(nrow = nreps,ncol=K)
    for(i in 1:nreps){
      for(k in 1:K){
        temp = t.test(p_hat_sieve[k,out$rep_info[[i]]$groups==1,i],p_hat_sieve[k,out$rep_info[[i]]$groups==2,i])
        sd.naive.rnasieve[i,k] = temp$stderr
        cover.naive.rnasieve[i,k] = (dif[k]>temp$conf.int[1])&(dif[k]<temp$conf.int[2])
      }
    }

    diff_hat_array_my = matrix(nrow=nreps,ncol=K)
    diff_hat_array_se = matrix(nrow=nreps,ncol=K)
    for(r in 1:nreps){
      diff_hat_array_my[r,] = out$my_fit[[r]]$two_group_res$diff_hat
      diff_hat_array_se[r,] = out$my_fit[[r]]$two_group_res$diff_hat_se
    }
    temp = abs(diff_hat_array_my - rep(1,nreps)%*%t(dif))/diff_hat_array_se
    cover.asy.cv = temp<1.96


    cc = c(mean(cover.naive0),
           mean(cover.naive.music),
           mean(cover.naive.cibersort),
           mean(cover.naive.rnasieve,na.rm=TRUE),
           mean(cover.naive),
           mean(cover.asy.cv,na.rm=TRUE))

    res = cbind(res,cc)


  }
}



rownames(res)=(c('t+truep','t+music','t+cibersort','t+rnasieve','t+mea.err+weight','asy.weight.hc3.cv'))
colnames(res) = cases
res










