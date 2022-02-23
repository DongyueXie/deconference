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

dirichlet.scales = c(20)
#cases = c("all_diff")
cases = c("null")
rmses = c()

#method_list = c('my','my_unweighted','music','cibersort','rnasieve','ols')
method_list = c('MEAD','MuSiC','CIBERSORT','RNA-Sieve','OLS')
#method_list = c('mea.err','my.QP')

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
    if('MEAD'%in%method_list){
      out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_",case,'.rds',sep=''))
      p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_my,p_true_array))
      meth.order = c(meth.order,'MEAD')
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

    if('my.QP'%in%method_list){
      out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_QP_",case,'.rds',sep=''))
      p_hat_array_my = array(dim = c(K,n_bulk,n_rep))
      p_true_array = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_my[,,r] = out$my_fit[[r]]$p_hat
        p_true_array[,,r] = out$rep_info[[r]]$b
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_my,p_true_array))
      meth.order = c(meth.order,'my.QP')
    }


    # music
    if('MuSiC'%in%method_list){
      out = readRDS(paste("output/manuscript/real/music/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
      p_hat_array_music = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        c.order = match(celltypes,colnames(out$music_fit[[r]]$Est.prop.weighted))
        p_hat_array_music[,,r] = t(out$music_fit[[r]]$Est.prop.weighted)[c.order,]
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_music,p_true_array))
      meth.order = c(meth.order,'MuSiC')
    }



    # cibersort
    if('CIBERSORT'%in%method_list){
      out = readRDS(paste("output/manuscript/real/cibersort/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
      p_hat_array_cibersort = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        c.order = match(celltypes,colnames(out$cibersort_fit[[r]][,1:K]))
        p_hat_array_cibersort[,,r] = t(out$cibersort_fit[[r]][,1:K])[c.order,]
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_cibersort,p_true_array))
      meth.order = c(meth.order,'CIBERSORT')
    }



    # rna-sieve
    if('RNA-Sieve'%in%method_list){
      p_hat_sieve = np$load(paste("output/manuscript/real/rnasieve/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'_p_hat.npy',sep=''))
      for(r in 1:n_rep){
        p_hat_sieve[,,r] = p_hat_sieve[match(celltypes,celltypes_sieve),,r]
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_sieve,p_true_array))
      meth.order = c(meth.order,'RNA-Sieve')
    }


    # ols
    if('OLS'%in%method_list){
      out = readRDS(paste("output/manuscript/real/ols/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
      p_hat_array_ols = array(dim = c(K,n_bulk,n_rep))
      for(r in 1:n_rep){
        p_hat_array_ols[,,r] = out$ols_fit[[r]]$p_hat
      }
      rmses = cbind(rmses,get_rmse_array_cell_type(p_hat_array_ols,p_true_array))
      meth.order = c(meth.order,'OLS')


    }


  }
}

rownames(rmses) = celltypes
colnames(rmses) = meth.order
round(rmses,3)
library(ggplot2)
library(reshape2)
# datax = data.frame(rmse = c(rmses),cell.type = rep(celltypes,5),method = rep(c('MEAD','MuSiC','CIBERSORT','RNA-Sieve','OLS'),each = 6))
#
# ggplot(datax, aes(x=cell.type, y=rmse, group=method)) +
#   geom_point(aes(shape=method,color = method))+ theme(panel.grid.minor = element_blank())

rmses
plot.meth.order = c('OLS','MuSiC','CIBERSORT','RNA-Sieve','MEAD')
df = melt(rmses[,match(meth.order,plot.meth.order)])
colnames(df) = c('Cell.type', 'Method','RMSE')
plot1 = ggplot(df, aes(x = Method, y = Cell.type, fill = RMSE)) +
  scale_fill_gradient2(
  low = 'steelblue', high = "red", mid = 'white', midpoint  = quantile(df$RMSE,0.4), limit = range(df$RMSE))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 11),axis.text.y = element_text(size = 11),axis.title.x = element_blank(),axis.title.y = element_blank())

# ggplot(df, aes(x = method, y = cell.type, fill = rmse)) +scale_fill_gradient(
#   low = 'white', high = "red", limit = range(df$rmse))+
#   geom_tile()
#
# ggplot(df, aes(x = method, y = cell.type, fill = rmse)) +scale_fill_gradient(
#   low = 'white', high = "blue", limit = range(df$rmse))+
#   geom_tile()
# plot(rmses[,1],ylab='RMSE',ylim=range(rmses),xaxt = 'n',xlab='cell types')
# axis(1, at=1:K, labels=celltypes)
# lines(rmses[,2],type='p',pch=2)
# lines(rmses[,3],type='p',pch=3)
# lines(rmses[,4],type='p',pch = 4)
# lines(rmses[,5],type='p',pch = 5)
# legend('topright',c('mea.err','music','cibersort','rnasieve','ols'),pch=c(1,2,3,4,5))


### coverage of p

dirichlet.scales = c(5)
cases = c("all_diff")
method_list = c('my','rnasieve')

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
      out = readRDS(paste("output/manuscript/real/my/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr03_",case,'.rds',sep=''))
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
dirichlet.scales = c(20)
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
    diff_hat_my = matrix(nrow=n_rep,ncol=K)
    for(i in 1:n_rep){
      diff_hat_my[i,] = rowMeans(p_hat_array_my[,out$rep_info[[i]]$groups==1,i]) - rowMeans(p_hat_array_my[,out$rep_info[[i]]$groups==2,i])
    }

    # music
    out = readRDS(paste("output/manuscript/real/music/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
    p_hat_array_music = array(dim = c(K,n_bulk,n_rep))
    for(r in 1:n_rep){
      c.order = match(celltypes,colnames(out$music_fit[[r]]$Est.prop.weighted))
      p_hat_array_music[,,r] = t(out$music_fit[[r]]$Est.prop.weighted)[c.order,]
    }
    diff_hat_music = matrix(nrow=n_rep,ncol=K)
    for(i in 1:n_rep){
      diff_hat_music[i,] = rowMeans(p_hat_array_music[,out$rep_info[[i]]$groups==1,i]) - rowMeans(p_hat_array_music[,out$rep_info[[i]]$groups==2,i])
    }

    # cibersort
    out = readRDS(paste("output/manuscript/real/cibersort/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'.rds',sep=''))
    p_hat_array_cibersort = array(dim = c(K,n_bulk,n_rep))
    for(r in 1:n_rep){
      c.order = match(celltypes,colnames(out$cibersort_fit[[r]][,1:K]))
      p_hat_array_cibersort[,,r] = t(out$cibersort_fit[[r]][,1:K])[c.order,]
    }
    diff_hat_cibersort = matrix(nrow=n_rep,ncol=K)
    for(i in 1:n_rep){
      diff_hat_cibersort[i,] = rowMeans(p_hat_array_cibersort[,out$rep_info[[i]]$groups==1,i]) - rowMeans(p_hat_array_cibersort[,out$rep_info[[i]]$groups==2,i])
    }

    # rna-sieve
    if(add_rnasieve){
      p_hat_sieve = np$load(paste("output/manuscript/real/rnasieve/add_bulk_bias/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_",case,'_p_hat.npy',sep=''))
      for(r in 1:n_rep){
        p_hat_sieve[,,r] = p_hat_sieve[match(celltypes,celltypes_sieve),,r]
      }

      diff_hat_sieve = matrix(nrow=n_rep,ncol=K)
      for(i in 1:n_rep){
        diff_hat_sieve[i,] = rowMeans(p_hat_sieve[,out$rep_info[[i]]$groups==1,i]) - rowMeans(p_hat_sieve[,out$rep_info[[i]]$groups==2,i])
      }
    }else{
      diff_hat_sieve = NA
    }


    diff_hat_true_p = matrix(nrow=n_rep,ncol=K)
    for(i in 1:n_rep){
      diff_hat_true_p[i,] = rowMeans(p_true_array[,out$rep_info[[i]]$groups==1,i]) - rowMeans(p_true_array[,out$rep_info[[i]]$groups==2,i])
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
      cover.naive.rnasieve = NA
    }


    diff_true = matrix(dif,nrow=n_rep,ncol=K,byrow = T)


    diff_res = cbind(sqrt(colMeans(diff_hat_true_p - diff_true)^2),
                       sqrt(colMeans(diff_hat_music- diff_true)^2),
                       sqrt(colMeans(diff_hat_cibersort - diff_true)^2),
                       sqrt(colMeans(diff_hat_sieve - diff_true)^2),
                       sqrt(colMeans(diff_hat_my - diff_true)^2))

    rownames(diff_res) = celltypes
    colnames(diff_res) = c('True.p','MuSiC','CIBERSORT','RNA-Sieve','MEAD')


    #sqrt(colMeans(diff_hat_sieve - diff_true)^2)
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
colnames(res_mat) = c('True.p','MuSiC','CIBERSORT','RNA-Sieve','MEAD')

datax = data.frame(Coverage = c(res_mat),Cell.type = rep(celltypes,5),Method = rep(c('True.p','MuSiC','CIBERSORT','RNA-Sieve','MEAD'),each=6))
datax$Method = factor(datax$Method,levels = unique(datax$Method))
ggplot(datax, aes(x=Method, y=Coverage, group=Cell.type)) +
  geom_hline(yintercept = 0.95,linetype = "dashed")+geom_point(aes(shape=Cell.type,color = Cell.type))+
  theme(panel.grid.minor = element_blank(),axis.text.x = element_text(size = 11),axis.text.y = element_text(size = 11))



# plot(res_mat[,5],ylab='Coverage',ylim=range(res_mat),xaxt = 'n',xlab='cell types')
# abline(h=0.95,lty=2)
# axis(1, at=1:K, labels=celltypes)
# lines(res_mat[,2],type='p',pch=2)
# lines(res_mat[,3],type='p',pch=3)
# lines(res_mat[,4],type='p',pch = 4)
# lines(res_mat[,1],type='p',pch = 5)
# legend('bottomleft',c('mea.err','music','cibersort','rnasieve','truep'),pch=c(1,2,3,4,5))
#

# get diff hat



round(diff_res,3)

datax = data.frame(rmse = c(diff_res),cell.type = rep(celltypes,5),method = rep(c('True.p','MuSiC','CIBERSORT','RNA-Sieve','MEAD'),each=6))
ggplot(datax, aes(x=cell.type, y=rmse, group=method)) +geom_point(aes(shape=method,color = method))+ theme(panel.grid.minor = element_blank())

datax = data.frame(RMSE = c(diff_res),Cell.type = rep(celltypes,5),Method = rep(c('True.p','MuSiC','CIBERSORT','RNA-Sieve','MEAD'),each=6))
datax$Method = factor(datax$Method,levels = unique(datax$Method))


 plot2 =ggplot(datax, aes(x = Method, y = Cell.type, fill = RMSE)) +scale_fill_gradient2(
  low = 'steelblue', high = "red", mid = 'white', midpoint  = quantile(datax$RMSE,0.9), limit = range(datax$RMSE))+
  geom_tile()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=11),axis.text.y = element_text(size = 11),axis.title.x = element_blank(),axis.title.y = element_blank())

 ggplot(datax, aes(x = method, y = cell.type, fill = rmse)) +scale_fill_gradient2(
   low = 'steelblue', high = "red", mid = 'white', midpoint = mean(datax$rmse), limit = range(datax$rmse))+
   geom_tile()

ggplot(datax, aes(x = method, y = cell.type, fill = rmse)) +scale_fill_gradient(
  low = 'white', high = "red", limit = range(datax$rmse))+
  geom_tile()

ggplot(datax, aes(x = method, y = cell.type, fill = rmse)) +scale_fill_gradient(
  low = 'white', high = "steelblue", limit = range(datax$rmse))+
  geom_tile()
## put rmses and diff_rmse together
library(gridExtra)
grid.arrange(plot1, plot2, ncol=2)


#### why RNA-sieve biased??



u_hat = p_hat_sieve[6,,]
dim(u_hat)

u = p_true_array[6,,]

par(mfrow=c(2,3))

plot(u,u_hat,pch='.',xlab='true p of Uneur, 86 bulk and 100 replicates',ylab='RNA-sieve estimated',ylim=c(0,1),xlim=c(0,1))
abline(a=0,b=1)

plot(u[1:43,],u_hat[1:43,],pch='.',xlab='true p of Uneur, 43 bulk from group 1 and 100 replicates',ylab='RNA-sieve estimated',ylim=c(0,1),xlim=c(0,1))
abline(a=0,b=1)

plot(u[44:86,],u_hat[44:86,],pch='.',xlab='true p of Uneur, 43 bulk from group 2 and 100 replicates',ylab='RNA-sieve estimated',ylim=c(0,1),xlim=c(0,1))
abline(a=0,b=1)


sqrt(mean(rowMeans((u-u_hat)^2)))

plot(rowMeans(u),ylim=c(0,0.5),xlab='bulks, first half in group 1, second half in group 2',ylab = 'averaged p-Uneur over 100 replicates')
lines(rowMeans(u_hat),type='p',col=4)
legend('topright',c('use true p','use rna-sieve p_hat'),pch= c(1,1),col=c(1,4))

u_diff = colMeans(u[1:43,]) - colMeans(u[44:86,])
u_diff_hat = colMeans(u_hat[1:43,]) - colMeans(u_hat[44:86,])


plot(u_diff,ylim=c(0,0.4),xlab='replicates',ylab='diff hat of two groups')
lines(u_diff_hat,type='p',col=4)
abline(h = 0.15,col='grey80')
legend('topright',c('use true p','rna sieve','true diff'),pch=c(1,1,NA),lty=c(NA,NA,1),col=c(1,4,'grey80'))




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




