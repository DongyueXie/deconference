plot_simu = function(results){

  G_list = unlist(lapply(results,function(z){z$simu_param$Ng}))

  mse_adj = c()
  mse_unadj = c()

  bias_adj = c()
  bias_unadj = c()

  se_adj = c()
  se_unadj=c()

  coverage_adj = c()
  coverage_unadj = c()

  bias_diff_adj = c()
  bias_diff_unadj = c()
  coverage_diff_adj = c()
  coverage_diff_unadj = c()

  for(i in 1:length(results)){
    mse_adj = rbind(mse_adj,results[[i]]$mse_adj)
    mse_unadj = rbind(mse_unadj,results[[i]]$mse_unadj)

    bias_adj = rbind(bias_adj,results[[i]]$mean_est_adj - results[[i]]$p)
    bias_unadj = rbind(bias_unadj,results[[i]]$mean_est_unadj - results[[i]]$p)

    se_adj = rbind(se_adj,results[[i]]$se_est_adj)
    se_unadj = rbind(se_unadj,results[[i]]$se_est_unadj)

    coverage_adj = rbind(coverage_adj,results[[i]]$covergae_adj)
    coverage_unadj = rbind(coverage_unadj,results[[i]]$covergae_unadj)

    if(!is.null(results[[i]]$covergae_diff_adj)){
      bias_diff_adj = rbind(bias_diff_adj,apply(results[[i]]$diff_adj,2,mean)-(results[[i]]$p-results[[i]]$p2))
      bias_diff_unadj = rbind(bias_diff_unadj,apply(results[[i]]$diff_unadj,2,mean)-(results[[i]]$p-results[[i]]$p2))
      coverage_diff_adj = rbind(coverage_diff_adj,results[[i]]$covergae_diff_adj)
      coverage_diff_unadj = rbind(coverage_diff_unadj,results[[i]]$covergae_diff_unadj)
    }

  }

  # plot mse

  nplot = length(results[[1]]$p)
  par(mfrow=c(1,nplot))

  for(i in 1:nplot){
    plot(G_list,mse_adj[,i],type='l',col=2,ylim=range(c(mse_adj[,i],mse_unadj[,i])),
         xlab='G',ylab='mse',main=paste('cell type',i))
    lines(G_list,mse_unadj[,i],col=4)
    legend('topright',c('adjusted','unajusted'),lty=c(1,1),col=c(2,4))
  }
  title("Mean Squared Error",  line = -1, outer = TRUE)
  #2. bias

  par(mfrow=c(1,nplot))

  for(i in 1:nplot){
    plot(G_list,bias_adj[,i],type='l',col=2,ylim=range(c(bias_adj[,i],bias_unadj[,i])),
         xlab='G',ylab='bias',main=paste('cell type',i))
    lines(G_list,bias_unadj[,i],col=4)
    abline(h=0,lty=2)

  }

  #legend('right',c('adjusted','unajusted'),lty=c(1,1),col=c(2,4))
  title("Bias",  line = -1, outer = TRUE)

  #3. standard error
  par(mfrow=c(1,nplot))

  for(i in 1:nplot){
    plot(G_list,se_adj[,i],type='l',col=2,ylim=range(c(se_adj[,i],se_unadj[,i])),
         xlab='G',ylab='standard error',main=paste('cell type',i))
    lines(G_list,se_unadj[,i],col=4)
    #legend('topright',c('adjusted','unajusted'),lty=c(1,1),col=c(2,4))
  }
  title("Standard Error",  line = -1, outer = TRUE)

  #4. Coverage

  par(mfrow=c(1,nplot))

  for(i in 1:nplot){

    plot(G_list,coverage_adj[,i],type='l',col=2,ylim=c(0,1),
         xlab='G',ylab='coverage',main=paste('cell type',i))
    lines(G_list,coverage_unadj[,i],col=4)
    abline(h=0.95,lty=2)
    #legend('bottomleft',c('adjusted','unajusted'),lty=c(1,1),col=c(2,4))

  }
  title("Coverage",  line = -1, outer = TRUE)

  #4. diff, Coverage of diff

  if(!is.null(results[[1]]$covergae_diff_adj)){

    p = results[[1]]$p
    p2 = results[[1]]$p2

    par(mfrow=c(1,nplot))


    for(i in 1:nplot){

      plot(G_list,bias_diff_adj[,i],type='l',col=2,ylim=range(c(bias_diff_adj[,i],0,bias_diff_unadj[,i])),
           xlab='G',ylab='coverage',main=paste('cell type',i))
      lines(G_list,bias_diff_unadj[,i],col=4)
      abline(h=0,lty=2)
      #legend('bottomleft',c('adjusted','unajusted'),lty=c(1,1),col=c(2,4))

    }
    title("bias of diff",  line = -1, outer = TRUE)

    for(i in 1:nplot){

      plot(G_list,coverage_diff_adj[,i],type='l',col=2,ylim=c(0,1),
           xlab='G',ylab='coverage',main=paste('cell type',i))
      lines(G_list,coverage_diff_unadj[,i],col=4)
      abline(h=0.95,lty=2)
      #legend('bottomleft',c('adjusted','unajusted'),lty=c(1,1),col=c(2,4))

    }
    title("Coverage of diff",  line = -1, outer = TRUE)


  }




}
