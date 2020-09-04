


ploter = function(results){

  par(mfrow=c(2,2))

  plot(log(results$mse_adj),type='p',ylim = range(log(c(results$mse_adj,results$mse_unadj))),col=2,ylab='log(mse)',xlab='coefs',main='mse')
  lines(log(results$mse_unadj),type='p',pch=2,col=4)


  plot(results$covergae_adj,ylim = range(c(results$covergae_adj,results$covergae_unadj,results$covergae_unadj_cv,1)),
       col=2,type='b',ylab='coverage',xlab='coefs',main='coverage of p')
  lines(results$covergae_unadj,type='b',pch=2,col=4)
  lines(results$covergae_unadj_cv,type='b',pch=2,col=3)
  abline(h=0.95,lty=2)

  plot(results$covergae_diff_adj,ylim = range(c(results$covergae_diff_adj,results$covergae_diff_unadj,results$covergae_diff_unadj_cv,1)),
       col=2,type='b',ylab='coverage',xlab='coefs',main='coverage of difference')
  lines(results$covergae_diff_unadj,type='b',pch=2,col=4)
  lines(results$covergae_diff_unadj_cv,type='b',pch=2,col=3)
  abline(h=0.95,lty=2)

  p_order = order(abs(results$p-results$p2))

  plot(results$power_adj[p_order],ylim = range(c(results$power_adj,results$power_unadj,results$power_unadj_cv)),
       col=2,ylab="power",xlab='',main='power')
  lines(results$power_unadj[p_order],type='p',pch=2,col=4)
  lines(results$power_unadj_cv[p_order],type='p',pch=2,col=3)

  legend('topleft',c('adjusted','unadj_hc0',"unadj_lm"),col=c(2,4,3),pch=c(1,2,2))

  par(mfrow=c(1,1))

}



ploter_multibulk = function(results){

  par(mfrow=c(2,2))

  plot(log(results$mse_adj),type='p',ylim = range(log(c(results$mse_adj,results$mse_unadj))),col=2,ylab='log(mse)',xlab='coefs',main='mse')
  lines(log(results$mse_unadj),type='p',pch=2,col=4)


  plot(results$covergae_adj,ylim = range(c(results$covergae_adj,results$covergae_unadj,results$covergae_unadj_cv,1)),
       col=2,type='b',ylab='coverage',xlab='coefs',main='coverage of p')
  lines(results$covergae_unadj,type='b',pch=2,col=4)
  lines(results$covergae_unadj_cv,type='b',pch=2,col=3)
  abline(h=0.95,lty=2)

  plot(results$covergae_diff_adj,ylim = range(c(results$covergae_diff_adj,results$covergae_diff_unadj,results$covergae_diff_unadj_cv,1)),
       col=2,type='b',ylab='coverage',xlab='coefs',main='coverage of difference')
  lines(results$covergae_diff_unadj,type='b',pch=2,col=4)
  lines(results$covergae_diff_unadj_cv,type='b',pch=2,col=3)
  abline(h=0.95,lty=2)

  p_order = order(abs(results$p-results$p2))

  plot(results$power_adj[p_order],ylim = range(c(results$power_adj,results$power_unadj,results$power_unadj_cv)),
       col=2,ylab="power",xlab='',main='power')
  lines(results$power_unadj[p_order],type='p',pch=2,col=4)
  lines(results$power_unadj_cv[p_order],type='p',pch=2,col=3)

  legend('topleft',c('adjusted','unadj_hc0',"unadj_lm"),col=c(2,4,3),pch=c(1,2,2))

  par(mfrow=c(1,1))

}


