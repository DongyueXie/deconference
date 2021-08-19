rmse = function(x,y){sqrt(mean((x-y)^2))}


summary_neuron = function(out,out.music=NULL){
  rmse_ols = c()
  rmse_err = c()
  rmse_music = c()
  coverage = c()
  median_std = c()
  wald= list()

  for(i in 1:length(out)){
    rmse_ols[i]=rmse(out[[i]]$fit.ols$beta_hat,out[[i]]$input$b)
    rmse_err[i] = rmse(out[[i]]$fit.err.hc0$beta_hat,out[[i]]$input$b)
    if(!is.null(out.music)){
      rmse_music[i] = rmse(out.music[[i]],out[[i]]$input$b)
    }

    waldi = list()
    waldi[[1]] = (out[[i]]$fit.ols$beta_hat-out[[i]]$input$b)/out[[i]]$fit.ols$ols.out$beta_se
    waldi[[2]] = (out[[i]]$fit.ols$beta_hat-out[[i]]$input$b)/out[[i]]$fit.ols$sand.out$beta_se
    waldi[[3]] = (out[[i]]$fit.ols$beta_hat-out[[i]]$input$b)/out[[i]]$fit.ols$sand.out.hc3$beta_se
    waldi = c(waldi,lapply(2:5,function(j){(out[[i]][[j]]$beta_hat-out[[i]]$input$b)/out[[i]][[j]]$beta_se}))
    wald[[i]] = waldi
    coverage = rbind(coverage,unlist(lapply(waldi,function(z){mean(abs(z)<=1.96,na.rm = T)})))
    median_std = rbind(median_std,c(median(c(out[[i]]$fit.ols$ols.out$beta_se)),
                                    median(c(out[[i]]$fit.ols$sand.out$beta_se)),
                                    median(c(out[[i]]$fit.ols$sand.out.hc3$beta_se)),
                                    unlist(lapply(2:5,function(j){median(c(out[[i]][[j]]$beta_se),na.rm = T)}))))
  }
  colnames(coverage) = c('ols.cv','ols.hc0','ols.hc3','err.hc0','err.hc3','err.cor.hc0','err.cor.hc3')
  colnames(median_std)  = c('ols.cv','ols.hc0','ols.hc3','err.hc0','err.hc3','err.cor.hc0','err.cor.hc3')

  return(list(rmse_ols=rmse_ols,
  rmse_err=rmse_err,
  rmse_music=rmse_music,

  coverage=coverage,
  median_std=median_std))

}

