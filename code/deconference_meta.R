#'@title Q-statistics - (n-1) in Meta-analysis
#'@param sigma2 random effect variance
#'@param x effects of studies, a vector
#'@param v variance of x, a vector
#'@return Q value - (n-1)

Q_stats = function(sigma2,x,v){
  sum(1/(sigma2+v)*(x-(sum(x/(v+sigma2)))/(sum(1/(sigma2+v))))^2) - (length(x)-1)
}

#'@title PM meta-analysis random effecr estimator
#'@param x effects of studies, a vector
#'@param v variance of x, a vector
#'@param ub upper bound of random effect variance sigma^2
#'@return estimated random effect variance sigma^2, mu_hat and var(mu_hat)
PMmeta = function(x,v,ub=0.25,sigma2=NULL,meta_var='plug_in'){

  l = length(x)
  # do not use v=0
  rm.v = which(v==0)


  if(length(rm.v)<(l-1)){
    if(length(rm.v)!=0){
      x = x[-rm.v]
      v = v[-rm.v]
    }

    if(is.null(sigma2)){
      Q0 = Q_stats(0,x,v)
      Qub = Q_stats(ub,x,v)
      if(Q0>0){
        if(Qub<0){
          sigma2 = uniroot(Q_stats,c(0,ub),x=x,v=v)$root
        }else{
          sigma2 = ub
        }
      }else{
        sigma2 = 0
      }
    }


    w = 1/(v+sigma2)
    w = w/sum(w)
    mu_hat = sum(w*x)
    #browser()
    if(meta_var=='adjust'){
      var_mu_hat = sum(w^2/(1-w)*(x-mu_hat)^2)
    }
    if(meta_var=='plug_in'){
      var_mu_hat = 1/sum(1/(sigma2+v))
    }

    return(list(sigma2=sigma2,mu_hat=mu_hat,var_mu_hat=var_mu_hat))
  }else{

    if(is.null(sigma2)){
      return(list(sigma2=0,mu_hat=mean(x),var_mu_hat=mean(v)/l))
    }else{
      return(list(sigma2=sigma2,mu_hat=mean(x),var_mu_hat=mean(v+sigma2)/l))
    }


  }


}
