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
PMmeta = function(x,v,ub=0.25){

  l = length(x)
  # do not use v=0
  rm.v = which(v==0|is.na(v))

  if(length(rm.v)<(l-1)){
    if(length(rm.v)!=0){
      x = x[-rm.v]
      v = v[-rm.v]
    }


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

  }else{
    sigma2 = 0
  }

  sigma2

}



meta_analysis = function(x,v,sigma2,meta_var='plug_in'){

  l = length(x)
  # do not use v=0, and v=NA
  rm.v = which(v==0|is.na(v))


  if(length(rm.v)<(l-1)){
    if(length(rm.v)!=0){
      x = x[-rm.v]
      v = v[-rm.v]
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

    return(list(mu_hat=mu_hat,var_mu_hat=var_mu_hat))
  }else{

    return(list(mu_hat=mean(x,na.rm=TRUE),var_mu_hat=mean(v+sigma2,na.rm=TRUE)/l))

  }


}
