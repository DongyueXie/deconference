##########################
#########################
# meta mode: all
##########################
##########################


#'@title Q-statistics - (n-1) in Meta-analysis
#'@param sigma2 random effect variance
#'@param x effects of studies, a vector
#'@param v variance of x, a vector
#'@return Q value - (n-1)

Q_stats_vector = function(sigma2,x,v){
  sum(1/(sigma2+v)*(x-(sum(x/(v+sigma2)))/(sum(1/(sigma2+v))))^2) - (length(x)-1)
}


#'@title PM meta-analysis random effecr estimator
#'@param x effects of studies, a vector
#'@param v variance of x, a vector
#'@param ub upper bound of random effect variance sigma^2
#'@return estimated random effect variance sigma^2, mu_hat and var(mu_hat)
PMmeta_vector = function(x,v,ub=var(x,na.rm = TRUE),tol = .Machine$double.eps^0.25){

  l = length(x)
  # do not use v=0
  rm.v = which(v==0|is.na(v))

  if(length(rm.v)<(l-1)){
    if(length(rm.v)!=0){
      x = x[-rm.v]
      v = v[-rm.v]
    }


    Q0 = Q_stats_vector(0,x,v)
    Qub = Q_stats_vector(ub,x,v)
    if(Q0>0){
      if(Qub<0){
        sigma2 = uniroot(Q_stats_vector,c(0,ub),x=x,v=v)$root
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

    if(length(rm.v)!=0){
      w.temp = rep(0,l)
      w.temp[-rm.v] = w
    }else{
      w.temp = w
    }

    #print(length(w.temp))

    return(list(mu_hat=mu_hat,var_mu_hat=var_mu_hat,w=w.temp))
  }else{
    return(list(mu_hat=mean(x,na.rm=TRUE),var_mu_hat=mean(v+sigma2,na.rm=TRUE)/l,w=rep(1/l,l)))
  }

}





##########################
#########################
# meta mode: one
############################
##########################

#'@param X array, G bu K by I
#'@param V array.

Q_stats_array = function(sigma2,X,V){

  GKI = dim(X)
  G = GKI[1]
  K = GKI[2]
  I = GKI[3]
  Q = 0
  N = prod(GKI)
  df_lost = G*K
  for(g in 1:G){
    for(k in 1:K){
      x = X[g,k,]
      v = V[g,k,]

      l = I
      # do not use v=0
      rm.v = which(v==0|is.na(v))
      if(length(rm.v)!=0){
        x = x[-rm.v]
        v = v[-rm.v]
        N = N - length(rm.v)
      }

      if(length(x)!=0){
        Q = Q + sum(1/(sigma2+v)*(x-(sum(x/(v+sigma2)))/(sum(1/(sigma2+v))))^2)
      }else{
        df_lost = df_lost - 1
      }
    }
  }

  Q - (N-df_lost)

}


PMmeta_array = function(X,V,ub=var(as.vector(X),na.rm = TRUE),tol=.Machine$double.eps^0.25){

  Q0 = Q_stats_array(0,X,V)
  Qub = Q_stats_array(ub,X,V)

  if(Q0>0){
    if(Qub<0){
      sigma2 = uniroot(Q_stats_array,c(0,ub),X=X,V=V,tol=tol)$root
    }else{
      sigma2 = ub
    }
  }else{
    sigma2 = 0
  }
  sigma2
}



# test PMmeta_one

# G = 100
# K = 4
# I = 6
# GKI = G*K*I
#
# V = array(runif(GKI,0,0.5),dim=c(G,K,I))
#
# M = matrix(rnorm(GKI),ncol=K)
#
# X = array(dim = c(G,K,I))
#
# for(i in 1:I){
#   X[,,i] = matrix(rnorm(G*K,as.vector(M),sqrt(as.vector(V[,,i])+0.5)),ncol=K)
# }
#
# PMmeta_array(X,V,100)
# PMmeta_matrix(X[,1,],V[,1,],100)
# PMmeta_matrix(X[1,,],V[1,,],100)
#
# sigma_hat = matrix(nrow=G,ncol=K)
# for(g in 1:G){
#   for(k in 1:K){
#     sigma_hat[g,k] = PMmeta_vector(X[g,k,],V[g,k,],100)
#   }
# }
# mean(sigma_hat)






####################################
####################################
##### meta mode' by cell type / gene ######
####################################
###################################

#'@param X matrix G by I
#'@param V array.
#'
Q_stats_matrix = function(sigma2,X,V){

  GI = dim(X)
  G = GI[1]
  I = GI[2]
  Q = 0
  N = prod(GI)
  df_lost = G
  for(g in 1:G){

    x = X[g,]
    v = V[g,]

    l = I
    # do not use v=0
    rm.v = which(v==0|is.na(v))
    if(length(rm.v)!=0){
      x = x[-rm.v]
      v = v[-rm.v]
      N = N - length(rm.v)
    }

    if(length(x)!=0){
      Q = Q + sum(1/(sigma2+v)*(x-(sum(x/(v+sigma2)))/(sum(1/(sigma2+v))))^2)
    }else{
      df_lost = df_lost - 1
    }

  }

  Q - (N-df_lost)

}

PMmeta_matrix = function(X,V,ub=var(as.vector(X),na.rm = TRUE),tol = .Machine$double.eps^0.25){

  Q0 = Q_stats_matrix(0,X,V)
  Qub = Q_stats_matrix(ub,X,V)

  if(Q0>0){
    if(Qub<0){
      sigma2 = uniroot(Q_stats_matrix,c(0,ub),X=X,V=V)$root
    }else{
      sigma2 = ub
    }
  }else{
    sigma2 = 0
  }
  sigma2
}



calc_cov = function(x1,x2,w1,w2,mu1,mu2){

  temp = w1*w2*(x1-mu1)*(x2-mu2)

  out = sum(temp,na.rm = TRUE)

  out


}

