
#'@param y response
#'@param X covariates
#'@param b coefficietns
#'@param V error in variable covariance matrix: 1. a K by K matrix, assume the same for each x; 2. a n by K matrix,
#'@param w weights
#'@param type type of sandwich estimator


sandwich_var = function(y,X,V,w=NULL,type='hc',correction=TRUE,a=ncol(X)+4,res.var = FALSE){

  n = nrow(X)


  if(is.null(w)){
    #temp = solve(t(X)%*%X-n*V)
    w = rep(1,n)
  }

    Xw = X*sqrt(w)
    yw = cbind(y*sqrt(w))
    if(nrow(V)==n){
      Vgw = V*w
      Vw = diag(c(colSums(Vgw)))
    }else{
      Vw = V*sum(w)
    }

    A = t(Xw)%*%Xw

  if(correction){

    #browser()

    if(res.var){

      M = crossprod(cbind(yw,Xw))
      B = diag(c(1/colSums(yw),1/diag(Vw)))
      #browser()
      lambda = min(eigen(M%*%B)$values)

    }else{

      m_xy = t(Xw)%*%yw
      lambda = min(eigen((A-m_xy%*%t(m_xy)/sum(yw^2))%*%diag(c(1/diag(Vw))))$values)

    }

    if(lambda>(1+1/n)){
      A = A -(1-a/n)*Vw
    }else{
      A = A - (lambda-1/n-a/n)*Vw
    }
  }else{
    A = A-Vw
  }

  temp = solve(A)

  b_hat = temp%*%t(Xw)%*%yw
  res = c(y - X%*%b_hat)




  #browser()

  if(type=='hc'){
    asyV = temp%*%(t(X)%*%diag(w^2*res^2)%*%X+calc_remainSand(X,y,V,res,w,b_hat))%*%temp
  }
  if(type=='hc2'){
    h = diag(X%*%temp%*%t(X)%*%diag(w))
    #print(h)
    asyV = temp%*%(t(X)%*%diag(w^2*res^2/(1-h))%*%X+calc_remainSand(X,y,V,res,w,b_hat))%*%temp
  }
  if(type=='hc3'){
    h = diag(X%*%temp%*%t(X)%*%diag(w))
    asyV = temp%*%(t(X)%*%diag(w^2*res^2/(1-h)^2)%*%X+calc_remainSand(X,y,V,res,w,b_hat))%*%temp
  }

  return(list(asyV=asyV,b_hat=b_hat))

}

calc_remainSand = function(X,y,V,res,w,b_hat){

  n = nrow(X)
  out = 0
  temp_b = b_hat%*%t(b_hat)
  if(nrow(V)==n){
    for(i in 1:n){
      temp = X[i,]%*%t(b_hat)%*%diag(c(V[i,]))
      out = out + w[i]^2*res[i]*(temp+t(temp)) + w[i]^2*diag(c(V[i,]))%*%temp_b%*%diag(c(V[i,]))
   }
  }else{

    for(i in 1:n){
      temp = X[i,]%*%t(b_hat)%*%V
      out = out + w[i]^2*res[i]*(temp+t(temp)) + w[i]^2*V%*%temp_b%*%V
    }

  }

  out


}


jackknife_var = function(y,X,V,w=NULL){

  n = nrow(X)
  p = ncol(X)

  b_jack = matrix(nrow=n,ncol=p)

  if(is.null(w)){
    for(j in 1:n){
    temp = solve(t(X[-j,])%*%X[-j,] - (n-1)*V)
    b_jack[j,] = temp%*%t(X[-j,])%*%y[-j]
   }
  }else{
    for(j in 1:n){
      temp = solve(t(X[-j,])%*%diag(w[-j])%*%X[-j,] - sum(w[-j])*V)
      b_jack[j,] = temp%*%t(X[-j,])%*%diag(w[-j])%*%y[-j]
    }
  }


  bb = b_jack - rep(1,n)%*%t(apply(b_jack,2,mean))
  covb_jack = t(bb)%*%bb*(n-1)/n
  covb_jack

}
