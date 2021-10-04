

library(Matrix)
library(Rfast)



#'@param X gene by individual matrix
#'@param method testing or thresholding
#'@param alpha FDR level

get_cor_pairs2 = function(X,alpha=0.05,method='testing'){

  p = nrow(X)
  n = ncol(X)



  X.center = scale(t(X),center=TRUE,scale=FALSE)
  S = cova(X.center,center = TRUE)

  # calc S2
  S2 = 0
  for(k in 1:n){
    S2 = S2+(tcrossprod(X.center[k,]))^2
  }


  #browser()
  if(method=='thresholding'){

    Threshmat = 2/n*sqrt(log(p)*(S2+(2-n)*S^2))
    cor.idx = which(abs(S)>=Threshmat,arr.ind = T)
    cor.idx = cor.idx[(cor.idx[,1]!=cor.idx[,2]),]


  }

  if(method=='testing'){

    # for each cell type, get p-value
    Tmat = S*(n-1)/sqrt(S2+(2-n)*S^2)
    P = 2*(1-pnorm(abs(Tmat[lower.tri(Tmat)])))
    P.order = sort(P,decreasing = TRUE)



    nt = length(P.order)
    P.adj = P.order*(p^2-p)/2/(nt:1)

    for(t in 1:nt){
      if(P.adj[t]<=alpha){
        break
      }
    }

    bp = 2*(1-pnorm(sqrt(4*log(p)-2*log(log(p)))))
    if(P.order[t]<bp){
      thresh = 2*(1-pnorm(sqrt(4*log(p))))
    }else{
      thresh = P.order[t]
    }

    P.rej = c()
    P.rej[P>thresh] = 0
    P.rej[P<=thresh] = 1
    P.mat = Matrix(0,nrow=p,ncol=p,sparse = T)
    P.mat[lower.tri(P.mat)] = P.rej
    cor.idx = which(P.mat!=0,arr.ind = T)
    cor.idx = rbind(cor.idx,cbind(cor.idx[,2],cor.idx[,1]))

  }


  if(length(cor.idx)==0){
    cor.idx = NULL
  }

  cor.idx

}




get_cor_pairs = function(X_array,alpha=0.05){

  p = dim(X_array)[1]
  n = dim(X_array)[3]
  K = dim(X_array)[2]

  # for each cell type, get p-value
  P = 0
  for(k in 1:K){
    P = P + 1/K*tan((0.5-get_cor_p(t(X_array[,k,])))*pi)
  }
  P = 0.5 - atan(P)/pi

  P.order = sort(P,decreasing = TRUE)



  nt = length(P.order)
  P.adj = P.order*(p^2-p)/2/(nt:1)

  for(t in 1:nt){
    if(P.adj[t]<=alpha){
      break
    }
  }

  bp = 2*(1-pnorm(sqrt(4*log(p)-2*log(log(p)))))
  if(P.order[t]<bp){
    thresh = 2*(1-pnorm(sqrt(4*log(p))))
  }else{
    thresh = P.order[t]
  }

  P.rej = c()
  P.rej[P>thresh] = 0
  P.rej[P<=thresh] = 1
  P.mat = Matrix(0,nrow=p,ncol=p,sparse = T)
  P.mat[lower.tri(P.mat)] = P.rej
  cor.idx = which(P.mat!=0,arr.ind = T)
  cor.idx = rbind(cor.idx,cbind(cor.idx[,2],cor.idx[,1]))
  cor.idx
}



# returns a vector of the lower tri mat

get_cor_p = function(X){

  n = nrow(X)
  p = ncol(X)

  X.center = scale(X,center=TRUE,scale=FALSE)
  S = cova(X.center,center = TRUE)

  # calc S2

  S2 = 0
  for(k in 1:n){
    S2 = S2+(tcrossprod(X.center[k,]))^2
  }

  # calc T statistics

  Tmat = S*(n-1)/sqrt(S2+(2-n)*S^2)

  2*(1-pnorm(abs(Tmat[lower.tri(Tmat)])))

  # diag(Tmat) = 0
  # Tmat[upper.tri(Tmat)] = 0
  # Tmat[lower.tri(Tmat)] = 2*(1-pnorm(Tmat[lower.tri(Tmat)] ))
  # Tmat = Matrix::Matrix(Tmat,sparse = T)



}


