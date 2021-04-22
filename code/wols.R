#'@param W variance, either a length G vector or a G*n_bulk matrix
#'@param b beta for evaluating weights
#'@param X reference matrix
#'@param Y bulk data
#'@param Vg variance
#'@param w.mode "equal", "res", "res+ref_var", "ref_var", "res+bvb", "bvb"

wols = function(X,Y,
                Vg=NULL,
                b=NULL,
                w.mode = "res",
                adj.v = FALSE,
                nu = 1e-4,
                nu.q = 0.05,
                marker = NULL,
                scale.y = FALSE,
                X_array = NULL,
                hc.type = 'hc0',
                use.weight.for.var=TRUE
){
  if(!is.null(marker)){
    marker.x.idx = match(marker,rownames(X))
    marker.x.idx = marker.x.idx[complete.cases(marker.x.idx)]
    X = X[marker.x.idx,]
    marker.y.idx = match(marker,rownames(Y))
    marker.y.idx = marker.y.idx[complete.cases(marker.y.idx)]
    Y = Y[marker.y.idx,]
    marker.vg.idx = match(marker,rownames(Vg))
    marker.vg.idx = marker.vg.idx[complete.cases(marker.vg.idx)]
    Vg = Vg[marker.vg.idx,]
    if(!is.null(X_array)){
      marker.array.idx = match(marker,rownames(X_array))
      marker.array.idx = marker.array.idx[complete.cases(marker.array.idx)]
      X_array = X_array[marker.array.idx,,]
    }
  }

  ## remove rows of X that are all 0

  rm0 = which(rowSums(X)==0)
  if(length(rm0)>0){
    X = X[-rm0,]
    Y = Y[-rm0,]
    Vg = Vg[-rm0,]
    X_array = X_array[-rm0,,]
  }


  n = nrow(X)
  K = ncol(X)
  ni = ncol(Y)

  if(scale.y){
    Y = apply(Y,2,function(z){z/sum(z)*100})
  }




  # if(is.null(b)){
  #   # then use b as
  #   b = c(rowMeans(solve(t(X)%*%X)%*%t(X)%*%Y))
  # }

  # Sigma_sample = c()
  # for(i in 1:dim(X_array)[3]){
  #   Sigma_sample = cbind(Sigma_sample,(X_array[,,i] - X)%*%b)
  # }
  # w = diag.cov(t(Sigma_sample))
  # W = w%*%t(rep(1,ni))

  # if(is.null(dim(b))){
  #     b = b%*%t(rep(1,ni))
  #  }
  #
  # r = Y - X%*%b

  beta_hat = matrix(nrow=ni,ncol=K)
  beta_hat_se = matrix(nrow=ni,ncol=K)
  weight = matrix(nrow=n,ncol=ni)

  var.res = matrix(nrow=n,ncol=ni)
  var.ref = matrix(nrow=n,ncol=ni)


  for(j in 1:ni){

    if(sum(Y[,j]==0)>0){
      idx = which(Y[,j]!=0)
    }else{
      idx = 1:n
    }

    yj = Y[idx,j]
    Xj = X[idx,]
    Vgj = Vg[idx,]
    X_arrayj = X_array[idx,,]

    if(is.null(b)){
      # use ols estimate
      bj = solve(t(Xj)%*%Xj)%*%t(Xj)%*%yj
      bj = pmax(bj,0)
    }else{
      bj = b[,j]
    }


    Sigma_sample = c()
    for(i in 1:dim(X_arrayj)[3]){
      Sigma_sample = cbind(Sigma_sample,(X_arrayj[,,i] - Xj)%*%bj)
    }
    w = diag.cov(t(Sigma_sample))
    r = c(yj - Xj%*%bj)

    var.res[idx,j] = r^2
    var.ref[idx,j] = w

    if(w.mode=="equal"){
      wj = rep(1,length(idx))
    }

    if(w.mode=="res"){
      if(is.null(nu)){
        nu = quantile(r^2,nu.q,na.rm=TRUE)
      }
      wj = 1/(nu + r^2)
    }
    if(w.mode=="res+ref_var"){
      if(is.null(nu)){
        nu = quantile(r^2+w,nu.q,na.rm=TRUE)
      }
      wj = 1/(nu + r^2+w)
    }
    if(w.mode=="ref_var"){
      if(is.null(nu)){
        nu = quantile(w,nu.q,na.rm=TRUE)
      }
      wj = 1/(nu +w)
    }



    if(w.mode=="res+bvb"){
      noise_var = apply(Vgj,1,function(v){
        if(length(v)==K^2){
          v = matrix(v,ncol = K)
        }else if(length(v)==K){
          v = diag(v)
        }
        t(bj)%*%v%*%bj
      })
      if(is.null(nu)){
        nu = quantile(r^2 + c(noise_var),nu.q,na.rm=TRUE)
      }
      wj = 1/(nu + r^2 + c(noise_var))
    }
    if(w.mode=="bvb"){
      noise_var = apply(Vgj,1,function(v){
        if(length(v)==K^2){
          v = matrix(v,ncol = K)
        }else if(length(v)==K){
          v = diag(v)
        }
        t(bj)%*%v%*%bj
      })
      if(is.null(nu)){
        nu = quantile(c(noise_var),nu.q,na.rm=TRUE)
      }
      wj = 1/(nu + c(noise_var))
    }

    wj = wj/sum(wj)
    weight[idx,j] = wj

    #w.temp = 1/(W[j,]+nu)
    #w.temp = w.temp/sum(w.temp)
    #print(cor(wj,w.temp))

    Xw = Xj*sqrt(wj)
    yw = yj*sqrt(wj)

    A = t(Xw)%*%Xw
    if(adj.v){
      Vw = Vgj*wj
      d_Vg = ncol(Vgj)
      if(d_Vg==K^2){
        V = matrix(c(colSums(Vw)),ncol = K)
      }else if(d_Vg==K){
        V = diag(c(colSums(Vw)))
      }else{
        stop('check dimension of Vg')
      }
    }else{
      V = 0
    }

    A_inv = solve(A-V)
    bhat = A_inv%*%t(Xw)%*%yw
    bhat = pmax(bhat,0)

    ## calculate covariance matrix

    if(adj.v){

      if(d_Vg==K^2){
        Vbi = t(apply(Vw,1,function(z){v = matrix(z,ncol=K);v%*%bhat}))
      }
      if(d_Vg==K){
        Vbi = (Vw)%*%diag(c(bhat))
      }

    }else{
      Vbi = 0
    }



    #Hi = t(t(X%*%A_inv%*%t(X))*w)
    hi = rowSums((Xj%*%A_inv)*Xj)*wj
    ri = yj - Xj%*%bhat

    if(use.weight.for.var){

      if(hc.type == 'hc0'){
        Sigma_j = crossprod(c(ri)*wj*Xj+Vbi)
      }else if(hc.type == 'hc2'){
        Sigma_j = crossprod(c(ri)/sqrt(1-pmax(pmin(hi,1-1/n),0))*wj*Xj+Vbi)
      }else if(hc.type == 'hc3'){
        Sigma_j = crossprod(c(ri)/(1-pmax(pmin(hi,1-1/n),0))*wj*Xj+Vbi)
      }

      covb = A_inv%*%Sigma_j%*%A_inv

    }else{


      A = t(Xj)%*%Xj
      if(adj.v){
        d_Vg = ncol(Vgj)
        if(d_Vg==K^2){
          V = matrix(c(colSums(Vgj)),ncol = K)
        }else if(d_Vg==K){
          V = diag(c(colSums(Vgj)))
        }else{
          stop('check dimension of Vg')
        }
      }else{
        V = 0
      }

      A_inv = solve(A-V)


      if(adj.v){

        if(d_Vg==K^2){
          Vbi = t(apply(Vgj,1,function(z){v = matrix(z,ncol=K);v%*%bhat}))
        }
        if(d_Vg==K){
          Vbi = (Vgj)%*%diag(c(bhat))
        }

      }else{
        Vbi = 0
      }


      if(hc.type == 'hc0'){
        Sigma_j = crossprod(c(ri)*Xj+Vbi)
      }else if(hc.type == 'hc2'){
        Sigma_j = crossprod(c(ri)/sqrt(1-pmax(pmin(hi,1-1/n),0))*Xj+Vbi)
      }else if(hc.type == 'hc3'){
        Sigma_j = crossprod(c(ri)/(1-pmax(pmin(hi,1-1/n),0))*Xj+Vbi)
      }

      covb = A_inv%*%Sigma_j%*%A_inv

    }





    J = J_sum2one(bhat,K)


    asyV = (J)%*%covb%*%t(J)

    beta_hat_se[j,] = sqrt(diag(asyV))
    bhat = bhat/sum(bhat)
    beta_hat[j,] = bhat
  }

  rownames(weight) = rownames(X)
  rownames(var.res) = rownames(X)
  rownames(var.ref) = rownames(X)

  list(beta_hat=beta_hat,
       beta_hat_se=beta_hat_se,
       weight=weight,
       X=X,
       Y=Y,
       Vg=Vg,
       X_array=X_array,
       var.res=var.res,
       var.ref=var.ref)
}
