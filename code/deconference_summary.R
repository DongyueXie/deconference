

two_group_test = function(out, groups){

  if(is.factor(groups)){
    group_name = levels(groups)
  }else{
    group_name = levels(as.factor(groups))
  }

  group1_idx = which(groups==group_name[1])
  a = c()
  a[group1_idx] = 1/length(group1_idx)

  group2_idx = which(groups==group_name[2])
  a[group2_idx] = -1/length(group2_idx)

  K = ncol(out$input$X)

  diff_group = rowMeans(out$beta_hat[,group1_idx,drop=F]) - rowMeans(out$beta_hat[,group2_idx,drop=F])

  # N_indi = length(group1_idx) + length(group2_idx)

  V_tilde = 0

  idx = c(group1_idx,group2_idx)
  for(i in idx){
    for(j in idx){
      V_tilde = V_tilde + a[i]*a[j]*out$cov_beta_hat[((i-1)*K+1):(i*K),((j-1)*K+1):(j*K)]
    }
  }

  z_score = diff_group/sqrt(diag(V_tilde))

  p_value = (1-pnorm(abs(z_score)))*2

  return(list(p_value = p_value,z_score=z_score,V_tilde=V_tilde,
              diff_group=diff_group,diff_se = sqrt(diag(V_tilde)),a=a))

}



build_ci = function(out,alpha=0.05){

  ci_l = pmax(out$beta_hat - qnorm(1-alpha/2)*out$beta_se,0)
  ci_r = pmin(out$beta_hat + qnorm(1-alpha/2)*out$beta_se,1)

  ci = list(ci_l = ci_l,ci_r=ci_r)

  ci

}
