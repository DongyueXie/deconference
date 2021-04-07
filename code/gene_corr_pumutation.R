

##

nreps = 1000

X_xin <- readRDS("data/pancreas/X_xin.rds")
X_array = X_xin$X_array
X_array = X_array[, , -c(4, 7)]

N = dim(X_array)[1]
K = dim(X_array)[2]
NI = dim(X_array)[3]







library(coop)

P = list()

for(k in 1:K){
  print(paste('working on cell type', k))

  gene_mat_k = X_array[,k,]
  cor_k = abs(pcor(t(gene_mat_k)))

  I_k = 0
  for(r in 1:nreps){
    if(r%%100==0){
      print(paste('done',r,'out of',nreps))
    }
    gene_mat_k_perm = apply(gene_mat_k,1,sample,size=NI)
    cor_k_rand = pcor(gene_mat_k_perm)
    I_k = I_k + (abs(cor_k_rand)>cor_k)
  }
  P[[k]] = I_k/nreps
}


Pmin = matrix(nrow=N,ncol=N)
for(i in 1:N){
  for(j in (i+1):N){
    Pmin[i,j] = min(c((P[[1]])[i,j],(P[[2]])[i,j],(P[[3]])[i,j]))
  }
}


Pave = matrix(nrow=N,ncol=N)
for(i in 1:N){
  for(j in (i+1):N){
    Pave[i,j] = mean(c((P[[1]])[i,j],(P[[2]])[i,j],(P[[3]])[i,j]))
  }
}


# remove some genes

rm.gene = which(rowSums(rowSums(X_array,dims=2))==0)

for(k in 1:K){
  gene_mat_k = X_array[,k,]
  rm.gene = c(rm.gene,which(rowSums(gene_mat_k)==0))
}


P.min = Pmin[-rm.gene,-rm.gene]
X_array = X_array[-rm.gene,,]
P.adjusted = matrix(p.adjust(P.min,method = 'BH'),nrow=dim(P.min)[1])

out = corr_matrix_prune(P.adjusted,n_var = 100,decreasing = TRUE)

simu_result(X_array[out,,],c(0.2,0.3,0.5))

rm.idx = which(P.adjusted<0.9,arr.ind = T)

# will remove all genes

rm.idx = unique(c(rm.idx))

##
P.adjusted[is.na(P.adjusted)] = 0
P.adjusted = (P.adjusted + t(P.adjusted))
diag(P.adjusted) = rep(1,dim(P.adjusted)[1])
Reject.mat = 1*(P.adjusted<0.05)
out = corr_matrix_prune(Reject.mat,n_var = 100)
simu_result(X_array[out,,],c(0.2,0.3,0.5))


out = corr_matrix_prune(Reject.mat,n_var = 1000)
simu_result(X_array[out,,],c(0.2,0.3,0.5))


##

P.ave = Pave[-rm.gene,-rm.gene]
P.ave.adjusted = matrix(p.ave.adjust(P.ave,method = 'BH'),nrow=dim(P.ave)[1])

P.ave.adjusted[is.na(P.ave.adjusted)] = 0
P.ave.adjusted = (P.ave.adjusted + t(P.ave.adjusted))
diag(P.ave.adjusted) = rep(1,dim(P.ave.adjusted)[1])
Reject.mat.ave = 1*(P.ave.adjusted<0.05)
hist(rowSums(Reject.mat.ave),breaks = 1000)



P = matrix(nrow=N,ncol = N)

for(i in 1:N){

  if(i%%10==0){
    print(paste('looking at gene',i))
  }

  xi = X_array[i,,]
  for(j in (i+1):N){

    print(j)


    xj = X_array[j,,]

    cor_ij = c()
    for(k in 1:K){
      cor_ij[k] = cor(xi[k,],xj[k,])
    }
    cor_ij = max(cor_ij)

    cor_ij_random = rep(0,nreps)
    for(r in 1:nreps){

      cor_ij_random_k = c()
      for(k in 1:K){
        cor_ij_random_k[k] = cor(xi[k,],sample(xj[k,],NI))
      }
      cor_ij_random[r] = max(cor_ij_random_k)

    }

    P[i,j] = sum(abs(cor_ij_random)>abs(cor_ij))/nreps
  }
}
