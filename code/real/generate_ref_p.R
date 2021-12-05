

# This function generate ref index and bulk proportions p

generate_ref_p = function(dirichlet.scales,
                          cases,
                          n_bulks,
                          n_refs,
                          n_indi = 97,
                          n_rep=100,
                          K = 6,
                          celltypes = NULL,
                          p1.prop = 0.5,
                          p_null=c(0.3,0.2,0.15,0.15,0.1,0.1),
                          p1_diff=c(0.15,0.15,0.1,0.1,0.2,0.3),
                          p2_diff=c(0.1,0.1,0.2,0.3,0.15,0.15)){
  library(reticulate)
  library(gtools)
  #n_bulk = n_indi - n_ref

  for(aa in dirichlet.scales){
    for(case in cases){

      if(case=='null'){
        p1 = p_null
        p2 = p_null
      }else if(case=='all_diff'){
        p1 = p1_diff
        p2 = p2_diff
      }

      for(n_bulk in n_bulks){
        n.p1.temp =  round(n_bulk*p1.prop)
        for(n_ref in n_refs){



          print(paste('Running:',case,'dirichlet.scale=',aa,"n_bulk=",n_bulk,'n_ref=',n_ref))

          #out = readRDS(paste("output/manuscript/real/neuron_ref",n_ref,"_rep",n_rep,"_dirichlet",aa,"_corfdr005_",case,'.rds',sep=''))

          ref_idx_mat=matrix(nrow=n_rep,ncol=n_ref)
          bulk_idx_mat=matrix(nrow=n_rep,ncol=n_bulk)
          bulk_p_array=array(dim=c(K,n_bulk,n_rep))
          groups_mat = matrix(nrow=n_rep,ncol=n_bulk)
          for(i in 1:n_rep){
            ref.idx = sample(1:n_indi,n_ref)
            ref_idx_mat[i,] = ref.idx
            bulk_idx_mat[i,] = sample((1:n_indi)[-ref.idx],n_bulk,replace = ifelse(n_bulk>(n_indi-n_ref),TRUE,FALSE))
            groups = c(rep(1,n.p1.temp),rep(2,n_bulk-n.p1.temp))
            if(aa!=0){
              b = cbind(t(rdirichlet(n.p1.temp,p1*aa)),t(rdirichlet(n_bulk-n.p1.temp,p2*aa)))
            }else{
              b = cbind(p1%*%t(rep(1,n.p1.temp)),p2%*%t(rep(1,n_bulk-n.p1.temp)))
            }

            rownames(b) = celltypes
            colnames(b) = 1:n_bulk
            #perm.temp = sample(1:n_bulk)
            #groups = groups[perm.temp]
            #b = b[,perm.temp]
            bulk_p_array[,,i] = b
            groups_mat[i,] = groups
          }


          saveRDS(ref_idx_mat,paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_ref_idx.rds',sep=''))
          saveRDS(bulk_p_array,paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_bulk_p.rds',sep=''))
          saveRDS(groups_mat,paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_groups_mat.rds',sep=''))
          saveRDS(bulk_idx_mat,paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_bulk_idx_mat.rds',sep=''))

          ref_idx_mat = r_to_py(ref_idx_mat)
          py_save_object(ref_idx_mat,filename=paste("data/neuron/real_manu/rnasieve/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_",case,'_ref_idx',sep=''))

          bulk_p_array = r_to_py(bulk_p_array)
          py_save_object(bulk_p_array,filename=paste("data/neuron/real_manu/rnasieve/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_",case,'_bulk_p',sep=''))

          groups_mat = r_to_py(groups_mat)
          py_save_object(groups_mat,filename=paste("data/neuron/real_manu/rnasieve/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_",case,'_groups_mat',sep=''))

          bulk_idx_mat = r_to_py(bulk_idx_mat)
          py_save_object(bulk_idx_mat,filename=paste("data/neuron/real_manu/rnasieve/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_corfdr005_",case,'_bulk_idx_mat',sep=''))

        }
      }
    }
  }
}

indis_ref = readRDS('data/neuron/indis_ref_12400by6by97.rds')
celltypes = dimnames(indis_ref)[[2]]
set.seed(12345)
generate_ref_p(c(5,20,50),c("null","all_diff"),celltypes = celltypes,n_bulks = c(86,500),n_refs=11)
