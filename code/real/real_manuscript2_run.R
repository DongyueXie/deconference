

#################2021/12/03##########################
#####################################################
#####################################################
# try different cor fdr level

source('code/real/real_manuscript2.R')
indis_ref = readRDS('data/neuron/indis_ref_12400by6by97.rds')


n_rep = 100
dirichlet.scales = c(5)
n_bulks = c(86)
n_refs = 11
cases = c('null','all_diff')
cor_fdrs = c(0.1,0.3,0.5)
for(n_bulk in n_bulks){
  for(case in cases){
    for(aa in dirichlet.scales){
      for(n_ref in n_refs){



        ref_idx_mat = readRDS(paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_ref_idx.rds',sep=''))
        bulk_p_array = readRDS(paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_bulk_p.rds',sep=''))
        groups_mat = readRDS(paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_groups_mat.rds',sep=''))
        bulk_idx_mat = readRDS(paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_bulk_idx_mat.rds',sep=''))

        for(cor_fdr in cor_fdrs){
          print(paste('Running:',case,'dirichlet.scale=',aa,"n_bulk=",n_bulk,'n_ref=',n_ref,'cor_fdr',cor_fdr))

          if(cor_fdr==0.1){
            R01 = readRDS("data/neuron/R01_neuron_G12400_alpha01.rds")
            R01.dist = as.dist(1-R01)
            clusters01 = cluster::pam(R01.dist,10,pamonce=5)
            folds01 = clusters01$clustering
            table(folds01)
            rm(R01.dist)
          }
          if(cor_fdr==0.3){
            R01 = readRDS("data/neuron/R01_neuron_G12400_alpha03.rds")
            R01.dist = as.dist(1-R01)
            clusters01 = cluster::pam(R01.dist,10,pamonce=5)
            folds01 = clusters01$clustering
            table(folds01)
            rm(R01.dist)
          }
          if(cor_fdr==0.5){
            R01 = readRDS("data/neuron/R01_neuron_G12400_alpha05.rds")
            R01.dist = as.dist(1-R01)
            clusters01 = cluster::pam(R01.dist,10,pamonce=5)
            folds01 = clusters01$clustering
            table(folds01)
            rm(R01.dist)
          }


          neuron_simu_study(indis_ref,
                            ref_idx_mat,
                            bulk_idx_mat,
                            bulk_p_array,
                            groups_mat,
                            case,
                            R01 = R01,
                            folds = folds01,
                            cor_fdr = cor_fdr,
                            add_bulk_bias=TRUE,
                            verbose = TRUE,
                            method_list = c('my'))

        }

      }
    }
  }
}

##########################################
##########################################




source('code/real/real_manuscript2.R')
indis_ref = readRDS('data/neuron/indis_ref_12400by6by97.rds')
R01 = readRDS("data/neuron/R01_neuron_G12400_alpha005.rds")



#case = c('null')

R01.dist = as.dist(1-R01)
clusters01 = cluster::pam(R01.dist,10,pamonce=5)
folds01 = clusters01$clustering
table(folds01)
rm(R01.dist)

n_rep = 100
dirichlet.scales = c(5,20)
n_bulks = c(86)
n_refs = 11
cases = c('null','all_diff')
for(n_bulk in n_bulks){
  for(case in cases){
    for(aa in dirichlet.scales){
      for(n_ref in n_refs){

        print(paste('Running:',case,'dirichlet.scale=',aa,"n_bulk=",n_bulk,'n_ref=',n_ref))

        ref_idx_mat = readRDS(paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_ref_idx.rds',sep=''))
        bulk_p_array = readRDS(paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_bulk_p.rds',sep=''))
        groups_mat = readRDS(paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_groups_mat.rds',sep=''))
        bulk_idx_mat = readRDS(paste("data/neuron/real_manu/input/neuron_ref",n_ref,"_rep",n_rep,"_bulk",n_bulk,"_dirichlet",aa,"_12400by97_",case,'_bulk_idx_mat.rds',sep=''))

        neuron_simu_study(indis_ref,
                          ref_idx_mat,
                          bulk_idx_mat,
                          bulk_p_array,
                          groups_mat,
                          case,
                          R01 = R01,
                          folds = folds01,
                          add_bulk_bias=TRUE,
                          verbose = TRUE,
                          method_list = c('my'))

      }
    }
  }
}



#
# R01 = sparseMatrix(i=corGene_idx_lower_cpm_alpha03[,1],j = corGene_idx_lower_cpm_alpha03[,2],dims=c(20293,20293))
# idx = match(gene_name_12400,gene_name_20293)
# R01=R01[idx,idx]
# R01 = R01 + t(R01)
# diag(R01) = 1
# rownames(R01) = gene_name_12400
# colnames(R01) = gene_name_12400
# R01[1:5,1:5]
# saveRDS(R01,"data/neuron/R01_neuron_G12400_alpha03.rds")
#
