
source('code/real/real_manuscript.R')
source('code/deconference_main.R')
# For neuron dataset,
indis_ref = readRDS('data/neuron/indis_ref_12400by6by97.rds')
day30sce = readRDS('data/neuron/day30sce.rds')
R01 = readRDS("data/neuron/R01_neuron_G12400_alpha005.rds")

clusters01 = pam(R01.dist,10,pamonce=5)
folds01 = clusters01$clustering
table(folds01)

p1 = c(0.1,0.1,0.15,0.15,0.2,0.3)
p2 = c(0.1,0.15,0.25,0.3,0.1,0.1)

real_out = neuron_simu_study(indis_ref,day30sce,
                             p1,p2,dirichlet=TRUE,p1.prop=0.5,
                             bulk_lib_size = 500,
                             n_ref=11,n_rep=10,seed=12345,
                             verbose=TRUE,
                             R01=R01)
saveRDS(real_out,file = 'output/manuscript/real/neuron_ref11_rep10_fdr005.rds')




################# do not save data, my, music, cibersort##################
###################10/25/2021############################################
source('code/real/real_manuscript.R')
source('code/deconference_main.R')
source('code/CIBERSORT.R')
# For neuron dataset,
indis_ref = readRDS('data/neuron/indis_ref_12400by6by97.rds')
#day30sce = readRDS('data/neuron/day30sce.rds')
R01 = readRDS("data/neuron/R01_neuron_G12400_alpha005.rds")

R01.dist = as.dist(1-R01)
clusters01 = cluster::pam(R01.dist,10,pamonce=5)
folds01 = clusters01$clustering
table(folds01)

rm(R01.dist)

p1 = c(0.1,0.1,0.15,0.15,0.2,0.3)
p2 = c(0.1,0.15,0.25,0.3,0.1,0.1)
dirichlet.scale = 10
n_ref=11
n_rep=100

real_out = neuron_simu_study(indis_ref,
                             p1=p1,
                             p2=p2,
                             dirichlet=TRUE,
                             p1.prop=0.5,
                             dirichlet.scale = dirichlet.scale,
                             bulk_lib_size = 500,
                             n_ref=n_ref,
                             n_rep=n_rep,
                             seed=12345,
                             verbose=TRUE,
                             R01=R01,
                             folds=folds01,
                             save_data=FALSE,
                             file.name = paste("neuron_ref",n_ref,"_rep",n_rep,"_dirichlet",dirichlet.scale,"_corfdr005",sep=''))
#saveRDS(real_out,file = 'output/manuscript/real/neuron_ref11_rep100_fdr005.rds')


################# do not save data, my, music, cibersort, more cases##################
###################10/26/2021############################################
source('code/real/real_manuscript.R')
source('code/deconference_main.R')
source('code/CIBERSORT.R')
# For neuron dataset,
indis_ref = readRDS('data/neuron/indis_ref_12400by6by97.rds')
#day30sce = readRDS('data/neuron/day30sce.rds')
R01 = readRDS("data/neuron/R01_neuron_G12400_alpha005.rds")

R01.dist = as.dist(1-R01)
clusters01 = cluster::pam(R01.dist,10,pamonce=5)
folds01 = clusters01$clustering
table(folds01)

rm(R01.dist)

n_ref=11
n_rep=100

dirichlet.scales = c(5,10)
cases = c("null","all_diff")

for(aa in dirichlet.scales){
  for(case in cases){

    if(case=='null'){
      p1 = c(0.3,0.2,0.15,0.15,0.1,0.1)
      p2 = c(0.3,0.2,0.15,0.15,0.1,0.1)
    }else if(case=='all_diff'){
      p1 = c(0.15,0.15,0.1,0.1,0.2,0.3)
      p2 = c(0.1,0.1,0.2,0.3,0.15,0.15)
    }

    print(paste('Running:',case,'dirichlet.scale=',aa))

    real_out = neuron_simu_study(indis_ref,
                                 p1=p1,
                                 p2=p2,
                                 dirichlet=TRUE,
                                 p1.prop=0.5,
                                 dirichlet.scale = aa,
                                 bulk_lib_size = 500,
                                 n_ref=n_ref,
                                 n_rep=n_rep,
                                 seed=12345,
                                 verbose=TRUE,
                                 R01=R01,
                                 folds=folds01,
                                 save_data=FALSE,
                                 file.name = paste("neuron_ref",n_ref,"_rep",n_rep,"_dirichlet",aa,"_corfdr005_",case,sep=''))


  }
}





