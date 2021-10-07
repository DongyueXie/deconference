
source('code/real/real_manuscript.R')
source('code/deconference_main.R')
# For neuron dataset,
indis_ref = readRDS('data/neuron/indis_ref_12400by6by97.rds')
day30sce = readRDS('data/neuron/day30sce.rds')
R01 = readRDS("data/neuron/R01_neuron_G12400_alpha005.rds")


p1 = c(0.1,0.1,0.15,0.15,0.2,0.3)
p2 = c(0.1,0.15,0.25,0.3,0.1,0.1)

real_out = neuron_simu_study(indis_ref,day30sce,
                             p1,p2,dirichlet=TRUE,p1.prop=0.5,
                             bulk_lib_size = 500,
                             n_ref=11,n_rep=10,seed=12345,
                             verbose=TRUE,
                             R01=R01)
saveRDS(real_out,file = 'output/manuscript/real/neuron_ref11_rep10_fdr005.rds')
