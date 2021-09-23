

# For neuron dataset,
indis_ref = readRDS('data/neuron/indis_ref_12400by6by97.rds')
day30sce = readRDS('data/neuron/day30sce.rds')

individuals = dimnames(indis_ref)[[3]]

neuron_simu_study = function(indis_ref,day30sce,p1,p2,n_ref=11,n_rep=1,seed=12345){
  gene_to_use = dimnames(indis_ref)[[1]]
  individuals = dimnames(indis_ref)[[3]]
  celltypes = dimnames(indis_ref)[[2]]

  set.seed(seed)

  for(r in 1:n_rep){
    # randomly draw
  }

}
