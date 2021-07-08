

day30bulk = readRDS("data/day30bulk.rds")
dim(day30bulk)
day30bulk[1:5,1:5]
hist(colSums(day30bulk),breaks = 100)

# filter out genes
gene_to_use = apply(day30bulk, MARGIN = 1, FUN = function( row, exp_th, min_sam ) {
  sum( row > exp_th ) >= min_sam
}, exp_th = 0.1, min_sam = 60 )
sum(gene_to_use)

day30bulk = day30bulk[gene_to_use,]

# transform to cpm
day30bulk = apply(day30bulk,2,function(z){z/sum(z)*1e6})


source('code/simulation/get_cor_pairs.R')

cor.idx = get_cor_pairs2(day30bulk,alpha=0.5,method='testing')

cor.idx.log = get_cor_pairs2(log(day30bulk+0.01),alpha=0.5,method='testing')


