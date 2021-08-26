

day30bulk = readRDS("data/neuron/day30bulk.rds")
dim(day30bulk)
day30bulk[1:5,1:5]
hist(colSums(day30bulk),breaks = 100)

# filter out genes
# a gene that has expression in at least 1/3 individuals.
gene_to_use = apply(day30bulk, MARGIN = 1, FUN = function( row, exp_th, min_sam ) {
  sum( row > exp_th ) >= min_sam
}, exp_th = 0.1, min_sam = 60 )
sum(gene_to_use)

day30bulk = day30bulk[gene_to_use,]

# transform to cpm
day30bulk = apply(day30bulk,2,function(z){z/sum(z)*1e6})


source('code/simulation/get_cor_pairs.R')

cor.idx = get_cor_pairs2(day30bulk,alpha=0.01,method='testing')

cor.idx.log = get_cor_pairs2(log(day30bulk+0.01),alpha=0.5,method='testing')


saveRDS('data/neuron/corGene_idx_lower_cpm_alpha001.rds')

#########################################
#########################################
#########################################

day30bulk = readRDS("data/neuron/day30bulk.rds")
gene_name_12400 = readRDS('data/neuron/gene_name_12400.rds')
bulk = day30bulk[match(gene_name_12400,rownames(day30bulk)),]
dim(bulk)
bulk = apply(bulk,2,function(z){z/sum(z)*1e6})
X = bulk
p = nrow(X)
n = ncol(X)

X.center = scale(t(X),center=TRUE,scale=FALSE)
S = Rfast::cova(X.center,center = TRUE)

# calc S2
S2 = 0
for(k in 1:n){
  S2 = S2+(tcrossprod(X.center[k,]))^2
}
Tmat = S*(n-1)/sqrt(S2+(2-n)*S^2)
P = 2*(1-pnorm(abs(Tmat[lower.tri(Tmat)])))
P.order = sort(P,decreasing = TRUE)

nt = length(P.order)
P.adj = P.order*(p^2-p)/2/(nt:1)
alpha =0.005
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

saveRDS(cor.idx,file='data/neuron/gene12400_cor_idx_alpha0005.rds')
