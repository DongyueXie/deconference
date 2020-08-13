library(MuSiC)
library(xbioc)
XinT2Deset <- readRDS("data/MuSiC/XinT2Deset.rds")
XinT2D.construct.full = bulk_construct(XinT2Deset, clusters = 'cellType', samples = 'SubjectName')
XinT2D.construct.full$prop.real = relative.ab(XinT2D.construct.full$num.real, by.col = FALSE)
XinT2D.construct.full$prop.real

#GSE50244bulkeset = readRDS("data/MuSiC/GSE50244bulkeset.rds")
#gene_names = rownames(exprs(XinT2Deset))[which(rownames(exprs(XinT2Deset))%in%rownames(exprs(GSE50244bulkeset)))]
#w = rowSums(exprs(GSE50244bulkeset))
#w = w[match(gene_names,rownames(exprs(GSE50244bulkeset)))]

#gene_idx = match(gene_names,rownames(exprs(XinT2Deset)))
subject_name = levels(XinT2Deset$SubjectName)
ref_indi_name = subject_name[1:8]
ref_samples = which(XinT2Deset$SubjectName%in%(ref_indi_name))
Y = exprs(XinT2Deset)[,ref_samples]
cell_type_idx = XinT2Deset$cellType[ref_samples]
indi_idx = XinT2Deset$sampleID[ref_samples]
y = exprs(XinT2D.construct.full$Bulk.counts)[,-(which(XinT2D.construct.full$Bulk.counts$SubjectName%in%ref_indi_name))]

source('~/deconference/code/deconference_main.R')




datax = set_data_decon(y,Y,cell_type_idx = cell_type_idx,indi_idx = indi_idx,w=rep(1,nrow(Y)),gene_thresh = 0.1)
dim(datax$y)
out_universal_hc3 = deconference(datax,est_pop_var = TRUE,meta_mode = 'universal',correction = T)
save(out_universal_hc3,file = 'output/Xin_universal_hc3.RData')

out_universal_hc2 = deconference(datax,est_pop_var = TRUE,meta_mode = 'universal',hc.type = 'hc2',correction = T)
save(out_universal_hc2,file = 'output/Xin_universal_hc2.RData')

out_universal_hc0 = deconference(datax,est_pop_var = TRUE,meta_mode = 'universal',hc.type = 'hc0',correction = T)
save(out_universal_hc0,file = 'output/Xin_universal_hc0.RData')

out_bycelltype_hc3 = deconference(datax,est_pop_var = TRUE,meta_mode = 'by_celltype',correction = T)
save(out_bycelltype_hc3,file = 'output/Xin_bycelltype_hc3.RData')

# out_bycelltype_hc2 = deconference(datax,est_pop_var = TRUE,meta_mode = 'by_celltype',hc.type = 'hc2',correction = T)
# save(out_bycelltype_hc2,file = 'output/Xin_bycelltype_hc2.RData')

out_bygene_hc3 = deconference(datax,est_pop_var = TRUE,meta_mode = 'by_gene',hc.type = 'hc3',correction = T)
save(out_bygene_hc3,file = 'output/Xin_bygene_hc3.RData')

out_bygene_hc2 = deconference(datax,est_pop_var = TRUE,meta_mode = 'by_gene',hc.type = 'hc2',correction = T)
save(out_bygene_hc2,file = 'output/Xin_bygene_hc2.RData')

out_bygene_hc3 = deconference(datax,est_pop_var = TRUE,meta_mode = 'by_gene',hc.type = 'hc3',correction = T)
save(out_bygene_hc3,file = 'output/Xin_bygene_hc3.RData')

out_local_hc2 = deconference(datax,est_pop_var = TRUE,meta_mode = 'local',hc.type = 'hc2',correction = T)
save(out_local_hc2,file = 'output/Xin_local_hc2.RData')

out_local_hc3 = deconference(datax,est_pop_var = TRUE,meta_mode = 'local',hc.type = 'hc3',correction = T)
save(out_local_hc3,file = 'output/Xin_local_hc3.RData')

out_sampleM_hc2 = deconference(datax,est_pop_var = FALSE,correction = T,hc.type = 'hc2')
save(out_sampleM_hc2,file = 'output/Xin_sampleM_hc2.RData')

out_unadj = unadjusted_lm(out_sampleM_hc3$input$y,out_sampleM_hc3$input$X,out_sampleM_hc3$input$w)
save(out_unadj,file = 'output/Xin_sampleM_unadj.RData')



source('code/simu_func1.R')

ref = out_universal$input$X
Sigma = out_universal$input$Sigma
b = c(0.05,0.05,0.1,0.8)
b2 = c(0.05,0.15,0.3,0.5)
Ng = 3e3

set.seed(12345)
out_simu = simu_study(ref,Ng,b,
                      ref_type='multi_sc',
                      nreps = 100,
                      sc_lib_size = 20,
                      printevery = 10,same_indi = F,tau2 = Sigma*nrow(ref)^2/Ng^2,sigma2 = Sigma*nrow(ref)^2/Ng^2,
                      tau2known=F,sigma2known = F,est_pop_var  = T,correction = T,
                      b2=b2,weight = 'default',hc.type = 'hc0',n_indi = 10,meta_mode = 'universal')

