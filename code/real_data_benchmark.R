
XinT2Deset <- readRDS("data/MuSiC/XinT2Deset.rds")
XinT2D.construct.full = bulk_construct(XinT2Deset, clusters = 'cellType', samples = 'SubjectName')
XinT2D.construct.full


subject_name = unique(XinT2Deset$SubjectName)
ref_indi_name = subject_name[1:6]
ref_samples = which(XinT2Deset$SubjectName%in%(ref_indi_name))
Y = exprs(XinT2Deset)[,ref_samples]
cell_type_idx = XinT2Deset$cellType[ref_samples]
indi_idx = XinT2Deset$sampleID[ref_samples]
y = exprs(XinT2D.construct.full$Bulk.counts)[,-(which(XinT2D.construct.full$Bulk.counts$SubjectName%in%ref_indi_name))]

datax = set_data_decon(y,Y,cell_type_idx = cell_type_idx,indi_idx = indi_idx)
out = deconference(datax,est_sigma2 = TRUE)
