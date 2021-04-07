rmse = function(x,y){sqrt(mean((x-y)^2))}

music_prop = function(bulk.eset, sc.eset, markers = NULL, clusters, samples,
                      select.ct = NULL, cell_size = NULL, ct.cov = FALSE, verbose = TRUE,
                      iter.max = 1000, nu = 1e-04, eps = 0.01, centered = FALSE,
                      normalize = FALSE, ...)
{
  bulk.gene = rownames(bulk.eset)[rowMeans(exprs(bulk.eset)) !=
                                    0]
  bulk.eset = bulk.eset[bulk.gene, , drop = FALSE]
  if (is.null(markers)) {
    sc.markers = bulk.gene
  }
  else {
    sc.markers = intersect(bulk.gene, unlist(markers))
  }
  sc.basis = music_basis(sc.eset, non.zero = TRUE, markers = sc.markers,
                         clusters = clusters, samples = samples, select.ct = select.ct,
                         cell_size = cell_size, ct.cov = ct.cov, verbose = verbose)
  cm.gene = intersect(rownames(sc.basis$Disgn.mtx), bulk.gene)
  if (is.null(markers)) {
    if (length(cm.gene) < 0.2 * min(length(bulk.gene), nrow(sc.eset)))
      stop("Too few common genes!")
  }
  else {
    if (length(cm.gene) < 0.2 * length(unlist(markers)))
      stop("Too few common genes!")
  }
  if (verbose) {
    message(paste("Used", length(cm.gene), "common genes..."))
  }
  m.sc = match(cm.gene, rownames(sc.basis$Disgn.mtx))
  m.bulk = match(cm.gene, bulk.gene)
  D1 = sc.basis$Disgn.mtx[m.sc, ]
  M.S = colMeans(sc.basis$S, na.rm = T)
  if (!is.null(cell_size)) {
    if (!is.data.frame(cell_size)) {
      stop("cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes")
    }
    else if (sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))) {
      stop("Cell type names in cell_size must match clusters")
    }
    else if (any(is.na(as.numeric(cell_size[, 2])))) {
      stop("Cell sizes should all be numeric")
    }
    my_ms_names <- names(M.S)
    cell_size <- cell_size[my_ms_names %in% cell_size[, 1],
                           ]
    M.S <- cell_size[match(my_ms_names, cell_size[, 1]),
                     ]
    M.S <- M.S[, 2]
    names(M.S) <- my_ms_names
  }
  Yjg = relative.ab(exprs(bulk.eset)[m.bulk, ])
  N.bulk = ncol(bulk.eset)
  if (ct.cov) {
    Sigma.ct = sc.basis$Sigma.ct[, m.sc]
    if (sum(Yjg[, i] == 0) > 0) {
      D1.temp = D1[Yjg[, i] != 0, ]
      Yjg.temp = Yjg[Yjg[, i] != 0, i]
      Sigma.ct.temp = Sigma.ct[, Yjg[, i] != 0]
      if (verbose)
        message(paste(colnames(Yjg)[i], "has common genes",
                      sum(Yjg[, i] != 0), "..."))
    }
    else {
      D1.temp = D1
      Yjg.temp = Yjg[, i]
      Sigma.ct.temp = Sigma.ct
      if (verbose)
        message(paste(colnames(Yjg)[i], "has common genes",
                      sum(Yjg[, i] != 0), "..."))
    }
    lm.D1.weighted = music.iter.ct(Yjg.temp, D1.temp, M.S,
                                   Sigma.ct.temp, iter.max = iter.max, nu = nu, eps = eps,
                                   centered = centered, normalize = normalize)
    Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
    Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
    weight.gene.temp = rep(NA, nrow(Yjg))
    weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
    Weight.gene = cbind(Weight.gene, weight.gene.temp)
    r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
    Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
  }
  else {
    Sigma = sc.basis$Sigma[m.sc, ]
    valid.ct = (colSums(is.na(Sigma)) == 0) & (colSums(is.na(D1)) ==
                                                 0) & (!is.na(M.S))
    if (sum(valid.ct) <= 1) {
      stop("Not enough valid cell type!")
    }
    if (verbose) {
      message(paste("Used", sum(valid.ct), "cell types in deconvolution..."))
    }
    D1 = D1[, valid.ct]
    M.S = M.S[valid.ct]
    D1 = D1/M.S[1]
    M.S = M.S/M.S[1]

    Sigma = Sigma[, valid.ct]
    Est.prop.allgene = NULL
    Est.prop.weighted = NULL
    Weight.gene = NULL
    r.squared.full = NULL
    Var.prop = NULL
    for (i in 1:N.bulk) {
      if (sum(Yjg[, i] == 0) > 0) {
        D1.temp = D1[Yjg[, i] != 0, ]
        Yjg.temp = Yjg[Yjg[, i] != 0, i]
        Sigma.temp = Sigma[Yjg[, i] != 0, ]
        if (verbose)
          message(paste(colnames(Yjg)[i], "has common genes",
                        sum(Yjg[, i] != 0), "..."))
      }
      else {
        D1.temp = D1
        Yjg.temp = Yjg[, i]
        Sigma.temp = Sigma
        if (verbose)
          message(paste(colnames(Yjg)[i], "has common genes",
                        sum(Yjg[, i] != 0), "..."))
      }
      lm.D1.weighted = music.iter(Yjg.temp, D1.temp, M.S,
                                  Sigma.temp, iter.max = iter.max, nu = nu, eps = eps,
                                  centered = centered, normalize = normalize)
      Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
      Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
      weight.gene.temp = rep(NA, nrow(Yjg))
      weight.gene.temp[Yjg[, i] != 0] = lm.D1.weighted$weight.gene
      Weight.gene = cbind(Weight.gene, weight.gene.temp)
      r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
      Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
    }
  }
  colnames(Est.prop.weighted) = colnames(D1)
  rownames(Est.prop.weighted) = colnames(Yjg)
  colnames(Est.prop.allgene) = colnames(D1)
  rownames(Est.prop.allgene) = colnames(Yjg)
  names(r.squared.full) = colnames(Yjg)
  colnames(Weight.gene) = colnames(Yjg)
  rownames(Weight.gene) = cm.gene
  colnames(Var.prop) = colnames(D1)
  rownames(Var.prop) = colnames(Yjg)
  return(list(Est.prop.weighted = Est.prop.weighted, Est.prop.allgene = Est.prop.allgene,
              Weight.gene = Weight.gene, r.squared.full = r.squared.full,
              Var.prop = Var.prop,X=D1,S=M.S))
}





summary_real = function(out,true_b,groups){

  # compare fitted vs true
  b_hat = t(out$beta_hat)
  rownames(b_hat) = rownames(true_b)
  colnames(b_hat) = colnames(true_b)
  print(knitr::kable(
    list(t0 = round(b_hat,2), t1 = round(true_b,2)),
    caption = 'estimated vs true beta',
    booktabs = TRUE, valign = 't'))


  ci_length = t(round(out$beta_se,2))
  colnames(ci_length) = colnames(true_b)
  rownames(ci_length) = rownames(true_b)
  print(knitr::kable(ci_length,caption = 'beta hat se'))

  # summary statistics
  rmse = sqrt(mean((t(out$beta_hat)-true_b)^2))
  mad = mean(abs(t(out$beta_hat)-true_b))
  r2 = cor(c(t(out$beta_hat)),c(true_b))
  ss = rbind(c(rmse,mad,r2))
  colnames(ss) = c('rmse','mad','r2')
  rownames(ss) = NULL

   print(knitr::kable(ss,caption='error measure'))

  # coverage, length of CI
  ci = build_ci(out)
  coverage = c(true_b>=t(ci$ci_l))&(true_b<=t(ci$ci_r))

  print(knitr::kable(coverage,caption='coverage'))

  # two group test
  group_name = levels(as.factor(groups))
  g1 = which(groups==group_name[1])
  g2 = which(groups==group_name[2])
  test2 = two_group_test(out,groups)

  test_out = rbind(round(colMeans(true_b[g1,]) - colMeans(true_b[g2,]),3),
                   round(test2$diff_group,3),
                   round(sqrt(diag(test2$V_tilde)),3),
                   round(test2$p_value,3))
  rownames(test_out) = c('true diff','estimated diff','se','p value')
  print(knitr::kable(test_out,caption= ' two group test'))

}

get_coverage_beta = function(out,K,nbulk){
  cov_adj = matrix(out$coverage_adj,ncol=nbulk)
  cov_unadj = matrix(out$coverage_unadj,ncol=nbulk)
  cov_unadj_hc3 = matrix(out$coverage_unadj_hc3,ncol=nbulk)
  cov_unadj_cv = matrix(out$coverage_unadj_cv,ncol=nbulk)
  rname = paste('cell',1:K,sep = '')
  cname = paste('bulk',1:nbulk,sep = '')
  rownames(cov_adj) = rname
  colnames(cov_adj) = cname
  rownames(cov_unadj) = rname
  colnames(cov_unadj) = cname
  rownames(cov_unadj_hc3) = rname
  colnames(cov_unadj_hc3) = cname
  rownames(cov_unadj_cv) = rname
  colnames(cov_unadj_cv) = cname

  print(knitr::kable(cov_adj,caption = 'coverage adjusted'))
  print(knitr::kable(cov_unadj,caption = 'coverage unadjusted hc0'))
  print(knitr::kable(cov_unadj_hc3,caption = 'coverage unadjusted hc3'))
  print(knitr::kable(cov_unadj_cv,caption = 'coverage unadjusted lm'))
}

get_se_beta = function(out,K,nbulk){
  se_adj = round(matrix(colMeans(out$se_adj),ncol=nbulk),3)
  se_unadj = round(matrix(colMeans(out$se_unadj),ncol=nbulk),3)
  se_unadj_hc3 = round(matrix(colMeans(out$se_unadj_hc3),ncol=nbulk),3)
  se_unadj_cv = round(matrix(colMeans(out$se_unadj_cv),ncol=nbulk),3)
  rname = paste('cell',1:K,sep = '')
  cname = paste('bulk',1:nbulk,sep = '')
  rownames(se_adj) = rname
  colnames(se_adj) = cname
  rownames(se_unadj) = rname
  colnames(se_unadj) = cname
  rownames(se_unadj_hc3) = rname
  colnames(se_unadj_hc3) = cname
  rownames(se_unadj_cv) = rname
  colnames(se_unadj_cv) = cname

  print(knitr::kable(se_adj,caption = 'mean se adjusted'))
  print(knitr::kable(se_unadj,caption = 'mean se unadjusted hc0'))
  print(knitr::kable(se_unadj_hc3,caption = 'mean se unadjusted hc3'))
  print(knitr::kable(se_unadj_cv,caption = 'mean se unadjusted lm'))
}

get_coverage_diff = function(out,K,bulk){
  cov_diff = cbind(out$coverage_diff_adj,out$coverage_diff_unadj,out$coverage_diff_unadj_cv,out$coverage_diff_unadj_hc3)
  rname = paste('cell',1:K,sep = '')
  cname = c('adj','unadj-hc0','unadj-hc3','unadj-lm')
  rownames(cov_diff) = rname
  colnames(cov_diff) = cname
  print(knitr::kable(cov_diff,caption = 'coverage of group diff'))
}

get_se_diff = function(out,K,nbulk){
  se_diff = cbind(colMeans(out$diff_adj_se),colMeans(out$diff_unadj_se),colMeans(out$diff_unadj_se_hc3),colMeans(out$diff_unadj_se_cv))
  se_diff = round(se_diff,3)
  rname = paste('cell',1:K,sep = '')
  cname = c('adj','unadj-hc0','unadj-hc3','unadj-lm')
  rownames(se_diff) = rname
  colnames(se_diff) = cname
  print(knitr::kable(se_diff,caption = 'mean se of group diff'))
}


# random select a half from each dataset as reference
# data.obj is a list of single cell dataset
# return 1. bulk data; 2. the rest as reference; 3. the bulk as reference; 4. all as reference
create_realsimu_data = function(refs,cell_types,proph,propd,ncells = 1000,indis_ref = NULL){
  G = nrow(refs[[1]])
  bulk_counts = matrix(0,nrow=G,ncol=50)
  bulk_prop = c()
  all_indis = c()
  condition = c()
  indi_counter = 0
  ref1 = list()
  ref2 = list()
  ref_same = list()
  for(d in 1:length(refs)){

    # random sample a half individuals
    data.obj = refs[[d]]

    ## these indis are for bulk data
    if(is.null(indis_ref)){
      c.idx = which(data.obj$cell_type%in%cell_types)
      temp = table(data.obj$individual[c.idx],data.obj$cell_type[c.idx])
      temp = temp[,match(cell_types,colnames(temp))]
      valid_indi = which(rowSums(temp==0)==0)
      indis = sample(rownames(temp)[valid_indi],ceiling(length(valid_indi)/2/2)*2)
    }else{
      indis = indis_ref[[d]]
    }
    all_indis = c(all_indis,indis)
    condition = c(condition,c(rep('health',length(indis)/2),rep('disease',length(indis)/2)))
    #print(indis)
    # first half health
    # second half disease
    for(i in 1:(length(indis)/2)){
      idx.temp = which(data.obj$individual==indis[i]&data.obj$cell_type%in%cell_types)
      cell_counts = table(data.obj$cell_type[idx.temp])
      if(is.null(ncells)){
        ncells = round(max(cell_counts/proph))
      }
      bcounts = 0
      for(j in 1:length(cell_types)){
        idx_ij = which(data.obj$individual==indis[i]&data.obj$cell_type==cell_types[j])
        #print(idx_ij)
        #print(data.obj[,sample(idx_ij,ncells*prop[j],replace = TRUE)])
        bcounts = bcounts + rowSums(counts(data.obj[,sample(idx_ij,ncells*proph[j],replace = TRUE)]))
      }
      bulk_prop = rbind(bulk_prop,proph)
      indi_counter = indi_counter + 1
      bulk_counts[,indi_counter] = bcounts

    }


    ## disease
    for(i in (length(indis)/2+1):(length(indis))){
      idx.temp = which(data.obj$individual==indis[i]&data.obj$cell_type%in%cell_types)
      cell_counts = table(data.obj$cell_type[idx.temp])
      if(is.null(ncells)){
        ncells = round(max(cell_counts/propd))
      }
      bcounts = 0
      for(j in 1:length(cell_types)){
        idx_ij = which(data.obj$individual==indis[i]&data.obj$cell_type==cell_types[j])
        #print(idx_ij)
        #print(data.obj[,sample(idx_ij,ncells*prop[j],replace = TRUE)])
        bcounts = bcounts + rowSums(counts(data.obj[,sample(idx_ij,ncells*propd[j],replace = TRUE)]))
      }
      bulk_prop = rbind(bulk_prop,propd)
      indi_counter = indi_counter + 1
      bulk_counts[,indi_counter] = bcounts
    }


    ## the rest indis are for reference data

    ### the same dataset for reference
    ref_same_idx = which(data.obj$individual%in%indis&data.obj$cell_type%in%cell_types)
    ref_same[[d]] = data.obj[,ref_same_idx]
    ### the first reference data uses half the rest data
    rest_indis = unique(data.obj$individual)[-match(indis,unique(data.obj$individual))]
    ref1_indis = sample(rest_indis,length(rest_indis)/2)
    ref1_idx = which(data.obj$individual%in%ref1_indis&data.obj$cell_type%in%cell_types)
    ref1[[d]] = data.obj[,ref1_idx]

    ### the second reference data uses all of the half dataset
    ref2_idx = which(data.obj$individual%in%rest_indis&data.obj$cell_type%in%cell_types)
    ref2[[d]] = data.obj[,ref2_idx]
  }

  #browser()
  bulk_counts = bulk_counts[,which(colSums(bulk_counts)!=0)]
  colnames(bulk_counts) = all_indis
  rownames(bulk_counts) = rownames(data.obj)
  #rownames(bulk_ncell) = indis
  rownames(bulk_prop) = all_indis
  colnames(bulk_prop) = cell_types
  #ii = which(data.obj$individual%in%indis&data.obj$cell_type%in%cell_types)
  #bulk_prop = table(data.obj$individual[ii],data.obj$cell_type[ii])
  #bulk_prop = bulk_prop/rowSums(bulk_prop)
  coldata = DataFrame(individual = all_indis,condition=condition)
  #browser()
  out = SingleCellExperiment(assays = list(counts = as.matrix(bulk_counts)),
                             colData = coldata)
  return(list(bulk = out,bulk_prop = bulk_prop,group = condition,ref1=ref1,ref2=ref2,ref_same=ref_same))
}

create_bulk_abitrary_prop = function(refs,cell_types,proph,propd,ncells = 1000,seed=12345){
  set.seed(seed)
  G = nrow(refs[[1]])
  bulk_counts = matrix(0,nrow=G,ncol=50)
  bulk_prop = c()
  all_indis = c()
  condition = c()
  indi_counter = 0

  for(d in 1:length(refs)){

    # random sample a half individuals
    data.obj = refs[[d]]

    ## these indis are for bulk data

      c.idx = which(data.obj$cell_type%in%cell_types)
      temp = table(data.obj$individual[c.idx],data.obj$cell_type[c.idx])
      temp = temp[,match(cell_types,colnames(temp))]
      valid_indi = which(rowSums(temp==0)==0)
      indis = rownames(temp)[valid_indi]

    all_indis = c(all_indis,indis)
    condition = c(condition,c(rep('health',length(indis)/2),rep('disease',length(indis)/2)))
    #print(indis)
    # first half health
    # second half disease
    for(i in 1:(length(indis)/2)){
      idx.temp = which(data.obj$individual==indis[i]&data.obj$cell_type%in%cell_types)
      cell_counts = table(data.obj$cell_type[idx.temp])
      if(is.null(ncells)){
        ncells = round(max(cell_counts/proph))
      }
      bcounts = 0
      for(j in 1:length(cell_types)){
        idx_ij = which(data.obj$individual==indis[i]&data.obj$cell_type==cell_types[j])
        #print(idx_ij)
        #print(data.obj[,sample(idx_ij,ncells*prop[j],replace = TRUE)])
        bcounts = bcounts + rowSums(counts(data.obj[,sample(idx_ij,ncells*proph[j],replace = TRUE)]))
      }
      bulk_prop = rbind(bulk_prop,proph)
      indi_counter = indi_counter + 1
      bulk_counts[,indi_counter] = bcounts

    }


    ## disease
    for(i in (length(indis)/2+1):(length(indis))){
      idx.temp = which(data.obj$individual==indis[i]&data.obj$cell_type%in%cell_types)
      cell_counts = table(data.obj$cell_type[idx.temp])
      if(is.null(ncells)){
        ncells = round(max(cell_counts/propd))
      }
      bcounts = 0
      for(j in 1:length(cell_types)){
        idx_ij = which(data.obj$individual==indis[i]&data.obj$cell_type==cell_types[j])
        #print(idx_ij)
        #print(data.obj[,sample(idx_ij,ncells*prop[j],replace = TRUE)])
        bcounts = bcounts + rowSums(counts(data.obj[,sample(idx_ij,ncells*propd[j],replace = TRUE)]))
      }
      bulk_prop = rbind(bulk_prop,propd)
      indi_counter = indi_counter + 1
      bulk_counts[,indi_counter] = bcounts
    }
  }

  #browser()
  bulk_counts = bulk_counts[,which(colSums(bulk_counts)!=0)]
  colnames(bulk_counts) = all_indis
  rownames(bulk_counts) = rownames(data.obj)
  #rownames(bulk_ncell) = indis
  rownames(bulk_prop) = all_indis
  colnames(bulk_prop) = cell_types
  #ii = which(data.obj$individual%in%indis&data.obj$cell_type%in%cell_types)
  #bulk_prop = table(data.obj$individual[ii],data.obj$cell_type[ii])
  #bulk_prop = bulk_prop/rowSums(bulk_prop)
  coldata = DataFrame(individual = all_indis,condition=condition)
  #browser()
  out = SingleCellExperiment(assays = list(counts = as.matrix(bulk_counts)),
                             colData = coldata)
  return(list(bulk = out,bulk_prop = bulk_prop,group = condition))
}

sce_to_es = function(sce){
  #bulk.data = ExpressionSet(assayData = counts(enge_bulk$bulk))
  pheno_Data = data.frame(cell_type = sce$cell_type,individual = sce$individual)
  rownames(pheno_Data) = colnames(sce)
  ref.data = ExpressionSet(assayData = counts(sce),phenoData = AnnotatedDataFrame(pheno_Data))
  ref.data
}

merge_ref = function(refs){

  ref1_es = refs[[1]]
  colData(ref1_es) = DataFrame(cell_type = ref1_es$cell_type,individual = ref1_es$individual)
  rowData(ref1_es) = DataFrame(feature_symbol = rowData(ref1_es)$feature_symbol)
  #colnames(ref1_es) = ref1_es$cell_type
  reducedDims(ref1_es) = NULL
  assays(ref1_es)$logcounts = NULL
  if(length(refs)>1){
    for(i in 2:length(datax$ref1)){
      temp = refs[[i]]
      colData(temp) = DataFrame(cell_type = temp$cell_type,individual = temp$individual)
      rowData(temp) = DataFrame(feature_symbol = rowData(temp)$feature_symbol)
      #colnames(temp) = temp$cell_type
      reducedDims(temp) = NULL
      if(!is.null(assays(temp)$logcounts)){
        assays(temp)$logcounts = NULL
      }
      ref1_es = cbind(ref1_es,temp)
    }
  }

  ref1_es

}



# obj, a lm fitted model
# Omega, correlation materix
vcovHCC = function(obj,Omega=NULL,hc.type='hc3'){
  library(Rfast)
  X = as.matrix(obj$model[,-1])
  if(is.null(Omega)){
    A = mat.mult(mat.mult(solve(Crossprod(X,X)),t(X)),as.matrix(obj$residuals))
    V = Tcrossprod(A,A)
  }else{
    A = mat.mult(solve(Crossprod(X,X)),t(X))
    d = abs(obj$residuals)
    if(hc.type=='hc3'){
      d = d/(1-influence(obj)$hat)
    }
    Sigma = t(t(Omega*d)*d)
    V = mat.mult(mat.mult(A,Sigma),t(A))
  }
  V
}
