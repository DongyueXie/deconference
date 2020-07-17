source('code/deconference_setdata.R')
source('code/deconference_meta.R')
source('code/deconference_estfunc.R')
source('code/unadjusted_lm.R')

#'@title deconvolution inference
#'@description There are three cases: 1. bulk reference data, then input Y (gene by cell type);
#'2. single-cell reference data from the same individual, then input Y (gene by cell) and cell_type_idx indicator;
#'3. single-cell reference data from multiple-subject, then input Y (gene by cell), cell_type_idx and individual index
#'@param y a vector of bulk sample; or matrix whose columns are bulk samples
#'@param Y a reference count matrix, bulk ref: gene by cell type; single cell ref: gene by cells
#'@param ref_type bulk, sc, multi_sc
#'@param cell_type_idx a numerical vector indicates cell types of Y if using single cell reference.
#'@param indi_idx a vector indicates individuals
#'@param tau2 varaince of gene expression across single cells, gene by cell type matrix; if NULL, will be estimated.
#'@param sigma2 variance of gene expression across individuals, gene by cell type matrix.
#'@param w gene weights, length(w) = number of genes. will be adjusted to sum to G.
#'@param x_estimator separate or aggregate
#'@param est_sigma2 whether estimate sigma^2, the variance of X across individuals
#'@param meta_var 'plug_in', or 'adjust'
#'@param adjust whether adjust for uncertainty in reference matrix
#'@param alpha significance level
#'@param hc.type hc3, hc2, or hc0
#'@param correction whether perform fuller's small sample correction.
#'@param a alpha in the Fuller's small sample correction
#'@param eps adjust of zero variane if a gene has no expression observed in one cell type
#'@param gene_thresh remove genes that appear in less than number of  cells
#'@return a list from estimation_func


deconference = function(data.obj,
                        marker_gene = NULL,
                        hc.type = 'hc3',
                        x_estimator = 'separate',
                        est_sigma2 = FALSE,
                        meta_var = 'plug_in',
                        a=0, correction=FALSE,eps=0){

  ref_type = data.obj$ref_type
  w = data.obj$w
  y = data.obj$y
  Y = data.obj$Y

  #browser()

  if(ref_type=='bulk'){

    design.mat = bulkRef_proc(Y)

  }else{

    cell_type_idx = data.obj$cell_type_idx
    tau2 = data.obj$tau2

    if(ref_type=='sc'){

      design.mat = scRef1_proc(Y,cell_type_idx,estimator=x_estimator,tau2=tau2)
    }else if(ref_type=='multi_sc'){
      indi_idx = data.obj$indi_idx
      sigma2 = data.obj$sigma2

      # multiple individual single cell reference samples, estimate sigma^2
      design.mat = scRef_multi_proc(Y,cell_type_idx,indi_idx,estimator=x_estimator,tau2=tau2,
                                    sigma2=sigma2,est_sigma2 = est_sigma2,eps=eps,meta_var=meta_var)
    }else{
      stop("unspported reference type")
    }

  }

  #browser()

  out = estimation_func(y=y,X=design.mat$X,Vg=design.mat$Vg,design.mat$Sigma,marker_gene=marker_gene,
                        w=w,hc.type=hc.type,a=a,correction=correction)
  return(out)
}
