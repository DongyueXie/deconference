#'@title select variables such that the matrix is close to diagonal mat.
#'@param n_var how many variables to keep.
corr_matrix_prune = function(X,n_var=500,decreasing=FALSE){
  X[which(is.na(X))] = 1
  total_cor = rowSums(abs((X)))
  selected = order(total_cor,decreasing=decreasing)[1:n_var]
  selected
}
