# 1. Sample a sketching matrix Mat \in R^{m,p}
# m: Sketching Dimension
# p: Original Dimension of X
# orthog: Toggle to return orthogonal sketching matrix

get_sketch_mat <- function(m, p, orthog = FALSE){
  # browser()
  if(isFALSE(orthog)){
    Mat <- matrix(rnorm(m*p), nrow = ml, ncol = p)
  }else{
    Mat <- t(qr.Q(qr(t(Mat))))
  }
  return(Mat)
}
