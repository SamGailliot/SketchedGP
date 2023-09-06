# Sample a sketching matrix Mat \in R^{m,p}
# m: Sketching Dimension
# p: Original Dimension of X
# orthog: Toggle to return orthogonal sketching matrix

get_sketch_mat <- function(m, p, orthog = FALSE){
  # browser()
    Mat <- matrix(rnorm(m*p), nrow = m, ncol = p)
    if(isTRUE(orthog)){
    Mat <- t(qr.Q(qr(t(Mat))))
  }
  return(Mat)
}

# Calculate exponential kernel
exp_kernel <- function(Mat, lam){
  eps <- sqrt(.Machine$double.eps)
  Dists <- plgp::distance(Mat)
  Kern <- exp(-lam * sqrt(Dists)) + diag(eps, nrow(Mat))
  return(Kern)
}
