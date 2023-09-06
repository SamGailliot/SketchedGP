# Sample a sketching matrix Mat \in R^{m,p}
# m: Sketching Dimension
# p: Original Dimension of X
# orthog: Toggle to return orthogonal sketching matrix

get_sketch_mat <- function(m, p, orthog = FALSE){
    Mat <- matrix(rnorm(m*p), nrow = m, ncol = p)
    if(isTRUE(orthog)){
    Mat <- t(qr.Q(qr(t(Mat))))
  }
  return(Mat)
}

# Calculate exponential kernel
exp_kernel <- function(M1, M2 = NULL, lam){
  if(is.null(M2)){
    eps <- sqrt(.Machine$double.eps)
    Dists <- plgp::distance(M1)
    Kern <- exp(-lam * sqrt(Dists)) + diag(eps, nrow(M1))
  return(Kern)
  }
  else{
    if (ncol(M1) != ncol(M2)){
      stop("col dim mismatch for M1 & M2. Please ensure data have same dimension")
    }
    Dists <- plgp::distance(M1, M2)
    Kern <- exp(-lam * sqrt(Dists))
    return(Kern)
  }
}

# Fitting the GP?
