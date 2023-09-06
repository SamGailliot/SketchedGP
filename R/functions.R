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

# Helper function for get_lambda
# Get the pieces for the ratio of the different lambda_i posteriors
get_log_fs <- function(lambdas, Psi.XX, y){
  log.fs <- numeric(length = length(lambdas))
  n <- nrow(Psi.XX)
  for(ii in 1:length(log.fs)){
    K1 <- squared_exp_kernel(Psi.XX, lambdas[ii])
    G <- K1 + diag(1, n, n)
    log.fs[ii] <- -0.5 * log(det(G)) - n/2 * log(colSums(y * (solve(G) %*% y))) + log(dgamma(lambdas[ii], 1, 1))
  }
  return(log.fs)
}

# Helper function for get_lambda
get_ws <- function(log.fs){
  ws <- exp(log.fs - max(log.fs))
  Ws <- ws / sum(ws)
  return(Ws)
}

# Given a compression matrix \Psi construct lambda
get_lambda <- function(y, Psi.XX, dmin, dmax, k = 100, a = 1, b = 1){
  lambdas <- runif(k, min = 3/dmax, max = 3/dmin)

  # For each lambda we obtain f(lambda_i | y, Psi)
  log.fs <- get_log_fs(lambdas = lambdas, Psi.XX = Psi.XX, y = y)
  ws <- get_ws(log.fs)

  lambda <- sample(lambdas, size = 1, prob = ws)
  return(lambda)
}


# Given y, X, m & K return the matrix of loo densities
# get_probs <- function(y, X, m, K){
#   # Initialize constants and storage
#   n <- nrow(X); p <- ncol(X)
#   Psis <- vector("list", length = K)
#   lambdas <- numeric(length = K)
#   snrs <- numeric(length = K)
#   Pmat <- matrix(0, nrow = n, ncol = K)
#
#   for(ii in 1:K){
#     Psi <- get_sketch_mat(m, p)
#     sketched_X <- t(Psi %*% t(X))
#
#
#   }
# }

# What happens next?

get_swiss_roll <- function(n, p, tau){
  x.mat <- matrix(data = NA, nrow = n, ncol = p)
  y.vec <- numeric(length = n)
  for(ii in 1:n){
    d <- rnorm(p, 0, sd = tau)
    t <- runif(1, (3*pi)/2, (9*pi)/2)
    h <- runif(1, 0, 3)

    # Sample Noise s
    #d <- rnorm(p, 0, sd = 0.05)

    x <- numeric(length = length(d))
    x[1] <- t * cos(t); x[2] <- h; x[3] <- t * sin(t)
    x.mat[ii, ] <- x + d

    y.vec[ii] <- sin(5 * pi * t) + h^2 + rnorm(1, 0, 0.02)
  }
  out <- list("X" = x.mat, "y" = y.vec)
  return(out)
}
