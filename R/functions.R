# Sample a sketching matrix Mat \in R^{m,p}
# m: Sketching Dimension
# p: Original Dimension of X
# orthog: Toggle to return orthogonal sketching matrix

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

# Sample length-scale parameter theta
get_theta <- function(y, sketched_X, dmin, dmax, N, SNR){
  # browse()
  n <- nrow(sketched_X)
  Identity_Mat <- diag(1, n,n)
  # Sample thetas
  thetas <- runif(N, min = 3/dmax, max = 3/dmin)

  # Obtain the log densities f(theta_i|y, Pn, psi^2)
  log.fs <- numeric(length = N)

  for(ii in 1:N){
    C <- exp_kernel(sketched_X, lam = thetas[ii])
    G <- (SNR*C) + Identity_Mat
    Ginv <- chol2inv(chol(G))
    logdetG <- log(det(G))
    log.fs[ii] <- -0.5 * logdetG - n/2 * log(colSums(y * (Ginv %*% y))) # multiply by uniform prior over lambda
  }

  # Get probability weights
  ws <- exp(log.fs - max(log.fs))
  Ws <- ws / sum(ws)

  # Samps theta according to weights
  theta <- sample(thetas, size = 1, prob = Ws)

  return(list("thetas" = thetas, "ws" = log.fs, "theta" = theta))
  # return(theta)
}

# Given y, X, m & K return the matrix of loo densities
get_probs <- function(y, X, m, K, SNRs = c(0.1, 0.5, 1, 2), N.thetas = 50){
  # Initialize constants and storage
  n <- nrow(X); p <- ncol(X)

  P.list <- vector("list", length = K)
  thetas <- numeric(length = K)

  dist.mat <- plgp::distance(X)
  dmax <- max(dist.mat)
  dmin <- min( dist.mat[dist.mat != 0] )
  # maxs <- numeric(length = K)
  # mins <- numeric(length = K)
  Pmat <- matrix(0, nrow = n, ncol = K)

  # Construct K Ps
  for(ii in 1:K){
    P.list[[ii]] <- get_sketch_mat(m, p, orthog = FALSE)
  }

  for(ii in 1:K){
    sketched_X <- t(P.list[[ii]] %*% t(X))

    # Now for each member of the snrs list we
    for(jj in 1:length(SNRs)){
      # For each SNR, P pair we sample an appropriate theta
      get_theta(y, sketched_X = sketched_X, dmin = dmin, dmax = dmax, N = N.thetas, SNR = SNRs[jj])
      # Fill in one column of the Prob Mat each time
    }

  }
}

# What happens next?
