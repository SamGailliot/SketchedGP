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

get_NIS_scores <- function(X, y, dn = NULL, DOPAR = TRUE){

  if(isTRUE(DOPAR)){
    n_cores <- detectCores()
    registerDoParallel(cores = n_cores)
    n <- nrow(X)
    p <- ncol(X)
    if(is.null(dn)){dn <- floor(n^(1/5) + 2)}
    SSPs <- numeric(p)

    SSPs <- foreach(jj = 1:p, .combine = c) %dopar% {
      spline.mat <- bs(X[,jj], degree = dn)
      fit.lm.bs <- lm(y ~ spline.mat)

      (1/n) * sum(fit.lm.bs$fitted.values^2)
    }
    return(SSPs)
  }else{
    n <- nrow(X)
    p <- ncol(X)
    if(is.null(dn)){dn <- floor(n^(1/5) + 2)}
    SSPs <- numeric(p)

    for(jj in 1:p){
      spline.mat <- bs(X[,jj], degree = dn)
      fit.lm.bs <- lm(y ~ spline.mat)
      SSPs[jj] <- (1/n) * sum(fit.lm.bs$fitted.values^2)
    }

    return(SSPs)
  }
}

get_sketch_mat <- function(m, p, orthog = FALSE){
    Mat <- matrix(rnorm(m*p), nrow = m, ncol = p)
    if(isTRUE(orthog)){
    Mat <- t(qr.Q(qr(t(Mat))))
  }
  return(Mat)
}

# Calculate exponential kernel
exp_kernel <- function(M1, M2 = NULL, theta){
  if(is.null(M2)){
    eps <- sqrt(.Machine$double.eps)
    Dists <- plgp::distance(M1)
    Kern <- exp(-theta * sqrt(Dists)) + diag(eps, nrow(M1))
    #Kern <- exp(sqrt(Dists)/theta) + diag(eps, nrow(M1))
    #Kern <- exp(-theta * Dists) + diag(eps, nrow(M1))
  return(Kern)
  }
  else{
    if (ncol(M1) != ncol(M2)){
      stop("col dim mismatch for M1 & M2. Please ensure data have same dimension")
    }
    Dists <- plgp::distance(M1, M2)
    Kern <- exp(-theta * sqrt(Dists))
    #Kern <- exp(sqrt(Dists)/theta)
    #Kern <- exp(-theta * Dists)
    return(Kern)
  }
}

# Sample length-scale parameter theta
get_theta <- function(y, sketched.X, dmin, dmax, N, SNR){
  # browser()
  n <- nrow(sketched.X)
  Identity_Mat <- diag(1, n,n)
  # Sample thetas
  thetas <- runif(N, min = 3/dmax, max = 3/dmin)

  # Obtain the log densities f(theta_i|y, Pn, psi^2)
  log.fs <- numeric(length = N)

  # for(ii in 1:N){
  #   C <- exp_kernel(sketched.X, theta = thetas[ii])
  #   G <- (SNR*C) + Identity_Mat + diag(sqrt(.Machine$double.eps), nrow(C))
  #   Ginv <- chol2inv(chol(G))
  #   logdetG <- det(G, log = TRUE)
  #   # log.fs[ii] <- -0.5 * logdetG - n/2 * log(colSums(y * (Ginv %*% y))) # multiply by uniform prior over thetas
  #   log.fs[ii] <- det(chol(G), log = TRUE) - n/2 * log(colSums(y * (Ginv %*% y))) # multiply by uniform prior over thetas
  # }

  log.fs <- foreach(ii = 1:N, .combine = c) %dopar%{
    C <- exp_kernel(sketched.X, theta = thetas[ii])
    G <- (SNR*C) + Identity_Mat + diag(sqrt(.Machine$double.eps), nrow(C))
    Ginv <- chol2inv(chol(G))
    logdetG <- det(G, log = TRUE)
    # log.fs[ii] <- -0.5 * logdetG - n/2 * log(colSums(y * (Ginv %*% y))) # multiply by uniform prior over thetas

    det(chol(G), log = TRUE) - n/2 * log(colSums(y * (Ginv %*% y))) # multiply by uniform prior over thetas
    }

  # Get probability weights
  ws <- exp(log.fs - max(log.fs))
  Ws <- ws / sum(ws)

  # Sample theta according to weights
  theta <- sample(thetas, size = 1, prob = Ws)

  # return(list("thetas" = thetas, "ws" = log.fs, "theta" = theta))
  return(theta)
}

# Fill in the ii'th column of the loo densities matrix
get_model_densities <- function(y, sketched.X, theta, SNR){
  # browser()

  n <- nrow(sketched.X)
  n.old <- (n-1)
  p.out <- numeric(length = n)
  mu.out <- numeric(length = n)

  # Get the exponential kernel | theta
  C <- exp_kernel(sketched.X, theta = theta)

  # for i = 1:n get the loo density
  # for(ii in 1:n){
  p.out <- foreach(ii = 1:n, .combine = c) %dopar%{
    # Covariance between non loo locs
    C.loo <- C[-ii, -ii]

    # Covariance between loo and non-loo locs
    C.pred <- C[-ii, ii]

    # non-loo response
    y.loo <- as.matrix(y[-ii], nrow = (n-1), ncol = 1)

    G <- SNR*C.loo + diag( rep(1, (n-1)) )

    G.inv <- chol2inv(chol(G))

    mu.pred = SNR*as.numeric(t(C.pred) %*% G.inv %*% y.loo)
    b1 <-  (t(y.loo) %*% G.inv %*% y.loo) / 2

    sigsq.pred <- as.numeric((2*b1/n.old) * (1 + SNR - (SNR^2)*(t(C.pred) %*% G.inv %*% C.pred)) )

    y.pred <- y[ii]

    T.pred <- (y.pred - mu.pred)/sqrt(sigsq.pred)

    dt(T.pred, n)
    # mu.out[ii] <- mu.pred
  }
  #return(list("dens" = p.out, "mu" = mu.out))
  return(p.out)
}

# Given y, X, m & K return the matrix of loo densities
get_densities <- function(y, X, m, K, SNRs = c(0.1, 0.5, 1, 2), N.thetas = 50){
  # Initialize constants and storage
  # browser()
  n_cores <- detectCores()
  registerDoParallel(cores = 8)
  n <- nrow(X); p <- ncol(X); N.snrs <- length(SNRs)

  P.list <- vector("list", length = K)
  thetas <- matrix(0, nrow = K, ncol = N.snrs) # K by N.snrs matrix which holds all length scales

  dist.mat <- plgp::distance(X)
  dmax <- max(dist.mat)
  dmin <- min( dist.mat[dist.mat != 0] )
  # maxs <- numeric(length = K)
  # mins <- numeric(length = K)
  Density.Mat <- matrix(0, nrow = n, ncol = K*N.snrs)

  # Construct K Ps
  for(ii in 1:K){
    P.list[[ii]] <- get_sketch_mat(m, p, orthog = FALSE)
  }

  #density.mat.col <- 1
  # for(ii in 1:K){
  Density.Mat <- foreach(ii = 1:K, .combine = cbind)%dopar%{
    Density.Mat <- matrix(data = 0, nrow = (n+3), ncol = N.snrs)
    thetas <- matrix(data = 0, nrow = 1, ncol = N.snrs)
    sketched.X <- t(P.list[[ii]] %*% t(X))

    # Now for each member of the snrs list we
    for(jj in 1:length(SNRs)){
      # For each SNR, P pair we sample an appropriate theta
      #thetas[ii, jj] <- get_theta(y, sketched.X = sketched.X, dmin = dmin, dmax = dmax, N = N.thetas, SNR = SNRs[jj])
      theta <- get_theta(y, sketched.X = sketched.X, dmin = dmin, dmax = dmax, N = N.thetas, SNR = SNRs[jj])
      # Fill in one column of the Density Mat each time
      C <- exp_kernel(sketched.X, theta = thetas[jj])

      # Calculate all loo densities.
      # The top three rows are the theta values and then the repeated ii, jj values
      # this is a hacky solution to the problem of tracking and saving both
      # the densities and the theta values in the foreach loop...

      Density.Mat[,jj] <-c(theta,ii,jj, get_model_densities(y, sketched.X, theta = theta, SNRs[jj]))
      #density.mat.col <- density.mat.col + 1
    }
    Density.Mat
  }
  Thetas <- matrix(Density.Mat[1, ], nrow = K, ncol = N.snrs, byrow = TRUE)
  Density.Mat <- Density.Mat[-c(1,2,3), ]
  return(list("Density.Mat" = Density.Mat, "Thetas" = Thetas, "P.list" = P.list))
}

# Given y, X, m & K return the matrix of loo densities
get_densities_parallel <- function(y, X, m, K, SNRs = c(0.1, 0.5, 1, 2), N.thetas = 50){
  # Initialize constants and storage
  # browser()
  n <- nrow(X); p <- ncol(X); N.snrs <- length(SNRs)

  P.list <- vector("list", length = K)
  thetas <- matrix(0, nrow = K, ncol = N.snrs) # K by N.snrs matrix which holds all length scales

  dist.mat <- plgp::distance(X)
  dmax <- max(dist.mat)
  dmin <- min( dist.mat[dist.mat != 0] )
  # maxs <- numeric(length = K)
  # mins <- numeric(length = K)
  Density.Mat <- matrix(0, nrow = n, ncol = K*N.snrs)

  # Construct K Ps
  for(ii in 1:K){
    P.list[[ii]] <- get_sketch_mat(m, p, orthog = FALSE)
  }

  density.mat.col <- 1
  for(ii in 1:K){
    sketched.X <- t(P.list[[ii]] %*% t(X))

    # Now for each member of the snrs list we
    for(jj in 1:length(SNRs)){
      # For each SNR, P pair we sample an appropriate theta
      thetas[ii, jj] <- get_theta(y, sketched.X = sketched.X, dmin = dmin, dmax = dmax, N = N.thetas, SNR = SNRs[jj])

      # Fill in one column of the Density Mat each time
      C <- exp_kernel(sketched.X, theta = thetas[ii, jj])

      # Calculate
      Density.Mat[, density.mat.col] <- get_model_densities(y, sketched.X, theta = thetas[ii, jj], SNRs[jj])

      density.mat.col <- density.mat.col + 1
    }
  }


  return(Density.Mat)
}

# What happens next?
