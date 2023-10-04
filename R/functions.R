# Sample a sketching matrix Mat \in R^{m,p}
# m: Sketching Dimension
# p: Original Dimension of X
# orthog: Toggle to return orthogonal sketching matrix

get_torus <- function(n, p, tau){
  P = p - 3 # Leftover dimensions to fill with noise
  tor <- geozoo::torus(p = 3, n = n, radius = c(2,1)) # Generate torus
  X <- tor$points + matrix(rnorm(n*3,0, tau), n, 3)
  X <- cbind(X, matrix(rnorm(P*n, 0, tau), n, P))
  y <- scale(X[,2]^2 - sin(5*pi*X[,3])) + rnorm(n, 0, 0.1)
  out <- list("X" = X, "y" = y)
  return(out)
}

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
    Mat <- matrix(rnorm(m*p, 0, 1), nrow = m, ncol = p)
    if(isTRUE(orthog)){
    Mat <- t(qr.Q(qr(t(Mat))))
  }
  return(Mat)
}

# Calculate exponential kernel
exp_kernel <- function(M1, M2 = NULL, theta){
  # browser()
  if(is.null(M2)){
    eps <- sqrt(.Machine$double.eps)
    Dists <- sqrt(plgp::distance(M1))
    #Dists <- plgp::distance(M1)
    Kern <- exp(-theta * Dists) + diag(eps, nrow(M1))
    #Kern <- exp(sqrt(Dists)/theta) + diag(eps, nrow(M1))
    #Kern <- exp(-theta * Dists) + diag(eps, nrow(M1))
  return(Kern)
  }
  else{
    if (ncol(M1) != ncol(M2)){
      stop("col dim mismatch for M1 & M2. Please ensure data have same dimension")
    }
    Dists <- sqrt(plgp::distance(M1, M2))
    #Dists <- plgp::distance(M1, M2)
    Kern <- exp(-theta * Dists)
    #Kern <- exp(sqrt(Dists)/theta)
    #Kern <- exp(-theta * Dists)
    return(Kern)
  }
}

# Sample length-scale parameter theta
get_theta <- function(y, sketched.X, dmin, dmax, N, SNR){
  # browser()
  # SNR = 1
  n <- nrow(sketched.X)
  Identity_Mat <- diag(1, n,n)
  # Sample thetas
  thetas <- runif(N, min = 3/dmax, max = 3/dmin)

  # Obtain the log densities f(theta_i|y, Pn, psi^2)
  log.fs <- numeric(length = N)

  # for(ii in 1:N){
  #   C <- exp_kernel(sketched.X, theta = thetas[ii])
  #   G <- (SNR*C) + Identity_Mat + diag(sqrt(.Machine$double.eps), nrow(C))
  #   G.chol <- chol(G)
  #   Ginv <- chol2inv(G.chol)
  #   # logdetG <- det(G, log = TRUE)
  #   logdetG <- 2*sum(log(diag(G.chol)))
  #
  #   log.fs[ii] <- -0.5 * logdetG - n/2 * log(colSums(y * (Ginv %*% y))) # multiply by uniform prior over thetas
  #   #log.fs[ii] <- det(chol(G), log = TRUE) - n/2 * log(colSums(y * (Ginv %*% y))) # multiply by uniform prior over thetas
  # }

  log.fs <- foreach(ii = 1:N, .combine = c) %dopar%{
    C <- exp_kernel(sketched.X, theta = thetas[ii])
    G <- (SNR*C) + Identity_Mat + diag(sqrt(.Machine$double.eps), nrow(C))
    Gchol <- chol(G)
    Ginv <- chol2inv(chol(G))
    #logdetG <- det(G, log = TRUE)
    logdetG <- 2*sum(log(diag(Gchol))) # Compare with SVDs
    # log.fs[ii] <- -0.5 * logdetG - n/2 * log(colSums(y * (Ginv %*% y))) # multiply by uniform prior over thetas

    #det(chol(G), log = TRUE) - n/2 * log(colSums(y * (Ginv %*% y))) # multiply by uniform prior over thetas
    #det(chol(G), log = TRUE) - n/2 * log(colSums(y * (Ginv %*% y)))# + log(dgamma(thetas[ii], 1, 1))
    -0.5 * logdetG - n/2 * log(colSums(y * (Ginv %*% y)))
    }

  # Get probability weights
  ws <- exp(log.fs - max(log.fs))
  Ws <- ws / sum(ws)

  # Sample theta according to weights
  theta <- sample(thetas, size = 1, prob = Ws)
  # theta <- thetas[which.max(Ws)]
  # return(list("thetas" = thetas, "ws" = log.fs, "theta" = theta))
  return(theta)
}

# Fill in the ii'th column of the loo densities matrix
get_model_densities <- function(y, sketched.X, theta, SNR){
  #browser()

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

    # G <- SNR*C.loo + diag( rep(1, (n-1)) )
    #G.inv <- chol2inv(chol(G))

    G.inv <- chol2inv(chol(SNR*C.loo + diag( rep(1, (n-1)) )))

    mu.pred = SNR*as.numeric(t(C.pred) %*% G.inv %*% y.loo)
    b1 <-  (t(y.loo) %*% G.inv %*% y.loo) / 2

    sigsq.pred <- as.numeric((2*b1/n) * (1 + SNR - (SNR^2)*(t(C.pred) %*% G.inv %*% C.pred)) )

    y.pred <- y[ii]

    T.pred <- (y.pred - mu.pred)/sqrt(sigsq.pred)

    # dt(T.pred, n)
    #mu.pred
    c(dt(T.pred, n), mu.pred)
  }
  #return(list("dens" = p.out, "mu" = mu.out))
  return(p.out)
}

# Fill in the ii'th column of the loo densities matrix
get_model_preds <- function(y, sketched.X, theta, SNR){
  #browser()

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

    # G <- SNR*C.loo + diag( rep(1, (n-1)) )
    #G.inv <- chol2inv(chol(G))

    G.inv <- chol2inv(chol(SNR*C.loo + diag( rep(1, (n-1)) )))

    mu.pred = SNR*as.numeric(t(C.pred) %*% G.inv %*% y.loo)
    b1 <-  (t(y.loo) %*% G.inv %*% y.loo) / 2

    sigsq.pred <- as.numeric((2*b1/n) * (1 + SNR - (SNR^2)*(t(C.pred) %*% G.inv %*% C.pred)) )

    y.pred <- y[ii]

    T.pred <- (y.pred - mu.pred)/sqrt(sigsq.pred)

    c(dt(T.pred, n), mu.pred)
  }
  #return(list("dens" = p.out, "mu" = mu.out))
  return(p.out)
}

# get_model_preds_matrix <- function(X, Xstar, y, models.mat, Sketch.Mat.List){
#
# }


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
    P.list[[ii]] <- get_sketch_mat(m, p, orthog = TRUE)
  }

  density.mat.col <- 1
  #for(ii in 1:K){
  Density.Mat <- foreach(ii = 1:K, .combine = cbind)%dopar%{
    Density.Mat <- matrix(data = 0, nrow = (n+3), ncol = N.snrs)
    #thetas <- matrix(data = 0, nrow = 1, ncol = N.snrs)
    sketched.X <- t(P.list[[ii]] %*% t(X))

    # Now for each member of the snrs list we
    for(jj in 1:length(SNRs)){
      # For each SNR, P pair we sample an appropriate theta
      #thetas[ii, jj] <- get_theta(y, sketched.X = sketched.X, dmin = dmin, dmax = dmax, N = N.thetas, SNR = SNRs[jj])
      theta <- get_theta(y, sketched.X = sketched.X, dmin = dmin, dmax = dmax, N = N.thetas, SNR = SNRs[jj])
      # Fill in one column of the Density Mat each time
      C <- exp_kernel(sketched.X, theta = theta)

      # Calculate all loo densities.
      # The top three rows are the theta values and then the repeated ii, jj values
      # this is a hacky solution to the problem of tracking and saving both
      # the densities and the theta values in the foreach loop...

      Density.Mat[,jj] <-c(theta,ii,jj, get_model_densities(y, sketched.X, theta = theta, SNRs[jj]))
      #Density.Mat[,density.mat.col] <-c(theta,ii,jj, get_model_densities(y, sketched.X, theta = theta, SNRs[jj]))
      #density.mat.col <- density.mat.col + 1
    }
    Density.Mat
  }
  Thetas <- matrix(Density.Mat[1, ], nrow = K, ncol = N.snrs, byrow = TRUE)
  Density.Mat <- Density.Mat[-c(1,2,3), ]
  return(list("Density.Mat" = Density.Mat, "Thetas" = Thetas, "P.list" = P.list))
}

get_preds_parallel <- function(Xstar, X, y, Sketch.Mat.List, thetas, SNRs, stack.weights, Mins, Maxs, sketch.mat.num){
  # browser()
  if(is.null(dim(Xstar))){Xstar <- matrix(Xstar, nrow = 1)}

  weights.inds <- which(stack.weights != 0)
  K = length(thetas)
  n = nrow(X)
  nstar = nrow(Xstar)
  I <- diag(1, nrow = n, ncol = n)
  Istar <- diag(1, nrow = nstar, ncol = nstar)
  y <- as.matrix(y, nrow = length(y), ncol = 1)

  # SIGS <- vector("list", length(w.inds))
  sigs <- matrix(data = 0, nrow = length(weights.inds), ncol = nstar)
  mus <- matrix(data = 0, nrow = length(weights.inds), ncol = nstar)

  #jj = 1
  out.loop <- foreach(ii = weights.inds) %dopar% {
  #for(ii in weights.inds){
    Sketched.X <- t(Sketch.Mat.List[[sketch.mat.num[ii]]] %*% t(X))
    Sketched.Xstar <- t(Sketch.Mat.List[[sketch.mat.num[ii]]] %*% t(Xstar))
    m = ncol(Sketched.X)
    # How do I restandardize...
    mins <- apply(Sketched.X, 2, min)
    maxs <- apply(Sketched.X, 2, max)
    Sketched.X <- ( Sketched.X - matrix(mins, nrow = n, ncol = m, byrow = TRUE) ) / (matrix(maxs, nrow = n, ncol = m, byrow = TRUE) - matrix(mins, nrow = n, ncol = m, byrow = TRUE))
    # Sketched.X <- ( Sketched.X - matrix( Mins[[ii]] , nrow = nrow(Sketched.X), ncol = ncol(Sketched.X), byrow = TRUE ) ) / ( matrix( Maxs[[ii]] , nrow = nrow(Sketched.X), ncol = ncol(Sketched.X), byrow = TRUE ) - matrix( Mins[[ii]] , nrow = nrow(Sketched.X), ncol = ncol(Sketched.X), byrow = TRUE ) )
    #Psi.XXstar <- ( Psi.XXstar - matrix( Mins[[ii]] , nrow = nrow(Psi.XXstar), ncol = ncol(Psi.XXstar), byrow = TRUE ) ) / ( matrix( Maxs[[ii]] , nrow = nrow(Psi.XXstar), ncol = ncol(Psi.XXstar), byrow = TRUE ) - matrix( Mins[[ii]] , nrow = nrow(Psi.XXstar), ncol = ncol(Psi.XXstar), byrow = TRUE ) )
    Sketched.Xstar <- ( Sketched.Xstar - matrix( apply(Sketched.Xstar, 2, min) , nrow = nrow(Sketched.Xstar), ncol = ncol(Sketched.Xstar), byrow = TRUE ) ) / ( matrix( apply(Sketched.Xstar, 2, max), nrow = nrow(Sketched.Xstar), ncol = ncol(Sketched.Xstar), byrow = TRUE ) - matrix( apply(Sketched.Xstar, 2, min) , nrow = nrow(Sketched.Xstar), ncol = ncol(Sketched.Xstar), byrow = TRUE ) )
    K1 <- exp_kernel(Sketched.X, theta = thetas[ii])
    # K.1.pred <- exp(-lams[ii] *cdist(Psi.XX, Psi.XXstar)); K.pred.1 <- t(K.1.pred)
    # K.1.pred <- exp(-lams[ii] *cdist(X, Xstar)); K.pred.1 <- t(K.1.pred)
    K.1.pred <- exp_kernel(Sketched.X, Sketched.Xstar, theta = thetas[ii]); K.pred.1 <- t(K.1.pred)
    # K.pred <- squared_exp_kernel(Psi.XXstar, lambda = lams[ii])
    K.pred <- exp_kernel(Sketched.Xstar, theta = thetas[ii])
    # K.pred <- exp_kernel(Xstar, theta = thetas[ii])

    SNR = SNRs[ii]
    G <- SNR*K1 + I
    G.chol <- chol(G)
    G.inv <- chol2inv(G.chol)

    mu.pred = as.numeric(SNR*K.pred.1 %*% G.inv %*% y)
    b1 <-  (t(y) %*% G.inv %*% y) / 2
    SIG.pred <- diag(as.numeric((2*b1/n)) * (Istar + SNR*K.pred - (SNR^2)*K.pred.1 %*% G.inv %*% K.1.pred))

    # mus[jj,] <- mu.pred
    # sigs[[jj,]] <- SIG.pred
    # jj = jj + 1
    list("mus" = mu.pred, "SIGS" = SIG.pred)
  }
  for(jj in 1:length(weights.inds)){
    mus[jj, ] <- out.loop[[jj]]$mus
    sigs[jj,] <- out.loop[[jj]]$SIGS
  }
  return(list("mus" = mus, "sigs" = sigs))
}

sketched_GP <- function(y, X, y.star, X.star, m = 60, K, SNRs, n.theta, prediction = FALSE){
  n <- nrow(X); p <- ncol(X)
  n.snr <- length(SNRs)
  n.models <- K*n.snr
  n_cores <- detectCores()
  registerDoParallel(cores = 8)

  # Sample K p x m sketching matrices
  Sketch.Mat.List <- list()
  for(ii in 1:K){
    Sketch.Mat.List[[ii]] <- get_sketch_mat(m, p, orthog = FALSE)
  }

  # I want some kind of matrix that tells me where to look for things, or their values
  models.mat <- matrix(0, nrow = 3, ncol = n.models)
  models.mat[1, ] <- rep(1:K, each = n.snr)
  models.mat[2, ] <- rep(SNRs, K)

  dist.mat = plgp::distance(X)
  dmax <- max(dist.mat)
  dmin <- min( dist.mat[dist.mat != 0] )

  Sketched.X.List <- list()
  mins.mat <- matrix(0, nrow = K, ncol = m)
  maxs.mat <- matrix(0, nrow = K, ncol = m)

  print("Sampling thetas")
  for(ii in 1:K){
  sketched.X <- X%*%t(Sketch.Mat.List[[ii]])
  #Psi.XX.l <- t(Psil %*% t(X))
  mins <- apply(sketched.X, 2, min)
  maxs <- apply(sketched.X, 2, max)
  sketched.X <- ( sketched.X - matrix(mins, nrow = n, ncol = m, byrow = TRUE) ) / (matrix(maxs, nrow = n, ncol = m, byrow = TRUE) - matrix(mins, nrow = n, ncol = m, byrow = TRUE))

  # Normalize the sketched X's to the unit cube before fitting GP
  mins.mat[ii, ] <- mins
  maxs.mat[ii, ] <- maxs
  Sketched.X.List[[ii]] <- sketched.X
    for(jj in 1:n.snr){
      # For each sketched.X, snr pair sample an appropriate theta
      ind <- (ii - 1)*n.snr + jj
      print(ind)
      snr = models.mat[2, ind]
      models.mat[3, ind] <- get_theta(y = y, sketched.X = Sketched.X.List[[ii]], dmin = dmin, dmax = dmax, N = n.theta, SNR = snr)

      # Test how changing the theta parameter changes the results
      # models.mat[3, ind] <- 0.1
    }
  }

  # Now we have all of the model parameters.
  # So we can get the model densities in parallel.
  # browser()
  print("Getting stacking densities")
  out.mat <- matrix(data = NA, nrow = 2*n, ncol = n.models)
  Densities.Mat <- matrix(data = NA, nrow = n, ncol = n.models)
  Preds.Mat <- matrix(data = NA, nrow = n, ncol = n.models)

  for(ii in 1:n.models){
    # browser()
    print(ii)
    n.mat = models.mat[1, ii]
    snr = models.mat[2, ii]
    theta = models.mat[3, ii]

    #sketched.X <- X%*%t(Sketch.Mat.List[[n.mat]])
    sketched.X <- Sketched.X.List[[n.mat]]
    out.mat[,ii] <- get_model_densities(y, sketched.X, theta, snr)
    Densities.Mat[,ii] <- out.mat[seq(from = 1, to = 2*n, by = 2), ii]
    Preds.Mat[,ii] <- out.mat[seq(from = 2, to = 2*n, by = 2), ii]
  }

  print("Solving stacking problem")
  # Solve for the stacking weights
  stack.weights <- stacking_weights(lpd_point = log(Densities.Mat))

  print("Fitting on X")
  fit.out <- get_preds_parallel(Xstar = X, X = X,
                                y = y,
                                Sketch.Mat.List = Sketch.Mat.List,
                                thetas = models.mat[3, ],
                                SNR = models.mat[2, ],
                                stack.weights = stack.weights,
                                Mins = mins.mat, Maxs= maxs.mat,
                                sketch.mat.num = models.mat[1, ])
  fit.mus <- fit.out$mus
  fit.sigs <- fit.out$sigs

  if(isTRUE(prediction)){
    print("Predicting on X*")
    preds.out <- get_preds_parallel(Xstar = X.star, X = X,
                       y = y,
                       Sketch.Mat.List = Sketch.Mat.List,
                       thetas = models.mat[3, ],
                       SNR = models.mat[2, ],
                       stack.weights = stack.weights,
                       Mins = mins.mat, Maxs= maxs.mat,
                       sketch.mat.num = models.mat[1, ])
    preds.mus <- preds.out$mus
    preds.sigs <- preds.out$sigs
    return(list("Densities.Mat" = Densities.Mat, "loo.Preds.mat"= Preds.Mat, "stack.weights" = stack.weights, "Models" = models.mat, "mu.pred" = preds.mus, "sig.pred" = preds.sigs, "mu.fit" = fit.mus, "sig.fit" = fit.sigs))
  }


  return(list("Densities.Mat" = Densities.Mat, "loo.Preds.mat"= Preds.Mat, "stack.weights" = stack.weights, "Models" = models.mat, "mu.fit" = fit.mus, "sig.fit" = fit.sigs))
}

# Stacking Weights -- cite this
stacking_weights <-function(lpd_point,
                            optim_method = "BFGS",
                            optim_control = list()) {

  stopifnot(is.matrix(lpd_point))
  N <- nrow(lpd_point)
  K <- ncol(lpd_point)
  if (K < 2) {
    stop("At least two models are required for stacking weights.")
  }

  exp_lpd_point <- exp(lpd_point)
  negative_log_score_loo <- function(w) {
    # objective function: log score
    stopifnot(length(w) == K - 1)
    w_full <- c(w, 1 - sum(w))
    sum <- 0
    for (i in 1:N) {
      sum <- sum + log(exp(lpd_point[i, ]) %*% w_full)
    }
    return(-as.numeric(sum))
  }

  gradient <- function(w) {
    # gradient of the objective function
    stopifnot(length(w) == K - 1)
    w_full <- c(w, 1 - sum(w))
    grad <- rep(0, K - 1)
    for (k in 1:(K - 1)) {
      for (i in 1:N) {
        grad[k] <- grad[k] +
          (exp_lpd_point[i, k] - exp_lpd_point[i, K]) / (exp_lpd_point[i,]  %*% w_full)
      }
    }
    return(-grad)
  }

  ui <- rbind(rep(-1, K - 1), diag(K - 1))  # K-1 simplex constraint matrix
  ci <- c(-1, rep(0, K - 1))
  w <- constrOptim(
    theta = rep(1 / K, K - 1),
    f = negative_log_score_loo,
    grad = gradient,
    ui = ui,
    ci = ci,
    method = optim_method,
    control = optim_control
  )$par

  wts <- structure(
    c(w, 1 - sum(w)),
    names = paste0("model", 1:K),
    class = c("stacking_weights")
  )

  return(wts)
}

get_preds <- function(Xstar, X, y, P.list, thetas, stacking.weights, SNRs){
  #browser()
  if(is.null(dim(Xstar))){Xstar <- matrix(Xstar, nrow = 1)}
  w.inds <- which(stacking.weights != 0)

  K = nrow(thetas)
  n = nrow(X)
  nstar = nrow(Xstar)
  n.snrs = length(SNRs)

  I <- diag(1, nrow = n, ncol = n)
  Istar <- diag(1, nrow = nstar, ncol = nstar)

  y <- as.matrix(y, nrow = length(y), ncol = 1)

  SIGS <- vector("list", length(w.inds))
  mus <- matrix(data = 0, nrow = length(w.inds), ncol = nstar)

  out.mat <- foreach(ii = 1:K, .combine = rbind) %dopar% {
  sketched.X <- t(P.list[[ii]] %*% t(X))
  sketched.Xstar <- t(P.list[[ii]] %*% t(Xstar))
  sigs.mat <- matrix(0, nrow = n.snrs, ncol = nstar)
  mus.mat <- matrix(0, nrow = n.snrs, ncol = nstar)
    for(jj in 1:n.snrs){
      C.old <- exp_kernel(sketched.X, theta = thetas[ii, jj])
      C.new.old <- exp_kernel(M1 = sketched.Xstar, M2 = sketched.X, theta = thetas[ii, jj])
      C.new <- exp_kernel(sketched.Xstar, theta = thetas[ii, jj])

      G.inv <- chol2inv(chol(SNRs[jj]*C.old + I))
      mu.pred <- SNRs[jj] * C.new.old %*% G.inv %*% y
      mus.mat[jj, ] <- mu.pred
      b1 = as.numeric((t(y) %*% G.inv %*% y) / 2)
      SIG.pred = diag((2*b1 / n) * (Istar + (SNRs[jj]* C.new) - (SNRs[jj]^2)*(C.new.old %*% (G.inv %*% t(C.new.old))) ))
      sigs.mat[jj, ] <- SIG.pred
      out.mat <- cbind(mus.mat, sigs.mat)
    }
  out.mat
  }

  mus.mat <- out.mat[,1:nstar]
  sigs.mat <- out.mat[,(nstar+1):(2*nstar)]

  return(list("mus" = mus.mat, "Sigs" = sigs.mat))
}


# What happens next?
