# Simulations for paper
library(BART)
library(glmnet)
library(tidyverse)
library(fields)
library(tseries)
library(foreach)
library(doParallel)
library(plgp)
library(splines)
library(randomForest)
library(GpGp)
library(GauPro)

source("R/functions.R")
# Compete with the following:
# 1. BART
# 2. LASSO
# 3. Random Forest
# 4. Regular GP
# 5. Neural Network???

# Repeat for n = 100, 500
#   for torus, swiss roll
#     for noise = (0.01, 0.03, 0.05, 0.1)
#       for

# Read in appropriate data set
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/TORI01.rds")

inds.train <- 1:100
inds.test <- 101:200
PREDS.skgp <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
PREDS.bart <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
PREDS.cbart <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
PREDS.lasso <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
PREDS.rf <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
PREDS.crf <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
PREDS.gp <-  matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))

SIGS.skgp <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
SIGS.gp <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))

LOWER.bart <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
LOWER.cbart <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
LOWER.rf <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
LOWER.crf <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))

UPPER.bart <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
UPPER.cbart <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
UPPER.rf <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))
UPPER.crf <- matrix(data = NA, nrow = length(DATA), ncol = length(inds.test))

for(ii in 1:length(DATA)){
  print(paste0("#################### ", ii, " ####################"))
  # Pull out data
  DD <- DATA[[ii]]

  X <- DD$X[inds.train,]; Xstar <- DD$X[inds.test,]
  y <- DD$y[inds.train]; ystar <- DD$y[inds.test]

  # Center the data
  X <- scale(X, center = TRUE, scale = FALSE); Xstar <- scale(Xstar, center = TRUE, scale = FALSE)
  y <- scale(y, center = TRUE, scale = FALSE); ystar <- scale(ystar, center = TRUE, scale = FALSE)

  # Get screening scores
  screening.scores <- get_NIS_scores(X, y)

  # Screen down to 1000
  n.keep = 2000
  covars.save <- order(screening.scores, decreasing = TRUE)[1:n.keep]

  # Get the screened X's
  X.s <- X[,covars.save]
  Xstar.s <- Xstar[,covars.save]

  # Compress the screened X's for some of the methods
  P1 <- get_sketch_mat(m = 60, p = n.keep)
  Sketched.X.s <- X.s%*%t(P1)
  Sketched.Xstar.s <- Xstar.s%*% t(P1)

  # Fit the different methods
  # Sketched GP
  SkGP <- sketched_GP(y, X.s, ystar, Xstar.s, m = 60, K = 20,
                    n.theta = 10, n.snrs = 25, prediction = TRUE,
                    snr.method = "sample", snr.max = 10,
                    stacking.method = "LOO",n.folds = 10)
  stack.weights <- SkGP$stack.weights
  mu.pred <- t(SkGP$mu.pred)
  sig.pred <- t(SkGP$sig.pred)

  skgp.preds <- rowSums(mu.pred*matrix(data = unname(stack.weights), nrow = nrow(mu.pred), ncol = ncol(mu.pred), byrow = TRUE))
  skgp.sigs <- sqrt(rowSums(sig.pred*matrix(data = unname(stack.weights), nrow = nrow(sig.pred), ncol = ncol(sig.pred), byrow = TRUE)))

  # BART
  print("BART")
  BART <- mc.wbart(X.s, y, Xstar.s, ndpost = 1000, nskip = 1000, mc.cores = 6, seed = 113)
  bart.preds <- BART$yhat.test.mean
  bart.quantiles <- apply(t(BART$yhat.test), 1, quantile, probs = c(0.025, 0.975))
  # How to get coverage?

  # Compressed BART
  print("CBART")
  CBART <- mc.wbart(Sketched.X.s, y, Sketched.Xstar.s, ndpost = 1000, nskip = 1000, mc.cores = 6, seed = 113)
  cbart.preds <- CBART$yhat.test.mean
  cbart.quantiles <- apply(t(CBART$yhat.test), 1, quantile, probs = c(0.025, 0.975))
  # How to get coverage?

  # LASSO
  print("LASSO")
  LASSO <- cv.glmnet(x=X.s, y=y, alpha = 1)
  l.lasso.min <- LASSO$lambda.min
  lasso.model <- glmnet(x=X.s, y=y,
                      alpha  = 1,
                      lambda = l.lasso.min)

  lasso.preds <- predict(lasso.model, Xstar.s )

  # RANDOM FOREST
  print("RF")
  RF <- randomForest(x=X.s, y=y, ytest = ystar, ntree = 500, mtry = sqrt(ncol(X.s)))
  rf.preds.obj <- predict(RF, Xstar.s, predict.all = TRUE)
  rf.preds <- rf.preds.obj$aggregate
  rf.quantiles <- apply(rf.preds.obj$individual, 1, quantile, probs = c(0.025, 0.975))

  # Compressed RANDOM FOREST
  print("CRF")
  CRF <- randomForest(x=Sketched.X.s, y=y, ntree = 500, mtry = sqrt(ncol(Sketched.X.s)))
  crf.preds.obj <- predict(CRF, Sketched.Xstar.s, predict.all = TRUE)
  crf.preds <- crf.preds.obj$aggregate
  crf.quantiles <- apply(crf.preds.obj$individual, 1, quantile, probs = c(0.025, 0.975))

  # Regular GP
  # No compression but same GP parameters as before
  # Exponential kernel with SNR = 1? and sample thetas...
  print("GP")
  dist.mat = plgp::distance(X.s)
  dmax <- max(dist.mat)
  dmin <- min( dist.mat[dist.mat != 0] )

  theta <- get_theta(y, X.s, dmin, dmax, 100, 1)
  K1 <- exp_kernel(X.s, theta = theta)
  K.1.pred <- exp_kernel(X.s, Xstar.s, theta = theta); K.pred.1 <- t(K.1.pred)
  K.pred <- exp_kernel(Xstar.s, theta = theta)

  SNR = 1; n <- nrow(X.s); nstar <- nrow(Xstar.s)
  I <- diag(rep(1, n)); Istar <- diag(rep(1, nstar))
  G <- SNR*K1 + I
  G.chol <- chol(G)
  G.inv <- chol2inv(G.chol)

  gp.preds = as.numeric(SNR*K.pred.1 %*% G.inv %*% y)
  b1 <-  (t(y) %*% G.inv %*% y) / 2
  gp.sigs <- diag(as.numeric((2*b1/n)) * (Istar + SNR*K.pred - (SNR^2)*K.pred.1 %*% G.inv %*% K.1.pred))


  # Save everything out
  PREDS.skgp[ii, ] <- skgp.preds
  PREDS.bart[ii, ] <- bart.preds
  PREDS.cbart[ii, ] <- cbart.preds
  PREDS.lasso[ii, ] <- lasso.preds
  PREDS.rf[ii, ] <- rf.preds
  PREDS.crf[ii, ] <- crf.preds
  PREDS.gp[ii, ] <- gp.preds

  SIGS.skgp[ii, ] <- skgp.sigs
  SIGS.gp[ii, ] <- gp.sigs

  LOWER.bart[ii, ] <- bart.quantiles[1, ]
  LOWER.cbart[ii, ] <- cbart.quantiles[1, ]
  LOWER.rf[ii, ] <- rf.quantiles[1, ]
  LOWER.crf[ii, ] <- crf.quantiles[1, ]

  UPPER.bart[ii, ] <- bart.quantiles[2, ]
  UPPER.cbart[ii, ] <- cbart.quantiles[2, ]
  UPPER.rf[ii, ] <- rf.quantiles[2, ]
  UPPER.crf[ii, ] <- crf.quantiles[2, ]
}

# Save results out
RESULTS <- list("PREDS.skgp" = PREDS.skgp, "PREDS.gp" = PREDS.gp,"PREDS.bart" = PREDS.bart, "PREDS.cbart" = PREDS.cbart,
                "PREDS.lasso" = PREDS.lasso, "PREDS.rf" = PREDS.rf, "PREDS.crf" = PREDS.crf,
                "SIGS.skgp" = SIGS.skgp, "SIGS.gp" = SIGS.gp,
                "UPPER.bart" = UPPER.bart, "UPPER.cbart" = UPPER.cbart, "UPPER.rf" = UPPER.rf, "UPPER.crf" = UPPER.crf,
                "LOWER.bart" = LOWER.bart, "LOWER.cbart" = LOWER.cbart, "LOWER.rf" = LOWER.rf, "LOWER.crf" = LOWER.crf,
                "inds.test" = inds.test, "inds.train" = inds.train)

# Make sure you have the right name...
saveRDS(RESULTS, "Results/Sims_TORI01_p2k.rds")

################################################################################
# Extras #######################################################################
################################################################################
Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_ROLLS005_p10K.rds")
# DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/TORI003.rds")
attach(Sims)
NN <- length(DATA)
pmse.skgp <- numeric(length = NN)
pmse.bart <- numeric(length = NN)
pmse.cbart <- numeric(length = NN)
pmse.lasso <- numeric(length = NN)
pmse.rf <- numeric(length = NN)
pmse.crf <- numeric(length = NN)
pmse.gp <- numeric(length = NN)

coverage.skgp <- numeric(length = NN)
coverage.bart <- numeric(length = NN)
coverage.cbart <- numeric(length = NN)
coverage.rf <- numeric(length = NN)
coverage.crf <- numeric(length = NN)
coverage.gp <- numeric(length = NN)

par(mfrow = c(1, 3))
for(ii in 1:NN){
  y <- scale(DATA[[ii]]$y[inds.test], center = TRUE, scale = FALSE)

  pmse.skgp[ii] <- mean((PREDS.skgp[ii, ] - y)^2)
  pmse.bart[ii] <- mean((PREDS.bart[ii, ] - y)^2)
  pmse.cbart[ii] <- mean((PREDS.cbart[ii, ] - y)^2)
  pmse.lasso[ii] <- mean((PREDS.lasso[ii, ] - y)^2)
  pmse.rf[ii] <- mean((PREDS.rf[ii, ] - y)^2)
  pmse.crf[ii] <- mean((PREDS.crf[ii, ] - y)^2)
  pmse.gp[ii] <- mean((PREDS.gp[ii, ] - y)^2)

  coverage.skgp[ii] <- sum((y > (PREDS.skgp[ii, ] - 2*SIGS.skgp[ii, ])) & (y < (PREDS.skgp[ii, ] + 2*SIGS.skgp[ii, ]))) / length(y)
  coverage.gp[ii] <- sum((y > (PREDS.gp[ii, ] - 2*SIGS.gp[ii, ])) & (y < (PREDS.gp[ii, ] + 2*SIGS.gp[ii, ]))) / length(y)
  coverage.bart[ii] <- sum((y > LOWER.bart[ii]) & (y < UPPER.bart[ii])) / length(y)
  coverage.cbart[ii] <- sum((y > LOWER.cbart[ii]) & (y < UPPER.cbart[ii])) / length(y)
  coverage.rf[ii] <- sum((y > LOWER.rf[ii]) & (y < UPPER.rf[ii])) / length(y)
  coverage.crf[ii] <- sum((y > LOWER.crf[ii]) & (y < UPPER.crf[ii])) / length(y)
}

boxplot(cbind(pmse.skgp, pmse.gp, pmse.bart, pmse.cbart, pmse.lasso, pmse.rf, pmse.crf),
        labels = c("SkGP", "BART", "CBART", "LASSO", "RF", "CRF"))
boxplot(cbind(coverage.skgp, coverage.gp, coverage.bart, coverage.cbart, coverage.rf, coverage.crf))

ii = 1
Xstar <- scale(DATA[[ii]]$X[inds.test,], center = TRUE, scale = FALSE)
ystar <- scale(DATA[[ii]]$y[inds.test], center = TRUE, scale = FALSE)
ord = order(Xstar[,2])
plot(Xstar[ord,2], ystar[ord])
lines(Xstar[ord,2], PREDS.skgp[ii, ord], col = 2, lwd = 2)
lines(Xstar[ord,2], PREDS.gp[ii, ord], col = 3, lwd = 2)
lines(Xstar[ord,2], PREDS.bart[ii, ord], col = 4, lwd = 2)
lines(Xstar[ord,2], PREDS.cbart[ii, ord], col = 5, lwd = 2)
lines(Xstar[ord,2], PREDS.lasso[ii, ord], col = 6, lwd = 2)
lines(Xstar[ord,2], PREDS.rf[ii, ord], col = 7, lwd = 2)
lines(Xstar[ord,2], PREDS.crf[ii, ord], col = 8, lwd = 2)
legend("topleft", col = c(1,2,3,4,5,6,7,8), pch = rep(19, 8),
       legend = c("Truth", "SkGP", "GP", "BART", "CBART","LASSO", "RF", "CRF"))


################################################################################
# Messing Around

ii = 1
DD <- DATA[[ii]]

X <- DD$X[inds.train,]; Xstar <- DD$X[inds.test,]
y <- DD$y[inds.train]; ystar <- DD$y[inds.test]

# X <- DD$X
# y <- DD$y
# Center the data
X <- scale(X, center = TRUE, scale = FALSE); Xstar <- scale(Xstar, center = TRUE, scale = FALSE)
y <- scale(y, center = TRUE, scale = FALSE); ystar <- scale(ystar, center = TRUE, scale = FALSE)

screening.scores <- get_NIS_scores(X, y)

# Screen down to 1000
n.keep = 2500
covars.save <- order(screening.scores, decreasing = TRUE)[1:n.keep]

# Get the screened X's
X.s <- X[,covars.save]
Xstar.s <- Xstar[,covars.save]

library(GPvecchia)
help("vecchia_specify")

vech.spec <- GPvecchia::vecchia_specify(X.s, m = 30)

GPvecchia::vecchia_prediction(vecchia.approx = vech.spec, z,covparms = c(1,2,0.5),nuggets = 1)

help("createU")
U <- GPvecchia::createU(vech.spec, covparms = c(1, 1 , 1), nuggets = 0, covmodel = "")
V <- as.matrix(Matrix(U2V(U), sparse = FALSE))


## function to reverse-order a matrix
revMat=function(mat) mat[nrow(mat):1,ncol(mat):1,drop=FALSE]

######  compute V for posterior inference   #######

U2V=function(U.obj){
  ### when changing this function make sure it returns a dtCMatrix!
  ### Otherwise solve in the parent function will be very slow

  U.y=U.obj$U[U.obj$latent,]

  if(U.obj$cond.yz=='zy') {

    V.ord=revMat(U.y[,U.obj$latent,drop=FALSE])

  } else if(U.obj$ord.pred!='obspred'){

    W=Matrix::tcrossprod(U.y)
    W.rev=revMat(W)

    if(U.obj$ic0){
      V.ord=Matrix::t(ichol(W.rev))
    } else {
      V.ord=Matrix::t(Matrix::chol(W.rev))
    }

    V.ord = methods::as(V.ord, 'dtCMatrix')


  } else {  # for obspred ordering

    last.obs=max(which(!U.obj$latent))
    latents.before=sum(U.obj$latent[1:last.obs])
    latents.after=sum(U.obj$latent[-(1:last.obs)])

    # pred columns are unchanged
    V.pr=revMat(U.y[,(last.obs+1):ncol(U.y),drop=FALSE])

    # have to compute cholesky for obs block
    U.oo=U.y[1:latents.before,1:last.obs]
    A=Matrix::tcrossprod(U.oo)
    A.rev=revMat(A)
    if(U.obj$ic0){ V.oor=Matrix::t(ichol(A.rev))
    }     else   V.oor=Matrix::t(Matrix::chol(A.rev))

    # combine the blocks into one matrix
    zeromat.sparse=Matrix::sparseMatrix(c(),c(),dims=c(latents.after,latents.before))
    V.or=rbind(zeromat.sparse,V.oor)

    V.ord=methods::as(cbind(V.pr,V.or),'dtCMatrix')
  }
  return(V.ord)
}
