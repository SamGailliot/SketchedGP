# Real Data for Real

# Required Libraries
library(tidyverse)
library(fields)
library(tseries)
library(foreach)
library(doParallel)
library(plgp)
library(splines)
library(dplyr)
library(zoo)
library(BART)
library(Matrix)
library(intRinsic)
library(mcmc)
source("R/functions.R")

# Functions
vectorizor <- function(data){
  # Create a matrix which will hold the data
  # Make sure that you check the correct dimensions
  data.mat <- matrix(NA,  nrow = dim(data)[1], ncol = dim(data)[2]*dim(data)[3])
  # Loop over n or t or whatever
  for(ii in 1:dim(data)[1]){
    # Vectorize the appropriate slice
    # c(t(.)) does it by row since c(.) is by column
    data.mat[ii, ] <- c(t(data[ii,,]))
  }
  return(data.mat)
}

make_screened_image_plots <- function(X, dt = 1, inds, covars.save){
  image.plot(t(matrix(X[1, inds], 128, 128)))
  KEEPERS <- inds %in% covars.save
  X.KEPT <- X[1, inds]
  X.KEPT[!KEEPERS] <- NA
  image.plot(t(matrix(X.KEPT, 128, 128)))
}

SkGP_Unpack <- function(SkGP.fit, y.train, y.test){
  stack.weights <- SkGP.fit$stack.weights
  mu.pred <- t(SkGP.fit$mu.pred)
  sig.pred <- t(SkGP.fit$sig.pred)
  mu.fit <- t(SkGP.fit$mu.fit)
  sig.fit <- t(SkGP.fit$sig.fit)

  stack.mod.pred <- rowSums(mu.pred*matrix(data = unname(stack.weights), nrow = nrow(mu.pred), ncol = ncol(mu.pred), byrow = TRUE))
  stack.sd.pred <-  sqrt(rowSums(sig.pred*matrix(data = unname(stack.weights)^2, nrow = nrow(sig.pred), ncol = ncol(sig.pred), byrow = TRUE)))

  stack.mod.fit <- rowSums(mu.fit*matrix(data = unname(stack.weights), nrow = nrow(mu.fit), ncol = ncol(mu.fit), byrow = TRUE))
  stack.sd.fit <-  sqrt(sqrt(rowSums(sig.fit*matrix(data = unname(stack.weights)^2, nrow = nrow(sig.fit), ncol = ncol(sig.fit), byrow = TRUE))))

  n.train <- length(y.train)
  n.test <- length(y.test)

  cvg.fit <- sum((y.train > (stack.mod.fit - 1.96*stack.sd.fit)) & (y.train < (stack.mod.fit + 1.96*stack.sd.fit))) / length(y.train)
  cvg.pred <- sum((y.test > (stack.mod.pred - 1.96*stack.sd.pred)) & (y.test < (stack.mod.pred + 1.96*stack.sd.pred))) / length(y.test)

  length.fit <- mean((stack.mod.fit + 1.96*stack.sd.fit) - (stack.mod.fit - 1.96*stack.sd.fit))
  length.pred <- mean(2*(1.96*stack.sd.pred))

  mse <- mean((y.train - stack.mod.fit)^2)
  pmse <- mean((y.test - stack.mod.pred)^2)

  return(list("yhat" = stack.mod.fit, "yhat.pred" = stack.mod.pred,
              "CI" = stack.sd.fit, "CI.pred" = stack.sd.pred,
              "cvg" = cvg.fit, "cvg.pred" = cvg.pred,
              "avg.CI.length" = length.fit, "avg.CI.length.pred" = length.pred,
              "mse" = mse, "pmse" = pmse))
}

SNR_Search_Debugger <- function(X, y,n.theta, n.snr, snr.max, m){
  #set.seed(69420)
  m = m; N = 100; n = nrow(X)
  sketched.X <- X%*%t(get_sketch_mat(m, ncol(X)))

  mins <- apply(sketched.X, 2, min)
  maxs <- apply(sketched.X, 2, max)
  sketched.X <- ( sketched.X - matrix(mins, nrow = n, ncol = m, byrow = TRUE) ) / (matrix(maxs, nrow = n, ncol = m, byrow = TRUE) - matrix(mins, nrow = n, ncol = m, byrow = TRUE))

  dist.mat = plgp::distance(X)
  dmax <- max(dist.mat)
  dmin <- min( dist.mat[dist.mat != 0] )

# dist.mat <- plgp::distance(sketched.X)
# dmax <- max(dist.mat)
# dmin <- min(dist.mat[dist.mat != 0])

  n <- nrow(sketched.X)
  Identity_Mat <- diag(1, n,n)

  # The double loop below is what I need for the (theta, snr) sampling procedure
  #n.theta <- 10
  #n.snr <- 25
  thetas <- seq(3/dmax, 3/dmin, length.out = n.theta)
  SNR = c(0.01, 0.1, 0.5, seq(1, snr.max, length.out  = n.snr))
  # SNR = c(0.001, 0.01, 0.1, 0.5, 10^seq(1, log(snr.max, 10), length.out = n.snr))

  logdets <- numeric(length(SNR))
  quadratic.exp <- numeric(length(SNR))
  log.outs <- matrix(NA, length(thetas), length(SNR))
  logdets.out <- matrix(NA, length(thetas), length(SNR))
  quad.out <- matrix(NA, length(thetas), length(SNR))

  for(ii in 1:length(thetas)){
  print(ii)
  log.out <- foreach(jj = 1:length(SNR), .combine = c)%dopar%{
    print(jj)
    C <- exp_kernel(sketched.X, theta = thetas[ii])
    G <- (SNR[jj]*C) + Identity_Mat + diag(sqrt(.Machine$double.eps), nrow(C))
    Gchol <- chol(G)
    Ginv <- chol2inv(chol(G))
    logdets[jj] <- 2*sum(log(diag(Gchol)))
    quadratic.exp[jj] <- log(colSums(y * (Ginv %*% y)))
    return(logdets[jj] + (n)*quadratic.exp[jj])
  }
  log.outs[ii, ] <- log.out
  logdets.out[ii, ] <- logdets
  quad.out[ii, ]  <- quadratic.exp
  }

  min.ind <- which(log.outs == min(log.outs), arr.ind = TRUE)

  par(mfrow = c(1,1))
  plot(SNR, log.outs[1, ], "l", ylim = range(log.outs), xlab = "SNR", ylab = "Objective Function",
     main = paste0("Best Theta = ", round(thetas[min.ind[1]], 5), "; Best SNR = ", round(SNR[min.ind[2]], 3)))
  for(ii in 2:length(thetas)){
    lines(SNR, log.outs[ii, ], col = ii)
  }
  lines(SNR, log.outs[min.ind[1],], col = "red", lwd = 2)
  abline(v = SNR[min.ind[2]], col = "red", lty = 2, lwd = 2)

  contour(SNR,thetas, t(log.outs), xlab = "SNR", ylab = "Theta")
}


################################################################################
# Load in Data #################################################################
################################################################################
setwd("/Users/samuelgailliot/Desktop/Research/SketchedGP/")
load("~/Desktop/Research/SketchedGP/Data/pm_images/image_data.rda")
# Load in the response readings
pmFRM <- read_csv("~/Desktop/Research/SketchedGP/Data/pm_images/pmFRM.csv")
dates <- ymd(pmFRM$date)

BLUE <- vectorizor(imgs_a[,,,1]); which.nonzero1 <- which(colSums(BLUE)!=0); which.zero1 <- which(colSums(BLUE)==0)
BLUE <- BLUE[,which.nonzero1]
GREEN <- vectorizor(imgs_a[,,,2]); which.nonzero2 <- which(colSums(GREEN)!=0); which.zero2 <- which(colSums(GREEN)==0)
GREEN <- GREEN[,which.nonzero2]
RED <- vectorizor(imgs_a[,,,3]); which.nonzero3 <- which(colSums(RED)!=0); which.zero3 <- which(colSums(RED)==0)
RED <- RED[,which.nonzero3]
INFRA <- vectorizor(imgs_a[,,,4]); which.nonzero4 <- which(colSums(INFRA)!=0); which.zero4 <- which(colSums(INFRA)==0)
INFRA <- INFRA[,which.nonzero4]
# DIFF <- vectorizor(imgs_a[,,,3] - imgs_a[,,,4]); which.nonzero5 <- which(colSums(DIFF)!=0)

# Row normalize the channels
row_norm <- function(X){
  maxs <- apply(X, 1, max)
  mins <- apply(X, 1, min)
  n <- nrow(X); p <- ncol(X)

  MAXS <- matrix(maxs, nrow = n, ncol = p, byrow = FALSE)
  MINS <- matrix(mins, nrow = n, ncol = p, byrow = FALSE)

  X <- (X - MINS) / (MAXS - MINS)
  return(X)
}

row_stand <- function(X){
  means <- apply(X, 1, mean)
  sds <- apply(X, 1, sd)
  n <- nrow(X); p <- ncol(X)
  MEANS <- matrix(means, nrow = n, ncol = p, byrow = FALSE)
  SDS <- matrix(sds, nrow = n, ncol = p, byrow = FALSE)
  X <- (X - MEANS)/SDS
  return(X)
}

# BLUE  <- row_norm(BLUE)
# GREEN <- row_norm(GREEN)
# RED   <- row_norm(RED)
# INFRA  <- row_norm(INFRA)

BLUE  <- row_stand(BLUE)
GREEN <- row_stand(GREEN)
RED   <- row_stand(RED)
INFRA  <- row_stand(INFRA)

# X <- t(scale(t(cbind(INFRA, RED, BLUE, GREEN))))
X <- cbind(INFRA, RED, BLUE, GREEN)
# X[is.nan(X)] <- 0
y <- pmFRM$pmFRM
dates <- ymd(pmFRM$date)

# library(intRinsic)
# hid <- intRinsic::Hidalgo(INFRA)
#
# DistMat.INFRA <- plgp::distance(INFRA)
# DISTMat.X <- plgp::distance(X)
#
#
# summary(twonn(dist_mat = DISTMat.X, alpha = 0.99))
# summary(twonn(dist_mat = DistMat.INFRA, alpha = 0.99))
################################################################################
# Analysis #####################################################################
################################################################################
# Split the Data
set.seed(01131995)
n.train <- 900
n.test <- 434
inds.train <- 1:n.train
inds.test <- (n.train+1):(n.train + n.test)
n.tot <- n.test + n.train

inds.tot <- (1:length(dates))[which(dates >= "2019-01-01")]

inds.test <- seq(4, length(inds.tot), 4)
inds.train <- 1:length(inds.tot)[-inds.test]
dates.tot <- dates[inds.tot]
dates.test <- dates[inds.test]
dates.train <- dates[inds.train]

X <- X[inds.tot, ]
y.orig <- y[inds.tot]
y <- scale(log(y.orig))

X.train <- X[inds.train, ]; y.train <- y[inds.train];
X.test <- X[inds.test, ]; y.test <- y[inds.test];

# Center and scale the data?
################################################################################
# Screening ####################################################################
################################################################################
quant.thresh = 0.99
screening.scores <- get_NIS_scores(X.train, y.train)

#Permutation Test from NIS paper
#Doesnt do anything on the full data
perm <- sample(1:nrow(X.train), size = nrow(X.train), replace = FALSE)
screening.scores.perm <- get_NIS_scores(X.train[perm, ], y.train)

par(mfrow = c(1,1))
plot(1:ncol(X.train),screening.scores)
points(1:ncol(X.train), screening.scores, col = ifelse(screening.scores < quantile(screening.scores.perm, quant.thresh), "black", "red"))
abline(h = quantile(screening.scores.perm, quant.thresh), col = "red", lty = 2, lwd = 3)
points(1:ncol(X.train), screening.scores.perm, col = "blue")

# Create an average image that is based on all 4 wavelengths. Uses the screening
# Scores as weights.
# DD <- dim(BLUE)[2]
# inds1 <- 1:(DD)
# inds2 <- (DD+1):(2*DD)
# inds3 <- (2*DD+1):(3*DD)
# inds4 <- (3*DD+1):(4*DD)
# P1 <- rbind(Matrix(diag(screening.scores[inds1])),
#             Matrix(diag(screening.scores[inds2])),
#             Matrix(diag(screening.scores[inds3])),
#             Matrix(diag(screening.scores[inds4])))
#
# P1 <- P1/colSums(P1)
# Xavg.test <- as.matrix(X[inds.test, ] %*% P1)
# Xavg.train <- as.matrix(X[inds.train,] %*% P1)
#
# # Do another round of screening?
# quant.thresh <- 0.99
# screening.scores.avg <- get_NIS_scores(Xavg.train, y.train)
# screening.scores.perm.avg <- get_NIS_scores(Xavg.train[perm, ], y.train)
# # covars.save.avg <- order(screening.scores.avg, decreasing = TRUE)[1:8267]
# covars.save.avg <- which(screening.scores.avg >= quantile(screening.scores.perm.avg, quant.thresh))

# par(mfrow = c(1,1))
# plot(1:ncol(Xavg.train),screening.scores.avg)
# points(1:ncol(Xavg.train), screening.scores.avg, col = ifelse(screening.scores.avg < quantile(screening.scores.perm.avg, quant.thresh), "black", "red"))
# abline(h = quantile(screening.scores.perm, quant.thresh), col = "red", lty = 2, lwd = 3)
# points(1:ncol(Xavg.train), screening.scores.perm.avg, col = "blue")

# How many covariates/locations to keep in original data?
# n.keep = 33068
# covars.save <- order(screening.scores, decreasing = TRUE)[1:n.keep]

quant.thresh <- 0.99
covars.save <- which(screening.scores >= quantile(screening.scores.perm, quant.thresh))


# Plots
# par(mfrow = c(4,2))
# make_screened_image_plots(X, dt = 1, inds = inds1, covars.save = covars.save)
# make_screened_image_plots(X, dt = 1, inds = inds2, covars.save = covars.save)
# make_screened_image_plots(X, dt = 1, inds = inds3, covars.save = covars.save)
# make_screened_image_plots(X, dt = 1, inds = inds4, covars.save = covars.save)

# make_screened_image_plots(Xavg.test, dt = 1, inds = inds1, covars.save = covars.save.avg)
# make_screened_image_plots(X, dt = 1, inds = inds2, covars.save = covars.save)

# Screen the Xs
# Xavg.test <- Xavg.test[,covars.save.avg]
# Xavg.train <- Xavg.train[,covars.save.avg]
X.test <- X.test[,covars.save]
X.train <- X.train[, covars.save]
################################################################################
# Fit my method on both sets! ##################################################
################################################################################
# Find good SNR values        ##################################################
################################################################################
# par(mfrow = c(1, 2))
# SNR_Search_Debugger(Xavg.train, y.train, n.theta = 20, n.snr = 75, snr.max = 100)
SNR_Search_Debugger(X.train,  y.train, n.theta = 10, n.snr = 100, snr.max = 10^4, m = 60)
# SNR_Search_Debugger(INFRA,  y.train, n.theta = 3, n.snr = 50, snr.max = 10^6)

# SkGP.fit.avg <- sketched_GP(y.train, Xavg.train, y.test, Xavg.test,
#                             m = 100, K = 5, SNRs = c(10000 ,20000, 30000),
#                             n.theta = 50, n.snrs = 25, prediction = TRUE, snr.method = "set",
#                             stacking.method = "kfold", n.folds = 20,
#                             snr.max = 40000)
# res.avg <- SkGP_Unpack(SkGP.fit.avg, y.train, y.test)

###### SKGP ######
SkGP.fit.full <- sketched_GP(y.train, X.train, y.test, X.test,
                   m = 100, K = 15, SNRs = c(100000,200000),
                   n.theta = 100, n.snrs = 50, prediction = TRUE, snr.method = "set",
                   stacking.method = "kfold", n.folds = 20,
                   snr.max = 40000)
res.full <- SkGP_Unpack(SkGP.fit.full, y.train, y.test)
#

# best.res <- list("SkGP.results" = res.full,"BART.results" = BART, "y" = y.test, "dates" = dates.test)
# saveRDS(best.res, "Results/best_results.rds")

###### BART ######
BART <- mc.wbart(X.train, y.train, X.test, ndpost = 1000, nskip = 1000, mc.cores = 6, seed = 113)
bart.preds <- BART$yhat.test.mean
bart.quantiles <- apply(t(BART$yhat.test), 1, quantile, probs = c(0.025, 0.975))
cvg.bart <- sum( (y.test > (bart.quantiles[1, ]) ) & ( y.test < (bart.quantiles[2,]) )) / length(y.test)

# BART.avg <- mc.wbart(Xavg.train, y.train, Xavg.test, ndpost = 1000, nskip = 1000, mc.cores = 6, seed = 113)
# bart.avg.preds <- BART.avg$yhat.test.mean
# bart.avg.quantiles <- apply(t(BART.avg$yhat.test), 1, quantile, probs = c(0.025, 0.975))
# cvg.bart.avg <- sum( (y.test > (bart.avg.quantiles[1, ]) ) & ( y.test < (bart.avg.quantiles[2,]) )) / length(y.test)

# print("RF")
# RF <- randomForest(x=X.train, y=y.train, ntree = 1000, mtry = sqrt(ncol(X.train)))
# rf.preds.obj <- predict(RF, X.test, predict.all = TRUE)
# rf.preds <- rf.preds.obj$aggregate
# rf.quantiles <- apply(rf.preds.obj$individual, 1, quantile, probs = c(0.025, 0.975))



P <- get_sketch_mat(m = 60, p = ncol(X.train))
Sketched.X.train <- X.train%*%t(P)
Sketched.X.test <- X.test%*% t(P)

###### CBART ######
CBART <- mc.wbart(Sketched.X.train, y.train, Sketched.X.test, ndpost = 1000, nskip = 1000, mc.cores = 6, seed = 113)
cbart.preds <- CBART$yhat.test.mean
cbart.quantiles <- apply(t(CBART$yhat.test), 1, quantile, probs = c(0.025, 0.975))
cvg.cbart <- sum( (y.test > (cbart.quantiles[1, ]) ) & ( y.test < (cbart.quantiles[2,]) )) / length(y.test)

CRF <- randomForest::randomForest(x=Sketched.X.train, y=y.train, ntree = 1000, mtry = sqrt(ncol(Sketched.X.train)))
crf.preds.obj <- predict(CRF, Sketched.X.test, predict.all = TRUE)
crf.preds <- crf.preds.obj$aggregate
crf.quantiles <- apply(crf.preds.obj$individual, 1, quantile, probs = c(0.025, 0.975))


pmse.bart <- mean((bart.preds - y.test)^2)
mse.bart <- mean((BART$yhat.train.mean - y.train)^2)
#pmse.bart.avg <- mean((bart.avg.preds - y.test)^2)
#mse.bart.avg <- mean((bart.avg.preds - y.test)^2)
pmse.cbart <- mean((cbart.preds - y.test)^2)
pmse.crf <- mean((crf.preds - y.test)^2)
#mse.cbard <- mean((cbart))

pmse.bart
pmse.cbart
pmse.crf
#pmse.bart.avg
pmse.full <- res.full$pmse; pmse.full
pmse.full
mse.full <- res.full$mse; mse.full
#pmse.avg <- res.avg$pmse; pmse.avg
#mse.avg <- res.avg$mse; mse.avg

#res.avg$pmse

res.full$CI.pred
res.full$cvg.predsum( (y.test > (bart.quantiles[1, ]) ) & ( y.test < (bart.quantiles[2,]) )) / length(y.test)
sum( (y.test > (cbart.quantiles[1, ]) ) & ( y.test < (cbart.quantiles[2,]) )) / length(y.test)
sum( (y.test > (crf.quantiles[1, ]) ) & ( y.test < (crf.quantiles[2,]) )) / length(y.test)

min.preds <- apply(SkGP.fit.full$mu.pred, MARGIN = 2, FUN = min)
max.preds <- apply(SkGP.fit.full$mu.pred, MARGIN = 2, FUN = max)
sum( (y.test > (min.preds) ) & ( y.test < (max.preds) )) / length(y.test)

par(mfrow = c(1,1))
plot(dates.train, y.train, main = "Fit to Training Data")
legend("topleft", legend = c(paste0("BART: MSPE = ", round(mse.bart, 3)),
                             paste0("SkGp: ", round(mse.full, 3))),
       col = c("darkblue", "maroon"), lty = c(1,1), lwd = c(2,2)
      )
lines(dates.train, BART$yhat.train.mean, col = "darkblue", lwd = 2)
lines(dates.train, res.full$yhat, col = "maroon", lwd = 2)
#lines(dates.train, res.avg$yhat, col = "lightskyblue", lwd = 2)

plot(dates.test, y.test, main = "Full Test Data", pch = 19)
legend("topleft", legend = c(paste0("SkGp: ", round(pmse.full, 3))),
       col = c("maroon"), lty = c(1), lwd = c(2)
)
#lines(dates.test, bart.preds, col = "darkblue", lwd = 2, lty = 1)
lines(dates.test, res.full$yhat.pred, col = "maroon", lwd = 2, lty = 1)
lines(dates.test, res.full$yhat.pred - 1.96*res.full$CI.pred, col = "darkred", lwd = 1.5, lty = 2)
lines(dates.test, res.full$yhat.pred + 1.96*res.full$CI.pred, col = "darkred", lwd = 1.5, lty = 2)
#lines(dates.test, cbart.preds)

plot(dates.test, y.test, main = "Full Test Data", pch = 19)
legend("topleft", legend = c(paste0("SkGp: ", round(pmse.full, 3))),
       col = c("maroon"), lty = c(1), lwd = c(2)
)
#lines(dates.test, bart.preds, col = "darkblue", lwd = 2, lty = 1)
lines(dates.test, res.full$yhat.pred, col = "maroon", lwd = 2, lty = 1)
lines(dates.test, bart.preds, col = 1, lwd = 2, lty = 1)
lines(dates.test, cbart.preds, col = 2, lwd = 2, lty = 1)
lines(dates.test, crf.preds, col = 3, lwd = 2, lty = 1)


library(RColorBrewer)

cols <- brewer.pal(8, "Dark2")
pdf(file = "Plots/Real_Data_forcast_Preds.pdf", width = 30, height = 30)
par(mfrow = c(3, 1))
par(mar=c(6, 6, 3, 6))
# SkGP
plot(dates.test, y.test, main = "SkGP", pch = 19,
     ylab = "PM", xlab = NA, cex.lab = 2, cex.axis = 2, cex.main = 2)
lines(dates.test, res.full$yhat.pred, col = cols[1], lwd = 2, lty = 1)
lines(dates.test, res.full$yhat.pred - 1.96*res.full$CI.pred , col = cols[1], lwd = 1, lty = 2)
lines(dates.test, res.full$yhat.pred + 1.96*res.full$CI.pred , col = cols[1], lwd = 1, lty = 2)


plot(dates.test, y.test, main = "BART", pch = 19,
     ylab = "PM", xlab = NA, cex.lab = 2, cex.axis = 2, cex.main = 2)
lines(dates.test, bart.preds, col = cols[3], lwd = 2, lty = 1)
lines(dates.test, bart.quantiles[1, ], col = cols[3], lwd = 1, lty = 2)
lines(dates.test, bart.quantiles[2, ], col = cols[3], lwd = 1, lty = 2)

plot(dates.test, y.test, main = "SkBART", pch = 19,
     ylab = "PM", xlab = NA, cex.lab = 2, cex.axis = 2, cex.main = 2)
lines(dates.test, cbart.preds, col = cols[4], lwd = 2, lty = 1)
lines(dates.test, cbart.quantiles[1, ], col = cols[4], lwd = 1, lty = 2)
lines(dates.test, cbart.quantiles[2, ], col = cols[4], lwd = 1, lty = 2)

# plot(dates.test, y.test, main = "SkRF", pch = 19,
#      ylab = "PM", xlab = "Date", cex.lab = 2, cex.axis = 2, cex.main = 2)
# lines(dates.test, crf.preds, col = cols[6], lwd = 2, lty = 1)
# lines(dates.test, crf.quantiles[1, ], col = cols[6], lwd = 1, lty = 2)
# lines(dates.test, crf.quantiles[2, ], col = cols[6], lwd = 1, lty = 2)
dev.off()


par(mfrow = c(1,1))
plot(dates.test, y.test)
for(ii in 1:20){
lines(dates.test, SkGP.fit.full$mu.pred[ii, ])
}


min.preds <- apply(SkGP.fit.full$mu.pred, MARGIN = 2, FUN = min)
max.preds <- apply(SkGP.fit.full$mu.pred, MARGIN = 2, FUN = max)
#pdf(file = "Plots/Real_Data_SkGP_Preds.pdf", width = 30, height = 10)

#dev.off()`

plot(dates.test, y.test, main = "Full Test Data", pch = 19)
# plot(dates.test, y.test, main = "Pixel Wise Predictive Averaged Data", pch = 19)
# legend("topleft", legend = c(paste0("BART: ", round(pmse.bart.avg, 3)),
#                              paste0("SkGp: ", round(pmse.avg, 3))),
#        col = c("darkblue", "maroon"), lty = c(1,1), lwd = c(2,2)
# )
# #lines(dates.test, bart.avg.preds, col = "darkblue", lwd = 2, lty = 1)
# lines(dates.test, res.avg$yhat.pred, col = "maroon", lwd = 2, lty = 1)
# lines(dates.test, res.avg$yhat.pred - 1.96*res.avg$CI.pred, col = "darkred", lwd = 1, lty = 2)
# lines(dates.test, res.avg$yhat.pred + 1.96*res.avg$CI.pred, col = "darkred", lwd = 1, lty = 2)
################################################################################
# Plots for Paper ##############################################################
################################################################################
setwd("/Users/samuelgailliot/Desktop/Research/SketchedGP/")
load("~/Desktop/Research/SketchedGP/Data/pm_images/image_data.rda")
# Load in the response readings
pmFRM <- read_csv("~/Desktop/Research/SketchedGP/Data/pm_images/pmFRM.csv")
dates <- ymd(pmFRM$date)

# Read in the data
# Remove zeros
BLUE.full <- vectorizor(imgs_a[,,,1]); which.nonzero1 <- which(colSums(BLUE.full)!=0); which.zero1 <- which(colSums(BLUE.full)==0)
BLUE <- BLUE.full[,which.nonzero1]
GREEN.full <- vectorizor(imgs_a[,,,2]); which.nonzero2 <- which(colSums(GREEN.full)!=0); which.zero2 <- which(colSums(GREEN.full)==0)
GREEN <- GREEN.full[,which.nonzero2]
RED.full <- vectorizor(imgs_a[,,,3]); which.nonzero3 <- which(colSums(RED.full)!=0); which.zero3 <- which(colSums(RED.full)==0)
RED <- RED.full[,which.nonzero3]
INFRA.full <- vectorizor(imgs_a[,,,4]); which.nonzero4 <- which(colSums(INFRA.full)!=0); which.zero4 <- which(colSums(INFRA.full)==0)
INFRA <- INFRA.full[,which.nonzero4]

row_stand <- function(X){
  means <- apply(X, 1, mean)
  sds <- apply(X, 1, sd)
  n <- nrow(X); p <- ncol(X)
  MEANS <- matrix(means, nrow = n, ncol = p, byrow = FALSE)
  SDS <- matrix(sds, nrow = n, ncol = p, byrow = FALSE)
  X <- (X - MEANS)/SDS
  return(X)
}

# Standardize images
BLUE  <- row_stand(BLUE)
GREEN <- row_stand(GREEN)
RED   <- row_stand(RED)
INFRA  <- row_stand(INFRA)

# Put the standardized data together
XX <- cbind(BLUE, GREEN, RED, INFRA)
yy <- pmFRM$pmFRM
dates <- ymd(pmFRM$date)
inds.tot <- (1:length(dates))[which(dates >= "2019-01-01")]
inds.test <- seq(3, length(inds.tot), 3)
inds.train <- 1:length(inds.tot)[-inds.test]
dates.tot <- dates[inds.tot]
dates.test <- dates[inds.test]
dates.train <- dates[inds.train]
XX <- XX[inds.tot, ]
yy.orig <- yy[inds.tot]
yy <- scale(log(yy.orig))
XX.train <- XX[inds.train, ]; yy.train <- yy[inds.train];
XX.test <- XX[inds.test, ]; yy.test <- yy[inds.test];

quant.thresh = 0.999
screening.scores <- get_NIS_scores(XX.train, yy.train)

#Permutation Test from NIS paper
#Doesnt do anything on the full data
perm <- sample(1:nrow(XX.train), size = nrow(XX.train), replace = FALSE)
screening.scores.perm <- get_NIS_scores(XX.train[perm, ], yy.train)

par(mfrow = c(1,1))
plot(1:ncol(XX.train),screening.scores)
points(1:ncol(XX.train), screening.scores, col = ifelse(screening.scores < quantile(screening.scores.perm, quant.thresh), "black", "red"))
abline(h = quantile(screening.scores.perm, quant.thresh), col = "red", lty = 2, lwd = 3)
points(1:ncol(XX.train), screening.scores.perm, col = "blue")

covars.save <- which(screening.scores >= quantile(screening.scores.perm, quant.thresh))
keepers <- 1:ncol(XX) %in% covars.save

keepers.blue <- keepers[1:8267]
keepers.green <- keepers[(8267+1):(2*8267)]
keepers.red <- keepers[(2*8267+1):(3*8267)]
keepers.infra <- keepers[(3*8267+1):(4*8267)]

BLUE.kept <- BLUE.full
BLUE.kept[,which.nonzero1][,!keepers.blue] <- NA

GREEN.kept <- GREEN.full
GREEN.kept[,which.nonzero2][,!keepers.green] <- NA

RED.kept <- RED.full
RED.kept[,which.nonzero3][,!keepers.red] <- NA

INFRA.kept <- INFRA.full
INFRA.kept[,which.nonzero4][,!keepers.infra] <- NA

WIDTH = 10
HEIGHT = 10

day_ind <- 809
par(mfrow = c(4,2))
image.plot(t(matrix(BLUE.full[day_ind,], 128, 128)),
           xaxt = 'n', yaxt = 'n', main = "Blue")
image.plot(t(matrix(BLUE.kept[day_ind,], 128, 128)))

image.plot(t(matrix(GREEN.full[day_ind,], 128, 128)))
image.plot(t(matrix(GREEN.kept[day_ind,], 128, 128)))

image.plot(t(matrix(RED.full[day_ind,], 128, 128)))
image.plot(t(matrix(RED.kept[day_ind,], 128, 128)))

image.plot(t(matrix(INFRA.full[day_ind,], 128, 128)))
image.plot(t(matrix(INFRA.kept[day_ind,], 128, 128)))

day_ind <- 1
par(mfrow = c(1,1))
pdf(file = "Plots/BLUE_Full.pdf", width = WIDTH, height = HEIGHT)
image.plot(t(matrix(BLUE.full[day_ind,], 128, 128)),
           xaxt = 'n', yaxt = 'n', main = "Blue",cex.main = 3)
dev.off()
pdf(file = "Plots/GREEN_Full.pdf", width = WIDTH, height = HEIGHT)
image.plot(t(matrix(GREEN.full[day_ind,], 128, 128)),
           xaxt = 'n', yaxt = 'n', main = "Green",cex.main = 3)
dev.off()
pdf(file = "Plots/RED_Full.pdf", width = WIDTH, height = HEIGHT)
image.plot(t(matrix(RED.full[day_ind,], 128, 128)),
           xaxt = 'n', yaxt = 'n', main = "Red",cex.main = 3)
dev.off()
pdf(file = "Plots/INFRA_Full.pdf", width = WIDTH, height = HEIGHT)
image.plot(t(matrix(INFRA.full[day_ind,], 128, 128)),
           xaxt = 'n', yaxt = 'n', main = "Near Infrared",cex.main = 3)
dev.off()

day_ind <- 809
par(mfrow = c(1,1))
pdf(file = "Plots/BLUE_screen.pdf", width = WIDTH, height = HEIGHT)
image.plot(t(matrix(BLUE.kept[day_ind,], 128, 128)),
           xaxt = 'n', yaxt = 'n', main = "Blue",cex.main = 3)
dev.off()
pdf(file = "Plots/GREEN_screen.pdf", width = WIDTH, height = HEIGHT)
image.plot(t(matrix(GREEN.kept[day_ind,], 128, 128)),
           xaxt = 'n', yaxt = 'n', main = "Green",cex.main = 3)
dev.off()
pdf(file = "Plots/RED_screen.pdf", width = WIDTH, height = HEIGHT)
image.plot(t(matrix(RED.kept[day_ind,], 128, 128)),
           xaxt = 'n', yaxt = 'n', main = "Red",cex.main = 3)
dev.off()
pdf(file = "Plots/INFRA_screen.pdf", width = WIDTH, height = HEIGHT)
image.plot(t(matrix(INFRA.kept[day_ind,], 128, 128)),
           xaxt = 'n', yaxt = 'n', main = "Near Infrared",cex.main = 3)
dev.off()

pdf(file = "Plots/pmfrm.pdf", width = 30, height = 10)
par(mfrow = c(2, 1))
par(mar=c(2, 6, 1.5, 0.5))
plot(dates.tot, pmFRM$pmFRM[inds.tot], "b", xaxt = "n",
     ylab = "Original Scale", xlab = NA,
     cex.axis = 2, cex.lab = 2)
plot(dates.tot, scale(log(pmFRM$pmFRM[inds.tot])), "b",
     ylab = "Logged-Centered Scale", xlab = "Date",
     cex.axis = 2, cex.lab = 2)
dev.off()

