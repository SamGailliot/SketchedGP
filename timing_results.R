# Timing Results
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

p.vec <- seq(0, 10000, 1000)[-1]
n.reps <- 2

#SkGP.loo.times   <- matrix(data = NA, nrow = 3, ncol = n.reps)
SkGP.10fold.times <- matrix(data = NA, nrow = length(p.vec), ncol = n.reps)
#SkGP.20fold.times <- matrix(data = NA, nrow = length(p.vec), ncol = n.reps)
#SkGP.5fold.times <- matrix(data = NA, nrow = length(p.vec), ncol = n.reps)
GP.times          <- matrix(data = NA, nrow = length(p.vec), ncol = n.reps)
BART.times        <- matrix(data = NA, nrow = length(p.vec), ncol = n.reps)
CBART.times       <- matrix(data = NA, nrow = length(p.vec), ncol = n.reps)
RF.times          <- matrix(data = NA, nrow = length(p.vec), ncol = n.reps)
CRF.times         <- matrix(data = NA, nrow = length(p.vec), ncol = n.reps)
LASSO.times       <- matrix(data = NA, nrow = length(p.vec), ncol = n.reps)


for(ii in 1:length(p.vec)){
  print(paste0("#################### n = ", p.vec[ii], " ####################"))
  for(jj in 1:n.reps){
    print(paste0("#################### ", jj, " ####################"))

    n.roll <- 200
    ROLL <- get_swiss_roll(n.roll, p.vec[ii], 0.01)
    y <- ROLL$y[1:100]
    Xstar.s <- ROLL$X[101:200, ]
    ystar <- ROLL$y[101:200]

    # Pull out data of correct size
    X.s <- ROLL$X[1:100, ]
    y <- ROLL$y[1:100]

    M <- ceiling(sqrt(p.vec[ii]))
    P1 <- get_sketch_mat(m = M, p = ncol(Xstar.s))
    Sketched.X.s <- X.s%*%t(P1)
    Sketched.Xstar.s <- Xstar.s%*% t(P1)

    print("SkGP 10 Fold")
    # SkGP timing
    # 10 FOLD
    t1 <- Sys.time()
    out <- sketched_GP(y, X.s, ystar, Xstar.s, m = M, K = 2,
                   SNRs = c(10), n.theta = 2, n.snrs = n.snrs,prediction = FALSE,
                   snr.method = "set", snr.max = 8750,
                   stacking.method = "kfold", n.folds = 10)
    t2 <- Sys.time()
    SkGP.10fold.times[ii, jj] <- difftime(t2, t1, units = "secs")

    # print("SkGP 20 Fold")
    # # SkGP timing
    # # 10 FOLD
    # t1 <- Sys.time()
    # out <- sketched_GP(y, X.s, ystar, Xstar.s, m = 60, K = 10,
    #                    SNRs = c(10), n.theta = 10, n.snrs = n.snrs,prediction = FALSE,
    #                    snr.method = "set", snr.max = 8750,
    #                    stacking.method = "kfold", n.folds = 20)
    # t2 <- Sys.time()
    # SkGP.20fold.times[ii, jj] <- difftime(t2, t1, units = "secs")

    # print("SkGP 5 Fold")
    # # SkGP timing
    # # 10 FOLD
    # t1 <- Sys.time()
    # out <- sketched_GP(y, X.s, ystar, Xstar.s, m = 60, K = 10,
    #                    SNRs = c(10), n.theta = 10, n.snrs = n.snrs,prediction = FALSE,
    #                    snr.method = "set", snr.max = 8750,
    #                    stacking.method = "kfold", n.folds = 5)
    # t2 <- Sys.time()
    # SkGP.5fold.times[ii, jj] <- difftime(t2, t1, units = "secs")

    # if(ii <= 3){
    # print("SkGP LOO")
    # # LOO
    # t1 <- Sys.time()
    # out <- sketched_GP(y, X.s, ystar, Xstar.s, m = 60, K = 10,
    #                    SNRs = c(10), n.theta = 10, n.snrs = n.snrs,prediction = FALSE,
    #                    snr.method = "set", snr.max = 8750,
    #                    stacking.method = "LOO", n.folds = 10)
    # t2 <- Sys.time()
    # SkGP.loo.times[ii, jj] <- difftime(t2, t1, units = "secs")
    # }

    print("BART")
    # Bart timing
    t1 <- Sys.time()
    #BART <- mc.wbart(X.s, y, Xstar.s, ndpost = 1000, nskip = 1000, mc.cores = 6, seed = 113)
    BART <- wbart(X.s, y, ndpost = 500, nskip = 100)
    t2 <- Sys.time()
    BART.times[ii, jj] <- difftime(t2, t1, units = "secs")

    print("CBART")
    # CBart timing
    t1 <- Sys.time()
    #CBART <- mc.wbart(Sketched.X.s, y, Sketched.Xstar.s, ndpost = 1000, nskip = 500, mc.cores = 6, seed = 113)
    CBART <- wbart(Sketched.X.s, y, ndpost = 500, nskip = 100)
    t2 <- Sys.time()
    CBART.times[ii, jj] <- difftime(t2, t1, units = "secs")

    print("RF")
    # RF timing
    t1 <- Sys.time()
    RF <- randomForest(x=X.s, y=y, ytest = ystar, ntree = 500)
    #rf.preds.obj <- predict(RF, Xstar.s, predict.all = TRUE)
    t2 <- Sys.time()
    RF.times[ii, jj] <- difftime(t2, t1, units = "secs")

    print("CRF")
    # CRF timing
    t1 <- Sys.time()
    CRF <- randomForest(x=Sketched.X.s, y=y, ntree = 500)
    #crf.preds.obj <- predict(CRF, Sketched.Xstar.s, predict.all = TRUE)
    t2 <- Sys.time()
    CRF.times[ii, jj] <- difftime(t2, t1, units = "secs")

    print("LASSO")
    # LASSO timing
    t1 <- Sys.time()
    LASSO <- cv.glmnet(x=X.s, y=y, alpha = 1)
    l.lasso.min <- LASSO$lambda.min
    lasso.model <- glmnet(x=X.s, y=y,
                          alpha  = 1,
                          lambda = l.lasso.min)

    lasso.preds <- predict(lasso.model, Xstar.s )
    t2 <- Sys.time()
    LASSO.times[ii, jj] <- difftime(t2, t1, units = "secs")

    # print("GP")
    # # GP timing
    # t1 <- Sys.time()
    # dist.mat = plgp::distance(X.s)
    # dmax <- max(dist.mat)
    # dmin <- min( dist.mat[dist.mat != 0] )
    #
    # theta <- get_theta(y, X.s, dmin, dmax, 100, 1)
    # K1 <- exp_kernel(X.s, theta = theta)
    # K.1.pred <- exp_kernel(X.s, Xstar.s, theta = theta); K.pred.1 <- t(K.1.pred)
    # K.pred <- exp_kernel(Xstar.s, theta = theta)
    #
    # SNR = 1; n <- nrow(X.s); nstar <- nrow(Xstar.s)
    # I <- diag(rep(1, n)); Istar <- diag(rep(1, nstar))
    # G <- SNR*K1 + I
    # G.chol <- chol(G)
    # G.inv <- chol2inv(G.chol)
    #
    # gp.preds = as.numeric(SNR*K.pred.1 %*% G.inv %*% y)
    # b1 <-  (t(y) %*% G.inv %*% y) / 2
    # gp.sigs <- diag(as.numeric((2*b1/n)) * (Istar + SNR*K.pred - (SNR^2)*K.pred.1 %*% G.inv %*% K.1.pred))
    # t2 <- Sys.time()
    # GP.times[ii, jj] <- difftime(t2, t1, units = "secs")
  }
}



LASSO.times
SkGP.10fold.times
#SkGP.loo.times
BART.times
CBART.times
RF.times
CRF.times
#GP.times


#lasso.means <- apply(LASSO.times, 1, mean)
#skgp.loo.means <- c(apply(SkGP.loo.times, 1, mean), NA, NA)
#skgp.5fold.means <- apply(SkGP.5fold.times, 1, mean)
skgp.10fold.means <- apply(SkGP.10fold.times, 1, mean)
#skgp.20fold.means <- apply(SkGP.20fold.times, 1, mean)
bart.means <- apply(BART.times, 1, mean)
cbart.means <- apply(CBART.times, 1, mean)
rf.means <- apply(RF.times, 1, mean)
crf.means <- apply(CRF.times, 1, mean)
#gp.means <- apply(GP.times, 1, mean)

# plot(1:5, skgp.loo.means, )
library(RColorBrewer)
yrange = range(c(skgp.10fold.means, bart.means, cbart.means, rf.means, crf.means))
par(mfrow = c(1,1))
xx.grid <- 1:length(p.vec)
cols <- brewer.pal(8, "Dark2")
plot(xx.grid, skgp.10fold.means, "b", col = cols[1], lwd = 3, pch = 1,
     ylab = "Time (sec)", xaxt = "n", xlab = "Number of Features (p)", cex.axis = 2, cex.lab = 2,
     ylim = yrange)
axis(1, at = xx.grid, labels = p.vec, cex.axis = 1.5)
legend = legend("topleft", legend = c("SkGP", "BART", "SkBART", "RF", "SkRF"),
                pch = c(19,3,4,5,6),
                col = c(cols[1], cols[3], cols[4], cols[5], cols[6]),
                cex = 1.5)
lines(xx.grid, bart.means, "b", col = cols[3], lwd = 3, pch = 3)
lines(xx.grid, cbart.means, "b", col = cols[4], lwd = 3, pch = 4)
lines(xx.grid, rf.means, "b", col = cols[5], lwd = 3, pch = 5)
lines(xx.grid, crf.means, "b", col = cols[6], lwd = 3, pch = 6)
lines(xx.grid, skgp.10fold.means, "b", col = cols[1], lwd = 3, pch = 19)


bart.means


