source("R/functions.R")
library(foreach)
library(tidyverse)
library(fields)
library(tseries)
library(foreach)
library(doParallel)
library(plgp)
library(splines)
# Figures for Slides
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/TORI001.rds")
d1 <- DATA[[1]]
XX <- scale(d1$X, center = TRUE, scale = FALSE)
yy <- scale(d1$y, center = TRUE, scale = FALSE)

inds.train <- 1:400
inds.test <- 401:600

p.keep <- 500

X <- XX[inds.train, 1:p.keep]; y <- yy[inds.train]
Xstar <- XX[inds.test, 1:p.keep]; ystar <- yy[inds.test]

# Figure: Regular GP regression on full swiss roll
dist.mat = plgp::distance(X)
dmax <- max(dist.mat)
dmin <- min( dist.mat[dist.mat != 0] )

theta <- get_theta(y, X, dmin, dmax, 10, 1)
#theta =
K1 <- exp_kernel(X, theta = theta)
K.1.pred <- exp_kernel(X, Xstar, theta = theta); K.pred.1 <- t(K.1.pred)
K.pred <- exp_kernel(Xstar, theta = theta)

SNR = 1; n <- nrow(X); nstar <- nrow(Xstar)
I <- diag(rep(1, n)); Istar <- diag(rep(1, nstar))
G <- SNR*K1 + I
G.chol <- chol(G)
G.inv <- chol2inv(G.chol)

gp.preds = as.numeric(SNR*K.pred.1 %*% G.inv %*% y)
b1 <-  (t(y) %*% G.inv %*% y) / 2
gp.sigs <- diag(as.numeric((2*b1/n)) * (Istar + SNR*K.pred - (SNR^2)*K.pred.1 %*% G.inv %*% K.1.pred))

gp.preds

par(mfrow = c(1,1))
pmse0 <- mean((ystar - gp.preds)^2)
ord = order(Xstar[,2])
plot(Xstar[ord,2], ystar[ord], ylab = "Response", xlab = "X2",
     main = paste0("Predictions with Regular GP: PMSE = ", round(pmse0, 3)))
lines(Xstar[ord,2], gp.preds[ord], col = "maroon", lwd = 2)

SkGP <- sketched_GP(y, X, ystar, Xstar, m = 60, K = 30,
                    n.theta = 100, n.snrs = 10, SNRs = c(8750), prediction = TRUE,
                    snr.method = "sample", snr.max = 8750)

stack.weights <- SkGP$stack.weights
mu.pred <- t(SkGP$mu.pred)
sig.pred <- t(SkGP$sig.pred)

skgp.preds <- rowSums(mu.pred*matrix(data = unname(stack.weights), nrow = nrow(mu.pred), ncol = ncol(mu.pred), byrow = TRUE))
skgp.sigs <- sqrt(rowSums(sig.pred*matrix(data = unname(stack.weights), nrow = nrow(sig.pred), ncol = ncol(sig.pred), byrow = TRUE)))

ord = order(Xstar[,2])
plot(Xstar[ord,2], ystar[ord])
lines(Xstar[ord,2], gp.preds[ord], col = "red")

lines(Xstar[ord,2], skgp.preds[ord], col = "blue")

colfunc <- colorRampPalette(c("red", "blue"))
plot.cols <-colfunc(5)

plot(SkGP$stack.weights)
par(mfrow = c(1,2))

pmse1 <- mean((ystar - mu.pred[,1])^2)
plot(Xstar[ord, 2], ystar[ord],
     ylab = "Response", xlab = "X2",
     main = paste0("Model 1: PMSE = ", round(pmse1, 3)))
lines(Xstar[ord, 2], mu.pred[ord, 1], col = plot.cols[1], lwd = 2)

pmse2 <- mean((ystar - mu.pred[,2])^2)
plot(Xstar[ord, 2], ystar[ord],
     ylab = NA, xlab = "X2",
     main = paste0("Model 2: PMSE = ", round(pmse2, 3)))
lines(Xstar[ord, 2], mu.pred[ord, 2], col = plot.cols[5], lwd = 2)


par(mfrow = c(1,1))
a <- round(stack.weights[1], 3)
b <- round(stack.weights[2], 3)
pmse3 <- mean((ystar - skgp.preds)^2)
plot(Xstar[ord, 2], ystar[ord], xlab = "X2", ylab = "Response",
     main = paste0("Stacked Model K = ", 30, " pscreen = 500",": PMSE = ", round(pmse3, 3)))
lines(Xstar[ord, 2], skgp.preds[ord], col = plot.cols[3], lwd = 3)
lines(Xstar[ord, 2], mu.pred[ord, 1], col = plot.cols[1], lwd = 1, lty = 2)
lines(Xstar[ord, 2], mu.pred[ord, 2], col = plot.cols[5], lwd = 2, lty = 2)
legend("topleft",
       legend = c("Stacked Model",
         bquote(paste( 'Model 1: w'[1]*' = ',.(a) )),
         bquote(paste( 'Model 2: w'[2]*' = ',.(b) ))),
       col = c(plot.cols[3], plot.cols[1], plot.cols[5]),
       lty = c(1,2,2), lwd = c(3,1,1)
)


mu.fit <- t(SkGP$mu.fit)
sig.fit <- t(SkGP$sig.fit)

skgp.preds.fit <- rowSums(mu.fit*matrix(data = unname(stack.weights), nrow = nrow(mu.fit), ncol = ncol(mu.fit), byrow = TRUE))
skgp.sigs.fit <- sqrt(rowSums(sig.fit*matrix(data = unname(stack.weights), nrow = nrow(sig.fit), ncol = ncol(sig.fit), byrow = TRUE)))

ord <- order(X[,2])
par(mfrow = c(1,1))
plot(X[ord, 2], y[ord])
lines(X[ord, 2], skgp.preds.fit[ord], col = plot.cols[5], lwd = 2)



lines(Xstar[ord, 2], mu.pred[ord, 2])
lines(Xstar[ord, 2], mu.pred[ord, 3])
lines(Xstar[ord, 2], mu.pred[ord, 4])
lines(Xstar[ord, 2], mu.pred[ord, 5])
lines(Xstar[ord, 2], skgp.preds[ord], col = "red")

dim(mu.preds)
