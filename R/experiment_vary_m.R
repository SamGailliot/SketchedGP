# Changing m experiment
library(tidyverse)
library(fields)
library(tseries)
library(foreach)
library(doParallel)
library(plgp)
library(splines)
# changing p_screen experiment
source("R/functions.R")
ROLLS003 <- readRDS("~/Desktop/Research/SketchedGP/Data/ROLLS001.rds")

n.screen = 2000
m.vals <- c(5, 10,20,30,40,50,100,seq(0, 5000, 500)[-1])
K = 30; n.theta = 4; n.snrs = 3

RESULTS <- list()

for(jj in 1:40){
  print(paste0("############ Roll ", jj, " ############"))
  ROLL <- ROLLS003[[jj]]
  y <- ROLL$y[1:100]
  X <- ROLL$X[1:100, ]

  ystar <- ROLL$y[301:400]
  Xstar <- ROLL$X[301:400, ]

  X <- scale(X, center = TRUE, scale = FALSE)
  y <- scale(y, center = TRUE, scale = FALSE)
  Xstar <- scale(Xstar, center = TRUE, scale = FALSE)
  ystar <- scale(ystar, center = TRUE, scale = FALSE)

  # Get screening scores
  screening.scores <- get_NIS_scores(X, y)

  PREDS <- matrix(data = NA, nrow = length(m.vals), ncol = length(ystar))
  SIGS <- matrix(data = NA, nrow = length(m.vals), ncol = length(ystar))
  pmses <- numeric(length = length(m.vals))

  for(ii in 1:length(m.vals)){
    print(paste0("############ Run ", ii, "; m = ", m.vals[ii], " ############"))
    # Save the first n.keep covariates with the highest scores
    n.keep = n.screen
    covars.save <- order(screening.scores, decreasing = TRUE)[1:n.keep]

    # Get the screened X's
    X.s <- X[,covars.save]
    Xstar.s <- Xstar[,covars.save]

    out <- sketched_GP(y, X.s, ystar, Xstar.s, m = m.vals[ii], K = 10,
                        n.theta = 10, n.snrs = 25, prediction = TRUE,
                        snr.method = "sample", snr.max = 10,
                        stacking.method = "LOO",n.folds = 10)


    stack.weights <- out$stack.weights
    mu.pred <- t(out$mu.pred)
    sig.pred <- t(out$sig.pred)
    # mu.fit <- t(out$mu.fit)
    # sig.fit <- t(out$sig.fit)

    PREDS[ii, ] <- rowSums(mu.pred*matrix(data = unname(stack.weights), nrow = nrow(mu.pred), ncol = ncol(mu.pred), byrow = TRUE))
    SIGS[ii, ] <-  sqrt(rowSums(sig.pred*matrix(data = unname(stack.weights), nrow = nrow(sig.pred), ncol = ncol(sig.pred), byrow = TRUE)))
    pmses[ii] <- sum((PREDS[ii, ] - ystar)^2)/length(ystar)

    # stack.mod.fit <- rowSums(mu.fit*matrix(data = unname(stack.weights), nrow = nrow(mu.fit), ncol = ncol(mu.fit), byrow = TRUE))
    # stack.sd.fit <-  sqrt(rowSums(sig.fit*matrix(data = unname(stack.weights), nrow = nrow(sig.fit), ncol = ncol(sig.fit), byrow = TRUE)))
  }
  RESULTS[[jj]] <- list("PREDS" = PREDS, "SIGS" = SIGS, "pmses" = pmses)
}

pmse.pred <- sum((PREDS[ii,] - ystar)^2)/length(ystar)
ord = order(Xstar[,2])
plot(Xstar[ord,2], ystar[ord], main = paste0("Fit to Test Data: PMSE =  ", round(pmse.pred, 3)))
lines(Xstar[ord,2], PREDS[ii,ord], col = "Maroon", lwd = 2)
lines(Xstar[ord,2], PREDS[ii,ord] + SIGS[ii,ord], col = "lightblue", lwd = 1)
lines(Xstar[ord,2], PREDS[ii,ord] - SIGS[ii,ord], col = "lightblue", lwd = 1)
#points(Xstar[ord,2], ystar[ord])
out$Models

colfunc <- colorRampPalette(c("red", "blue"))
plot.cols <-colfunc(length(m.vals))
#
plot(Xstar[ord, 2], ystar[ord])
for(ii in 1:length(m.vals)){
  lines(Xstar[ord, 2], PREDS[ii, ord], col = plot.cols[ii])
}

PMSES <- matrix(nrow = 40, ncol = length(m.vals))
for(ii in 1:40){
  PMSES[ii, ] <- RESULTS[[ii]]$pmses
}

boxplot(PMSES)
col1 <- "lightskyblue"
col2 <- "azure4"
boxplot(PMSES, col = col1, xaxt = "n", ylab = "PMSE", xlab = "Sketching Dimension", cex.label = 2, cex.axis = 1.5)
axis(1, at = seq(1,length(m.vals),1), labels = m.vals)

