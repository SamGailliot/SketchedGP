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

n.screen <- c(250, 500, seq(1000, 10000, 1000))
K = 10; n.theta = 25; n.snrs = 25

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

  PREDS <- matrix(data = NA, nrow = length(n.screen), ncol = length(ystar))
  SIGS <- matrix(data = NA, nrow = length(n.screen), ncol = length(ystar))
  pmses <- numeric(length = length(n.screen))

for(ii in 1:length(n.screen)){
  print(paste0("############ Run ", ii, "; nscreen = ", n.screen[ii], " ############"))
  # Save the first n.keep covariates with the highest scores
  n.keep = n.screen[ii]
  covars.save <- order(screening.scores, decreasing = TRUE)[1:n.keep]

  # Get the screened X's
  X.s <- X[,covars.save]
  Xstar.s <- Xstar[,covars.save]

  out <- sketched_GP(y, X.s, ystar, Xstar.s, m = 60, K = K,
                   SNRs = c(10000), n.theta = n.theta, n.snrs = n.snrs,prediction = TRUE,
                   snr.method = "set", snr.max = 8750,
                   stacking.method = "kfold", n.folds = 10)


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

# par(mfrow = c(1,2))
# ord = order(X[,2])
# pmse.fit <- sum((stack.mod.fit - y)^2)/length(y)
# plot(X[ord,2], y[ord], main = paste0("Fit to Training Data: PMSE =  ", round(pmse.fit, 3)))
# lines(X[ord,2], stack.mod.fit[ord], col = "Maroon", lwd = 2)
# lines(X[ord,2], stack.mod.fit[ord] + stack.sd.fit[ord], col = "lightblue", lwd = 1)
# lines(X[ord,2], stack.mod.fit[ord] - stack.sd.fit[ord], col = "lightblue", lwd = 1)
#points(X[ord,2], y[ord])

#saveRDS(RESULTS, file = "Results/change_pscreen_ROLL01.rds")

pmse.pred <- sum((PREDS[ii,] - ystar)^2)/length(ystar)
ord = order(Xstar[,2])
plot(Xstar[ord,2], ystar[ord], main = paste0("Fit to Test Data: PMSE =  ", round(pmse.pred, 3)))
lines(Xstar[ord,2], PREDS[ii,ord], col = "Maroon", lwd = 2)
lines(Xstar[ord,2], PREDS[ii,ord] + SIGS[ii,ord], col = "lightblue", lwd = 1)
lines(Xstar[ord,2], PREDS[ii,ord] - SIGS[ii,ord], col = "lightblue", lwd = 1)
#points(Xstar[ord,2], ystar[ord])
out$Models

colfunc <- colorRampPalette(c("red", "blue"))
plot.cols <-colfunc(length(n.screen))
#
plot(Xstar[ord, 2], ystar[ord])
for(ii in 1:length(n.screen)){
  lines(Xstar[ord, 2], PREDS[ii, ord], col = plot.cols[ii])
}

PMSES <- matrix(nrow = 40, ncol = 12)
for(ii in 1:40){
  PMSES[ii, ] <- RESULTS[[ii]]$pmses
}

boxplot(PMSES)
col1 <- "lightskyblue"
col2 <- "azure4"
boxplot(PMSES, col = col1, xaxt = "n", ylab = "PMSE", xlab = "Sketching Dimension")
axis(1, at = seq(1,length(n.screen),1), labels = n.screen)

################################################################################
# Plot All PMSEs ###############################################################
################################################################################
PMSES001 <- matrix(nrow = 40, ncol = 12)
PMSES003 <- matrix(nrow = 40, ncol = 12)
PMSES005 <- matrix(nrow = 40, ncol = 12)
PMSES01 <- matrix(nrow = 40, ncol = 12)

for(ii in 1:40){
  PMSES001[ii, ] <- change_pscreen_ROLL001[[ii]]$pmses
  PMSES003[ii, ] <- change_pscreen_ROLL003[[ii]]$pmses
  PMSES005[ii, ] <- change_pscreen_ROLL005[[ii]]$pmses
  PMSES01[ii, ] <- change_pscreen_ROLL01[[ii]]$pmses
}

n.screen <- c(250,500,seq(1000, 10000, 1000))/1000
n.screen01 <- c(250,500,seq(1000, 10000, 1000))/1000

labels <- c("0.25", "0.5", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10E3")
PMSES.FULL <- cbind(PMSES001, PMSES003, PMSES005, PMSES01)
n.screen.full <- c(n.screen, n.screen, n.screen, n.screen01)

boxplot(PMSES.FULL)
col1 <- "lightskyblue"
col2 <- "azure4"
boxplot(PMSES.FULL, col = col1, xaxt = "n", ylab = "PMSE", xlab = "Sketching Dimension", )
axis(1, at = seq(1,length(n.screen.full),1), labels = n.screen.full)
abline(v = 10.5)

options(scipen=1)
par(mfrow = c(1,4))
col1 <- "lightskyblue"
col2 <- "azure4"
boxplot(PMSES001, col = col1, xaxt = "n", ylab = "PMSE", xlab = "Number of Features", ylim = c(0.5, 9), main = "Noise = 0.01",
        cex.label = 2, cex.axis = 2)
axis(1, at = seq(1,length(n.screen),1), labels = labels)

col1 <- "lightskyblue"
col2 <- "azure4"
boxplot(PMSES003, col = col1, xaxt = "n", yaxt = "n", xlab = "Number of Features", ylim = c(0.5, 9), main = "0.03")
axis(1, at = seq(1,length(n.screen),1), labels = labels)

col1 <- "lightskyblue"
col2 <- "azure4"
boxplot(PMSES005, col = col1, xaxt = "n", yaxt = "n", ylab = NA, xlab = "Number of Features", ylim = c(0.5, 9), main = "0.05")
axis(1, at = seq(1,length(n.screen),1), labels = labels)

col1 <- "lightskyblue"
col2 <- "azure4"
boxplot(PMSES01, col = col1, xaxt = "n", yaxt = "n", ylab = NA, xlab = "Number of Features", ylim = c(0.5, 9), main = "0.1")
axis(1, at = seq(1,length(n.screen01),1), labels = labels)


