# Gonna need o update this for regular gp
get_pmse_coverage <- function(Sims, DATA){
  #browser()
  NN <- length(DATA)
  inds.test <- Sims$inds.test
  pmse.skgp <- numeric(length = NN)
  pmse.gp <- numeric(length = NN)
  pmse.bart <- numeric(length = NN)
  pmse.cbart <- numeric(length = NN)
  pmse.lasso <- numeric(length = NN)
  pmse.rf <- numeric(length = NN)
  pmse.crf <- numeric(length = NN)

  coverage.skgp <- numeric(length = NN)
  coverage.gp <- numeric(length = NN)
  coverage.bart <- numeric(length = NN)
  coverage.cbart <- numeric(length = NN)
  coverage.rf <- numeric(length = NN)
  coverage.crf <- numeric(length = NN)

  length.skgp <- numeric(length = NN)
  length.gp <- numeric(length = NN)
  length.bart <- numeric(length = NN)
  length.cbart <- numeric(length = NN)
  length.rf <- numeric(length = NN)
  length.crf <- numeric(length = NN)

  for(ii in 1:NN){
    y <- scale(DATA[[ii]]$y[inds.test], center = TRUE, scale = FALSE)

    pmse.skgp[ii] <- mean((Sims$PREDS.skgp[ii, ] - y)^2)
    pmse.gp[ii] <- mean((Sims$PREDS.gp[ii, ] - y)^2)
    pmse.bart[ii] <- mean((Sims$PREDS.bart[ii, ] - y)^2)
    pmse.cbart[ii] <- mean((Sims$PREDS.cbart[ii, ] - y)^2)
    pmse.lasso[ii] <- mean((Sims$PREDS.lasso[ii, ] - y)^2)
    pmse.rf[ii] <- mean((Sims$PREDS.rf[ii, ] - y)^2)
    pmse.crf[ii] <- mean((Sims$PREDS.crf[ii, ] - y)^2)

    coverage.skgp[ii] <- sum((y > (Sims$PREDS.skgp[ii, ] - 1.96*sqrt(Sims$SIGS.skgp[ii, ]))) & (y < (Sims$PREDS.skgp[ii, ] + 1.96*sqrt(Sims$SIGS.skgp[ii, ])))) / length(y)
    coverage.gp[ii] <- sum((y > (Sims$PREDS.gp[ii, ] - 1.96*sqrt(Sims$SIGS.gp[ii, ]))) & (y < (Sims$PREDS.gp[ii, ] + 1.96*sqrt(Sims$SIGS.gp[ii, ])))) / length(y)
    coverage.bart[ii] <- sum((y > Sims$LOWER.bart[ii]) & (y < Sims$UPPER.bart[ii])) / length(y)
    coverage.cbart[ii] <- sum((y > Sims$LOWER.cbart[ii]) & (y < Sims$UPPER.cbart[ii])) / length(y)
    coverage.rf[ii] <- sum((y > Sims$LOWER.rf[ii]) & (y < Sims$UPPER.rf[ii])) / length(y)
    coverage.crf[ii] <- sum((y > Sims$LOWER.crf[ii]) & (y < Sims$UPPER.crf[ii])) / length(y)

    length.skgp[ii] <- mean(2*1.96*sqrt(Sims$SIGS.skgp[ii, ]))
    length.gp[ii] <- mean(2*1.96*sqrt(Sims$SIGS.gp[ii, ]))
    length.bart[ii] <- mean(Sims$UPPER.bart[ii] - Sims$LOWER.bart[ii])
    length.cbart[ii] <- mean(Sims$UPPER.cbart[ii] - Sims$LOWER.cbart[ii])
    length.rf[ii] <- mean(Sims$UPPER.rf[ii] - Sims$LOWER.rf[ii])
    length.crf[ii] <- mean(Sims$UPPER.crf[ii] - Sims$LOWER.crf[ii])
  }
  return(list("pmse.skgp"=pmse.skgp, "pmse.gp"=pmse.gp, "pmse.bart"=pmse.bart, "pmse.cbart"=pmse.cbart, "pmse.lasso"=pmse.lasso,
              "pmse.rf"=pmse.rf, "pmse.crf"=pmse.crf,"coverage.skgp"=coverage.skgp,"coverage.gp"=coverage.gp,
              "coverage.bart"=coverage.bart, "coverage.cbart"=coverage.cbart, "coverage.rf"=coverage.rf,
              "coverage.crf"=coverage.crf,
              "length.skgp"=length.skgp,"length.gp"=length.gp,
              "length.bart"=length.bart, "length.cbart"=length.cbart, "length.rf"=length.rf,
              "length.crf"=length.crf))
}
################################################################################
# p = 2000 #####################################################################
################################################################################

# Read in results
Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_ROLLS001_p2k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/ROLLS001.rds")
out.R001 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_ROLLS003_p2k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/ROLLS003.rds")
out.R003 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_ROLLS005_p2k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/ROLLS005.rds")
out.R005 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_ROLLS01_p2k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/ROLLS01.rds")
out.R01 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_TORI001_p2k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/TORI001.rds")
out.T001 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_TORI005_p2k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/TORI005.rds")
out.T005 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_TORI003_p2k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/TORI003.rds")
out.T003 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_TORI01_p2k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/TORI01.rds")
out.T01 <- get_pmse_coverage(Sims, DATA)

# Write out and plot results

l.skgp  <- matrix(NA, nrow = 4, ncol = 4)
l.gp    <- matrix(NA, nrow = 4, ncol = 4)
l.bart  <- matrix(NA, nrow = 4, ncol = 4)
l.cbart <- matrix(NA, nrow = 4, ncol = 4)
l.rf    <- matrix(NA, nrow = 4, ncol = 4)
l.crf   <- matrix(NA, nrow = 4, ncol = 4)
###### Rolls 2K ######

WIDTH = 10
HEIGHT = 10
l.max = 10

##### R001 #####
OUT <- out.R001
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.R001.2 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)

{
pdf(file = "Plots/pmses_R001_2k.pdf", width = WIDTH, height = HEIGHT)
boxplot(pmses,
        names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
        col = "lightskyblue", cex.axis = 1.5)
dev.off()

pdf(file = "Plots/coverages_R001_2k.pdf", width = WIDTH, height = HEIGHT)
boxplot(cvgs,
        names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
        col = "lightskyblue", cex.axis = 1.5)
dev.off()

# pdf(file = "Plots/lengths_R001_2k.pdf", width = WIDTH, height = HEIGHT)
# barplot(colMeans(lengths),
#         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,12),
#         col = "lightskyblue")
# dev.off()
}

##### R003 #####
OUT <- out.R003
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.R003.2 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_R003_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_R003_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_R003_2k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

##### R005 #####
OUT <- out.R005
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.R005.2 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)

{
  pdf(file = "Plots/pmses_R005_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_R005_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_R005_2k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

##### R01 #####
OUT <- out.R01
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.R01.2 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_R01_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_R01_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_R01_2k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

l001 <- colMeans(lengths.R001.2)
l003 <- colMeans(lengths.R003.2)
l005 <- colMeans(lengths.R005.2)
l01  <- colMeans(lengths.R01.2)

l.skgp[1, ] <- c(l001[1], l003[1], l005[1], l01[1])
l.gp[1, ] <- c(l001[2], l003[2], l005[2], l01[2])
l.bart[1, ] <- c(l001[3], l003[3], l005[3], l01[3])
l.cbart[1, ] <- c(l001[4], l003[4], l005[4], l01[4])
l.rf[1, ] <- c(l001[5], l003[5], l005[5], l01[5])
l.crf[1, ] <- c(l001[6], l003[6], l005[6], l01[6])

# Tori 2K
# Write out and plot results

###### TORI 2K ######

WIDTH = 10
HEIGHT = 10
l.max = 4

##### T001 #####
OUT <- out.T001
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.T001.2 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_T001_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_T001_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_T001_2k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,12))
  # dev.off()
}

##### T003 #####
OUT <- out.T003
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.T003.2 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_T003_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_T003_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_T003_2k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

##### T005 #####
OUT <- out.T005
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.T005.2 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_T005_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_T005_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_T005_2k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

##### T01 #####
OUT <- out.T01
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.T01.2 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_T01_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_T01_2k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_T01_2k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

l001 <- colMeans(lengths.T001.2)
l003 <- colMeans(lengths.T003.2)
l005 <- colMeans(lengths.T005.2)
l01  <- colMeans(lengths.T01.2)

l.skgp[2, ] <- c(l001[1], l003[1], l005[1], l01[1])
l.gp[2, ] <- c(l001[2], l003[2], l005[2], l01[2])
l.bart[2, ] <- c(l001[3], l003[3], l005[3], l01[3])
l.cbart[2, ] <- c(l001[4], l003[4], l005[4], l01[4])
l.rf[2, ] <- c(l001[5], l003[5], l005[5], l01[5])
l.crf[2, ] <- c(l001[6], l003[6], l005[6], l01[6])
################################################################################
# p = 10000 ####################################################################
################################################################################

# Read in results
Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_ROLLS001_p10k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/ROLLS001.rds")
out.R001 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_ROLLS003_p10k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/ROLLS003.rds")
out.R003 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_ROLLS005_p10k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/ROLLS005.rds")
out.R005 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_ROLLS01_p10k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/ROLLS01.rds")
out.R01 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_TORI001_p10k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/TORI001.rds")
out.T001 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_TORI005_p10k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/TORI005.rds")
out.T005 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_TORI003_p10k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/TORI003.rds")
out.T003 <- get_pmse_coverage(Sims, DATA)

Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_TORI01_p10k.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/TORI01.rds")
out.T01 <- get_pmse_coverage(Sims, DATA)

# Write out and plot results

###### Rolls 2K ######

WIDTH = 10
HEIGHT = 10
l.max = 12

##### R001 #####
OUT <- out.R001
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.R001.10 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

X <- DATA[[1]]$X
plot(DATA[[1]]$X[,2], DATA[[1]]$y)
for(ii in 1:40){
  lines
}
colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_R001_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_R001_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_R001_10k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

##### R003 #####
OUT <- out.R003
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.R003.10 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_R003_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_R003_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_R003_10k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

##### R005 #####
OUT <- out.R005
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.R005.10 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_R005_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_R005_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_R005_10k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

##### R01 #####
OUT <- out.R01
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.R01.10 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_R01_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_R01_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_R01_10k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}


# Tori 2K
# Write out and plot results

l001 <- colMeans(lengths.R001.10)
l003 <- colMeans(lengths.R003.10)
l005 <- colMeans(lengths.R005.10)
l01  <- colMeans(lengths.R01.10)

l.skgp[3, ] <- c(l001[1], l003[1], l005[1], l01[1])
l.gp[3, ] <- c(l001[2], l003[2], l005[2], l01[2])
l.bart[3, ] <- c(l001[3], l003[3], l005[3], l01[3])
l.cbart[3, ] <- c(l001[4], l003[4], l005[4], l01[4])
l.rf[3, ] <- c(l001[5], l003[5], l005[5], l01[5])
l.crf[3, ] <- c(l001[6], l003[6], l005[6], l01[6])

###### TORI 2K ######

WIDTH = 10
HEIGHT = 10
l.max = 4

##### T001 #####
OUT <- out.T001
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.T001.10  <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_T001_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_T001_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_T001_10k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

##### T003 #####
OUT <- out.T003
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.T003.10 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_T003_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_T003_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_T003_10k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

##### T005 #####
OUT <- out.T005
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.T005.10 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_T005_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_T005_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_T005_10k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

##### T01 #####
OUT <- out.T01
pmses <- cbind(OUT$pmse.skgp, OUT$pmse.gp, OUT$pmse.bart, OUT$pmse.cbart, OUT$pmse.lasso, OUT$pmse.rf, OUT$pmse.crf)
cvgs <- cbind(OUT$coverage.skgp, OUT$coverage.gp, OUT$coverage.bart, OUT$coverage.cbart, OUT$coverage.rf, OUT$coverage.crf)
lengths.T01.10 <- cbind(OUT$length.skgp, OUT$length.gp, OUT$length.bart, OUT$length.cbart, OUT$length.rf, OUT$length.crf)

colMeans(pmses)
apply(pmses, 2, sd)
{
  pdf(file = "Plots/pmses_T01_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(pmses,
          names = c("SkGP", "GP","BART", "CBART", "LASSO", "RF", "CRF"), ylim = range(pmses), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  pdf(file = "Plots/coverages_T01_10k.pdf", width = WIDTH, height = HEIGHT)
  boxplot(cvgs,
          names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), ylim = c(0,1), main = NA,
          col = "lightskyblue", cex.axis = 1.5)
  dev.off()

  # pdf(file = "Plots/lengths_T01_10k.pdf", width = WIDTH, height = HEIGHT)
  # barplot(colMeans(lengths),
  #         names = c("SkGP", "GP","BART", "CBART", "RF", "CRF"), main = NA, horiz = TRUE, las = 2, xlim = c(0,l.max))
  # dev.off()
}

l001 <- colMeans(lengths.T001.10)
l003 <- colMeans(lengths.T003.10)
l005 <- colMeans(lengths.T005.10)
l01  <- colMeans(lengths.T01.10)

l.skgp[4, ] <- c(l001[1], l003[1], l005[1], l01[1])
l.gp[4, ] <- c(l001[2], l003[2], l005[2], l01[2])
l.bart[4, ] <- c(l001[3], l003[3], l005[3], l01[3])
l.cbart[4, ] <- c(l001[4], l003[4], l005[4], l01[4])
l.rf[4, ] <- c(l001[5], l003[5], l005[5], l01[5])
l.crf[4, ] <- c(l001[6], l003[6], l005[6], l01[6])

################################################################################
# Plotting Lengths #############################################################
################################################################################
xx.grid <- c(1,2,3,4)
xx.labs <- c(0.01, 0.03, 0.05, 0.1)
range(cbind(l.skgp, l.gp, l.bart, l.cbart, l.rf, l.crf))
library(RColorBrewer)
display.brewer.all()
cols <- brewer.pal(6, "Dark2")

pdf(file = "Plots/Swiss_Roll_CI_Length.pdf", width = WIDTH, height = HEIGHT)
par(mfrow = c(1,2))
par(mfrow = c(1, 2))
plot(xx.grid, l.skgp[1, ], "b", lwd = 2, pch = 1, col = cols[1], ylim = c(1, 12),
     xaxt = "n", xlab = "Noise", ylab = "Length of CI", main = "Swiss Roll: p = 2000",
     cex.lab = 1.5, cex.axis = 1.5)
axis(1, at = xx.grid, labels = xx.labs, cex.axis = 1.5)
legend("topleft", legend = c("SkGP", "GP", "BART", "CBART","RF", "CRF"), pch = c(1,2,3,4,5,6), col = cols)
lines(xx.grid, l.gp[1, ], "b", lwd = 2, pch = 2, col = cols[2])
lines(xx.grid, l.bart[1, ], "b", lwd = 2, pch = 3, col = cols[3])
lines(xx.grid, l.cbart[1, ], "b", lwd = 2, pch = 4, col = cols[4])
lines(xx.grid, l.rf[1, ], "b", lwd = 2, pch = 5, col = cols[5])
lines(xx.grid, l.crf[1, ], "b", lwd = 2, pch = 6, col = cols[6])

plot(xx.grid, l.skgp[3, ], "b", lwd = 2, pch = 1, col = cols[1], ylim = c(1, 12),
     xaxt = "n", xlab = "Noise", ylab = NA, main = "p = 10000",
     cex.lab = 1.5, cex.axis = 1.5)
axis(1, at = xx.grid, labels = xx.labs, cex.axis = 1.5)
lines(xx.grid, l.gp[3, ], "b", lwd = 2, pch = 2, col = cols[2])
lines(xx.grid, l.bart[3, ], "b", lwd = 2, pch = 3, col = cols[3])
lines(xx.grid, l.cbart[3, ], "b", lwd = 2, pch = 4, col = cols[4])
lines(xx.grid, l.rf[3, ], "b", lwd = 2, pch = 5, col = cols[5])
lines(xx.grid, l.crf[3, ], "b", lwd = 2, pch = 6, col = cols[6])
dev.off()

pdf(file = "Plots/Torus_CI_Length.pdf", width = WIDTH, height = HEIGHT)
par(mfrow = c(1,2))
plot(xx.grid, l.skgp[2, ], "b", lwd = 2, pch = 1, col = cols[1], ylim = c(1, 5),
     xaxt = "n", xlab = "Noise", ylab = "Length of CI", main = "Torus: p = 2000",
     cex.lab = 1.5, cex.axis = 1.5)
legend("topleft", legend = c("SkGP", "GP", "BART", "CBART","RF", "CRF"), pch = c(1,2,3,4,5,6), col = cols)
axis(1, at = xx.grid, labels = xx.labs, cex.axis = 1.5)
lines(xx.grid, l.gp[2, ], "b", lwd = 2, pch = 2, col = cols[2])
lines(xx.grid, l.bart[2, ], "b", lwd = 2, pch = 3, col = cols[3])
lines(xx.grid, l.cbart[2, ], "b", lwd = 2, pch = 4, col = cols[4])
lines(xx.grid, l.rf[2, ], "b", lwd = 2, pch = 5, col = cols[5])
lines(xx.grid, l.crf[2, ], "b", lwd = 2, pch = 6, col = cols[6])

plot(xx.grid, l.skgp[4, ], "b", lwd = 2, pch = 1, col = cols[1], ylim = c(1, 5),
     xaxt = "n",xlab = "Noise", ylab = NA, main = "p = 10000",
     cex.lab = 1.5, cex.axis = 1.5)
axis(1, at = xx.grid, labels = xx.labs, cex.axis = 1.5)
lines(xx.grid, l.gp[4, ], "b", lwd = 2, pch = 2, col = cols[2])
lines(xx.grid, l.bart[4, ], "b", lwd = 2, pch = 3, col = cols[3])
lines(xx.grid, l.cbart[4, ], "b", lwd = 2, pch = 4, col = cols[4])
lines(xx.grid, l.rf[4, ], "b", lwd = 2, pch = 5, col = cols[5])
lines(xx.grid, l.crf[4, ], "b", lwd = 2, pch = 6, col = cols[6])
dev.off()

################################################################################
Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_ROLLS001_p5K.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/ROLLS001.rds")
out.R001.p5K <- get_pmse_coverage(Sims, DATA)
par(mfrow = c(1,2))
pmses <- cbind(out.R001.p5K$pmse.skgp, out.R001.p5K$pmse.gp, out.R001.p5K$pmse.bart, out.R001.p5K$pmse.cbart, out.R001.p5K$pmse.lasso, out.R001.p5K$pmse.rf, out.R001.p5K$pmse.crf)
boxplot(pmses,
        names = c("SkGP", "GP","BART", "CBT", "L", "RF", "CRF"), ylim = range(pmses), main = "PMSE")
boxplot(cbind(out.R001.p5K$coverage.skgp, out.R001.p5K$coverage.gp, out.R001.p5K$coverage.bart, out.R001.p5K$coverage.cbart, out.R001.p5K$coverage.rf, out.R001.p5K$coverage.crf),
        names = c("SkGP", "GP","BART", "CBT", "RF", "CRF"), ylim = c(0,1), main = "Coverage")


Sims <- readRDS("~/Desktop/Research/SketchedGP/Results/Sims_TORI001_p5K.rds")
DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/TORI001.rds")
out.R001.p5K <- get_pmse_coverage(Sims, DATA)
par(mfrow = c(1,2))
pmses <- cbind(out.R001.p5K$pmse.skgp, out.R001.p5K$pmse.gp, out.R001.p5K$pmse.bart, out.R001.p5K$pmse.cbart, out.R001.p5K$pmse.lasso, out.R001.p5K$pmse.rf, out.R001.p5K$pmse.crf)
boxplot(pmses,
        names = c("SkGP", "GP","BART", "CBT", "L", "RF", "CRF"), ylim = range(pmses), main = "PMSE")
boxplot(cbind(out.R001.p5K$coverage.skgp, out.R001.p5K$coverage.gp, out.R001.p5K$coverage.bart, out.R001.p5K$coverage.cbart, out.R001.p5K$coverage.rf, out.R001.p5K$coverage.crf),
        names = c("SkGP", "GP","BART", "CBT", "RF", "CRF"), ylim = c(0,1), main = "Coverage")


for(ii in sample(1:NN, 2)){
Xstar <- scale(DATA[[ii]]$X[inds.test,], center = TRUE, scale = FALSE)
ystar <- scale(DATA[[ii]]$y[inds.test], center = TRUE, scale = FALSE)
ord = order(Xstar[,2])
plot(Xstar[ord,2], ystar[ord])
lines(Xstar[ord,2], PREDS.skgp[ii, ord], col = 2, lwd = 2)
lines(Xstar[ord,2], PREDS.bart[ii, ord], col = 3, lwd = 2)
lines(Xstar[ord,2], PREDS.cbart[ii, ord], col = 4, lwd = 2)
lines(Xstar[ord,2], PREDS.lasso[ii, ord], col = 5, lwd = 2)
lines(Xstar[ord,2], PREDS.rf[ii, ord], col = 6, lwd = 2)
lines(Xstar[ord,2], PREDS.crf[ii, ord], col = 7, lwd = 2)
lines(Xstar[ord,2], PREDS.skgp[ii, ord], col = 2, lwd = 2)
}

