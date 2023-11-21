DATA <- readRDS("~/Desktop/Research/SketchedGP/Data/ROLLS001.rds")
source("R/functions.R")

df <- DATA[[2]]

ds <- plgp::distance(DATA[[1]]$X)

P <- get_sketch_mat(60, 10000,orthog = FALSE)
pds <- plgp::distance(df$X%*%t(P))

par(mfrow = c(2,1))
hist(pds[lower.tri(pds)], breaks = 51, main = "Distances in Projected Space",
     probability = TRUE, xlab = NA)

hist(ds[lower.tri(ds)], breaks = 51, main = "Distances in Original Space",
     probability = TRUE)
