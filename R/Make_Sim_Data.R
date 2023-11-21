
# Create and save simulation data
set.seed(01131995)

n.sets <- 40
n.roll <- 600

source("R/functions.R")

ROLLS.001 <- list()
ROLLS.003 <- list()
ROLLS.005 <- list()
ROLLS.01  <- list()

TORI.001  <- list()
TORI.003 <- list()
TORI.005  <- list()
TORI.01   <- list()
for(ii in 1:n.sets){
  ROLLS.001[[ii]] <- get_swiss_roll(n.roll, 10000, 0.01)
  ROLLS.003[[ii]] <- get_swiss_roll(n.roll, 10000, 0.03)
  ROLLS.005[[ii]] <- get_swiss_roll(n.roll, 10000, 0.05)
  ROLLS.01[[ii]]  <- get_swiss_roll(n.roll, 10000, 0.1)
  TORI.001[[ii]]  <- get_torus(n.roll, 10000, 0.01)
  TORI.003[[ii]]  <- get_torus(n.roll, 10000, 0.03)
  TORI.005[[ii]]  <- get_torus(n.roll, 10000, 0.05)
  TORI.01[[ii]]   <- get_torus(n.roll, 10000, 0.1)
}

saveRDS(ROLLS.001, "Data/ROLLS001.rds")
saveRDS(ROLLS.003, "Data/ROLLS003.rds")
saveRDS(ROLLS.005, "Data/ROLLS005.rds")
saveRDS(ROLLS.01, "Data/ROLLS01.rds")
saveRDS(TORI.001, "Data/TORI001.rds")
saveRDS(TORI.003, "Data/TORI003.rds")
saveRDS(TORI.005, "Data/TORI005.rds")
saveRDS(TORI.01, "Data/TORI01.rds")
