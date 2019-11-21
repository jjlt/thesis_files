# Pre-text ----------------------------------------------------------------

  # Wind

# Headers -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')

# Steps -------------------------------------------------------------------

wind <- read.csv("wind_data.csv", na = c("#N/A", "VRB", "28s", "34s", "*"))
head(wind)
wind <- as.data.frame(wind)
wind <- na.omit(wind)
z_arg <- true_arg(wind$HourlyWindDirection)
z <- complex(modulus = wind$HourlyWindSpeed, argument = z_arg)
x_arg <- true_arg(wind$PrevHourlyWindDirection)
x <- complex(modulus = wind$PrevHourlyWindSpeed, argument = x_arg)
comb <- cbind(z, rep(1, length(z)), matrix(x))

set.seed(1)
tot <- seq(1, dim(comb)[1])
samp <- sample(tot, size = 200, replace = FALSE)

model_clm <- clm(comb[samp, 1], comb[samp, 2:3], tol = 0.001)
model_clm[[1]]

summary.clm(model_clm)
model_clm$coefficients