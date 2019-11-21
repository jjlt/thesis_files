# pre-text ----------------------------------------------------------------

# Look at the new time variates

# headers -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')
library(mvtnorm)

# steps -------------------------------------------------------------------

  # 1) Load new data, clean it, make it complex run clm

# code --------------------------------------------------------------------

wind <- read.csv("wind_data2.csv",
                 header = T,
                 na = c("#N/A", "VRB", "28s", "34s", "*"))
wind <- na.omit(as.data.frame(wind))
z <- complex(modulus = wind$HourlyWindSpeed, argument = true_arg(wind$HourlyWindDirection))
x <- complex(modulus = wind$PrevHourlyWindSpeed, argument = true_arg(wind$PrevHourlyWindDirection))
head(wind)
time1 <- complex(argument = wind$TimeOfDay * 2 * pi)
time2 <- complex(argument = wind$TimeOfYear * 2 * pi)
intercept <- rep(1, length(z))
data <- cbind(z, intercept, x, time1, time2)

n <- 500
set.seed(1)
samps <- sample(seq(1, length(z)), n)

testing <- data[samps, ]
model.clm.testing <- clm(testing[, 1], testing[, 2:5])
model.clm.testing
summary.clm(model.clm.testing)
