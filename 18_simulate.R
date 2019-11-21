# Pre-text ----------------------------------------------------------------

  # need to modify the data for wind to give 24hrs before hand, try with 
  # the newer data set, might have trouble validating

# Headers -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')

# Code --------------------------------------------------------------------

wind <- read.csv("wind_data2.csv",
                 header = T,
                 na = c("#N/A", "VRB", "28s", "34s", "*"))
wind <- na.omit(as.data.frame(wind))
z <- complex(modulus = wind$HourlyWindSpeed,
             argument = true_arg(wind$HourlyWindDirection))
x <- complex(modulus = wind$PrevHourlyWindSpeed,
             argument = true_arg(wind$PrevHourlyWindDirection))
head(wind)
time1 <- complex(argument = wind$TimeOfDay * 2 * pi)
time2 <- complex(argument = wind$TimeOfYear * 2 * pi)
intercept <- rep(1, length(z))
data <- cbind(z, intercept, x, time1, time2)
head(data, 100)
day_before <- data[1:(dim(data)[1] - 71), 1]
new_data <- cbind(data[72:dim(data)[1], ], day_before)
new_data <- new_data[, c(1,2,6)]

n <- 500
set.seed(1)
samps <- sample(seq(1, dim(new_data)[1]), n)
samped <- new_data[samps, ]
model.clm <- clm(samped[, 1], samped[, 2:3])
summary.clm(model.clm)

new_data2 <- cbind(data[72:dim(data)[1], ], day_before)
new_data2 <- new_data2[, c(1, 2, 6, 4, 5)]
samped2 <- new_data2[samps, ]
model.clm2 <- clm(samped2[, 1], samped2[, 2:5])
summary.clm(model.clm2)

#clm2SA <- sum(abs(samped2[, 1] - mean(samped[, 1])))
#clm2SR <- sum(abs(samped2[, 1] - model.clm2$coefficients %*% t(samped2[, 2:5])))

# might be worth looking into some sort of R squared