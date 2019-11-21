# pre-text ----------------------------------------------------------------

  # i want to do comparisons with lm and mvt tests, we'll see how these go

# headers -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')

# steps -------------------------------------------------------------------

  # 1) load the wind data, perhaps take a smaller cut of it
  # 2) make each of the models then compare the results

# code --------------------------------------------------------------------

# Lm taking abs method ----------------------------------------------------

wind <- read.csv("wind_data.csv",
                 header = T,
                 na = c("#N/A", "VRB", "28s", "34s", "*"))
wind <- na.omit(as.data.frame(wind))
training_predicting <- 500
predicting <- training_predicting / 5
training <- training_predicting - predicting
set.seed(2)
samps <- sample(seq(1, dim(wind)[1]), training_predicting)
training <- samps[1:training]
predicting <- samps[(length(training) + 1):training_predicting]
training_data <- wind[training, ]
predicting_data <- wind[predicting, ]
model.lm <- lm(HourlyWindSpeed ~ PrevHourlyWindSpeed, data = training_data)

# preds.lm <- predict(model.lm, newdata = predicting_data)
preds.lm2 <- cbind(
  rep(1, dim(predicting_data)[1]),
  predicting_data$PrevHourlyWindSpeed) %*% 
  model.lm$coefficients
# what if i did it the same method as clm?
#score.lm <- predicting_data$HourlyWindSpeed - preds.lm
#sum(abs(score.lm))
#sum(abs(resid(model.lm)))
score.lm2 <- predicting_data$HourlyWindSpeed - preds.lm2
sum(abs(score.lm2))

# Our CLM method ----------------------------------------------------------

z <- complex(modulus = wind$HourlyWindSpeed, argument = true_arg(wind$HourlyWindDirection))
x <- complex(modulus = wind$PrevHourlyWindSpeed, argument = true_arg(wind$PrevHourlyWindDirection))
complex_wind <- cbind(z, rep(1, length(z)), x)
rownames(complex_wind) <- NULL
complex_training_data <- complex_wind[training, ]
complex_predicting_data <- complex_wind[predicting, ]

model.clm <- clm(complex_training_data[, 1], complex_training_data[,2:3])
preds.clm <- complex_predicting_data[,2:3] %*% model.clm[1][[1]]
score.clm <- complex_predicting_data[, 1] - preds.clm
sum(abs(score.clm))

# Standard MLR method, lm cbind -------------------------------------------

mlm_wind <- cbind(cbind(Re(complex_wind[, 1]),
                  Im(complex_wind[, 1]),
                  Re(complex_wind[, 3]),
                  Im(complex_wind[, 3])))
mlm_wind_training <- mlm_wind[training, ]
mlm_wind_predicting <- mlm_wind[predicting, ]
model.mlm <- lm(
  cbind(mlm_wind_training[, 1], mlm_wind_training[, 2])
  ~ mlm_wind_training[, 3] + mlm_wind_training[, 4])
preds.mlm.real <- model.mlm$coefficients[2:3, 1] %*%
  t(mlm_wind_predicting[, 3:4]) + model.mlm$coefficients[1, 1]
preds.mlm.imag <- model.mlm$coefficients[2:3, 2] %*%
  t(mlm_wind_predicting[, 3:4]) + model.mlm$coefficients[1, 2]
score.mlm <- complex(real = preds.mlm.real - mlm_wind_predicting[, 1],
                     imaginary = preds.mlm.imag - mlm_wind_predicting[, 2])
score.mlm <- sum(abs(score.mlm))
score.mlm

# GLS method, symmetric covariance structure ------------------------------

mlm_wind <- rbind(mlm_wind_training, mlm_wind_predicting)
pairs_wind <- as.data.frame(mlm_wind)
colnames(pairs_wind) <- c("Wind", "Wind", "PrevWind", "PrevWind")

  # new_df will be stacking real parts on imaginary parts
  # new_df2 will be stacking real parts on imaginary parts only for the reponse
  # for the predictors, keeping both real and imaginary

pairs_wind$id <- seq(1, dim(pairs_wind)[1])
new_df <- rbind(pairs_wind[, c(1, 3, 5)], pairs_wind[, c(2, 4, 5)])
new_df <- new_df[order(new_df$id), ]
new_df$obsnum <- rep(c(1, 2), training_predicting)
library(nlme)
model.gls <- gls(Wind ~ PrevWind, data = new_df[(1:(2*length(training))), ],
                 correlation = corSymm(form = ~ obsnum | id), method = "ML")
new_df2 <- new_df[, c(1, 3:4)]
indices <- seq(1, dim(pairs_wind)[1])
new_df2[(2 * indices - 1), 4:5] <- pairs_wind[, 3:4]
new_df2[(2 * indices), 4:5] <- pairs_wind[, 3:4]
colnames(new_df2) <- c("Wind", "id", "obsnum", "ReWind", "ImWind")
model.gls2 <- gls(Wind ~ ReWind + ImWind, data = new_df2[(1:(2*length(training))), ],
                  correlation = corCompSymm(form = ~ obsnum | id), method = "ML")

new_df_pred <- new_df[(2 * length(training) + 1):(2 * training_predicting), ]
new_df2_pred <- new_df2[(2 * length(training) + 1):(2 * training_predicting), ]

score.gls <- new_df_pred$Wind - 
  cbind(rep(1, 2 * length(predicting)), new_df_pred$PrevWind) %*% 
  model.gls$coefficients
score.gls <- sum(abs(score.gls))
score.gls

score.gls2 <- new_df2_pred$Wind - 
  cbind(rep(1, 2 * length(predicting)), new_df2_pred$ReWind, new_df2_pred$ImWind) %*% 
  model.gls2$coefficients
score.gls2 <- sum(abs(score.gls2))
score.gls2
