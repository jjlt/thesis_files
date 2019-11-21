# Pre-text ---------------------------------------------------------------

# Checking the difference between this method and lm with modulus, hopefully
# it beats lm at determining a complex model

# Header -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')
library(mvtnorm)

# Testing ----------------------------------------------------------------

# 1) create the errors
# 2) create some x data
# 3) choose some beta
# 4) put them together into z
# 5) create ols estimate of beta
# 6) run fgls
# 7) inspect results

n <- 500
k <- 3
re_j <- 1
im_j <- -2

eps <- generate_eps(n, k, re_j, im_j)

x <- rnorm(2 * n, mean = 20, sd = 20) 
x <- complex(real = x[1:n], imaginary = x[(n + 1):(2 * n)])
x <- cbind(rep(1, n), x)

beta <- complex(real = c(1,2), imaginary = c(-2,3))

z <- x %*% beta + eps

beta_ols <- ols(z, x)
beta_gls <- fgls(z, x, beta_ols)[[1]]

abs_x <- abs(x)
lm_output <- lm(abs(z) ~ abs_x)

abs(beta_gls)
lm_output$coefficients
abs(beta)
  # first co-eff is NA because it is the same as the intercept
  # why is it so bad at doing intercepts?
