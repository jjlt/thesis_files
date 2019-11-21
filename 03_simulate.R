# Pre-text ---------------------------------------------------------------

  # Looking at if a larger x matrix works
  # Also uses largish covariance parameters

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
k <- 8
im_j <- -2
re_j <- 7

eps <- generate_eps(n, k, re_j, im_j)

x <- rnorm( 2 * 4 * n, mean = 50, sd = 10)
x <- complex( real = x[1:(4 * n)], imaginary = x[(4 * n + 1):( 2 * 4 * n)])
x <- matrix(x, ncol = 4, byrow = TRUE)
x <- cbind(rep(1, n), x)

re_beta <- c(1, 2, 3, 0, 5)
im_beta <- c(1, 2, 3, 4, 0)
beta <- complex(real= re_beta, imaginary = im_beta)

z <- x %*% beta + eps

beta_ols <- ols(z, x, method = 1)

output <- fgls(z, x, beta_ols)

beta_gls <- output[[1]]
gamma_gls <- output[[2]]
