# Pre-text ---------------------------------------------------------------

  # Want to try the t test, I think it should be t test, because you
  # are estimating the variance as well.

# Header -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')
library(mvtnorm)

# Code -------------------------------------------------------------------

k <- 5
re_j <- -2
im_j <- 1
n <- 5000

eps <- generate_eps(n, k, re_j, im_j)
x <- rnorm(2 * n, mean = 0, sd = 100)
x <- complex(real = x[1:n], imaginary = x[(n + 1):(2 * n)])
x <- cbind(rep(1, n), matrix(x))
beta <- complex(real = c(-2, 3), imaginary = c(1, -2))

z <- x %*% beta + eps

beta_ols <- ols(z, x)

output <- fgls(z, x, beta_ols)
gamma_gls <- output[[2]]
beta_gls <- output[[1]]
j <- complex(real = re_j, imaginary = im_j)
gamma <- matrix(c(k, j, Conj(j), k), nrow = 2, byrow = TRUE)
det(abs(gamma_gls - gamma))
beta - beta_gls
