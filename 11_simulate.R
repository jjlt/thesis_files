# Pre-text

  # I want to look into the lack of convergence of the gamma

# Header -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')
library(mvtnorm)

# Testing ----------------------------------------------------------------

set.seed(1)

n <- 2000
k <- runif(1, min = 3, max = 6)
im_j <- runif(1, min = -1, max = 1)
re_j <- runif(1, min = -2, max = 2)
j <- complex(real = re_j, imaginary = im_j)
paste0("This is the true value: ", j)
gamma <- matrix(c(k, j, Conj(j), k), nrow = 2, byrow= TRUE)
eps <- generate_eps(n, k, re_j, im_j)
re_x <- rnorm(n, mean = 10, sd = 30)
im_x <- rnorm(n, mean = 10, sd = 30) # rgamma(n, shape = 5, rate = 5)
x <- complex(real = re_x, imaginary = im_x)
x <- cbind(rep(1, n), x)
beta <- c(complex(real = c(1,5), imaginary = c(-2,3)))
z <- x %*% beta + eps

j_fgls1 <- (1 / n) * t(z - x %*% beta_ols) %*% (z - x %*% beta_ols)
k_fgls1 <- (1 / n) * t(Conj(z - x %*% beta_ols)) %*% (z - x %*% beta_ols)
gamma_fgls1 <- matrix(c(k_fgls1, j_fgls1, Conj(j_fgls1), k_fgls1),
                      nrow = 2, byrow = TRUE)
