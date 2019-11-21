# 

# Header -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')
library(mvtnorm)

# Testing ----------------------------------------------------------------

n <- 50
k <- 1
im_j <- 0
re_j <- 0
j <- complex(real = re_j, imaginary = im_j)

gamma <- matrix(c(k, j, Conj(j), k), nrow = 2, byrow= TRUE)
eps <- generate_eps(n, k, im_j, re_j)
re_x <- rnorm(n, mean = 10, sd = 30)
im_x <-rgamma(n, shape = 5, rate = 5)
x <- complex(real = re_x, imaginary = im_x)
x <- cbind(rep(1, n), x)
beta <- c(complex(real = c(1,5), imaginary = c(-2,3)))
z <- x %*% beta + eps

beta_ols <- ols(z, x)
output <- fgls(z, x, beta_ols)
beta_est <- output[[1]]
gamma_est <- output[[2]]

z_aug <- augment_z(z, modification = 1)
x_aug <- augment_x(x, modification = 1)
beta_aug <- augment_beta(beta_ols)
sumgam <- 0
for (i in 1:50) {
  sumgam <- (z_aug[(2 * i - 1):(2 * i)] -
               x_aug[(2 * i - 1):(2 * i),] %*% beta_aug) %*%
    t(Conj(z_aug[(2 * i - 1):(2 * i)] -
             x_aug[(2 * i - 1):(2 * i),] %*% beta_aug))
  + sumgam
}
sumgam
# should be (1 / n) * sumgam
# might need to make even more modifications to the z and beta etc,
# i think the statement before 2.2.5 doesn't make sense

k_est <- (1 / n) * sum(abs(z - x %*% beta_ols) ^ 2)
j_est <- (1 / n) * sum((z - x  %*% beta_ols) ^ 2)
gamma_jk <- matrix(
  c(k_est, j_est, Conj(j_est), k_est),
  nrow = 2,
  byrow = TRUE
)

