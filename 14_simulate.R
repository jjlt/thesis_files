# Pre-text ----------------------------------------------------------------

  # non-linear solve for the MLE

# Headers -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')
library(nleqslv)
library(mvtnorm)

# Steps -------------------------------------------------------------------

  # 1) Create some data
  # 2) Store the required constants for the non linear system
  # 3) Have to manually create the non linear system, find a way to automate
  # this

n <- 50
set.seed(1)
x <- rnorm(2 * n, sd = 20)
x <- complex(real = x[1:n], imaginary = x[(n + 1):(2 * n)])
x <- cbind(rep(1, n), matrix(x))
k <- 5
j <- complex(real = -3, imaginary = 2)

eps <- generate_eps(n, k, Re(j), Im(j))

beta <- complex(real = c(0, 4), imaginary = c(-2, 3))

z <- x %*% beta + eps

x_aug <- augment_x(x)
z_aug <- augment_z(z)

fn <- function(theta) {
  
  y <- numeric(7)
  gamma <- matrix(c(theta[5],
                    complex(real = theta[6], imaginary = theta[7]),
                    complex(real = theta[6], imaginary = - theta[7]),
                    theta[5]),
                    nrow = 2,
                    byrow = TRUE)
  result <- solve(t(Conj(x_aug)) %*%
                    kronecker(diag(1, n), solve(gamma)) %*%
                    x_aug) %*%
    (t(Conj(x_aug)) %*% kronecker(diag(1, n), solve(gamma)) %*% z_aug)
  y[1] <- (result[1] + result[3]) / 2
  y[2] <- (result[1] - result[3]) / 2
  y[3] <- (result[2] + result[4]) / 2
  y[4] <- (result[1] - result[3]) / 2
  y[5] <- (1 / n) * sum(abs(
    z - x %*% complex(
      real = c(theta[1], theta[3]),
      imaginary = c(theta[2], theta[4]))
    )^2)
  y[6] <- (1 / n) * Re(sum((
    z - x %*% complex(
      real = c(theta[1], theta[3]),
      imaginary = c(theta[2], theta[4]))
  )^2))
  y[7] <- (1 / n) * Im(sum((
    z - x %*% complex(
      real = c(theta[1], theta[3]),
      imaginary = c(theta[2], theta[4]))
  )^2))
  y
}

thetastart <- c(0.1, -1.9, 3.9, 2.9, 4.9, -2.9, 2)
fn(thetastart)
nleqslv(thetastart, fn, control = list(btol = 0.01))

# RIP
