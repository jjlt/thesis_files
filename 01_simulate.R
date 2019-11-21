# Pre-text ---------------------------------------------------------------

  # Initial sand box

# Header -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')
library(mvtnorm)

# Testing ----------------------------------------------------------------

n <- 50
k <- 2
re_j <- -1
im_j <- 1
tol <- 0.00000001
stopifnot(k^2 - re_j^2 - im_j^2 > 0)

sigma_chi <- matrix(
  c(
    (k + re_j) / 2,
    -(1/2) * im_j,
    -(1/2) * im_j,
    (k - re_j) / 2
  )
  , nrow = 2,
  byrow = TRUE
)

eps <- rmvnorm(n, mean = rep(0, 2), sigma = sigma_chi)
eps <- complex( real = eps[,1], imaginary = eps[,2])

x <- rnorm(2 * n, mean = 20, sd = 10)
x <- complex(real = x[1:n], imaginary = x[(n + 1):(2 * n)])
x <- cbind(
  rep(1,n),
  x
)

beta <- c(5, 2)

z <- x %*% beta + eps

z_aug <- augment_z(z)
x_aug <- augment_x(x)

beta_ols_conj <- ols(z, x, method = 2)
beta_ols <- ols(z, x, method = 1)

  # if we augment these, they just double, insert proof
  # but i think the augmented versions are need for the 
  # gls of beta

k_est1 <- (1 / n) * sum(abs(z - x %*% beta) ^ 2)
k_conj_est1 <- (1 / n) * sum(abs(z_aug - x_aug %*% beta_ols_conj) ^ 2)
# big hmm

j_est1 <- (1 / n) * sum((z - x  %*% beta) ^ 2)

gamma_est1 <- matrix(
  c(
    k_est1,
    j_est1,
    Conj(j_est1),
    k_est1
  ),
  nrow = 2,
  byrow = TRUE
)

beta_gls_est1 <- solve(Conj(t(x_aug)) %*% 
                         kronecker(diag(1, n), solve(gamma_est1)) %*% 
                         x_aug) %*%
  Conj(t(x_aug)) %*%
  kronecker(diag(1, n), solve(gamma_est1)) %*%
  z_aug

### BIG WARNING

sum((z - x  %*% beta_ols[1:2]) ^ 2)
sum((z_aug - x_aug  %*% beta_ols) ^ 2)

# this is because K and J should be calculated not with augmentation

# Try to find convergence with fgls -----------------------------------------
  
  # z and x are not augmented observation and data
  # beta is an ols estimate, whether that be using using transpose method
  # or transpose conjugate method.
  # try keep it very compartmentalised, so this function does as little as
  # possible.


output <- fgls(z,x,beta_ols)
beta_gls <- output[[1]]
gamma_gls <- output[[2]]
fgls(z,x,beta_ols_conj)
fgls(z,x,beta_ols)
# Doesn't matter, same convergence it seems

# Want to find the variance of the beta estimates
# xGx should be the fisher information in Fred's calculations
xGx<- t(x_aug) %*% kronecker(diag(1, n), solve(gamma_gls)) %*% Conj(x_aug)

# I'll try to figure what my term looks like, hopefully it is actually different

beta_aug <- rbind(as.matrix(beta_gls), as.matrix(Conj(beta_gls)))

xGx - xGx %*% Conj(beta_aug) %*% t(beta_aug) %*% xGx

# So this doesn't make too much sense

x_aug_1 <- x_aug[1:2,]

f_info_1 <- Conj(t(x_aug_1)) %*% solve(gamma_gls) %*% x_aug_1

# https://pdfs.semanticscholar.org/2250/d1af99672b9d415ee1bf1e2b7143527a2ecd.pdf
# http://ita.calit2.net/workshop/11/files/paper/paper_2379.pdf

# now i think i want to look at the distribution of the estimators and
# seeing how the fisher information calculations can help in making some 
# test statistics. 
