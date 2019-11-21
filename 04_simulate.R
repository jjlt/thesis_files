# Pre-text ---------------------------------------------------------------

  # Want to look at the convergence when the sample size increases
  # Do I have to sample at n = each multiple times?

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

n_vals <- c(50, 100, 250, 500, 750, 1000)
outcol <- length(n_vals)

k <- 8
im_j <- -2
re_j <- 7
resample <- 20
beta_gls_all <- matrix(list(c(NA, NA)), ncol = outcol, nrow = resample)
outcol <- 1
for( i in n_vals) {
  
  for(j in 1:resample) {
    
    eps <- generate_eps(i, k, re_j, im_j)
    
    x <- rnorm( 2 * i, mean = 50, sd = 10)
    x <- complex( real = x[1:i], imaginary = x[(i + 1):( 2 * i)])
    x <- as.matrix(x)
    x <- cbind(rep(1, i), x)
    
    re_beta <- c(1, 2)
    im_beta <- c(1, 2)
    beta <- complex(real= re_beta, imaginary = im_beta)
    
    z <- x %*% beta + eps
    
    beta_ols <- ols(z, x, method = 1)
    
    beta_gls_all[j, outcol][[1]] <- fgls(z, x, beta_ols)[[1]]
  }
  
  outcol <- outcol + 1
}

beta_0 <- matrix(NA, nrow = resample, ncol = length(n_vals))
beta_1 <- matrix(NA, nrow = resample, ncol = length(n_vals))
for(i in 1:ncol(beta_gls_all)) {
  
  for(j in 1:nrow(beta_gls_all)) {
    
    beta_0[j, i] <- beta_gls_all[j, i][[1]][1]
    beta_1[j, i] <- beta_gls_all[j, i][[1]][2]
  }
}

plot(
  beta[1],
  xlim = c(-1, 3),
  ylim = c(-1, 3),
  xlab = "Re(beta_0)",
  ylab = "Im(beta_0)",
  main = "Increasing accuracy of beta_0"
)
for(i in 1:ncol(beta_gls_all)) {
  
  for(j in 1:nrow(beta_gls_all)) {
    
    points(beta_0[j,i], col = i)
    
  }
}
plot(
  beta[2],
  xlim = c(1.95, 2.05),
  ylim = c(1.95, 2.05),
  xlab = "Re(beta_0)",
  ylab = "Im(beta_0)",
  main = "Increasing accuracy of beta_0"
)
for(i in 1:ncol(beta_gls_all)) {
  
  for(j in 1:nrow(beta_gls_all)) {
    
    points(beta_1[j,i], col = i)
  }
}

complex_covs <- apply(beta_0, 2,
                      function(x) t(x - mean(x)) %*% Conj((x - mean(x))))
psuedo_covs <- apply(beta_0, 2,
                     function(x) t(x - mean(x)) %*% (x - mean(x)))

plot(
  n_vals, Re(complex_covs),
  main = paste0("Sample variance of simulations of different sample sizes"),
  xlab = "Sample Size",
  ylab = "Sample Variance"
)
lines(n_vals, Re(complex_covs))

plot(
  n_vals, abs(psuedo_covs),
  main = paste0("Sample psuedo covariances of 
                simulations of different sample sizes"),
  xlab = "Sample Size",
  ylab = "Sample Psuedo Covariance"
)
lines(n_vals, abs(psuedo_covs))
# make a function that calcs the complex variance?