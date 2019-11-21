# Pre-text ---------------------------------------------------------------

  # Here I want to test the covariance estimator for the beta. In the past
  # the gamma as been shown to be consistent, so that should be fine right?
  # I will take the covariance of the betas, using both gamma estimate
  # and true gamma and changing the sample sizes, see which is better.

# Header -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')
library(mvtnorm)
library(ggplot2)
library(tidyr)

# Steps ------------------------------------------------------------------

  # 1) Loop through samplevals
  # 2) Random x, beta and gamma variables
  # 3) Calculate and store the differences in the covariance

# Code -------------------------------------------------------------------

samplevals <- c(50, 100, 200)
times <- 100
emptymat <- matrix(complex(real = NA, imaginary = NA), nrow = 4, ncol = 4)
cov_true_list <- matrix(list(emptymat), nrow = times, ncol = length(samplevals))
cov_gls_list <- matrix(list(emptymat), nrow = times, ncol = length(samplevals))
cov_btac_list <- matrix(list(emptymat), nrow = times, ncol = length(samplevals))
cov_bgac_list <- matrix(list(emptymat), nrow = times, ncol = length(samplevals))
coln <- 1
for (n in samplevals) {
  for (i in 1:times) {
    k <- runif(1, 4, 10)
    re_j <- runif(1, -2, 2)
    im_j <- runif(1, -2, 2)
    j <- complex(real = re_j, imaginary = im_j)
    gamma <- matrix(c(k, j, Conj(j), k), nrow = 2, byrow = TRUE)
    eps <- generate_eps(n, k, re_j, im_j)
    
    x <- rnorm(2 * n, mean = 20, sd = 100)
    x <- complex(real = x[1:n], imaginary = x[(n + 1):(2 * n)])
    x <- cbind(rep(1, n), matrix(x))
    beta <- rnorm(4, sd = 10)
    beta <- complex(real = beta[1:2], imaginary = beta[3:4])
    beta_aug <- augment_beta(beta)
    
    z <- x %*% beta + eps
    beta_ols <- ols(z, x)
    output <- fgls(z, x, beta_ols)
    gamma_gls <- output[[2]]
    beta_gls <- output[[1]]
    beta_gls_aug <- augment_beta(beta_gls)
    
    gamma_true_2n <- kronecker(diag(1,n), gamma)
    x_aug <- augment_x(x)
    z_aug <- augment_z(z)
    beta_true_aug <- solve(t(Conj(x_aug)) %*% gamma_true_2n %*% x_aug) %*%
      (t(Conj(x_aug)) %*% gamma_true_2n %*% z_aug)
    
    # beta_aug
    # beta_true_aug
    # beta_gls_aug
    
    cov_true <- solve(t(Conj(x_aug)) %*% solve(gamma_true_2n) %*% x_aug)
    cov_true_list[i, coln][[1]] <- cov_true
    cov_gls <- solve(t(Conj(x_aug)) %*% 
                       solve(kronecker(diag(1, n), gamma_gls)) %*% 
                       x_aug)
    cov_gls_list[i, coln][[1]] <- cov_gls
    # round(cov_true, 3)
    # round(cov_gls, 3)
    # feelsweirdman monkaShake
    
    btac <- beta_true_aug - beta_aug
    mu_btac <- apply(btac, 2, mean)
    cov_btac <- (btac - mu_btac) %*% t(Conj(btac - mu_btac))
    cov_btac_list[i, coln][[1]] <- cov_btac
    # round(cov_btac, 3)
    
    bgac <- beta_gls_aug - beta_aug
    mu_bgac <- apply(bgac, 2, mean)
    cov_bgac <- (bgac - mu_bgac) %*% t(Conj(bgac - mu_bgac))
    cov_bgac_list[i, coln][[1]] <- cov_bgac
    # round(cov_bgac, 4)
    print(n)
    print(i)
    
  }
  coln <- coln + 1
}

diff_list_mat <- function(list_mat1, list_mat2) {
  
  stopifnot(dim(list_mat1) == dim(list_mat2))
  fnorm_difference <- matrix(
    NA,
    ncol = dim(list_mat1)[2],
    nrow =dim(list_mat1)[1]
  )
  for (j in 1:dim(list_mat1)[2]) {
    for(i in 1:dim(list_mat1)[1]) {
      fnorm_difference[i, j] <- fnorm(list_mat1[i, j][[1]] -
                                        list_mat2[i, j][[1]])
    }
  }
  return(fnorm_difference)
}

true_true <- diff_list_mat(cov_btac_list, cov_true_list)
colnames(true_true) <- samplevals
true_true <-as.data.frame(true_true)
true_true <- gather(true_true, samplevals, fnorm_of_difference)
true_true$samplevals <- as.numeric(samplevals)

gls_true <- diff_list_mat(cov_bgac_list, cov_true_list)
colnames(gls_true) <- samplevals
gls_true <-as.data.frame(gls_true)
gls_true <- gather(gls_true, samplevals, fnorm_of_difference)
gls_true$samplevals <- as.numeric(samplevals)

gls_gls <- diff_list_mat(cov_bgac_list, cov_gls_list)
colnames(gls_gls) <- samplevals
gls_gls <-as.data.frame(gls_gls)
gls_gls <- gather(gls_gls, samplevals, fnorm_of_difference)
gls_gls$samplevals <- as.numeric(samplevals)

true_true_plot <- ggplot(data = true_true) +
  aes(x = samplevals, group = samplevals, y = fnorm_of_difference)
true_true_plot + geom_boxplot()
true_true_plot + geom_violin()

gls_true_plot <- ggplot(data = gls_true) +
  aes(x = samplevals, group = samplevals, y = fnorm_of_difference)
gls_true_plot + geom_boxplot()
gls_true_plot + geom_violin()

gls_gls_plot <- ggplot(data = gls_gls) +
  aes(x = samplevals, group = samplevals, y = fnorm_of_difference)
gls_gls_plot + geom_boxplot()
gls_gls_plot + geom_violin()

true_true$cov <- "TT"
gls_true$cov <- "GT"
gls_gls$cov <- "GG"

all_cov <- rbind(true_true, gls_true, gls_gls)

p <- ggplot(all_cov) +
  aes(y = fnorm_of_difference, x = factor(samplevals), colour = cov)
p + geom_boxplot()
p + geom_violin()

# doesn't seem to be much difference between all these methods at all levels