# Pre-text ----------------------------------------------------------------

  # Need to test the distribution, so the fact that we are estimating gamma
  # and not using the distribution of gamma is going to be a problem. But 
  # we might just try and skip that and say it is a consistent estimator
  # and judging by the sample size of the wind data, it shouldn't be too big
  # of a problem.
  # We want to do significance tests and maybe regions.

# Headers -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')
library(mvtnorm)
library(conics)

# Steps -------------------------------------------------------------------

  # 1) Generate some data, with some beta's at first
  # 2) Look at the estimates for the betas
  # 3) Figure out the theory on what the test statistic should be, perhaps on
  # each component at first with some adjustment?

# Code -------------------------------------------------------------------

set.seed(1)
n <- 500
k <- runif(1, min = 3, max = 6)
im_j <- runif(1, min = -1, max = 1)
re_j <- runif(1, min = -2, max = 2)
j <- complex(real = re_j, imaginary = im_j)
# paste0("This is the true value: ", j)
gamma <- matrix(c(k, j, Conj(j), k), nrow = 2, byrow= TRUE)
eps <- generate_eps(n, k, re_j, im_j)
re_x <- rnorm(n, mean = 10, sd = 30)
im_x <- rnorm(n, mean = 10, sd = 30) # rgamma(n, shape = 5, rate = 5)
x <- complex(real = re_x, imaginary = im_x)
x <- cbind(rep(1, n), x)
beta <- c(complex(real = c(1,5), imaginary = c(-2,3)))
z <- x %*% beta + eps

beta_ols <- ols(z, x)
output <- fgls(z, x, beta_ols)

gamma_est <- output[[2]]
beta_gls <- output[[1]]
beta_gls_aug <- augment_beta(beta_gls)

x_aug <- augment_x(x ,modification = 1)
f_info <- t(Conj(x_aug)) %*% kronecker(diag(1, n), gamma_est) %*% x_aug
V <- svd(f_info)
V_half <- V$u %*% diag(sqrt(V$d)) %*% t(Conj(V$v))
M <- generate_M(2)

M %*% V_half %*% beta_gls_aug
# create a CLM function that works similarly
# okay so the intercept is weaker

# want to look at how an F-Test of sorts will work
# i don't think that will work, because it requires
# maybe could be a likelihood ratio test suggested in ducharme

# ASR <- sum(abs(z - x %*% beta_gls))

# how can you do join then?
# if A and B are true, A union B

# confidence regions should be similar to that in multivariate
# 
# v<-c(64/15,-32/15,64/15,(16/15)*(2*4.7-8*1.5),
#      (16/15)*(2*1.5-8*4.7),(16/15)*(4*1.5^2-2*1.5*4.7+4*4.7^2)
#      -qchisq(0.9,2))
# conicPlot(v, xlab = "x1", ylab = "x2",ylim = c(0,8),xlim = c(0,8))
# find a way to pass complex numbers and the gamma estimate transformed

n <- 500
real_eps <- rnorm(n, mean = 1, sd = 1)
real_beta <- c(2, 3)
real_x <- rnorm(n, sd = 50)
real_x <- cbind(rep(1, n), matrix(real_x))
real_y <- real_x %*% real_beta + real_eps
real.model <- lm(real_y ~ real_x)
summary.lm(real.model)
# statistics not as big
real_beta_mle <- solve(t(real_x) %*% real_x) %*% (t(real_x) %*% real_y)
real_beta_mle - real.model$coefficients[c(1,3)]
# pretty much right

real_sig_mle <- (1 / n) * sum((real_y - real_x %*% real_beta_mle)^2)
real_sig_mle


real_beta_cov <- real_sig_mle * solve(t(real_x) %*% real_x)
udv <- svd(real_beta_cov)
real_beta_cov_half <- udv$u %*% diag(sqrt(udv$d)) %*% t(udv$v)
solve(real_beta_cov_half) %*% real_beta_mle

real_eps2 <- rgamma(n , shape = 0.5)
real_y2 <- log((real_x)^2) %*% real_beta + rnorm(n, mean = - 10, sd = 5)  + real_eps2

real.model2 <- lm(real_y2 ~ real_x)
summary(real.model2)

# let's try some insig data for the complex case
# there is a lot that can do wrong
  # not complex normal errors
  # the wrong relationship
  # mix of both
n <- 500
re_eps_bad <- rgamma(n, shape = 3, rate = 2) + 2
im_eps_bad <- runif(n, min = -5, max = -2)
eps_bad <- complex(real = re_eps_bad, imaginary = im_eps_bad)
z_bad1 <- x %*% beta + eps_bad

bad1_out <- clm(z_bad1, x)
z_scores <- abs(Re(bad1_out))
1 - pnorm(z_scores)
1 - pt(z_scores, df = n - 1)
# does this mean the intercept is insignificant?

re_eps_bad2 <- runif(n, min = - 50, max = 50)
im_eps_bad2 <- runif(n, min = - 50, max = 50)
eps_bad2 <- complex(real = re_eps_bad, imaginary = im_eps_bad)
z_bad2 <- x %*% beta + eps_bad2
bad2_out <- clm(z_bad2, x)

x_bad1 <- complex(real = log(Re(x)^2 + exp(1)) , imaginary = log(Im(x)^2 + 1))
x_bad1 <- matrix(x_bad1, nrow = n)
z_bad3 <- x_bad1 %*% c(0.1, 0.1) + eps_bad/100
bad3_out <- clm(z_bad3, x)

###############

z_bad4 <- rnorm(n, mean = 10, sd = 10)
x_bad2 <- cbind(rep(1, n), seq(1, n) + rnorm(n, sd = 1))
summary(bad4_lm <- lm(z_bad4 ~ x_bad2))
bad4_out <- clm(z_bad4, x_bad2)




######

system.time({out <- clm(z, x, tol = 0.001)})
out
beta0_pvals <- out[c(1, 3), ]
- 2 * log(beta0_pvals) # yea nah m8
# so it does dent the score, but not but much hmm, not noisy enough huh

# need clm to create it's own intercept