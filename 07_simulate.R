# Pre-text ---------------------------------------------------------------

  # Want to show that as the samples sizes increase, the variance of the 
  # residuals decrease, do quite a larger simulation for this, the code
  # will be similar to that of 04_simulate.R, but I'll also try to use
  # ggplot2

# Header -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')
library(mvtnorm)
library(ggplot2)
library(tidyr)
library(latex2exp)

# Method  -----------------------------------------------------------------

  # (1) Create a vector of how much I want the sample sizes to be
  # (2) Choose how much you want each sample to be tested. Temporarily,
  # keep the tests at some low value, before the actual execution.
  # (3) Create the loops that will create the data, using all the 
  # functions and yadda yadda. Will make all the data generation,
  # and the choosing beta, k, j all random
  # (4) Store the created data and create residuals. For everything,
  # the beta, k, j, reals and imaginaries.
  # (5) Find out how to plot the damn thing with bar graphs and stuff.
  # (6) Jazz it up. Figureout how I am going to plot all of it.
  # (7) Ramp up the testing for each sample size and save the results.

# Code --------------------------------------------------------------------

set.seed(1)

samplevals <- c(2000)
tests <- 10

empty_vec <- c(NA, NA)
empty_mat <- matrix(NA, nrow = 2, ncol = 2)
beta_res <- matrix(list(empty_vec), nrow = tests, ncol = length(samplevals))
gamma_res <- matrix(list(empty_mat), nrow = tests, ncol = length(samplevals))

coln <- 1
for (i in samplevals) {
  for (j in 1:tests) {
    k <- runif(1,min = 3, max = 5)
    re_j <- runif(1, min = -2, max = 2)
    im_j <- runif(1, min = -2, max = 2)
    true_j <- complex(real = re_j, imaginary = im_j)
    print(paste0("The true value: ", true_j ))
    eps <- generate_eps(i, k, re_j, im_j)
    beta <- runif(4, min = -5, max = 5)
    beta <- complex(real = beta[1:2], imaginary = beta[3:4])
    x <- rnorm(2 * i, mean = 0, sd = 20)
    x <- complex(real = x[1:i], imaginary = x[(i + 1):(2 * i)])
    x <- cbind(rep(1, i), matrix(x))
    z <- x %*% beta + eps
    
    beta_ols <- ols(z, x, method = 1)
    output <- fgls(z, x, beta_ols)
    b_res <- beta - output[[1]]
    gamma <- matrix(
      c(
        k,
        complex(real = re_j, imaginary = im_j),
        complex(real = re_j, imaginary =  - im_j), 
        k)
      , nrow = 2, byrow = TRUE)
    # print(paste0("This is the true: ", gamma))
    # print(paste0("This is the est: ", output[[2]]))
    g_res <- gamma - output[[2]]
    # print(i)
    print(j)
    beta_res[j, coln][[1]] <- b_res
    gamma_res[j, coln][[1]] <- g_res
  }
  coln <- coln + 1
}

# read the wiki on FLGS
# check how it works for he real case

beta0_res <- matrix(sapply(beta_res, '[[', 1), nrow = tests, byrow = FALSE)
beta1_res <- matrix(sapply(beta_res, '[[', 2), nrow = tests, byrow = FALSE)
k_res <- Re(matrix(sapply(gamma_res, '[[', 1), nrow = tests, byrow = FALSE))
j_res <- matrix(sapply(gamma_res, '[[', 3), nrow = tests, byrow = FALSE)
re_j_res <- abs(Re(j_res))
im_j_res <- abs(Im(j_res))
colnames(k_res) <- samplevals
k_res <- as.data.frame(abs(k_res))
colnames(re_j_res) <- samplevals
re_j_res <- as.data.frame(abs(re_j_res))
colnames(im_j_res) <- samplevals
im_j_res <- as.data.frame(abs(im_j_res))

k_res_long <- gather(k_res, Sample_size, Residual)
k_res_long$Parameter <- "k"
re_j_res_long <- gather(re_j_res, Sample_size, Residual)
re_j_res_long$Parameter <- "re_j"
im_j_res_long <- gather(im_j_res, Sample_size, Residual)
im_j_res_long$Parameter <- "im_j"

gamma_df <- rbind(k_res_long, re_j_res_long, im_j_res_long)
gamma_df$Sample_size <- as.numeric(gamma_df$Sample_size)
gamma_plot <- ggplot(gamma_df) +
  aes(x = factor(Sample_size),
      y = Residual,
      colour = Parameter,
      fill = Parameter
      ) +
  xlab('Sample Size') +
  ylab(TeX('| $\\widehat{\\theta} - \\theta$|')) +
  labs(title = TeX('Consistency in each of the components of $\\Gamma$')) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
gamma_plot + geom_boxplot(alpha = 0.3)
gamma_plot + geom_violin(position = "identity", alpha = 0.3)

re_b0_res <- abs(Re(beta0_res))
colnames(re_b0_res) <- samplevals
re_b0_res_long <- gather(as.data.frame(re_b0_res), Sample_size, Residual)
re_b0_res_long$Component <- "re"
im_b0_res <- abs(Im(beta0_res))
colnames(im_b0_res) <- samplevals
im_b0_res_long <- gather(as.data.frame(im_b0_res), Sample_size, Residual)
im_b0_res_long$Component <- "im"
b0_df <- rbind(im_b0_res_long, re_b0_res_long)
b0_df$Sample_size <- as.numeric(b0_df$Sample_size)
b0_plot <- ggplot(b0_df) +
  aes(x = factor(Sample_size),
      y = Residual,
      colour = Component,
      fill = Component
      ) +
  xlab('Sample Size') +
  ylab(TeX('| $\\widehat{\\beta_{0}} - \\beta_{0}$|')) +
  labs(title =
         TeX(
           'Consistency in the real and imaginary components of $\\beta_{0}$')
       ) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
b0_plot + geom_boxplot(alpha = 0.3)
b0_plot + geom_violin(alpha = 0.3, position = "identity")

re_b1_res <- abs(Re(beta1_res))
colnames(re_b1_res) <- samplevals
re_b1_res_long <- gather(as.data.frame(re_b1_res), Sample_size, Residual)
re_b1_res_long$Component <- "re"
im_b1_res <- abs(Im(beta1_res))
colnames(im_b1_res) <- samplevals
im_b1_res_long <- gather(as.data.frame(im_b1_res), Sample_size, Residual)
im_b1_res_long$Component <- "im"
b1_df <- rbind(im_b1_res_long, re_b1_res_long)
b1_df$Sample_size <- as.numeric(b1_df$Sample_size)
b1_plot <- ggplot(b1_df) +
  aes(x = factor(Sample_size),
      y = Residual,
      colour = Component,
      fill = Component
      ) +
  xlab('Sample Size') +
  ylab(TeX('| $\\widehat{\\beta_{1}} - \\beta_{1}$|')) +
  labs(title =
         TeX(
           'Consistency in the real and imaginary components of $\\beta_{1}$')
  ) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
b1_plot + geom_boxplot(alpha = 0.3)
b1_plot + geom_violin(alpha = 0.3, position = "identity")

