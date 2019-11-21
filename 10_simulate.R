# Pre-text ---------------------------------------------------------------

# Produce some violin plots of the real case.

# Header -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyr)
library(ggplot2)
library(latex2exp)

# Steps ------------------------------------------------------------------

  # 1) Have to generate some data. Create some real eps, create some beta
  # and x all in a loop, that goes through different sample sizes as well
  # 2) Inside, fit an lm, get the co-effs, plot em

# Code -------------------------------------------------------------------

samplevals <- c(50, 100, 200, 400, 800)
tests <- 1000
b0_diffs <- matrix(NA, nrow = tests, ncol = length(samplevals))
b1_diffs <- matrix(NA, nrow = tests, ncol = length(samplevals))
sig_diffs <- matrix(NA, nrow = tests, ncol = length(samplevals))
coln <- 1
for (n in samplevals) {
  for(i in 1:tests) {
    sigma_true <- runif(1, min = 0, max = 1)^2
    eps <- rnorm(n, sd = sigma_true)
    beta <- runif(2, min = -10, max = 10)
    x <- rnorm(n, sd= 50)
    x <- cbind(rep(1, n), matrix(x))
    y <- x %*% beta + eps
    
    output <- lm(y ~ x[,2])
    b0_diffs[i, coln] <- abs(beta[1] - output$coefficients[1])
    b1_diffs[i, coln] <- abs(beta[2] - output$coefficients[2])
    sig_diffs[i, coln] <- abs(sigma_true - sigma(output))
  }
  coln <- coln + 1
  
}

colnames(b0_diffs) <- samplevals
b0_diffs <- as.data.frame(b0_diffs)
b0_diffs <- gather(b0_diffs, samplevals, abs_estimate_true_diff)
b0_diffs$samplevals <- as.numeric(b0_diffs$samplevals)
b0_plot <- ggplot(data = b0_diffs) +
  aes(
    x = samplevals,
    group = samplevals,
    y = abs_estimate_true_diff,
    colours = 'cyan'
    ) + xlab('Sample Size') + 
  ylab(TeX('| $\\widehat{\\beta_{0}} - \\beta_{0}$|')) +
  labs(title = 'Consistency property at work') +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
b0_plot + geom_boxplot(alpha = 0.3, fill = 'cyan')
b0_plot + geom_violin(alpha = 0.3, fill = 'cyan')

colnames(b1_diffs) <- samplevals
b1_diffs <- as.data.frame(b1_diffs)
b1_diffs <- gather(b1_diffs, samplevals, abs_estimate_true_diff)
b1_diffs$samplevals <- as.numeric(b1_diffs$samplevals)
b1_plot <- ggplot(data = b1_diffs) +
  aes(
    x = samplevals,
    group = samplevals,
    y = abs_estimate_true_diff,
    colours = 'cyan'
  ) + xlab('Sample Size') + 
  ylab(TeX('| $\\widehat{\\beta_{1}} - \\beta_{1}$|')) +
  labs(title = 'Consistency property at work') +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
b1_plot + geom_boxplot(alpha = 0.3, fill = 'cyan')
b1_plot + geom_violin(alpha = 0.3, fill = 'cyan')

colnames(sig_diffs) <- samplevals
sig_diffs <- as.data.frame(sig_diffs)
sig_diffs <- gather(sig_diffs, samplevals, abs_estimate_true_diff)
sig_diffs$samplevals <- as.numeric(sig_diffs$samplevals)
sig_plot <- ggplot(data = sig_diffs) +
  aes(
    x = samplevals,
    group = samplevals,
    y = abs_estimate_true_diff,
    colours = 'cyan'
  ) + xlab('Sample Size') + 
  ylab(TeX('| $\\widehat{\\sigma^{2}} - \\sigma^{2}$|')) +
  labs(title = 'Consistency property at work') +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))
sig_plot + geom_boxplot(alpha = 0.3, fill = 'cyan')
sig_plot + geom_violin(alpha = 0.3, fill = 'cyan')
