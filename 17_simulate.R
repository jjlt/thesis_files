# Pre-text ---------------------------------------------------------------

  # Want to do a comparison of the models in estimating the true values of 
  # beta, there will also be graphical depictions

# Headers ----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')
library(nlme)

# Steps ------------------------------------------------------------------

  # 1) Generate data
  # 2) Obtain estimators with each of the methods
  # 3) Upscale to perform in a loop
  # 4) Graphically depict the differences

# Code -------------------------------------------------------------------

set.seed(1)
samplevals <- c(50, 100, 200, 400, 800)
times <- 1000
lm_diff_b0 <- matrix(NA, nrow = times, ncol = length(samplevals))
lm_diff_b1 <- lm_diff_b0
clm_diff_b0 <- lm_diff_b0
clm_diff_b1 <- lm_diff_b0
mlm_diff_b0 <- lm_diff_b0
mlm_diff_b1 <- lm_diff_b0
# gls_diff_b0 <- lm_diff_b0 
# gls_diff_b1 <- lm_diff_b0
j <- 1
for(n in samplevals) {
  for(i in 1:times) {
    k <- runif(1, min = 4, max = 10)
    re_j <- runif(1, min = -2, max = 2)
    im_j <- runif(1, min = -2, max = 2)
    eps <- generate_eps(n, k, re_j, im_j)
    
    beta <- rnorm(4, mean = 0, sd = 5)
    beta <- complex(real = beta[1:2], imaginary = beta[3:4])
    
    x <- rnorm(2 * n, mean = 20, sd = 30)
    x <- complex(real = x[1:n], imaginary = x[(n + 1):(2 * n)])
    X <- cbind(rep(1, n), matrix(x))
    
    z <- X %*% beta + eps
    df <- data.frame(cbind(Re(z),Im(z),Re(x),Im(x),seq(1, n)))
    colnames(df) <- c("z", "z", "x", "x","id")
    df <- rbind(df[, c(1,3,5)], df[, c(2,4,5)])
    df <- df[order(df$id), ]
    df$obsnum <- rep(c(1,2), n)
    
    temp_lm <- lm(abs(z) ~ abs(x))$coefficients
    temp_clm <- abs(clm(z, X)$coefficients)
    temp_mlm <- lm(cbind(Re(z),Im(z)) ~ Re(x) + Im(x))$coefficients
    temp_mlm <- apply(temp_mlm, 1, function(x) sqrt(sum(x^2)))
    #temp_gls <- gls(z ~ x, data = df,
         #corSymm(form = ~ obsnum | id), method = "ML")$coefficients
    b0 <- abs(beta)[1]
    b1 <- abs(beta)[2]
    lm_diff_b0[i, j] <- abs(abs(temp_lm[1]) - b0)
    lm_diff_b1[i, j] <- abs(abs(temp_lm[2]) - b1)
    clm_diff_b0[i, j] <- abs(abs(temp_clm[1]) - b0)
    clm_diff_b1[i, j] <- abs(abs(temp_clm[2]) - b1)
    mlm_diff_b0[i, j] <- abs(abs(temp_mlm[1]) - b0)
    mlm_diff_b1[i, j] <- abs(abs(mean(temp_mlm[2:3])) - b1)
    #gls_diff_b0[i, j] <- abs(abs(temp_gls[1]) - b0)
    #gls_diff_b1[i, j] <- abs(abs(temp_gls[2]) - b1)
  }
  j <- j + 1
}
apply(lm_diff_b0, 2, mean)
apply(lm_diff_b1, 2, mean)
apply(clm_diff_b0, 2, mean)
apply(clm_diff_b1, 2, mean)
apply(mlm_diff_b0, 2, mean)
apply(mlm_diff_b1, 2, mean)
#apply(gls_diff_b0, 2, mean)
#apply(gls_diff_b1, 2, mean)

library(ggplot2)
library(tidyr)

head(lm_diff_b0_long)
colnames(lm_diff_b0) <- samplevals
lm_diff_b0 <- as.data.frame(lm_diff_b0)
lm_diff_b0_long <- gather(lm_diff_b0, Sample_size, Difference)
lm_diff_b0_long$Method <- "lm"
colnames(clm_diff_b0) <- samplevals
clm_diff_b0 <- as.data.frame(clm_diff_b0)
clm_diff_b0_long <- gather(clm_diff_b0, Sample_size, Difference)
clm_diff_b0_long$Method <- "clm"
colnames(mlm_diff_b0) <- samplevals
mlm_diff_b0 <- as.data.frame(mlm_diff_b0)
mlm_diff_b0_long <- gather(mlm_diff_b0, Sample_size, Difference)
mlm_diff_b0_long$Method <- "mlm"

b0_df <- rbind(lm_diff_b0_long, mlm_diff_b0_long, clm_diff_b0_long)
b0_df$Sample_size <- as.numeric(b0_df$Sample_size)

b0_plot <- ggplot(b0_df) + 
  aes(x = factor(Sample_size),
      y = Difference,
      colour = Method,
      fill = Method
  ) + theme_bw() + xlab('Sample Size') +
  theme(plot.title = element_text(hjust = 0.5))
b0_plot + geom_boxplot(alpha = 0.3)
b0_plot + geom_violin(position = "identity", alpha = 0.3)

#####
b0_df <- rbind(lm_diff_b0_long)
#####
b0_df$Sample_size <- as.numeric(b0_df$Sample_size)

b0_plot <- ggplot(b0_df) + 
  aes(x = factor(Sample_size),
      y = Difference,
      colour = Method,
      fill = Method
  ) + theme_bw() + xlab('Sample Size') +
  theme(plot.title = element_text(hjust = 0.5))
b0_plot + geom_boxplot(alpha = 0.3)
b0_plot + geom_violin(position = "identity", alpha = 0.3)


colnames(lm_diff_b1) <- samplevals
lm_diff_b1 <- as.data.frame(lm_diff_b1)
lm_diff_b1_long <- gather(lm_diff_b1, Sample_size, Difference)
lm_diff_b1_long$Method <- "lm"
colnames(clm_diff_b1) <- samplevals
clm_diff_b1 <- as.data.frame(clm_diff_b1)
clm_diff_b1_long <- gather(clm_diff_b1, Sample_size, Difference)
clm_diff_b1_long$Method <- "clm"
colnames(mlm_diff_b1) <- samplevals
mlm_diff_b1 <- as.data.frame(mlm_diff_b1)
mlm_diff_b1_long <- gather(mlm_diff_b1, Sample_size, Difference)
mlm_diff_b1_long$Method <- "mlm"

##################################3

b1_df <- rbind(clm_diff_b1_long)
#################################
b1_df$Sample_size <- as.numeric(b1_df$Sample_size)

b1_plot <- ggplot(b1_df) + 
  aes(x = factor(Sample_size),
      y = Difference,
      colour = Method,
      fill = Method
  ) + theme_bw() + xlab('Sample Size') +
  theme(plot.title = element_text(hjust = 0.5))
b1_plot + geom_boxplot(alpha = 0.3)
b1_plot + geom_violin(position = "identity", alpha = 0.3)

