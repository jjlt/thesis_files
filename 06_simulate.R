# Pre-text ---------------------------------------------------------------

  # I think here I want to be able to simulate some xGx^1/2 (bhat - b)
  # and somehow show this is CN(0,1) then from there show that it is mvtn
  # or each component is normal or something.

# Header -----------------------------------------------------------------

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source('simulation_headers.R')
library(mvtnorm)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(RColorBrewer)

# Process ----------------------------------------------------------------

  # This is how I think I have to do this
  # Random some Gamma
  # Random some beta
  # Random some X
  # However keeping all the dimensions the same
  # Find the beta_gls
  # Store it somewhere
  # Plot it

  # Should find out how they do it for other methods.
  # What if you keep the same x?

# Well let's try? --------------------------------------------------------

set.seed(1)

reps <- 1000
stored <- matrix(rep(0, 4 * reps), nrow = reps)
M <- generate_M(2)
for(i in 1:reps) {
  n <- 1000
  k <- runif(1, min = 4, max = 10)
  re_j <- runif(1, min = -2, max = 2)
  im_j <- runif(1, min = -2, max = 2)
  eps <- generate_eps(n, k, re_j, im_j)
  re_beta <- runif(2, min = -5, max = 5)
  im_beta <- runif(2, min = -5, max = 5)
  beta <- complex(real = re_beta, imaginary = im_beta)
  x <- runif(2 * n, min = -50, max= 50)
  x <- complex(real = x[1:n], imaginary = x[(n + 1):(2 * n)])
  x <- as.matrix(x)
  x <- cbind(rep(complex(real = 1, imaginary = 1), n), x)
  z <- x %*% beta + eps
  j <- complex(real = re_j, imaginary = im_j)
  gamma <- matrix(c(k , j, Conj(j), k), nrow = 2, byrow =TRUE)
  
  beta_ols <- ols(z, x)
  output <- fgls(z, x, beta_ols)
  beta_gls <- output[[1]]
  gamma_gls <- output[[2]]
  x_aug <- augment_x(x, modification = 1)
  f_info <- t(Conj(x_aug)) %*% 
    kronecker(diag(1,n), solve(gamma_gls)) %*%
    x_aug
  Vsvd <- svd(f_info)
  V_half <- Vsvd$u %*% sqrt(diag(Vsvd$d)) %*% t(Conj(Vsvd$v))
  #temp <- M_half %*% 
  #  (augment_beta(beta_gls) - augment_beta(beta))
  
  # trans_beta_0 <- M %*% temp[c(1,3)]
  # trans_beta_1 <- M %*% temp[c(2,4)]
  # stored[i, ] <- rbind(trans_beta_0, trans_beta_1)
  stored[i, ] <- sqrt(2) * (M %*% V_half %*% 
                              (augment_beta(beta_gls) - augment_beta(beta)))
}

df <- data.frame(Re_b0 = Re(stored[, 1]),
                 Im_b0 = Re(stored[, 3]),
                 Re_b1 = Re(stored[, 2]),
                 Im_b1 = Re(stored[, 4])
)

bw <- 0.2
p <- ggplot(data = df)
p1 <- p + 
  geom_histogram(
    aes(
      x = Re_b0,
      y = ..density..
      ),
    fill = "#9999FF",
    colour = "#3333CC",
    binwidth = bw
    ) + 
  stat_function(
    fun = dnorm,
    size = 2,
    colour = "#0000CC"
  ) +
  theme(legend.position = "none")
p2 <- p + 
  geom_histogram(
    aes(
      x = Im_b0,
      y = ..density..
    ),
    fill = "#9999FF",
    colour = "#3333CC",
    binwidth = bw
  ) + 
  stat_function(
    fun = dnorm,
    size = 2,
    colour = "#0000CC"
  ) +
  theme(legend.position = "none")
p3 <- p + 
  geom_histogram(
    aes(
      x = Re_b1,
      y = ..density..
    ),
    fill = "#9999FF",
    colour = "#3333CC",
    binwidth = bw
  ) + 
  stat_function(
    fun = dnorm,
    size = 2,
    colour = "#0000CC"
  ) +
  theme(legend.position = "none")
p4 <- p + 
  geom_histogram(
    aes(
      x = Im_b1,
      y = ..density..
    ),
    fill = "#9999FF",
    colour = "#3333CC",
    binwidth = bw
  ) + 
  stat_function(
    fun = dnorm,
    size = 2,
    colour = "#0000CC"
  ) +
  theme(legend.position = "none")
grid.arrange(p1, p2, p3, p4, nrow = 2, 
 top = "Standard normal density overlaying histogram of simulated real and imaginary componenets of beta")
# TeX('Densities of the real and imaginary components of $\\widehat{\\beta}$'))

apply(df, 2, mean)
apply(df, 2, sd)

