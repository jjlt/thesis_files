# Pre-text ----------------------------------------------------------------

  # This will be the header file containing a lot of the reuseable
  # functions.

# Create the 2n by 2n M matrix function ----------------------------------

  # This is just useful to have around
  # The M matrix found in Ducharme

generate_M <- function(n) {
  
  M <- 1/2 * rbind(
    cbind(
      diag(n),
      diag(n)), 
    cbind(
      matrix(complex(imaginary = -diag(n)), nrow = n),
      matrix(complex(imaginary = diag(n)), nrow = n)
    )
  )
  return(M)
}

# Augment the data x -----------------------------------------------------

  # Description
    # Augments the data in two forms, modification = 1 gives the following
    #         |  x_1    0    |
    #         |   0    x_1*  |
    #         |   .     .    |
    # x_aug = |   .     .    |
    #         |   .     .    |
    #         |  x_n    0    |
    #         |   0    x_n*  |
    # modification = 2, gives this block format
    #         | x_data    0    |
    # x_aug = |                |
    #         |   0     x_data*|
  # Arguments
    # x - non-augmented data matrix

augment_x <- function(x, modification = 1) {
  
  l <- dim(x)[1]
  w <- dim(x)[2]
  x_aug <- matrix(rep(0, 4 * l * w), nrow = 2 * l)
  
  if (modification == 1) {
    
    indices <- 2 * seq(1, l) - 1
    x_aug[indices, 1:w] <- x
    x_aug[indices + 1, (w + 1):(2 * w)] <- Conj(x)
    
  } else if (modification == 2) {

    x_aug[(1:l), (1:w)] <- x
    x_aug[((l + 1):(2 * l)), ((w + 1):(2 * w))] <- Conj(x)
    
  } else {
    
    stop("Choose modification = 1 or 2")
    
  }
  return(x_aug)
}

# Augment the complex observations z -------------------------------------

  # Description
    # Makes a vector in the form
    # z_aug = (z_1, z_1*, ..., z_n, z_n*) or (z^{T}, z^{H})
  # Arguments
    # z - complex vector of observations
    # modificiation - 1) for first augment, 2) for second augment

augment_z <- function(z, modification = 1) {
  
  l <- length(z)
  z_aug <- rep(0, 2 * l)
  if (modification == 1) {
    indices <- 2 * seq(1, l) - 1
    z_aug[indices] <- z
    z_aug[indices + 1] <- Conj(z)
  } else if (modification == 2) {
    z <- matrix(z)
    z_aug = rbind(z, Conj(z))
  } else {
    
    stop("Choose modification = 1 or 2")
    
  }
  
  return(z_aug)
}

# OLS -------------------------------------------------------------------

  # Description
    # Generate an OLS estimate for the regression parameter beta. Uses
    # observations and data, z and x respectively. Two methods are
    # interesting, one using transpose of the data and the other using
    # the conjugate transpose.
  # Arguments
    # z - non-augmented observation vector.
    # x - non-augmented data vector (yet to test matrix).
    # method - method = 1 for transpose, method = 2 for the conjugate
    # transpose.

ols <- function(z, x, method = 1) {
  
  if (method == 1) {
    
    return(solve(t(x) %*% x) %*% t(x) %*% z)
    
  } else if (method == 2) {
    
    return(solve(Conj(t(x)) %*% x) %*% Conj(t(x)) %*% z)
    
  }
}

# FGLS ------------------------------------------------------------------

  # Description
    # Using an OLS estimator of the regression finds a gamma estimate
    # (error covariance), then uses that to find another beta estimate
    # (regression parameter) and so on until convergence with respect to
    # some tolerance or until a maximum number of iterations.
  # Arguments
    # z - non-augmented observation vector.
    # x - non-augmented data vector (yet to test matrix).
    # tol - tolerance for the convergence of the estimator.
    # max - maximum number of iterations before giving up converging.
  # Return
    # On no error, returns a list, [[1]] for converged beta estimate 
    # (regression parameters), [[2]] for the converged gamma matrix 
    # (error covariance matrix).
    # On error, "Did not converge ... iterations", returns estimates
    # of the about mentioned beta and gamma, but at the "max" iterations.

fgls <- function(z, x, beta, tol = 0.00000001, max = 10) {
  
  z_aug <- augment_z(z)
  x_aug <- augment_x(x)
  n <- length(z)
  
  for (i in (1:max)) {
    # print(beta)
    beta_prev <- beta
    k <- (1 / n) * sum(abs(z - x %*% beta) ^ 2)
    j <- (1 / n) * sum((z - x  %*% beta) ^ 2)
    # print(paste0("Calculated: ", j))
    # j <- complex(real = Im(j), imaginary = Re(j))
    # print(paste0("Manipulated: ", j))
    gamma <- matrix(
      c(k, j, Conj(j), k),
      nrow = 2,
      byrow = TRUE
    )
    # print(gamma)
    
    beta <- solve(Conj(t(x_aug)) %*% 
                    kronecker(
                      diag(1, n), 
                      solve(gamma)) %*% 
                    x_aug) %*%
      (Conj(t(x_aug)) %*%
      kronecker(
        diag(1, n), 
        solve(gamma)) %*%
      z_aug)
    beta <- beta[1:dim(x)[2]]
    # print(beta)
    if (sum(abs(beta_prev - beta)) < tol) {
      
      return(list(beta,gamma))
      
    }
  }
  
  return(list(beta,gamma))
  
  stop("Did not converge in ", max, " iterations")
  
}

# Generate n random errors -------------------------------------------------------

  # Description
    # Generates the covariance matrix for complex normal errors then for
    # however many samples you want. Requires some chosen covariance
    # parameters.

generate_eps <- function(n, k, re_j, im_j) {
  
  stopifnot(k^2 - re_j^2 - im_j^2 > 0)
  
  sigma_chi <- matrix(
    c(
      (k + re_j) / 2,
      (1 / 2) * im_j,
      (1 / 2) * im_j,
      (k - re_j) / 2
    )
    , nrow = 2,
    byrow = TRUE
  )
  eps <- rmvnorm(n, mean = rep(0, 2), sigma = sigma_chi)
  eps <- complex(real = eps[,1], imaginary = eps[,2])
  
  return(eps)
}

# Augment beta -----------------------------------------------------------

  # Description
    # I think this is only used for 06_simulate for now.
    # augmenting in a very simple way, beta_aug = (beta, Conj(beta))
  # Arguments
    # beta - the beta vector to be augmented

augment_beta <- function(beta) {
  beta <- as.matrix(beta)
  return(rbind(beta, Conj(beta)))
}

# Frobenius norm (complex matrix) ----------------------------------------

  # Description
    # Takes the mod of each of the elements and adds them, then square roots the
    # result. If you use a real matrix it will return the same thing as
    # norm(M, "f").
  # Arguments
    # complex_matrix - a matrix, does not have to be complex
fnorm <- function(complex_matrix) {
  return(sqrt(sum(Mod(complex_matrix)^2)))
}

# Complex Linear Model ----------------------------------------------------

  # Description
    # lm(), but complex, using all that we have done here
    # for now will give bonferroni adjusted confidence intervals
    # as well as intervals for the real and the complex
    # Think about what else do I want in here, god this thing is slow lol

clm <-function(z, x, tol = 0.00000001) {
  n <- length(z)
  beta_ols <- ols(z, x)
  output <- fgls(z, x, beta_ols, tol)
  beta_est <- output[[1]]
  gamma_est <- output[[2]]
  x_aug <- augment_x(x)
  f_info <- t(Conj(x_aug)) %*% kronecker(diag(1, n), gamma_est) %*% x_aug
  V <- svd(f_info)
  V_half <- V$u %*% diag(sqrt(V$d)) %*% t(Conj(V$v))
  M <- generate_M(dim(x)[2])
  beta_gls_aug <- augment_beta(beta_est)
  z_scores <- abs(Re(M %*% V_half %*% beta_gls_aug))
  pvals <- 1 - pnorm(z_scores)
  both_pvals <- beta_pvals(z_scores)
  return_list <- list(coefficients = beta_est,
                      modulus = Mod(beta_est),
                      argument = Arg(beta_est),
                      true_bearing = true_arg(Arg(beta_est), method = 2),
                      cov = gamma_est,
                      z_scores = z_scores,
                      re_im_sig = pvals,
                      beta_sig = both_pvals)
  return(return_list)
}


# CDF Bivariate Normal -------------------------------------------------

library(mvtnorm)
beta_pvals <- function(z_scores) {
  betas <- matrix(z_scores, ncol = 2)
  return(apply(betas, 1, function(x) 1 - pmvnorm(upper = x)))
}

# Summary for clm ------------------------------------------------------

summary.clm <- function(model) {
  print("        Summary output for Complex Linear Model        ")
  print("-------------------------------------------------------")
  print("       Complex Regression Co-efficient Estimates       ")
  print("-------------------------------------------------------")
  print("Parameter                     Estimate")
  p <- length(model[[1]])
  for (i in 0:(p - 1)) {
    print(paste0("Beta", i, "         ", model$coefficients[i + 1]))
  }
  print("-------------------------------------------------------")
  print("Complex Regression Co-efficient Estimates in Polar Form")
  print("-------------------------------------------------------")
  print("Parameter       Modulus                Argument")
  for (i in 0:(p - 1)) {
    print(paste0("Beta", i, "        ",
                 model$modulus[i + 1], "      ", model$argument[i + 1]))
  }
  print("-------------------------------------------------------")
  print("          Z-Scores of each component of beta           ")
  print("-------------------------------------------------------")
  print("Component                Score")
  for (i in 0:(p - 1)) {
    print(paste0("Real Beta", i, "         ", model$z_scores[i+ 1]))
    print(paste0("Imaginary Beta", i, "    ", model$z_scores[i + p + 1]))
  }
  print("-------------------------------------------------------")
  print("          p-values of each component of beta           ")
  print("-------------------------------------------------------")
  print("Component                p-value")
  for (i in 0:(p - 1)) {
    print(paste0("Real Beta", i, "         ", model$re_im_sig[i+ 1]))
    print(paste0("Imaginary Beta", i, "    ", model$re_im_sig[i + p + 1]))
  }
  print("-------------------------------------------------------")
  print("                p-values of each beta                  ")
  print("-------------------------------------------------------")
  print("Beta                p-value")
  for (i in 0:(p - 1)) {
    print(paste0("Beta", i, "         ", model$beta_sig[i + 1]))
  }
  print("-------------------------------------------------------")
  print("    The estimate for the gamma, covariance matrix      ")
  print("-------------------------------------------------------")
  print(model$cov)
}

# True Beraing to arg and back -------------------------------------------

true_arg <- function(vec, method = 1) {
  
  if (method == 1) {
    vec <- vec - 90
    vec <- vec - (vec >= 180)*360
    vec <- - 1 * vec
    vec <- vec *(pi / 180)
  } else if (method == 2) {
    vec <- vec*(180 / pi)
    vec <- - 1 * vec
    vec <- vec + (vec < - 90)*360
    vec <- vec + 90
  }
  return(vec)
}
