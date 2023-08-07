
# Function Contents -----------------------------------------------------------
# Internal:
#   logl: evaluates MVN log likelihood with zero mean
#   sample_g: conducts Metropolis Hastings sampling for nugget
#   sample_theta: conducts Metropolis Hastings sampling for theta
#   sample_w: conducts Elliptical Slice Sampling for w layer
#   sample_z: conducts Elliptical Slice Sampling for z layer

# Log Likelihood --------------------------------------------------------------

logl <- function(out_vec, in_dmat, g, theta, outer = TRUE, v, tau2 = FALSE,
                 mu = 0, scale = 1) {
  
  n <- length(out_vec)
  if (v == 999) {
    K <- scale * Exp2(in_dmat, 1, theta, g)
  } else K <- scale * Matern(in_dmat, 1, theta, g, v) 
  id <- invdet(K)
  quadterm <- t(out_vec - mu) %*% id$Mi %*% (out_vec - mu)
  
  if (outer) { # use profile log likelihood (with tau2 integrated out)
    logl <- (- n * 0.5) * log(quadterm) - 0.5 * id$ldet
  } else logl <- (- 0.5) * id$ldet - 0.5 * quadterm
  
  if (tau2) {
    tau2 <- c(quadterm) / n
  } else tau2 <- NULL
  
  return(list(logl = c(logl), tau2 = tau2))
}

# Log Likelihood SEPARABLE ----------------------------------------------------
# Only utilized in a one layer GP (i.e. outer = TRUE)

logl_sep <- function(out_vec, in_mat, g, theta, v, tau2 = FALSE) {
  
  n <- length(out_vec)
  if (v == 999) {
    K <- Exp2Sep(in_mat, in_mat, 1, theta, g)
  } else K <- MaternSep(in_mat, in_mat, 1, theta, g, v) 
  id <- invdet(K)
  quadterm <- t(out_vec) %*% id$Mi %*% out_vec
  
  # outer = TRUE, profile likelihood with tau2 integrated out
  logl <- (- n * 0.5) * log(quadterm) - 0.5 * id$ldet
  
  if (tau2) {
    tau2 <- c(quadterm) / n
  } else tau2 <- NULL
  
  return(list(logl = c(logl), tau2 = tau2))
}

# Sample G --------------------------------------------------------------------

sample_g <- function(out_vec, in_dmat, g_t, theta, alpha, beta, l, u, 
                     ll_prev = NULL, v) {
  
  # Propose value
  g_star <- runif(1, min = l * g_t / u, max = u * g_t / l)
  
  # Compute acceptance threshold
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl(out_vec, in_dmat, g_t, theta, outer = TRUE, v)$logl
  lpost_threshold <-  ll_prev + dgamma(g_t - eps, alpha, beta, log = TRUE) + 
    log(ru) - log(g_t) + log(g_star)
  
  ll_new <- logl(out_vec, in_dmat, g_star, theta, outer = TRUE, v)$logl
  
  # Accept or reject (lower bound of eps)
  new <- ll_new + dgamma(g_star - eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) {
    return(list(g = g_star, ll = ll_new))
  } else {
    return(list(g = g_t, ll = ll_prev))
  }
}

# Sample G SEPARABLE ----------------------------------------------------------
# Only used in one-layer GP (outer = TRUE, tau2 = FALSE only)

sample_g_sep <- function(out_vec, in_mat, g_t, theta, alpha, beta, l, u, 
                         ll_prev = NULL, v) {
  
  # Propose value
  g_star <- runif(1, min = l * g_t / u, max = u * g_t / l)
  
  # Compute acceptance threshold
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_sep(out_vec, in_mat, g_t, theta, v)$logl
  lpost_threshold <-  ll_prev + dgamma(g_t - eps, alpha, beta, log = TRUE) + 
    log(ru) - log(g_t) + log(g_star)
  
  ll_new <- logl_sep(out_vec, in_mat, g_star, theta, v)$logl
  
  # Accept or reject (lower bound of eps)
  new <- ll_new + dgamma(g_star - eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) {
    return(list(g = g_star, ll = ll_new))
  } else {
    return(list(g = g_t, ll = ll_prev))
  }
}

# Sample Theta ----------------------------------------------------------------

sample_theta <- function(out_vec, in_dmat, g, theta_t, alpha, beta, l, u, 
                         outer, ll_prev = NULL, v, tau2 = FALSE,
                         prior_mean = 0, scale = 1) {
  
  # Propose value
  theta_star <- runif(1, min = l * theta_t / u, max = u * theta_t / l)
  
  # Compute acceptance threshold
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl(out_vec, in_dmat, g, theta_t, outer, v, mu = prior_mean,
                    scale = scale)$logl
  
  lpost_threshold <- ll_prev + dgamma(theta_t - eps, alpha, beta, log = TRUE) + 
    log(ru) - log(theta_t) + log(theta_star)
  
  ll_new <- logl(out_vec, in_dmat, g, theta_star, outer, v, tau2 = tau2,
                 mu = prior_mean, scale = scale)

  # Accept or reject (lower bound of eps)
  new <- ll_new$logl + dgamma(theta_star - eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) {
    return(list(theta = theta_star, ll = ll_new$logl, tau2 = ll_new$tau2))
  } else {
    return(list(theta = theta_t, ll = ll_prev, tau2 = NULL))
  }
}

# Sample Theta SEPARABLE ------------------------------------------------------
# Only used in one-layer GP (outer = TRUE only)

sample_theta_sep <- function(out_vec, in_mat, g, theta_t, index = 1,
                             alpha, beta, l, u, ll_prev = NULL, v, tau2 = FALSE) {
  
  # Propose value
  theta_star <- runif(1, min = l * theta_t[index] / u, max = u * theta_t[index] / l)
  theta_t_updated <- theta_t
  theta_t_updated[index] <- theta_star
  
  # Compute acceptance threshold
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_sep(out_vec, in_mat, g, theta_t, v)$logl
  
  lpost_threshold <- ll_prev + dgamma(theta_t[index] - eps, alpha, beta, log = TRUE) + 
    log(ru) - log(theta_t[index]) + log(theta_star)
  
  ll_new <- logl_sep(out_vec, in_mat, g, theta_t_updated, v, tau2 = tau2)
  
  # Accept or reject (lower bound of eps)
  new <- ll_new$logl + dgamma(theta_star - eps, alpha, beta, log = TRUE)
  
  if (new > lpost_threshold) {
    return(list(theta = theta_star, ll = ll_new$logl, tau2 = ll_new$tau2))
  } else {
    return(list(theta = theta_t[index], ll = ll_prev, tau2 = NULL))
  }
}

# Elliptical Slice W ----------------------------------------------------------

sample_w <- function(out_vec, w_t, w_t_dmat, in_dmat, g, theta_y, theta_w,
                     ll_prev = NULL, v, 
                     prior_mean = matrix(0, nrow = nrow(w_t), ncol = ncol(w_t)),
                     scale = 1) {
  
  D <- ncol(w_t) # dimension of hidden layer

  if (is.null(ll_prev)) 
    ll_prev <- logl(out_vec, w_t_dmat, g, theta_y, outer = TRUE, v = v)$logl
  
  count <- vector(length = D)
  
  for (i in 1:D) { # separate sampling for each dimension of hidden layer
    
    # Draw from prior distribution
    if (v == 999) {
      w_prior <- mvtnorm::rmvnorm(1, mean = prior_mean[, i], 
                                  sigma = scale * Exp2(in_dmat, 1, theta_w[i], 0))
    } else {
      w_prior <- mvtnorm::rmvnorm(1, mean = prior_mean[, i], 
                                  sigma = scale * Matern(in_dmat, 1, theta_w[i], 0, v))
    }
    
    # Initialize a and bounds on a
    a <- runif(1, min = 0, max = 2 * pi)
    amin <- a - 2 * pi
    amax <- a
    
    # Compute acceptance threshold - based on all dimensions of previous w
    ru <- runif(1, min = 0, max = 1)
    ll_threshold <- ll_prev + log(ru)
    
    # Calculate proposed values, accept or reject, repeat if necessary
    accept <- FALSE
    count <- 0
    w_prev <- w_t[, i] # store for re-proposal
    
    while (accept == FALSE) {
      count <- count + 1
      
      # Calculate proposed values and new likelihood
      w_t[, i] <- w_prev * cos(a) + w_prior * sin(a)
      dw <- sq_dist(w_t)
      
      new_logl <- logl(out_vec, dw, g, theta_y, outer = TRUE, v = v)$logl
      
      # Accept or reject
      if (new_logl > ll_threshold) {
        ll_prev <- new_logl
        accept <- TRUE
      } else {
        # update the bounds on a and repeat
        if (a < 0) {
          amin <- a
        } else {
          amax <- a
        }
        a <- runif(1, amin, amax)
        if (count > 100) stop('reached maximum iterations of ESS')
      } # end of else statement
    } # end of while loop
  } # end of i for loop
  
  return(list(w = w_t, ll = ll_prev, dw = dw))
}

# Elliptical Slice Z ----------------------------------------------------------

sample_z <- function(out_mat, z_t, z_t_dmat, in_dmat, g, theta_w, theta_z,
                     ll_prev = NULL, v) {
  
  D <- ncol(z_t) # dimension of hidden layer
  
  if (is.null(ll_prev)) {
    ll_prev <- 0
    for (j in 1:D)
      ll_prev <- ll_prev + logl(out_mat[, j], z_t_dmat, g, theta_w[j], 
                                outer = FALSE, v = v)$logl
  }
  
  for (i in 1:D) { # separate sampling for each dimension of hidden layer
    
    # Draw from prior distribution
    if (v == 999) {
      z_prior <- mvtnorm::rmvnorm(1, sigma = Exp2(in_dmat, 1, theta_z[i], 0))
    } else {
      z_prior <- mvtnorm::rmvnorm(1, sigma = Matern(in_dmat, 1, theta_z[i], 0, v))
    }
    
    # Initialize a and bounds on a
    a <- runif(1, min = 0, max = 2 * pi)
    amin <- a - 2 * pi
    amax <- a
    
    # Compute acceptance threshold - based on all dimensions of previous z
    ru <- runif(1, min = 0, max = 1)
    ll_threshold <- ll_prev + log(ru)
    
    # Calculate proposed values, accept or reject, repeat if necessary
    accept <- FALSE
    count <- 0
    z_prev <- z_t[, i] # store for re-proposal
    
    while (accept == FALSE) {
      count <- count + 1
      
      # Calculate proposed values and new likelihood
      z_t[, i] <- z_prev * cos(a) + z_prior * sin(a)
      dz <- sq_dist(z_t)
      new_logl <- 0
      for (j in 1:D) new_logl <- new_logl + logl(out_mat[, j], dz, g, 
                                                 theta_w[j], outer = FALSE,
                                                 v = v)$logl
      
      # Accept or reject
      if (new_logl > ll_threshold) {
        ll_prev <- new_logl
        accept <- TRUE
      } else { 
        # update the bounds on a and repeat
        if (a < 0) {
          amin <- a
        } else {
          amax <- a
        }
        a <- runif(1, amin, amax)
        if (count > 100) stop('reached maximum iterations of ESS')
      } # end of else statement
    } # end of while loop
  } # end of i for loop
  
  return(list(z = z_t, ll = ll_prev, dz = dz))
}
