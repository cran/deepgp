
# Function Contents -----------------------------------------------------------
# Internal:
#   logl: evaluates MVN log likelihood with zero mean
#   sample_g: conducts Metropolis Hastings sampling for nugget
#   sample_theta: conducts Metropolis Hastings sampling for theta
#   sample_w: conducts Elliptical Slice Sampling for w layer
#   sample_z: conducts Elliptical Slice Sampling for z layer

# Log Likelihood Function -----------------------------------------------------
# Calculates log likelihood for multivariate normal distribution with zero mean

logl <- function(out_vec, in_dmat, g, theta, outer = TRUE) {
  
  K <- calc_K(in_dmat, theta, g)
  id <- invdet(K)
  
  if (outer) { # use profile log likelihood (with tau2 integrated out)
    n <- length(out_vec)
    logl <- (-n * 0.5) * log(t(out_vec) %*% id$Mi %*% (out_vec)) - 0.5 * id$ldet
  } else {
    logl <- (-0.5) * id$ldet - 0.5 * (t(out_vec) %*% id$Mi %*% (out_vec))
  }
  return(c(logl))
}

# Sample Nugget Function ------------------------------------------------------
# Completes one iteration of Metropolis Hastings sampling for nugget

sample_g <- function(out_vec, in_dmat, g_t, theta, alpha, beta, l, u, 
                     ll_prev = NULL) {
  
  # Propose value
  g_star <- runif(1, min = l * g_t / u, max = u * g_t / l)
  
  # Compute acceptance threshold
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl(out_vec, in_dmat, g_t, theta, outer = TRUE)
  lpost_threshold <-  ll_prev + dgamma(g_t, alpha, beta, log = TRUE) + 
    log(ru) - log(g_t) + log(g_star)
  
  ll_new <- logl(out_vec, in_dmat, g_star, theta, outer = TRUE)
  
  # Accept or reject
  if (ll_new + dgamma(g_star, alpha, beta, log = TRUE) > lpost_threshold) {
    return(list(g = g_star, ll = ll_new))
  } else{
    return(list(g = g_t, ll = ll_prev))
  }
}

# Sample Theta Function -------------------------------------------------------
# Completes one iteration of Metropolis Hastings sampling for length scale

sample_theta <- function(out_vec, in_dmat, g, theta_t, alpha, beta, l, u, 
                         outer, ll_prev = NULL) {
  
  # Propose value
  theta_star <- runif(1, min = l * theta_t / u, max = u * theta_t / l)
  
  # Compute acceptance threshold
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev))
    ll_prev <- logl(out_vec, in_dmat, g, theta_t, outer)
  lpost_threshold <- ll_prev + dgamma(theta_t, alpha, beta, log = TRUE) + 
    log(ru) - log(theta_t) + log(theta_star)
  
  ll_new <- logl(out_vec, in_dmat, g, theta_star, outer)
  
  # Accept or reject
  if (ll_new + dgamma(theta_star, alpha, beta, log = TRUE) > lpost_threshold) {
    return(list(theta = theta_star, ll = ll_new))
  } else{
    return(list(theta = theta_t, ll = ll_prev))
  }
}

# Elliptical Slice W Function -------------------------------------------------
# Completes one iteration of Elliptical Slice Sampling for a hidden layer

sample_w <- function(out_vec, w_t, w_t_dmat, in_dmat, g, theta_y, theta_w,
                     ll_prev = NULL) {
  
  D <- ncol(w_t) # dimension of hidden layer
  
  if (is.null(ll_prev)) 
    ll_prev <- logl(out_vec, w_t_dmat, g, theta_y, outer = TRUE)
  
  for (i in 1:D) { # separate sampling for each dimension of hidden layer
    
    # Draw from prior distribution
    w_prior <- rand_mvn(1, sigma = calc_K(in_dmat, theta_w[i]))
    
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
      new_logl <- logl(out_vec, dw, g, theta_y, outer = TRUE)
      
      # Accept or reject
      if (new_logl > ll_threshold) { # accept
        ll_prev <- new_logl
        accept <- TRUE
      } else { # reject
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

# Elliptical Slice Z Function -------------------------------------------------
# Completes one iteration of Elliptical Slice Sampling for a hidden layer

sample_z <- function(out_mat, z_t, z_t_dmat, in_dmat, g, theta_w, theta_z,
                     ll_prev = NULL) {
  
  D <- ncol(z_t) # dimension of hidden layer
  
  if (is.null(ll_prev)) {
    ll_prev <- 0
    for (j in 1:D)
      ll_prev <- ll_prev + logl(out_mat[, j], z_t_dmat, g, theta_w[j], 
                                outer = FALSE)
  }
  
  for (i in 1:D) { # separate sampling for each dimension of hidden layer
    
    # Draw from prior distribution
    z_prior <- rand_mvn(1, sigma = calc_K(in_dmat, theta_z[i]))
    
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
                                                 theta_w[j], outer = FALSE)
      
      # Accept or reject
      if (new_logl > ll_threshold) { # accept
        ll_prev <- new_logl
        accept <- TRUE
      } else { # reject
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
