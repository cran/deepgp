
# Function Contents -----------------------------------------------------------
# Internal:
#   logl_vec: calculates vecchia log likelihood
#   sample_g_vec: conducts Metropolis Hastings sampling for nugget
#   sample_theta_vec: conducts Metropolis Hastings sampling for theta
#   sample_w_vec: conducts Elliptical Slice Sampling for w layer
#   sample_z_vec: conducts Elliptical Slice Sampling for z layer

# Vecchia Log Likelihood-------------------------------------------------------

logl_vec <- function(out_vec, approx, g, theta, outer = TRUE, v, tau2 = FALSE) {
  
  n <- length(out_vec)
  out_vec_ord <- out_vec[approx$ord]
  U_mat <- create_U(approx, g, theta, v)
  Uty <- Matrix::crossprod(U_mat, out_vec_ord)
  ytUUty <- sum(Uty^2)
  logdet <- sum(log(Matrix::diag(U_mat)))
  
  if (outer) {
    logl <- logdet - (n * 0.5) * log(ytUUty)
  } else {
    logl <- logdet - 0.5 * ytUUty
  }
  
  if (tau2) {
    tau2 <- c(ytUUty) / n
  } else tau2 <- NULL
  
  return(list(logl = logl, tau2 = tau2))
}

# Sample G Vecchia ------------------------------------------------------------

sample_g_vec <- function(y, g_t, theta, alpha, beta, l, u, ll_prev = NULL, 
                         approx, v) {

  # Propose value
  g_star <- runif(1, min = l * g_t / u, max = u * g_t / l)

  # Compute acceptance threshold
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl_vec(y, approx, g_t, theta, outer = TRUE, v)$logl
  lpost_threshold <-  ll_prev + dgamma(g_t - eps, alpha, beta, log = TRUE) + 
                      log(ru) - log(g_t) + log(g_star)
 
  ll_new <- logl_vec(y, approx, g_star, theta, outer = TRUE, v)$logl
  
  # Accept or reject (lower bound of eps)
  new <- ll_new + dgamma(g_star - eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) {
    return(list(g = g_star, ll = ll_new))
  } else{
    return(list(g = g_t, ll = ll_prev))
  }
}

# Sample Theta Vecchia --------------------------------------------------------

sample_theta_vec <- function(y, g, theta_t, alpha, beta, l, u, outer, 
                              ll_prev = NULL, approx, v, tau2 = FALSE) {

  # Propose value
  theta_star <- runif(1, min = l * theta_t / u, max = u * theta_t / l)

  # Compute acceptance threshold
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev))
    ll_prev <- logl_vec(y, approx, g, theta_t, outer, v)$logl
  lpost_threshold <- ll_prev + dgamma(theta_t - eps, alpha, beta, log = TRUE) + 
                      log(ru) - log(theta_t) + log(theta_star)
  
  ll_new <- logl_vec(y, approx, g, theta_star, outer, v, tau2 = tau2)

  # Accept or reject (lower bound of eps)
  new <- ll_new$logl + dgamma(theta_star - eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) {
    return(list(theta = theta_star, ll = ll_new$logl, tau2 = ll_new$tau2))
  } else{
    return(list(theta = theta_t, ll = ll_prev, tau2 = NULL))
  }
}

# Elliptical Slice W Vecchia --------------------------------------------------

sample_w_vec <- function(y, w_approx, x_approx, g, theta_y, theta_w, 
                         ll_prev = NULL, v) {

  D <- ncol(w_approx$x_ord) # dimension of hidden layer

  if (is.null(ll_prev)) 
    ll_prev <- logl_vec(y, w_approx, g, theta_y, outer = TRUE, v)$logl

  for (i in 1:D) { # separate sampling for each dimension of hidden layer

    w_prior <- rand_mvn_vec(x_approx, theta_w[i], v)
    
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
    w_prev <- w_approx$x_ord[w_approx$rev_ord_obs, i] # store for re-proposal
    
    while (accept == FALSE) {
      count <- count + 1

      # Calculate proposed values and new likelihood
      w_proposal <- w_prev * cos(a) + w_prior * sin(a)
      
      # Incorporate proposal in vecchia approximation object
      w_approx <- update_obs_in_approx(w_approx, w_proposal, i)

      new_logl <- logl_vec(y, w_approx, g, theta_y, outer = TRUE, v)$logl

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

  return(list(w_approx = w_approx, ll = ll_prev))
}

# Elliptical Slice Z Vecchia --------------------------------------------------

sample_z_vec <- function(w, z_approx, x_approx, g, theta_w, theta_z, 
                         ll_prev = NULL, v) {
  
  D <- ncol(z_approx$x_ord) # dimension of hidden layer
  
  if (is.null(ll_prev)) {
    ll_prev <- 0
    for (j in 1:D)
      ll_prev <- ll_prev + logl_vec(w[, j], z_approx, g, theta_w[j], 
                                    outer = FALSE, v)$logl
  }
  
  for (i in 1:D) { # separate sampling for each dimension of hidden layer
    
    z_prior <- rand_mvn_vec(x_approx, theta_z[i], v)
    
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
    z_prev <- z_approx$x_ord[z_approx$rev_ord_obs, i] # store for re-proposal
    
    while (accept == FALSE) {
      count <- count + 1
      
      # Calculate proposed values and new likelihood
      z_proposal <- z_prev * cos(a) + z_prior * sin(a)
      
      # Incorporate proposal in vecchia approximation object
      z_approx <- update_obs_in_approx(z_approx, z_proposal, i)

      new_logl <- 0
      for (j in 1:D)
        new_logl <- new_logl + logl_vec(w[, j], z_approx, g, theta_w[j], 
                                        outer = FALSE, v)$logl
      
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
  
  return(list(z_approx = z_approx, ll = ll_prev))
}

