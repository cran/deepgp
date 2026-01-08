
# Function Contents -----------------------------------------------------------
# Internal:
#   logl_vec: calculates vecchia log likelihood
#   sample_g_vec: conducts Metropolis Hastings sampling for nugget
#   sample_theta_vec: conducts Metropolis Hastings sampling for theta
#   sample_w_vec: conducts elliptical slice ampling for w layer
#   sample_w_vec_grad: conducts gradient-enhanced elliptical slice sampling for w layer
#   sample_w_vec_mono: conducts monowarped elliptical slice sampling for w layer
#   sample_z_vec: conducts Elliptical Slice Sampling for z layer

# logl_vec --------------------------------------------------------------------

logl_vec <- function(y, approx, tau2 = 1, theta, g = 0, v, 
                     mu = 0, sep = FALSE, outer = FALSE) {

  # We do not need a grad_enhance argument, since the grad options are built
  # into approx (and create_U can handle both)

  n <- length(y)
  if (length(mu) > 1) mu_ordered <- mu[approx$ord] else mu_ordered <- mu
  y_ord <- y[approx$ord]

  U_mat <- create_U(approx, tau2 = tau2, theta = theta, g = g, v = v, sep = sep) 
  Uty <- Matrix::crossprod(U_mat, y_ord - mu_ordered)
  quadterm <- sum(Uty^2)
  logdet <- sum(log(Matrix::diag(U_mat)))
  
  if (outer) { # use profile log likelihood (with tau2 integrated out)
    ll <- logdet - (n*0.5)*log(quadterm)
  } else ll <- logdet - 0.5*quadterm
  
  return(list(ll = ll, tau2 = quadterm/n))
}

# sample_g_vec ----------------------------------------------------------------
# Outer layer only - always tau2 = 1, prior_mean = 0, outer = TRUE
# Note: any proposals below eps will be REJECTED

sample_g_vec <- function(y, approx, 
                         theta, g, v, 
                         alpha, beta, l, u, 
                         ll_prev = NULL,
                         sep = FALSE) {

  if (is.null(ll_prev)) {
    ll_prev <- logl_vec(y, approx = approx, tau2 = 1, theta = theta, g = g, v = v, 
                        sep = sep, outer = TRUE)$ll
  }

  # Propose value and compute acceptance threshold
  g_star <- runif(1, min = l*g/u, max = u*g/l)
  ru <- runif(1, min = 0, max = 1)  
  lpost_threshold <- ll_prev + dgamma(g - eps, alpha, beta, log = TRUE) + 
                        log(ru) - log(g) + log(g_star)
 
  # Calculate new likelihood
  ll_new <- logl_vec(y, approx = approx, tau2 = 1, theta = theta, g = g_star, v = v, 
                     sep = sep, outer = TRUE)
  
  # Accept or reject (lower bound of eps)
  new <- ll_new$ll + dgamma(g_star - eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) { # accept
    return(list(g = g_star, ll = ll_new$ll, tau2 = ll_new$tau2))
  } else { # reject
    return(list(g = g, ll = ll_prev, tau2 = NULL))
  }
}

# sample_theta_vec ------------------------------------------------------------
# Note: any proposals below eps will be REJECTED

sample_theta_vec <- function(y, approx, 
                             tau2, theta, g, v, 
                             alpha, beta, l, u, 
                             outer, ll_prev = NULL, prior_mean = 0,
                             sep = FALSE, index = 1) {

  if (!sep) index <- 1 # force index 1 since there is only one theta
  
  if (is.null(ll_prev)) {
    ll_prev <- logl_vec(y, approx = approx, tau2 = tau2, theta = theta, g = g, 
                        v = v, mu = prior_mean, sep = sep, outer = outer)$ll
  }
  
  # Propose value and compute acceptance threshold
  ru <- runif(1, min = 0, max = 1)
  theta_star <- theta
  theta_star[index] <- runif(1, min = l*theta[index]/u, max = u*theta[index]/l)

  lpost_threshold <- ll_prev + dgamma(theta[index] - eps, alpha, beta, log = TRUE) + 
                      log(ru) - log(theta[index]) + log(theta_star[index])
  
  # Calculate new likelihood
  ll_new <- logl_vec(y, approx = approx, tau2 = tau2, theta = theta_star, g = g, 
                     v = v, mu = prior_mean, sep = sep, outer = outer)

  # Accept or reject (lower bound of eps)
  new <- ll_new$ll + dgamma(theta_star[index] - eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) { # accept
    return(list(theta = theta_star[index], ll = ll_new$ll, tau2 = ll_new$tau2))
  } else { # reject
    return(list(theta = theta[index], ll = ll_prev, tau2 = NULL))
  }
}

# sample_w_vec ----------------------------------------------------------------

sample_w_vec <- function(y, w_approx, x_approx, tau2_w, theta_y, theta_w, g, v, 
                         ll_prev, prior_mean = NULL) {

  D <- ncol(w_approx$x_ord) # dimension of hidden layer
  if (length(tau2_w) == 1) tau2_w <- rep(tau2_w, times = D)
  if (length(theta_w) == 1) theta_w <- rep(theta_w, times = D)

  for (i in 1:D) { # separate sampling for each dimension of hidden layer
    # Check prior_mean
    if (is.null(prior_mean)) {
      pm <- rep(0, nrow(x_approx$x_ord))
    } else pm <- prior_mean[, i] # NOT ordered
    
    # Draw from prior distribution
    w_prior <- rand_mvn_vec(x_approx, tau2 = tau2_w[i], theta = theta_w[i], 
                            g = eps, v = v, prior_mean = pm)
    
    # Initialize a and bounds on a
    a <- runif(1, min = 0, max = 2*pi)
    amin <- a - 2*pi
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
      w_proposal <- w_prev*cos(a) + w_prior*sin(a)
      w_approx$x_ord[, i] <- w_proposal[w_approx$ord]

      ll_new <- logl_vec(y, approx = w_approx, tau2 = 1, theta = theta_y, g = g, 
                         v = v, outer = TRUE)

      # Accept or reject
      if (ll_new$ll > ll_threshold) {
        ll_prev <- ll_new$ll
        accept <- TRUE
      } else { 
        # update the bounds on a and repeat
        if (a < 0) {
          amin <- a
        } else {
          amax <- a
        }
        a <- runif(1, amin, amax)
        if (count > 100) stop("reached maximum iterations of ESS")
      } # end of else statement
    } # end of while loop
  } # end of i for loop

  return(list(w_approx = w_approx, ll = ll_new$ll, tau2_y = ll_new$tau2))
}

# sample_w_vec_grad -----------------------------------------------------------

sample_w_vec_grad <- function(y, w, dydx, w_approx, x_approx, tau2_w, theta_y, 
                              theta_w, g, v, ll_prev, prior_mean = NULL) {

  n <- length(y)
  D <- ncol(x_approx$x_ord) # dimension of x and w (forced to match)
  if (!is.matrix(w)) w <- as.matrix(w)
  if (length(tau2_w) == 1) tau2_w <- rep(tau2_w, times = D)
  if (length(theta_w) == 1) theta_w <- rep(theta_w, times = D)

  for (i in 1:D) { # separate sampling for each dimension of hidden layer
    # Check prior_mean
    if (is.null(prior_mean)) {
      pm <- rep(0, n*(D+1))
    } else pm <- prior_mean[, i]

    # Draw from prior distribution (includes gradients)
    w_prior <- rand_mvn_vec(x_approx, tau2 = tau2_w[i], theta = theta_w[i], 
                            g = eps, v = v, grad = TRUE, prior_mean = pm)
    
    # Initialize a and bounds on a
    a <- runif(1, min = 0, max = 2*pi)
    amin <- a - 2*pi
    amax <- a

    # Compute acceptance threshold - based on all dimensions of previous w
    ru <- runif(1, min = 0, max = 1)
    ll_threshold <- ll_prev + log(ru)

    # Calculate proposed values, accept or reject, repeat if necessary
    accept <- FALSE
    count <- 0
    w_prev <- w[, i] # store for re-proposal
    
    while (accept == FALSE) {
      count <- count + 1

      # Calculate proposed values and new likelihood
      w_proposal <- w_prev*cos(a) + w_prior*sin(a)
      w[, i] <- w_proposal
      w_approx$x_ord[, i] <- bind(w_proposal[w_approx$ord[1:n]], D)
      dydw <- get_dydw(w, dydx)
      y_all <- c(y, as.vector(dydw))
 
      ll_new <- logl_vec(y_all, approx = w_approx, tau2 = 1, theta = theta_y, g = g, 
                         v = v, outer = TRUE)

      # Accept or reject
      if (ll_new$ll > ll_threshold) {
        ll_prev <- ll_new$ll
        accept <- TRUE
      } else { 
        # update the bounds on a and repeat
        if (a < 0) {
          amin <- a
        } else {
          amax <- a
        }
        a <- runif(1, amin, amax)
        if (count > 100) stop("reached maximum iterations of ESS")
      } # end of else statement
    } # end of while loop
  } # end of i for loop

  return(list(y_all = y_all, w = w, w_approx = w_approx, ll = ll_new$ll, tau2_y = ll_new$tau2))
}

# sample_w_vec_mono -----------------------------------------------------------

sample_w_vec_mono <- function(y, w_approx, x, x_grid, w_grid, xdmat_grid,
                              tau2_w, theta_y, theta_w, g, v, 
                              ll_prev, prior_mean = NULL) {

  if (!is.matrix(w_grid)) w_grid <- as.matrix(w_grid)
  n <- length(y)
  D <- ncol(w_approx$x_ord) # dimension of x and hidden layer (forced to match)
  if (length(tau2_w) == 1) tau2_w <- rep(tau2_w, times = D)
  if (length(theta_w) == 1) theta_w <- rep(theta_w, times = D)

  for (i in 1:D) { # separate sampling for each dimension of hidden layer
    
    # Draw from prior distribution
    if (v == 999) {
      sigma <- Exp2(xdmat_grid, tau2 = tau2_w[i], theta = theta_w[i], g = eps)
    } else sigma <- Matern(xdmat_grid, tau2 = tau2_w[i], theta = theta_w[i], g = eps, v = v)
    w_grid_prior <- mvtnorm::rmvnorm(1, mean = prior_mean, sigma = sigma)
    
    # Initialize a and bounds on a
    a <- runif(1, min = 0, max = 2*pi)
    amin <- a - 2*pi
    amax <- a

    # Compute acceptance threshold - based on all dimensions of previous w
    ru <- runif(1, min = 0, max = 1)
    ll_threshold <- ll_prev + log(ru)

    # Calculate proposed values, accept or reject, repeat if necessary
    accept <- FALSE
    count <- 0
    w_grid_prev <- w_grid[, i]

    while (accept == FALSE) {
      count <- count + 1

      # Calculate proposed values and new likelihood
      w_grid[, i] <- w_grid_prev*cos(a) + w_grid_prior*sin(a)
      w_proposal <- monotransform(x[, i], x_grid, w_grid[, i])
      w_approx$x_ord[, i] <- w_proposal[w_approx$ord]
      
      ll_new <- logl_vec(y, approx = w_approx, tau2 = 1, theta = theta_y, g = g, 
                         v = v, outer = TRUE)

      # Accept or reject
      if (ll_new$ll > ll_threshold) {
        ll_prev <- ll_new$ll
        accept <- TRUE
      } else { 
        # update the bounds on a and repeat
        if (a < 0) {
          amin <- a
        } else {
          amax <- a
        }
        a <- runif(1, amin, amax)
        if (count > 100) stop("reached maximum iterations of ESS")
      } # end of else statement
    } # end of while loop
  } # end of i for loop

  return(list(w_approx = w_approx, w_grid = w_grid, ll = ll_new$ll, tau2_y = ll_new$tau2))
}

# sample_z_vec ----------------------------------------------------------------
# Always prior_mean = 0
# Likelihood calculation involves w which is always noise free

sample_z_vec <- function(w, z_approx, x_approx, tau2_w, tau2_z, theta_w, theta_z, v, 
                         ll_prev) {
  
  if (!is.matrix(w)) w <- as.matrix(w)
  D <- ncol(z_approx$x_ord) # dimension of hidden layer
  if (length(tau2_w) == 1) tau2_w <- rep(tau2_w, times = D)
  if (length(theta_w) == 1) theta_w <- rep(theta_w, times = D)
  if (length(tau2_z) == 1) tau2_z <- rep(tau2_z, times = D)
  if (length(theta_z) == 1) theta_z <- rep(theta_z, times = D)
  
  for (i in 1:D) { # separate sampling for each dimension of hidden layer
    
    z_prior <- rand_mvn_vec(x_approx, tau2 = tau2_z[i], theta = theta_z[i], g = eps, v = v)
    
    # Initialize a and bounds on a
    a <- runif(1, min = 0, max = 2*pi)
    amin <- a - 2*pi
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
      z_proposal <- z_prev*cos(a) + z_prior*sin(a)
      z_approx$x_ord[, i] <- z_proposal[z_approx$ord]

      ll_new <- 0
      for (j in 1:D) {
        ll_new <- ll_new + logl_vec(w[, j], approx = z_approx, tau2 = tau2_w[j], theta = theta_w[j], 
                                    g = eps, v = v, outer = FALSE)$ll
      }
      
      # Accept or reject
      if (ll_new > ll_threshold) {
        ll_prev <- ll_new
        accept <- TRUE
      } else {
        # update the bounds on a and repeat
        if (a < 0) {
          amin <- a
        } else {
          amax <- a
        }
        a <- runif(1, amin, amax)
        if (count > 100) stop("reached maximum iterations of ESS")
      } # end of else statement
    } # end of while loop
  } # end of i for loop
  
  return(list(z_approx = z_approx, ll = ll_prev))
}

