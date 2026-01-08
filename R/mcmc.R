
# Function Contents -----------------------------------------------------------
# Internal:
#   logl: evaluates MVN log likelihood (sep = TRUE and sep = FALSE options)
#   sample_g: conducts Metropolis Hastings sampling for nugget
#   sample_theta: conducts Metropolis Hastings sampling for theta
#   sample_w: conducts elliptical slice sampling for w layer
#   sample_w_grad: conducts gradient-enhanced elliptical slice sampling for w layer
#   sample_w_mono: conducts monowarped elliptical slice sampling for w layer
#   sample_z: conducts elliptical slice sampling for z layer

# logl ------------------------------------------------------------------------

logl <- function(y, xdmat = NULL, # used for sep = FALSE
                 x = NULL, # used for sep = TRUE and/or grad_enhance
                 tau2 = 1, theta, g = 0, v, mu = 0,
                 grad_enhance = FALSE,
                 sep = FALSE, outer = FALSE) {

  if (grad_enhance) {
    if (!is.matrix(x)) x <- as.matrix(x)
    d <- ncol(x)
    grad_indx <- rep(0:d, each = nrow(x)) # always the same gradient ordering
    x <- bind(x, d) # duplicate observations
    if (length(y) != nrow(x)) stop("dimension mismatch")
  }

  n <- length(y)
  if (outer & tau2 != 1) stop("outer = TRUE should always use tau2 = 1")
  if (v == 999) {
    if (grad_enhance) {
      if (sep) {
        K <- Exp2SepGrad(x, x, grad_indx, grad_indx, tau2, theta, g)
      } else K <- Exp2Grad(x, x, grad_indx, grad_indx, tau2, theta, g)
    } else {
      if (sep) {
        K <- Exp2Sep(x, x, tau2, theta, g)
      } else K <- Exp2(xdmat, tau2, theta, g)
    } 
  } else {
    if (sep) {
      K <- MaternSep(x, x, tau2, theta, g, v) 
    } else K <- Matern(xdmat, tau2, theta, g, v) 
  }
  id <- invdet(K)
  quadterm <- drop(t(y - mu) %*% id$Mi %*% (y - mu))
  
  if (outer) { # use profile log likelihood (with tau2 integrated out)
    ll <- (-n*0.5)*log(quadterm) - 0.5*id$ldet
  } else ll <- (-0.5)*id$ldet - 0.5*quadterm
  
  return(list(ll = ll, tau2 = quadterm/n))
}

# sample_g --------------------------------------------------------------------
# Outer layer only - always tau2 = 1, prior_mean = 0, outer = TRUE 
# Note: any proposals below eps will be REJECTED

sample_g <- function(y, xdmat = NULL, # used for sep = FALSE
                     x = NULL, # used for sep = TRUE or grad_enhance
                     theta, g, v, 
                     alpha, beta, l, u, 
                     ll_prev = NULL,
                     grad_enhance = FALSE,
                     sep = FALSE) {
  
  if (is.null(ll_prev)) {
    ll_prev <- logl(y, xdmat = xdmat, x = x, tau2 = 1, theta = theta, g = g, 
                    v = v, grad_enhance = grad_enhance, sep = sep, outer = TRUE)$ll
  }
  
  # Propose value and compute acceptance threshold
  g_star <- runif(1, min = l*g/u, max = u*g/l)
  ru <- runif(1, min = 0, max = 1)
  lpost_threshold <- ll_prev + dgamma(g - eps, alpha, beta, log = TRUE) + 
                        log(ru) - log(g) + log(g_star)
  
  # Calculate new likelihood
  ll_new <- logl(y, xdmat = xdmat, x = x, tau2 = 1, theta = theta, g = g_star, 
                 v = v, grad_enhance = grad_enhance, sep = sep, outer = TRUE)
  
  # Accept or reject (lower bound of eps)
  new <- ll_new$ll + dgamma(g_star - eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) { # accept
    return(list(g = g_star, ll = ll_new$ll, tau2 = ll_new$tau2))
  } else { # reject
    return(list(g = g, ll = ll_prev, tau2 = NULL))
  }
}

# sample_theta ----------------------------------------------------------------
# Note: any proposals below eps will be REJECTED

sample_theta <- function(y, xdmat = NULL, # used for sep = FALSE
                         x = NULL, # used for sep = TRUE and/or grad_enhance
                         tau2, theta, g, v,
                         alpha, beta, l, u, 
                         outer, ll_prev = NULL, prior_mean = 0,
                         grad_enhance = FALSE,
                         sep = FALSE, index = 1) {

  if (!sep) index <- 1 # force index 1 since there is only one theta
  
  if (is.null(ll_prev)) {
    ll_prev <- logl(y, xdmat = xdmat, x = x, tau2 = tau2, theta = theta, g = g, 
                    v = v, mu = prior_mean, grad_enhance = grad_enhance, sep = sep, 
                    outer = outer)$ll
  }

  # Propose value and compute acceptance threshold
  ru <- runif(1, min = 0, max = 1)
  theta_star <- theta
  theta_star[index] <- runif(1, min = l*theta[index]/u, max = u*theta[index]/l)
  
  lpost_threshold <- ll_prev + dgamma(theta[index] - eps, alpha, beta, log = TRUE) + 
                        log(ru) - log(theta[index]) + log(theta_star[index])
  
  # Calculate new likelihood
  ll_new <- logl(y, xdmat = xdmat, x = x, tau2 = tau2, theta = theta_star, g = g, 
                 v = v, mu = prior_mean, sep = sep, grad_enhance = grad_enhance,
                 outer = outer)

  # Accept or reject (lower bound of eps)
  new <- ll_new$ll + dgamma(theta_star[index] - eps, alpha, beta, log = TRUE)
  if (new > lpost_threshold) { # accept
    return(list(theta = theta_star[index], ll = ll_new$ll, tau2 = ll_new$tau2))
  } else { # reject
    return(list(theta = theta[index], ll = ll_prev, tau2 = NULL))
  }
}

# sample_w --------------------------------------------------------------------

sample_w <- function(y, w, xdmat, tau2_w, theta_y, theta_w, g, v, 
                     ll_prev, prior_mean = NULL) {
  
  if (!is.matrix(w)) w <- as.matrix(w)
  n <- length(y)
  D <- ncol(w) # dimension of hidden layer
  if (length(tau2_w) == 1) tau2_w <- rep(tau2_w, times = D)
  if (length(theta_w) == 1) theta_w <- rep(theta_w, times = D)
  
  for (i in 1:D) { # separate sampling for each dimension of hidden layer
    # Check prior_mean
    if (is.null(prior_mean)) {
      pm <- rep(0, n)
    } else pm <- prior_mean[, i]
    
    # Draw from prior distribution 
    if (v == 999) {
      sigma <- Exp2(xdmat, tau2 = tau2_w[i], theta = theta_w[i], g = eps)
    } else sigma <- Matern(xdmat, tau2 = tau2_w[i], theta = theta_w[i], g = eps, v = v)
    w_prior <- mvtnorm::rmvnorm(1, mean = pm, sigma = sigma)
    
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
      w[, i] <- w_prev*cos(a) + w_prior*sin(a)
      wdmat <- sq_dist(w)
      ll_new <- logl(y, xdmat = wdmat, tau2 = 1, theta = theta_y, g = g, v = v, 
                     outer = TRUE) 
      
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
  
  return(list(w = w, wdmat = wdmat, ll = ll_new$ll, tau2_y = ll_new$tau2))
}

# sample_w_grad ---------------------------------------------------------------

sample_w_grad <- function(y, dydx, w, x, tau2_w, theta_y, theta_w, g, v, 
                          ll_prev, prior_mean = NULL) { 
  
  if (!is.matrix(w)) w <- as.matrix(w)
  n <- nrow(x)
  D <- ncol(x) # dimension of x and w (forced to match)
  if (length(tau2_w) == 1) tau2_w <- rep(tau2_w, times = D)
  if (length(theta_w) == 1) theta_w <- rep(theta_w, times = D)
  
  x <- bind(x, D)
  grad_indx <- rep(0:D, each = n)
  
  for (i in 1:D) { # separate sampling for each dimension of hidden layer
    # Check prior_mean
    if (is.null(prior_mean)) {
      pm <- rep(0, n*(D+1))
    } else pm <- prior_mean[, i]

    # Draw from prior distribution (v = 999 only, includes gradients)
    sigma <- Exp2Grad(x, x, grad_indx, grad_indx, tau2 = tau2_w[i], 
                      theta = theta_w[i], g = eps)
    w_prior <- drop(mvtnorm::rmvnorm(1, mean = pm, sigma = sigma))
    
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
      w[, i] <- w_prev*cos(a) + w_prior*sin(a) # includes gradient
      dydw <- get_dydw(w, dydx)
      y_all <- c(y, as.vector(dydw))

      ll_new <- logl(y_all, x = w[1:n, ], tau2 = 1, theta = theta_y, g = g, v = v, 
                     grad_enhance = TRUE, outer = TRUE)      
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
        if (count > 300) stop("reached maximum iterations of ESS")
      } # end of else statement
    } # end of while loop
  } # end of i for loop
  
  return(list(y_all = y_all, w = w, ll = ll_new$ll, tau2_y = ll_new$tau2))
}

# sample_w_mono ---------------------------------------------------------------

sample_w_mono <- function(y, w, x, x_grid, w_grid, xdmat_grid, 
                          tau2_w, theta_y, theta_w, g, v,
                          ll_prev, prior_mean = NULL) { 
  
  if (!is.matrix(w)) w <- as.matrix(w)
  if (!is.matrix(w_grid)) w_grid <- as.matrix(w_grid)
  n <- length(y)
  D <- ncol(w) # dimension of x and hidden layer (forced to match)
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
    w_grid_prev <- w_grid[, i] # store for re-proposal
    
    while (accept == FALSE) {
      count <- count + 1
      
      # Calculate proposed values and new likelihood
      w_grid[, i] <- w_grid_prev*cos(a) + w_grid_prior*sin(a)
      w[, i] <- monotransform(x[, i], x_grid, w_grid[, i])
      wdmat <- sq_dist(w)
      ll_new <- logl(y, xdmat = wdmat, tau2 = 1, theta = theta_y, g = g, v = v,
                     outer = TRUE)
      
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
  
  return(list(w = w, wdmat = wdmat, w_grid = w_grid, 
              ll = ll_new$ll, tau2_y = ll_new$tau2))
}

# sample_z --------------------------------------------------------------------
# Always prior_mean = 0
# Likelihood calculation involves w which is always noise free

sample_z <- function(w, z, xdmat, tau2_w, tau2_z, theta_w, theta_z, v, ll_prev) {
  
  if (!is.matrix(w)) w <- as.matrix(w)
  if (!is.matrix(z)) z <- as.matrix(z)
  D <- ncol(z) # dimension of hidden layer
  if (length(tau2_w) == 1) tau2_w <- rep(tau2_w, times = D)
  if (length(theta_w) == 1) theta_w <- rep(theta_w, times = D)
  if (length(tau2_z) == 1) tau2_z <- rep(tau2_z, times = D)
  if (length(theta_z) == 1) theta_z <- rep(theta_z, times = D)
  
  for (i in 1:D) { # separate sampling for each dimension of hidden layer
    
    # Draw from prior distribution
    if (v == 999) {
      sigma <- Exp2(xdmat, tau2 = tau2_z[i], theta = theta_z[i], g = eps)
    } else sigma <- Matern(xdmat, tau2 = tau2_z[i], theta = theta_z[i], g = eps, v = v)
    z_prior <- mvtnorm::rmvnorm(1, sigma = sigma)
    
    # Initialize a and bounds on a
    a <- runif(1, min = 0, max = 2*pi)
    amin <- a - 2*pi
    amax <- a
    
    # Compute acceptance threshold - based on all dimensions of previous z
    ru <- runif(1, min = 0, max = 1)
    ll_threshold <- ll_prev + log(ru)
    
    # Calculate proposed values, accept or reject, repeat if necessary
    accept <- FALSE
    count <- 0
    z_prev <- z[, i] # store for re-proposal
    
    while (accept == FALSE) {
      count <- count + 1
      
      # Calculate proposed values and new likelihood
      z[, i] <- z_prev*cos(a) + z_prior*sin(a)
      zdmat <- sq_dist(z)
      ll_new <- 0
      for (j in 1:D) {
        ll_new <- ll_new + logl(w[, j], zdmat, tau2 = tau2_w[j], theta = theta_w[j], 
                                g = eps, v = v, sep = FALSE, outer = FALSE)$ll
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
  
  return(list(z = z, zdmat = zdmat, ll = ll_prev))
}
