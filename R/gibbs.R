
# Function Contents -----------------------------------------------------------
# Internal: 
#   gibbs_one_layer
#   gibbs_two_layer
#   gibbs_three_layer

# One layer Gibbs -------------------------------------------------------------

gibbs_one_layer <- function(x, y, nmcmc, verb, initial, true_g, settings, v) {
  
  dx <- sq_dist(x)
  g <- vector(length = nmcmc)
  if (is.null(true_g)) g[1] <- initial$g else g[1] <- true_g
  theta <- vector(length = nmcmc)
  theta[1] <- initial$theta
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll <- NULL
  
  for (j in 2:nmcmc) {
    
    if(verb) if(j %% 500 == 0) cat(j, '\n')
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g(y, dx, g[j - 1], theta[j - 1], alpha = settings$alpha$g, 
                       beta = settings$beta$g, l = settings$l, u = settings$u, 
                       ll_prev = ll, v = v)
      g[j] <- samp$g
      ll <- samp$ll
    } else g[j] <- true_g
    
    # Sample length scale (theta)
    samp <- sample_theta(y, dx, g[j], theta[j - 1], 
                         alpha = settings$alpha$theta,
                         beta = settings$beta$theta, l = settings$l, 
                         u = settings$u, outer = TRUE, ll_prev = ll, v = v, 
                         tau2 = TRUE)
    theta[j] <- samp$theta
    ll <- samp$ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
  } # end of j for loop
  
  return(list(g = g, theta = theta, tau2 = tau2))
}

# One layer Gibbs SEPARABLE ---------------------------------------------------

gibbs_one_layer_sep <- function(x, y, nmcmc, verb, initial, true_g, settings, v) {
  
  d <- ncol(x)
  g <- vector(length = nmcmc)
  if (is.null(true_g)) g[1] <- initial$g else g[1] <- true_g
  theta <- matrix(nrow = nmcmc, ncol = d)
  if (length(initial$theta) == 1) initial$theta <- rep(initial$theta, d)
  theta[1, ] <- initial$theta
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll <- NULL
  
  for (j in 2:nmcmc) {
    
    if(verb) if(j %% 500 == 0) cat(j, '\n')
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g_sep(y, x, g[j - 1], theta[j - 1, ], alpha = settings$alpha$g, 
                       beta = settings$beta$g, l = settings$l, u = settings$u, 
                       ll_prev = ll, v = v)
      g[j] <- samp$g
      ll <- samp$ll
    } else g[j] <- true_g
    
    # Sample length scale (theta)
    for (i in 1:d) {
      samp <- sample_theta_sep(y, x, g[j], theta[j - 1, ], index = i,
                               alpha = settings$alpha$theta,
                               beta = settings$beta$theta, l = settings$l, 
                               u = settings$u, ll_prev = ll, v = v,
                               tau2 = (i == d))
      theta[j, i] <- samp$theta
      ll <- samp$ll
      if (i == 1) { # update tau2 (repeat original value if nothing was accepted)
        if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
      } else { # only update tau2 if there was an acceptance
        if (!is.null(samp$tau2)) tau2[j] <- samp$tau2
      }
    }
  } # end of j for loop
  
  return(list(g = g, theta = theta, tau2 = tau2))
}

# Two layer Gibbs -------------------------------------------------------------

gibbs_two_layer <- function(x, y, nmcmc, D, verb, initial, true_g, 
                            settings, v) {
  
  dx <- sq_dist(x)
  dw <- sq_dist(initial$w)
  
  g <- vector(length = nmcmc)
  if (is.null(true_g)) g[1] <- initial$g else g[1] <- true_g
  theta_y <- vector(length = nmcmc)
  theta_y[1] <- initial$theta_y
  theta_w <- matrix(nrow = nmcmc, ncol = D)
  theta_w[1, ] <- initial$theta_w
  w <- list()
  w[[1]] <- initial$w
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll_outer <- NULL
  
  for (j in 2:nmcmc) {
    
    if(verb) if(j %% 500 == 0) cat(j, '\n')
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g(y, dw, g[j - 1], theta_y[j - 1], 
                       alpha = settings$alpha$g, beta = settings$beta$g, 
                       l = settings$l, u = settings$u, ll_prev = ll_outer, 
                       v = v)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else g[j] <- true_g
    
    # Sample outer length scale (theta_y)
    samp <- sample_theta(y, dw, g[j], theta_y[j - 1], 
                         alpha = settings$alpha$theta_y, 
                         beta = settings$beta$theta_y, l = settings$l, 
                         u = settings$u, outer = TRUE, ll_prev = ll_outer, 
                         v = v, tau2 = TRUE)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
    
    # Sample inner length scale (theta_w) - separately for each dimension
    for (i in 1:D) {
      samp <- sample_theta(w[[j - 1]][, i], dx, g = eps, theta_w[j - 1, i],
                           alpha = settings$alpha$theta_w, 
                           beta = settings$beta$theta_w, l = settings$l, 
                           u = settings$u, outer = FALSE, v = v)
      theta_w[j, i] <- samp$theta
    }
    
    # Sample hidden Gaussian layer (w)
    samp <- sample_w(y, w[[j - 1]], dw, dx, g[j], theta_y[j], theta_w[j, ], 
                     ll_prev = ll_outer, v = v)
    w[[j]] <- samp$w
    ll_outer <- samp$ll
    dw <- samp$dw
  } # end of j for loop
  
  return(list(g = g, theta_y = theta_y, theta_w = theta_w, w = w, tau2 = tau2))
}

# Three layer Gibbs -----------------------------------------------------------

gibbs_three_layer <- function(x, y, nmcmc, D, verb, initial, true_g, 
                              settings, v) {
  
  dx <- sq_dist(x)
  dz <- sq_dist(initial$z)
  dw <- sq_dist(initial$w)
  
  g <- vector(length = nmcmc)
  if (is.null(true_g)) g[1] <- initial$g else g[1] <- true_g
  theta_y <- vector(length = nmcmc)
  theta_y[1] <- initial$theta_y
  theta_w <- matrix(nrow = nmcmc, ncol = D)
  theta_w[1, ] <- initial$theta_w
  theta_z <- matrix(nrow = nmcmc, ncol = D)
  theta_z[1, ] <- initial$theta_z
  w <- list()
  w[[1]] <- initial$w
  z <- list()
  z[[1]] <- initial$z
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll_outer <- NULL
  
  for (j in 2:nmcmc) {
    
    if(verb) if(j %% 500 == 0) cat(j, '\n')
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g(y, dw, g[j - 1], theta_y[j - 1], 
                       alpha = settings$alpha$g, beta = settings$beta$g, 
                       l = settings$l, u = settings$u, ll_prev = ll_outer, 
                       v = v)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else g[j] <- true_g
    
    # Sample outer length scale (theta_y)
    samp <- sample_theta(y, dw, g[j], theta_y[j - 1], 
                         alpha = settings$alpha$theta_y,
                         beta = settings$beta$theta_y, l = settings$l, 
                         u = settings$u, outer = TRUE, ll_prev = ll_outer, 
                         v = v, tau2 = TRUE)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
    
    # Sample middle length scale (theta_w)
    ll_mid <- 0 # re-calculated each time since we have a new z
    for (i in 1:D) {
      samp <- sample_theta(w[[j - 1]][, i], dz, g = eps, theta_w[j - 1, i],
                           alpha = settings$alpha$theta_w, 
                           beta = settings$beta$theta_w, l = settings$l, 
                           u = settings$u, outer = FALSE, v = v)
      theta_w[j, i] <- samp$theta
      ll_mid <- ll_mid + samp$ll
    }
    
    # Sample inner length scale (theta_z)
    for (i in 1:D) {
      samp <- sample_theta(z[[j - 1]][, i], dx, g = eps, theta_z[j - 1, i], 
                           alpha = settings$alpha$theta_z, 
                           beta = settings$beta$theta_z, l = settings$l, 
                           u = settings$u, outer = FALSE, v = v)
      theta_z[j, i] <- samp$theta
    }
    
    # Sample inner hidden Gaussian layer (z)
    samp <- sample_z(w[[j - 1]], z[[j - 1]], dz, dx, g = eps, theta_w[j, ], 
                     theta_z[j, ], ll_prev = ll_mid, v = v)
    z[[j]] <- samp$z
    dz <- samp$dz
    
    # Sample middle hidden Gaussian layer (w)
    samp <- sample_w(y, w[[j - 1]], dw, dz, g[j], theta_y[j], theta_w[j, ], 
                     ll_prev = ll_outer, v = v)
    w[[j]] <- samp$w
    ll_outer <- samp$ll
    dw <- samp$dw
  } # end of j for loop
  
  return(list(g = g, theta_y = theta_y, theta_w = theta_w, theta_z = theta_z,
              w = w, z = z, tau2 = tau2))
}
