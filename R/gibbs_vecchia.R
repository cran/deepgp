
# Function Contents -----------------------------------------------------------
# Internal: 
#   gibbs_one_layer_vec
#   gibbs_two_layer_vec
#   gibbs_three_layer_vec

# One layer Gibbs with Vecchia ------------------------------------------------

gibbs_one_layer_vec <- function(x, y, nmcmc, verb, initial, true_g, settings, 
                                v, m, x_approx = NULL) {
  
  if (is.null(x_approx)) 
    x_approx <- create_approx(x, m)
  
  g <- vector(length = nmcmc)
  if (is.null(true_g)) g[1] <- initial$g else g[1] <- true_g
  theta <- vector(length = nmcmc)
  theta[1] <- initial$theta
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll <- NULL
  
  for (j in 2:nmcmc) {
    
    if(verb) if(j %% 500 == 0) cat(j, '\n')
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g_vec(y, g[j - 1], theta[j - 1], alpha = settings$alpha$g, 
                           beta = settings$beta$g, l = settings$l, 
                           u = settings$u, ll_prev = ll, approx = x_approx, 
                           v = v)
      g[j] <- samp$g
      ll <- samp$ll
    } else g[j] <- true_g
    
    # Sample length scale (theta)
    samp <- sample_theta_vec(y, g[j], theta[j - 1], 
                             alpha = settings$alpha$theta,
                             beta = settings$beta$theta, l = settings$l, 
                             u = settings$u, outer = TRUE, ll_prev = ll, 
                             approx = x_approx, v = v, tau2 = TRUE)
    theta[j] <- samp$theta
    ll <- samp$ll
    ll_store[j] <- ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
  } # end of j for loop
  
  return(list(g = g, theta = theta, tau2 = tau2, x_approx = x_approx,
              ll = ll_store))
}

# One layer Gibbs with Vecchia SEPARABLE ------------------------------------------------

gibbs_one_layer_vec_sep <- function(x, y, nmcmc, verb, initial, true_g, settings, 
                                v, m, x_approx = NULL) {
  
  d <- ncol(x)
  if (is.null(x_approx)) 
    x_approx <- create_approx(x, m)
  
  g <- vector(length = nmcmc)
  if (is.null(true_g)) g[1] <- initial$g else g[1] <- true_g
  theta <- matrix(nrow = nmcmc, ncol = d)
  if (length(initial$theta) == 1) initial$theta <- rep(initial$theta, d)
  theta[1, ] <- initial$theta
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll <- NULL
  
  for (j in 2:nmcmc) {
    
    if(verb) if(j %% 500 == 0) cat(j, '\n')
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g_vec(y, g[j - 1], theta[j - 1, ], alpha = settings$alpha$g, 
                           beta = settings$beta$g, l = settings$l, 
                           u = settings$u, ll_prev = ll, approx = x_approx, 
                           v = v, sep = TRUE)
      g[j] <- samp$g
      ll <- samp$ll
    } else g[j] <- true_g
    
    # Sample length scale (theta)
    for (i in 1:d) {
      samp <- sample_theta_vec_sep(y, g[j], theta[j - 1, ], index = i,
                                alpha = settings$alpha$theta,
                                beta = settings$beta$theta, l = settings$l, 
                                u = settings$u, ll_prev = ll, 
                                approx = x_approx, v = v, tau2 = TRUE)
      theta[j, i] <- samp$theta
      ll <- samp$ll
      ll_store[j] <- ll
      if (i == 1) { # update tau2 (repeat original value if nothing was accepted)
        if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
      } else { # only update tau2 if there was an acceptance
        if (!is.null(samp$tau2)) tau2[j] <- samp$tau2
      }
    }
  } # end of j for loop
  
  return(list(g = g, theta = theta, tau2 = tau2, x_approx = x_approx,
              ll = ll_store))
}

# Two layer Gibbs with Vecchia ------------------------------------------------

gibbs_two_layer_vec <- function(x, y, nmcmc, D, verb, initial, true_g, settings, 
                                v, m, x_approx = NULL, w_approx = NULL) {
  
  if (is.null(x_approx)) 
    x_approx <- create_approx(x, m)
  if (is.null(w_approx)) 
    w_approx <- create_approx(initial$w, m)
  
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
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll_outer <- NULL
  
  for (j in 2:nmcmc) {
   
    if(verb) if(j %% 500 == 0) cat(j, '\n')
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g_vec(y, g[j - 1], theta_y[j - 1], 
                           alpha = settings$alpha$g, beta = settings$beta$g, 
                           l = settings$l, u = settings$u, ll_prev = ll_outer, 
                           approx = w_approx, v = v)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else g[j] <- true_g
    
    # Sample outer length scale (theta_y)
    samp <- sample_theta_vec(y, g[j], theta_y[j - 1], 
                             alpha = settings$alpha$theta_y, 
                             beta = settings$beta$theta_y, l = settings$l, 
                             u = settings$u, outer = TRUE, ll_prev = ll_outer, 
                             approx = w_approx, v = v, tau2 = TRUE)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
    
    # Sample inner length scale (theta_w) - separately for each dimension
    for (i in 1:D) {
      if (settings$pmx) prior_mean <- x[, i] else prior_mean <- 0
      samp <- sample_theta_vec(w[[j - 1]][, i], g = eps, theta_w[j - 1, i],
                               alpha = settings$alpha$theta_w, 
                               beta = settings$beta$theta_w, l = settings$l, 
                               u = settings$u, outer = FALSE, 
                               approx = x_approx, v = v,
                               prior_mean = prior_mean,
                               scale = settings$inner_tau2)
      theta_w[j, i] <- samp$theta
    }
    
    # Sample hidden Gaussian layer (w)
    if (settings$pmx) prior_mean <- x else prior_mean = matrix(0, nrow(x), D)
    samp <- sample_w_vec(y, w_approx, x_approx, g[j], theta_y[j], theta_w[j, ], 
                         ll_prev = ll_outer, v = v, prior_mean = prior_mean,
                         scale = settings$inner_tau2)
    w_approx <- samp$w_approx
    w[[j]] <- w_approx$x_ord[w_approx$rev_ord_obs, , drop = FALSE]
    ll_outer <- samp$ll
    ll_store[j] <- ll_outer
  } # end of j for loop
  
  return(list(g = g, theta_y = theta_y, theta_w = theta_w, w = w, tau2 = tau2,
              w_approx = w_approx, x_approx = x_approx, ll = ll_store))
}

# Three layer Gibbs with Vecchia ----------------------------------------------

gibbs_three_layer_vec <- function(x, y, nmcmc, D, verb, initial, true_g, 
                                  settings, v, m, x_approx = NULL, 
                                  z_approx = NULL, w_approx = NULL) {
  
  if (is.null(x_approx)) 
    x_approx <- create_approx(x, m)
  if (is.null(z_approx))
    z_approx <- create_approx(initial$z, m)
  if (is.null(w_approx)) 
    w_approx <- create_approx(initial$w, m)
  
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
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll_outer <- NULL
  
  for (j in 2:nmcmc) {
    
    if(verb) if(j %% 500 == 0) cat(j, '\n')
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g_vec(y, g[j - 1], theta_y[j - 1], 
                           alpha = settings$alpha$g, beta = settings$beta$g, 
                           l = settings$l, u = settings$u, ll_prev = ll_outer, 
                           approx = w_approx, v = v)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else g[j] <- true_g
    
    # Sample outer length scale (theta_y)
    samp <- sample_theta_vec(y, g[j], theta_y[j - 1], 
                             alpha = settings$alpha$theta_y,
                             beta = settings$beta$theta_y, l = settings$l, 
                             u = settings$u, outer = TRUE, ll_prev = ll_outer, 
                             approx = w_approx, v = v, tau2 = TRUE)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
    
    # Sample middle length scale (theta_w)
    ll_mid <- 0 # re-calculated each time since we have a new z
    for (i in 1:D) {
      samp <- sample_theta_vec(w[[j - 1]][, i], g = eps, theta_w[j - 1, i],
                               alpha = settings$alpha$theta_w, 
                               beta = settings$beta$theta_w, l = settings$l, 
                               u = settings$u, outer = FALSE, 
                               approx = z_approx, v = v)
      theta_w[j, i] <- samp$theta
      ll_mid <- ll_mid + samp$ll
    }
    
    # Sample inner length scale (theta_z)
    for (i in 1:D) {
      samp <- sample_theta_vec(z[[j - 1]][, i], g = eps, theta_z[j - 1, i], 
                               alpha = settings$alpha$theta_z, 
                               beta = settings$beta$theta_z, l = settings$l, 
                               u = settings$u, outer = FALSE,
                               approx = x_approx, v = v)
      theta_z[j, i] <- samp$theta
    }
    
    # Sample inner hidden Gaussian layer (z)
    samp <- sample_z_vec(w[[j - 1]], z_approx, x_approx, g = eps, theta_w[j, ], 
                         theta_z[j, ], ll_prev = ll_mid, v = v)
    z_approx <- samp$z_approx
    z[[j]] <- z_approx$x_ord[z_approx$rev_ord_obs, , drop = FALSE]
    
    # Sample middle hidden Gaussian layer (w)
    samp <- sample_w_vec(y, w_approx, z_approx, g = g[j], theta_y[j], 
                         theta_w[j, ], ll_prev = ll_outer, v = v) 
    w_approx <- samp$w_approx
    w[[j]] <- w_approx$x_ord[w_approx$rev_ord_obs, , drop = FALSE]
    ll_outer <- samp$ll
    ll_store[j] <- ll_outer
  } # end of j for loop
  
  return(list(g = g, theta_y = theta_y, theta_w = theta_w, theta_z = theta_z,
              w = w, z = z, tau2 = tau2, w_approx = w_approx, z_approx = z_approx, 
              x_approx = x_approx, ll = ll_store))
}
