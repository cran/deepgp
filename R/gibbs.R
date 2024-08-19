
# Function Contents -----------------------------------------------------------
# Internal: 
#   gibbs_one_layer
#   gibbs_one_layer_sep
#   gibbs_two_layer
#   gibbs_two_layer_mono
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
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll <- NULL
  
  for (j in 2:nmcmc) {
    
    if (verb & (j %% 500 == 0)) cat(j, '\n')
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g(y, dx, g[j - 1], theta[j - 1], alpha = settings$alpha$g, 
                       beta = settings$beta$g, l = settings$l, u = settings$u, 
                       ll_prev = ll, v = v)
      g[j] <- samp$g
      ll <- samp$ll
    } else g[j] <- true_g
    
    # Sample lengthscale (theta)
    samp <- sample_theta(y, dx, g[j], theta[j - 1], 
                         alpha = settings$alpha$theta,
                         beta = settings$beta$theta, l = settings$l, 
                         u = settings$u, outer = TRUE, ll_prev = ll, v = v, 
                         calc_tau2 = TRUE)
    theta[j] <- samp$theta
    ll <- samp$ll
    ll_store[j] <- ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
  } # end of j for loop
  
  return(list(g = g, theta = theta, tau2 = tau2, ll = ll_store))
}

# One layer Gibbs SEPARABLE ---------------------------------------------------

gibbs_one_layer_sep <- function(x, y, nmcmc, verb, initial, true_g, settings, v) {
  
  d <- ncol(x)
  g <- vector(length = nmcmc)
  if (is.null(true_g)) g[1] <- initial$g else g[1] <- true_g
  theta <- matrix(nrow = nmcmc, ncol = d)
  theta[1, ] <- initial$theta
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll <- NULL
  
  for (j in 2:nmcmc) {
    
    if (verb & (j %% 500 == 0)) cat(j, '\n')
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g_sep(y, x, g[j - 1], theta[j - 1, ], alpha = settings$alpha$g, 
                       beta = settings$beta$g, l = settings$l, u = settings$u, 
                       ll_prev = ll, v = v)
      g[j] <- samp$g
      ll <- samp$ll
    } else g[j] <- true_g
    
    # Sample lengthscale (theta)
    theta_curr <- theta[j - 1, ]
    tau2[j] <- tau2[j - 1] # start repeating tau2 in case there are no acceptances
    for (i in 1:d) {
      samp <- sample_theta_sep(y, x, g[j], theta_curr, index = i,
                               alpha = settings$alpha$theta,
                               beta = settings$beta$theta, l = settings$l, 
                               u = settings$u, ll_prev = ll, v = v,
                               calc_tau2 = TRUE)
      theta_curr[i] <- samp$theta
      theta[j, i] <- samp$theta
      ll <- samp$ll
      ll_store[j] <- ll
      if (!is.null(samp$tau2)) tau2[j] <- samp$tau2
    }
  } # end of j for loop
  
  return(list(g = g, theta = theta, tau2 = tau2, ll = ll_store))
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
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll_outer <- NULL
  
  for (j in 2:nmcmc) {
    
    if (verb & (j %% 500 == 0)) cat(j, '\n')
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g(y, dw, g[j - 1], theta_y[j - 1], 
                       alpha = settings$alpha$g, beta = settings$beta$g, 
                       l = settings$l, u = settings$u, ll_prev = ll_outer, 
                       v = v)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else g[j] <- true_g
    
    # Sample outer lengthscale (theta_y) 
    samp <- sample_theta(y, dw, g[j], theta_y[j - 1], 
                         alpha = settings$alpha$theta_y, 
                         beta = settings$beta$theta_y, l = settings$l, 
                         u = settings$u, outer = TRUE, ll_prev = ll_outer, 
                         v = v, calc_tau2 = TRUE)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
    
    # Sample inner lengthscale (theta_w) - separately for each dimension
    for (i in 1:D) {
      if (settings$pmx) prior_mean <- x[, i] else prior_mean <- 0
      samp <- sample_theta(w[[j - 1]][, i], dx, 
                           g = eps, theta_w[j - 1, i],
                           alpha = settings$alpha$theta_w, 
                           beta = settings$beta$theta_w, l = settings$l, 
                           u = settings$u, outer = FALSE, v = v,
                           prior_mean = prior_mean,
                           scale = settings$inner_tau2)
      theta_w[j, i] <- samp$theta
    }
    
    # Sample hidden Gaussian layer (w)
    if (settings$pmx) prior_mean <- x else prior_mean <- NULL # defaults to zero
    samp <- sample_w(y, w[[j - 1]], dw, dx, 
                     g[j], theta_y[j], theta_w[j, ], 
                     ll_prev = ll_outer, v = v, prior_mean = prior_mean,
                     scale = settings$inner_tau2)
    w[[j]] <- samp$w
    ll_outer <- samp$ll
    ll_store[j] <- ll_outer
    dw <- samp$dw
  } # end of j for loop
  
  return(list(g = g, theta_y = theta_y, theta_w = theta_w, w = w, tau2 = tau2,
              ll = ll_store))
}

# Two layer Gibbs MONOTONE ----------------------------------------------------

gibbs_two_layer_mono <- function(x, y, x_grid, nmcmc, D, verb, initial, true_g, 
                                 settings, v) {

  ng <- nrow(x_grid)
  dx_grid <- list()
  for (i in 1:D) dx_grid[[i]] <- sq_dist(x_grid[, i])
  
  # Snap initial$w to grid
  w0 <- matrix(nrow = ng, ncol = D)
  for (i in 1:D) 
    w0[, i] <- fo_approx(x[, i], initial$w[, i], x_grid[, i]) # calculates index
  initial$w <- w0 
  
  g <- vector(length = nmcmc)
  if (is.null(true_g)) g[1] <- initial$g else g[1] <- true_g
  theta_y <- matrix(nrow = nmcmc, ncol = D)
  theta_y[1, ] <- initial$theta_y
  theta_w <- matrix(nrow = nmcmc, ncol = D)
  theta_w[1, ] <- initial$theta_w
  w_grid <- list()
  w_grid[[1]] <- initial$w
  tau2 <- vector(length = nmcmc)
  tau2[1] <- initial$tau2
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll_outer <- NULL

  grid_index <- fo_approx_init(x_grid, x)
  w_warp_curr <- monowarp_ref(x, x_grid, w_grid[[1]], grid_index)
  
  for (j in 2:nmcmc) {
    
    if (verb & (j %% 500 == 0)) cat(j, '\n')

    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g_sep(y, w_warp_curr, g[j - 1], theta_y[j - 1, ], 
                       alpha = settings$alpha$g, beta = settings$beta$g, 
                       l = settings$l, u = settings$u, ll_prev = ll_outer, 
                       v = v)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else g[j] <- true_g
    
    # Sample outer lengthscale (theta_y) 
    theta_curr <- theta_y[j - 1, ]
    tau2[j] <- tau2[j - 1] # start repeating tau2 in case there are no acceptances
    for (i in 1:D) {
      samp <- sample_theta_sep(y, w_warp_curr, g[j], theta_curr, index = i,
                               alpha = settings$alpha$theta_y, 
                               beta = settings$beta$theta_y, l = settings$l, 
                               u = settings$u, ll_prev = ll_outer, 
                               v = v, calc_tau2 = (i == D))
      theta_curr[i] <- samp$theta
      theta_y[j, i] <- samp$theta
      ll_outer <- samp$ll
      if (!is.null(samp$tau2)) tau2[j] <- samp$tau2
    }
    
    # Sample inner lengthscale (theta_w) - separately for each dimension
    for (i in 1:D) {
      if (settings$pmx) prior_mean <- x_grid[, i] else prior_mean <- 0
      samp <- sample_theta(w_grid[[j - 1]][, i], dx_grid[[i]], 
                           g = eps, theta_w[j - 1, i],
                           alpha = settings$alpha$theta_w, 
                           beta = settings$beta$theta_w, l = settings$l, 
                           u = settings$u, outer = FALSE, v = v,
                           prior_mean = prior_mean,
                           scale = settings$inner_tau2)
      theta_w[j, i] <- samp$theta
    }
    
    # Sample hidden Gaussian layer (w)
    if (settings$pmx) prior_mean <- x_grid else prior_mean <- NULL # defaults to zero
    samp <- sample_w_mono(y, w_grid[[j - 1]], w_warp_curr, x, x_grid, dx_grid,
                          grid_index, g[j], theta_y[j, ], theta_w[j, ], 
                          ll_prev = ll_outer, v = v, prior_mean = prior_mean,
                          scale = settings$inner_tau2)
    w_grid[[j]] <- samp$w_grid
    w_warp_curr <- samp$w_warp
    ll_outer <- samp$ll
    ll_store[j] <- ll_outer
  } # end of j for loop
  
  return(list(g = g, theta_y = theta_y, theta_w = theta_w, w_grid = w_grid, 
              tau2 = tau2, ll = ll_store))
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
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll_outer <- NULL
  
  for (j in 2:nmcmc) {
    
    if (verb & (j %% 500 == 0)) cat(j, '\n')
    
    # Sample nugget (g)
    if (is.null(true_g)) {
      samp <- sample_g(y, dw, g[j - 1], theta_y[j - 1], 
                       alpha = settings$alpha$g, beta = settings$beta$g, 
                       l = settings$l, u = settings$u, ll_prev = ll_outer, 
                       v = v)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else g[j] <- true_g
    
    # Sample outer lengthscale (theta_y)
    samp <- sample_theta(y, dw, g[j], theta_y[j - 1], 
                         alpha = settings$alpha$theta_y,
                         beta = settings$beta$theta_y, l = settings$l, 
                         u = settings$u, outer = TRUE, ll_prev = ll_outer, 
                         v = v, calc_tau2 = TRUE)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    if (is.null(samp$tau2)) tau2[j] <- tau2[j - 1] else tau2[j] <- samp$tau2
    
    # Sample middle lengthscale (theta_w)
    ll_mid <- 0 # re-calculated each time since we have a new z
    for (i in 1:D) {
      samp <- sample_theta(w[[j - 1]][, i], dz, g = eps, theta_w[j - 1, i],
                           alpha = settings$alpha$theta_w, 
                           beta = settings$beta$theta_w, l = settings$l, 
                           u = settings$u, outer = FALSE, v = v)
      theta_w[j, i] <- samp$theta
      ll_mid <- ll_mid + samp$ll
    }
    
    # Sample inner lengthscale (theta_z)
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
    ll_store[j] <- ll_outer
    dw <- samp$dw
  } # end of j for loop
  
  return(list(g = g, theta_y = theta_y, theta_w = theta_w, theta_z = theta_z,
              w = w, z = z, tau2 = tau2, ll = ll_store))
}
