
# Function Contents -----------------------------------------------------------
# Internal: 
#   gibbs_one_layer
#   gibbs_two_layer
#   gibbs_two_layer_grad
#   gibbs_two_layer_mono
#   gibbs_three_layer

# gibbs_one_layer -------------------------------------------------------------
# This one function handles sep = TRUE or FALSE and grad_enhance = TRUE or FALSE
# (things are much easier in one layer)

gibbs_one_layer <- function(x, y, dydx, nmcmc, verb, initial, true_g, settings, v) {
  
  d <- ncol(x)
  n <- nrow(x)
  
  grad_enhance <- !is.null(dydx)
  if (grad_enhance) { 
    y <- c(y, as.vector(dydx))
  }

  if (!settings$sep & !grad_enhance) {
    xdmat <- sq_dist(x)
  } else xdmat <- NULL

  if (is.null(true_g)) {
    g_store <- vector(length = nmcmc)
    g_store[1] <- initial$g
    g <- initial$g
  } else g <- true_g
  tau2 <- vector(length = nmcmc)
  tau2[1] <- NA
  if (settings$sep) {
    theta <- matrix(nrow = nmcmc, ncol = d)
    theta[1, ] <- initial$theta 
  } else {
    theta <- vector(length = nmcmc)
    theta[1] <- initial$theta
  }
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll <- NULL
  
  for (j in 2:nmcmc) {
    
    if (verb & (j %% 500 == 0)) cat(j, '\n')
    
    tau2[j] <- tau2[j-1] # start repeating tau2 in case there are no acceptances
    
    # Sample nugget (g) - only if true_g is not specified
    if (is.null(true_g)) {
      samp <- sample_g(y, xdmat = xdmat,
                       x = ifel(settings$sep | grad_enhance, x, NULL),
                       theta = ifel(settings$sep, theta[j-1, ], theta[j-1]), 
                       g = g, 
                       v = v,
                       alpha = settings$g$alpha, 
                       beta = settings$g$beta, 
                       l = settings$l, 
                       u = settings$u, 
                       ll_prev = ll,
                       grad_enhance = grad_enhance,
                       sep = settings$sep)
      g <- samp$g # use this value in later sampling
      g_store[j] <- g
      ll <- samp$ll
      if (!is.null(samp$tau2)) tau2[j] <- samp$tau2
    } 
    
    # Sample lengthscale (theta)
    if (settings$sep) {
      # Start by repeating previous samples, these will be overwritten if accepted
      theta[j, ] <- theta[j-1, ]
      for (i in 1:d) {
        samp <- sample_theta(y, x = x, 
                             tau2 = 1, # scale integrated out of outer layer 
                             theta = theta[j, ], # includes updated theta for j < i
                             g = g,
                             v = v,
                             alpha = settings$theta$alpha,
                             beta = settings$theta$beta, 
                             l = settings$l, 
                             u = settings$u, 
                             outer = TRUE,
                             ll_prev = ll,
                             grad_enhance = grad_enhance,
                             sep = TRUE,
                             index = i)
        theta[j, i] <- samp$theta
        ll <- samp$ll
        if (!is.null(samp$tau2)) tau2[j] <- samp$tau2
      }
    } else {
      samp <- sample_theta(y, xdmat = xdmat, 
                           x = ifel(grad_enhance, x, NULL),
                           tau2 = 1, # scale integrated out of outer layer 
                           theta = theta[j-1], 
                           g = g,
                           v = v,
                           alpha = settings$theta$alpha,
                           beta = settings$theta$beta, 
                           l = settings$l, 
                           u = settings$u, 
                           outer = TRUE,
                           ll_prev = ll,
                           grad_enhance = grad_enhance,
                           sep = FALSE)
      theta[j] <- samp$theta
      ll <- samp$ll
      if (!is.null(samp$tau2)) tau2[j] <- samp$tau2
    } # end of else statement
    ll_store[j] <- ll
  } # end of j for loop
  
  return(list(g = ifel(is.null(true_g), g_store, g), theta = theta, tau2 = tau2, 
              ll = ll_store))
}

# gibbs_two_layer -------------------------------------------------------------

gibbs_two_layer <- function(x, y, nmcmc, D, verb, initial, true_g, 
                            settings, v) {
  
  n <- length(y)
  xdmat <- sq_dist(x)
  
  if (is.null(true_g)) {
    g_store <- vector(length = nmcmc)
    g_store[1] <- initial$g
    g <- initial$g
  } else g <- true_g
  tau2_y <- vector(length = nmcmc)
  tau2_y[1] <- NA
  theta_y <- vector(length = nmcmc)
  theta_y[1] <- initial$theta_y
  theta_w <- matrix(nrow = nmcmc, ncol = D)
  theta_w[1, ] <- initial$theta_w
  w <- array(dim = c(nmcmc, n, D))
  w[1, , ] <- initial$w
  wdmat <- sq_dist(w[1, , ])
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA
  ll_outer <- NULL

  for (j in 2:nmcmc) {
    
    if (verb & (j %% 500 == 0)) cat(j, '\n')
    
    # Sample nugget (g) - only if true_g is not specified
    if (is.null(true_g)) {
      samp <- sample_g(y, xdmat = wdmat, 
                       theta = theta_y[j-1],
                       g = g,
                       v = v,
                       alpha = settings$g$alpha, 
                       beta = settings$g$beta, 
                       l = settings$l, 
                       u = settings$u, 
                       ll_prev = ll_outer)
      g <- samp$g
      g_store[j] <- g
      ll_outer <- samp$ll
    }
    
    # Sample outer lengthscale (theta_y) 
    samp <- sample_theta(y, xdmat = wdmat, 
                         tau2 = 1, # scale integrated out of outer layer 
                         theta = theta_y[j-1], 
                         g = g,
                         v = v,
                         alpha = settings$theta_y$alpha, 
                         beta = settings$theta_y$beta, 
                         l = settings$l, 
                         u = settings$u, 
                         outer = TRUE,
                         ll_prev = ll_outer)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll

    # Sample inner lengthscale (theta_w) - separately for each dimension
    # The inner likelihood needs to be updated every time since we have new W
    for (i in 1:D) {
      samp <- sample_theta(w[j-1, , i], xdmat = xdmat, 
                           tau2 = settings$tau2_w,
                           theta = theta_w[j-1, i],
                           g = eps, 
                           v = v,
                           alpha = settings$theta_w$alpha, 
                           beta = settings$theta_w$beta, 
                           l = settings$l, 
                           u = settings$u, 
                           outer = FALSE, 
                           ll_prev = NULL,
                           prior_mean = ifel(settings$pmx, x[, i], 0))
      theta_w[j, i] <- samp$theta
    }
    
    # Sample hidden Gaussian layer (w)
    samp <- sample_w(y, w = w[j-1, , ], 
                     xdmat = xdmat, 
                     tau2_w = settings$tau2_w,
                     theta_y = theta_y[j],
                     theta_w = theta_w[j, ],
                     g = g,
                     v = v,
                     ll_prev = ll_outer,
                     prior_mean = ifel(settings$pmx, x, NULL))
    w[j, , ] <- samp$w
    wdmat <- samp$wdmat
    ll_outer <- samp$ll
    ll_store[j] <- ll_outer
    tau2_y[j] <- samp$tau2_y # updated every iteration, since we always have a new W
  } # end of j for loop
  
  return(list(g = ifel(is.null(true_g), g_store, g), tau2_y = tau2_y, 
              theta_y = theta_y, theta_w = theta_w, w = w, ll = ll_store))
}

# gibbs_two_layer_grad --------------------------------------------------------

gibbs_two_layer_grad <- function(x, y, dydx, nmcmc, verb, initial, 
                                 true_g, settings, v) {
  
  n <- length(y)
  d <- ncol(x)
  xdmat <- sq_dist(x)

  if (is.null(true_g)) {
    g_store <- vector(length = nmcmc)
    g_store[1] <- initial$g
    g <- initial$g
  } else g <- true_g
  tau2_y <- vector(length = nmcmc)
  tau2_y[1] <- NA
  theta_y <- vector(length = nmcmc)
  theta_y[1] <- initial$theta_y
  theta_w <- matrix(nrow = nmcmc, ncol = d)
  theta_w[1, ] <- initial$theta_w
  w <- array(dim = c(nmcmc, n*(d+1), d)) # contains w and dwdx
  w[1, , ] <- initial$w

  dydw <- get_dydw(w[1, , ], dydx)
  y_all <- c(y, as.vector(dydw))
  
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- NA 
  ll_outer <- NULL

  for (j in 2:nmcmc) {

    if (verb & (j %% 500 == 0)) cat(j, '\n')
    
    # Sample nugget (g) - only if true_g is not specified
    if (is.null(true_g)) {
      samp <- sample_g(y_all, x = w[j-1, 1:n, ], 
                       theta = theta_y[j-1],
                       g = g,
                       v = v,
                       alpha = settings$g$alpha, 
                       beta = settings$g$beta, 
                       l = settings$l, 
                       u = settings$u, 
                       ll_prev = ll_outer,
                       grad_enhance = TRUE)
      g <- samp$g
      g_store[j] <- g
      ll_outer <- samp$ll
    }
    
    # Sample outer lengthscale (theta_y) 
    samp <- sample_theta(y_all, x = w[j-1, 1:n, ],
                         tau2 = 1, # scale integrated out of outer layer 
                         theta = theta_y[j-1], 
                         g = g,
                         v = v,
                         alpha = settings$theta_y$alpha, 
                         beta = settings$theta_y$beta, 
                         l = settings$l, 
                         u = settings$u, 
                         outer = TRUE,
                         ll_prev = ll_outer,
                         grad_enhance = TRUE)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll

    # Sample inner lengthscale (theta_w) - separately for each dimension
    for (i in 1:d) {
      samp <- sample_theta(w[j-1, , i], x = x, 
                           tau2 = settings$tau2_w,
                           theta = theta_w[j-1, i],
                           g = eps, 
                           v = v,
                           alpha = settings$theta_w$alpha, 
                           beta = settings$theta_w$beta, 
                           l = settings$l, 
                           u = settings$u, 
                           outer = FALSE, 
                           ll_prev = NULL, # needs to be recalculated    
                           prior_mean = ifel(settings$pmx, settings$w_prior_mean[, i], 0),
                           grad_enhance = TRUE)   
      theta_w[j, i] <- samp$theta
    }
    
    # Sample hidden Gaussian layer (w)
    samp <- sample_w_grad(y = y,
                          dydx = dydx,
                          w = w[j-1, , ],
                          x = x,
                          tau2_w = settings$tau2_w,
                          theta_y = theta_y[j],
                          theta_w = theta_w[j, ],
                          g = g,
                          v = v,
                          ll_prev = ll_outer,
                          prior_mean = ifel(settings$pmx, settings$w_prior_mean, NULL))
    y_all <- samp$y_all
    w[j, , ] <- samp$w
    ll_outer <- NULL # need to start fresh with our new w
    ll_store[j] <- samp$ll
    tau2_y[j] <- samp$tau2_y # updated every iteration, since we always have a new W
  } # end of j for loop
  
  return(list(g = ifel(is.null(true_g), g_store, g), tau2_y = tau2_y, 
              theta_y = theta_y, theta_w = theta_w, w = w, ll = ll_store))
}

# gibbs_two_layer_mono --------------------------------------------------------

gibbs_two_layer_mono <- function(x, y, x_grid, nmcmc, verb, initial, true_g, 
                                 settings, v) {
  
  n <- length(y)
  ng <- length(x_grid)
  d <- ncol(x)
  xdmat_grid <- sq_dist(x_grid)
  
  if (is.null(true_g)) {
    g_store <- vector(length = nmcmc)
    g_store[1] <- initial$g
    g <- initial$g
  } else g <- true_g
  tau2_y <- vector(length = nmcmc)
  tau2_y[1] <- NA
  theta_y <- vector(length = nmcmc)
  theta_y[1] <- initial$theta_y
  tau2_w <- matrix(nrow = nmcmc, ncol = d)
  tau2_w[1, ] <- rep(NA, d)
  theta_w <- matrix(nrow = nmcmc, ncol = d)
  theta_w[1, ] <- initial$theta_w
  w_grid <- array(dim = c(nmcmc, ng, d))
  for (i in 1:d) w_grid[1, , i] <- x_grid
  w <- array(dim = c(nmcmc, n, d))
  w[1, , ] <- monotransform(x, x_grid, w_grid[1, , ])
  wdmat <- sq_dist(w[1, , ])

  ll_outer <- logl(y, xdmat = wdmat, tau2 = 1, theta = initial$theta_y, g = g, 
                    v = v, outer = TRUE)$ll
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- ll_outer

  for (j in 2:nmcmc) {
    
    if (verb & (j %% 500 == 0)) cat(j, '\n')
    
    # Sample nugget (g) - only if true_g is not specified
    if (is.null(true_g)) {
      samp <- sample_g(y, xdmat = wdmat, 
                       theta = theta_y[j-1],
                       g = g,
                       v = v,
                       alpha = settings$g$alpha, 
                       beta = settings$g$beta, 
                       l = settings$l, 
                       u = settings$u, 
                       ll_prev = ll_outer)
      g <- samp$g
      g_store[j] <- g
      ll_outer <- samp$ll
    }
    
    # Sample outer lengthscale (theta_y) 
    samp <- sample_theta(y, xdmat = wdmat, 
                         tau2 = 1, # scale integrated out of outer layer 
                         theta = theta_y[j-1], 
                         g = g,
                         v = v,
                         alpha = settings$theta_y$alpha, 
                         beta = settings$theta_y$beta, 
                         l = settings$l, 
                         u = settings$u, 
                         outer = TRUE,
                         ll_prev = ll_outer)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll

    # Sample inner lengthscale (theta_w) - separately for each dimension
    # Acts on gridded, unwarped GP samples
    # AND save MLE estimates of tau2 on inner layer
    for (i in 1:d) {
      samp <- sample_theta(w_grid[j-1, , i], xdmat = xdmat_grid, 
                           tau2 = ifel(d == 1, settings$tau2_w, 1), # scale integrated out of inner layer when d > 1
                           theta = theta_w[j-1, i],
                           g = eps, 
                           v = v,
                           alpha = settings$theta_w$alpha, 
                           beta = settings$theta_w$beta, 
                           l = settings$l, 
                           u = settings$u, 
                           outer = ifel(d == 1, FALSE, TRUE), # in one dimension, fix tau2_w
                           ll_prev = NULL,
                           prior_mean = ifel(settings$pmx, x_grid, 0))
      theta_w[j, i] <- samp$theta
      if (d > 1) {
        if (is.null(samp$tau2)) { # we got a rejection but still need to calculate tau2
          tau2_w[j, i] <- logl(w_grid[j-1, , i], xdmat = xdmat_grid, tau2 = 1, 
                               theta = theta_w[j, i], g = g, v = v, 
                               mu = ifel(settings$pmx, x_grid, 0), outer = TRUE)$tau2
        } else tau2_w[j, i] <- samp$tau2
        if (tau2_w[j, i] < 0.01) tau2_w[j, i] <- 0.01
        if (tau2_w[j, i] > 10) tau2_w[j, i] <- 10
      }
    }

    # Sample hidden Gaussian layer (w)
    samp <- sample_w_mono(y, w = w[j-1, , ], x = x,
                          x_grid = x_grid,
                          w_grid = w_grid[j-1, , ], 
                          xdmat_grid = xdmat_grid, 
                          tau2_w = ifel(d == 1, settings$tau2_w, tau2_w[j, ]),
                          theta_y = theta_y[j],
                          theta_w = theta_w[j, ],
                          g = g,
                          v = v,
                          ll_prev = ll_outer,
                          prior_mean = ifel(settings$pmx, x_grid, rep(0, ng)))
    w[j, , ] <- samp$w
    wdmat <- samp$wdmat
    w_grid[j, , ] <- samp$w_grid
    ll_outer <- samp$ll
    ll_store[j] <- ll_outer
    tau2_y[j] <- samp$tau2_y
  } # end of j for loop
  
  return(list(g = ifel(is.null(true_g), g_store, g), tau2_y = tau2_y, 
              theta_y = theta_y, tau2_w = ifel(d == 1, NULL, tau2_w), 
              theta_w = theta_w, w = w, w_grid = w_grid,
              ll = ll_store))
}

# gibbs_three_layer -----------------------------------------------------------

gibbs_three_layer <- function(x, y, nmcmc, D, verb, initial, true_g, 
                              settings, v) {
  
  n <- length(y)
  xdmat <- sq_dist(x)
  zdmat <- sq_dist(initial$z)
  wdmat <- sq_dist(initial$w)
  
  if (is.null(true_g)) {
    g_store <- vector(length = nmcmc)
    g_store[1] <- initial$g 
    g <- initial$g
  } else g <- true_g
  tau2_y <- vector(length = nmcmc)
  tau2_y[1] <- NA
  theta_y <- vector(length = nmcmc)
  theta_y[1] <- initial$theta_y
  theta_w <- matrix(nrow = nmcmc, ncol = D)
  theta_w[1, ] <- initial$theta_w
  theta_z <- matrix(nrow = nmcmc, ncol = D)
  theta_z[1, ] <- initial$theta_z
  w <- array(dim = c(nmcmc, n, D))
  w[1, , ] <- initial$w
  z <- array(dim = c(nmcmc, n, D))
  z[1, , ] <- initial$z
  ll_outer <- logl(y, xdmat = wdmat, tau2 = 1, theta = initial$theta_y, g = g, 
                    v = v, outer = TRUE)$ll
  ll_store <- vector(length = nmcmc)
  ll_store[1] <- ll_outer
  
  for (j in 2:nmcmc) {
    
    if (verb & (j %% 500 == 0)) cat(j, '\n')
    
    # Sample nugget (g) - only if true_g is not specified
    if (is.null(true_g)) {
      samp <- sample_g(y, xdmat = wdmat, 
                       theta = theta_y[j-1], 
                       g = g,
                       v = v,
                       alpha = settings$g$alpha, 
                       beta = settings$g$beta, 
                       l = settings$l, 
                       u = settings$u, 
                       ll_prev = ll_outer)
      g <- samp$g
      g_store[j] <- g
      ll_outer <- samp$ll
    }
    
    # Sample outer lengthscale (theta_y)
    samp <- sample_theta(y, xdmat = wdmat, 
                         tau2 = 1, # scale integrated out of outer layer
                         theta = theta_y[j-1], 
                         g = g,
                         v = v,
                         alpha = settings$theta_y$alpha,
                         beta = settings$theta_y$beta, 
                         l = settings$l, 
                         u = settings$u, 
                         outer = TRUE, 
                         ll_prev = ll_outer)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll

    # Sample middle lengthscale (theta_w)
    ll_mid <- 0 # re-calculated each time since we have a new z
    for (i in 1:D) {
      samp <- sample_theta(w[j-1, , i], xdmat = zdmat, 
                           tau2 = settings$tau2_w,
                           theta = theta_w[j-1, i],
                           g = eps,
                           v = v,
                           alpha = settings$theta_w$alpha, 
                           beta = settings$theta_w$beta, 
                           l = settings$l, 
                           u = settings$u, 
                           outer = FALSE,
                           ll_prev = NULL)
      theta_w[j, i] <- samp$theta
      ll_mid <- ll_mid + samp$ll
    }
    
    # Sample inner lengthscale (theta_z)
    for (i in 1:D) {
      samp <- sample_theta(z[j-1, , i], xdmat = xdmat, 
                           tau2 = settings$tau2_z,
                           theta = theta_z[j-1, i], 
                           g = eps,  
                           v = v,
                           alpha = settings$theta_z$alpha, 
                           beta = settings$theta_z$beta, 
                           l = settings$l, 
                           u = settings$u, 
                           outer = FALSE,
                           ll_prev = NULL)
      theta_z[j, i] <- samp$theta
    }
    
    # Sample inner hidden Gaussian layer (z)
    samp <- sample_z(w[j-1, , ], z = z[j-1, , ], 
                     xdmat = xdmat, 
                     tau2_w = settings$tau2_w,
                     tau2_z = settings$tau2_z,
                     theta_w = theta_w[j, ],
                     theta_z = theta_z[j, ],
                     v = v,
                     ll_prev = ll_mid)
    z[j, , ] <- samp$z
    zdmat <- samp$zdmat
    
    # Sample middle hidden Gaussian layer (w)
    samp <- sample_w(y, w = w[j-1, , ], 
                     xdmat = zdmat, 
                     tau2_w = settings$tau2_w,
                     theta_y = theta_y[j],
                     theta_w = theta_w[j, ], 
                     g = g,
                     v = v,
                     ll_prev = ll_outer)
    w[j, , ] <- samp$w
    wdmat <- samp$wdmat
    ll_outer <- samp$ll
    ll_store[j] <- ll_outer
    tau2_y[j] <- samp$tau2_y
  } # end of j for loop
  
  return(list(g = ifel(is.null(true_g), g_store, g), tau2_y = tau2_y, 
              theta_y = theta_y, theta_w = theta_w, theta_z = theta_z, 
              w = w, z = z, ll = ll_store))
}
