
# Function Contents -----------------------------------------------------------
# External: (see documentation in predict.R file)
#   predict.gpvec
#   predict.dgp2vec
#   predict.dgp3vec
# Internal:
#   predict_vec (used in all three of the above)

# Predict One Layer Vecchia ---------------------------------------------------
#' @rdname predict
#' @export

predict.gpvec <- function(object, x_new, m = object$m, ordering_new = NULL,
                          lite = TRUE, return_all = FALSE, EI = FALSE,
                          entropy_limit = NULL, cores = 1, ...) {
  
  settings <- list(ordering_new = ordering_new, lite = lite, 
                   return_all = return_all, EI = EI, 
                   entropy_limit = entropy_limit, cores = cores)
  object <- predict_vec(object, x_new, m, settings, layers = 1)
  return(object)
}

# Predict Two Layer Vecchia ---------------------------------------------------
#' @rdname predict
#' @export
#' 
predict.dgp2vec <- function(object, x_new, m = object$m, ordering_new = NULL,
                            lite = TRUE, store_latent = FALSE, mean_map = TRUE, 
                            return_all = FALSE, EI = FALSE, entropy_limit = NULL,
                            cores = 1, ...) {

  settings <- list(ordering_new = ordering_new, lite = lite,
                   store_latent = store_latent, mean_map = mean_map,
                   return_all = return_all, EI = EI, 
                   entropy_limit = entropy_limit, cores = cores)
  object <- predict_vec(object, x_new, m, settings, layers = 2)
  return(object)
}
  
# Predict Three Layer Vecchia -------------------------------------------------
#' @rdname predict
#' @export
  
predict.dgp3vec <- function(object, x_new, m = object$m, ordering_new = NULL,
                            lite = TRUE, store_latent = FALSE, mean_map = TRUE,
                            return_all = FALSE, EI = FALSE, entropy_limit = NULL,
                            cores = 1, ...) {
  
  settings <- list(ordering_new = ordering_new, lite = lite,
                   store_latent = store_latent, mean_map = mean_map,
                   return_all = return_all, EI = EI, 
                   entropy_limit = entropy_limit, cores = cores)
  object <- predict_vec(object, x_new, m, settings, layers = 3)
  return(object)
}

# Predict Vecchia -------------------------------------------------------------

predict_vec <- function(object, x_new, m, settings, layers) {
  
  tic <- proc.time()[[3]]
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  monowarp <- (!is.null(object$x_grid))
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  object$m_pred <- m

  if (monowarp) {
    grid_index <- fo_approx_init(object$x_grid, object$x)
    grid_index_new <- fo_approx_init(object$x_grid, x_new)
  }

  if (layers >= 2) {
    if (monowarp) D <- ncol(object$w_grid[[1]]) else D <- ncol(object$w[[1]])
    if (!settings$mean_map) 
      stop("mean_map = FALSE is not available for `vecchia = TRUE` case")
  }

  if (layers == 1) {
    sep <- is.matrix(object$theta)
  } else sep <- monowarp # separable lengthscales only used in monotonic warpings
  
  if (settings$return_all & !settings$lite) 
    stop("return_all only offered when lite = TRUE")
  if (!is.null(settings$entropy_limit) & !is.numeric(settings$entropy_limit))
    stop("entropy_limit must be numeric")
  if (settings$EI) { # if no noise, use smallest observed value, else estimate f_min
    if (all(object$g <= 1e-6)) {
      f_min <- FALSE 
      y_min <- min(object$y)
    } else f_min <- TRUE
  } else f_min <- FALSE
  if (!is.null(settings$ordering_new)) {
    if (settings$lite) message("ordering_new is only relevant when lite = FALSE")
    test <- check_ordering(settings$ordering_new, nrow(x_new))
  }
  
  # If not monotonic, pre-calculate nearest neighbors for x
  if (!monowarp) {
    if (settings$lite) {
      NN_x_new <- FNN::get.knnx(object$x, x_new, m)$nn.index
      x_approx <-  NULL
    } else { 
      NN_x_new <- NULL
      x_approx <- add_pred_to_approx(object$x_approx, x_new, m, settings$ordering_new)
    }
  }
  
  # Initialize prior mean of 0 (these will only be changed if object$settings$pmx = TRUE)
  prior_mean_new <- 0
  prior_mean <- rep(0, length(object$y))
  prior_tau2 <- 1
  
  if (settings$cores == 1) { # run serial for loop
    
    mu_t <- matrix(nrow = n_new, ncol = object$nmcmc)
    if (settings$lite) {
      s2_sum <- rep(0, times = n_new)
      if (settings$return_all) s2_t <- matrix(nrow = n_new, ncol = object$nmcmc)
    } else sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    if (settings$EI) ei_sum <- rep(0, times = n_new)
    if (!is.null(settings$entropy_limit)) ent_sum <- rep(0, times = n_new)
    if (layers >= 2) {
      if (settings$store_latent) {
        w_new_list <- list()
        if (layers == 3) z_new_list <- list()
      }
    }
      
    for (t in 1:object$nmcmc) {
      
      if (layers == 3) {
        # 3 layers: map x_new to z_new (separately for each dimension)
        z_t <- object$z[[t]]
        z_new <- matrix(nrow = n_new, ncol = D)
        for (i in 1:D) { # mean_map = TRUE only
          k <- krig_vec(z_t[, i], object$theta_z[t, i], g = eps, tau2 = 1,
                        v = object$v, m = m, x = object$x, x_new = x_new,
                        NNarray_pred = NN_x_new) # lite = TRUE version only
          z_new[, i] <- k$mean
        } # end of i for loop
        if (settings$store_latent) z_new_list[[t]] <- z_new
      }
      
      if (layers >= 2) {
        # 2 layers: map x_new to w_new (separately for each dimension)
        # 3 layers: map z_new to w_new (separately for each dimension)
        if (monowarp) { # simply get the monowarped values at x_new locations
          w_grid_t <- object$w_grid[[t]]
          w_t <- monowarp_ref(object$x, object$x_grid, w_grid_t, grid_index)
          w_new <- monowarp_ref(x_new, object$x_grid, w_grid_t, grid_index_new)
          if (settings$lite) {
            w_approx <- NULL
          } else {
            w_approx <- update_obs_in_approx(object$w_approx, w_t)
            w_approx <- add_pred_to_approx(w_approx, w_new, m, settings$ordering_new)
          }
        } else {
          w_t <- object$w[[t]]
          w_new <- matrix(nrow = n_new, ncol = D)
          for (i in 1:D) {
            if (layers == 2) {
              if (object$settings$pmx) { # Optional prior mean of x
                prior_mean_new <- x_new[, i]
                prior_mean <- object$x[, i]
                prior_tau2 <- object$settings$inner_tau2
              }
            } 
            k <- krig_vec(w_t[, i], object$theta_w[t, i], g = eps, tau2 = prior_tau2,
                          v = object$v, m = m, x = ifel(layers == 2, object$x, z_t), 
                          x_new = ifel(layers == 2, x_new, z_new),
                          NNarray_pred = ifel(layers == 2, NN_x_new, NULL),
                          prior_mean = prior_mean, 
                          prior_mean_new = prior_mean_new) # lite = TRUE version only
            w_new[, i] <- k$mean
          } # end of i for loop
          if (settings$store_latent) w_new_list[[t]] <- w_new
          
          if (settings$lite) {
            w_approx <- NULL
          } else {
            w_approx <- update_obs_in_approx(object$w_approx, w_t)
            w_approx <- add_pred_to_approx(w_approx, w_new, m, settings$ordering_new)
          }
        } # end of monowarp else statement
      }
      
      # 1 layer: map x_new to y
      # 2 and 3 layers: map w_new to y
      if (sep) { # only occurs in one layer or monowarped two-layer
        k <- krig_vec(object$y, 
                      ifel(monowarp, object$theta_y[t, ], object$theta[t, ]), 
                      object$g[t], 
                      object$tau2[t],
                      s2 = settings$lite, sigma = !settings$lite, v = object$v, 
                      m = m, 
                      x = ifel(monowarp, w_t, object$x), 
                      x_new = ifel(monowarp, w_new, x_new), 
                      NNarray_pred = ifel(monowarp, NULL, NN_x_new),
                      approx = ifel(monowarp, w_approx, x_approx), 
                      sep = TRUE, f_min = f_min)
      } else {
        k <- krig_vec(object$y, ifel(layers == 1, object$theta[t], object$theta_y[t]), 
                      object$g[t], object$tau2[t],
                      s2 = settings$lite, sigma = !settings$lite, v = object$v, 
                      m = m,
                      x = ifel(layers == 1, object$x, w_t), 
                      x_new = ifel(layers == 1, x_new, w_new), 
                      NNarray_pred = ifel(layers == 1, NN_x_new, NULL),
                      approx = ifel(layers == 1, x_approx, w_approx),
                      f_min = f_min)
      }
      mu_t[, t] <- k$mean
      if (settings$lite) {
        s2_sum <- s2_sum + k$s2
        if (settings$return_all) s2_t[, t] <- k$s2
      } else sigma_sum <- sigma_sum + k$sigma
      if (settings$EI) {
        if (settings$lite) {
          sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
        ei_sum <- ei_sum + exp_improv(k$mean, sig2, ifel(f_min, k$f_min, y_min))
      }
      if (!is.null(settings$entropy_limit)) {
        if (settings$lite) {
          sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
        ent_sum <- ent_sum + calc_entropy(k$mean, sig2, settings$entropy_limit)
      }
    } # end of t for loop
    
  } else { # run in parallel using foreach
    
    iters <- 1:object$nmcmc
    chunks <- split(iters, sort(cut(iters, settings$cores, labels = FALSE)))
    if (settings$cores > detectCores()) 
      warning("cores is greater than available nodes")
    
    cl <- makeCluster(settings$cores)
    registerDoParallel(cl)
    
    thread <- NULL
    result <- foreach(thread = 1:settings$cores) %dopar% {
      out <- list()
      out$mu_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      if (settings$lite) {
        out$s2_sum <- rep(0, times = n_new)
        if (settings$return_all) 
          out$s2_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      } else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
      if (settings$EI) out$ei_sum <- rep(0, times = n_new)
      if (!is.null(settings$entropy_limit)) out$ent_sum <- rep(0, times = n_new)
      if (layers >= 2) {
        if (settings$store_latent) {
          out$w_new <- list()
          if (layers == 3) out$z_new <- list()
        }
      }
      
      # calculate predictions for each candidate MCMC iteration
      j <- 1
      for (t in chunks[[thread]]) {
        
        if (layers == 3) {
          # 3 layers: map x_new to z_new (separately for each dimension)
          z_t <- object$z[[t]]
          z_new <- matrix(nrow = n_new, ncol = D)
          for (i in 1:D) { # mean_map = TRUE only
            k <- krig_vec(z_t[, i], object$theta_z[t, i], g = eps, tau2 = 1,
                          v = object$v, m = m, x = object$x, x_new = x_new,
                          NNarray_pred = NN_x_new) # lite = TRUE version only
            z_new[, i] <- k$mean
          } # end of i for loop
          if (settings$store_latent) out$z_new[[j]] <- z_new
        }
        
        if (layers >= 2) {
          # 2 layers: map x_new to w_new (separately for each dimension)
          # 3 layers: map z_new to w_new (separately for each dimension)
          if (monowarp) { # simply get the monowarped values at x_new locations
            w_grid_t <- object$w_grid[[t]]
            w_t <- monowarp_ref(object$x, object$x_grid, w_grid_t, grid_index)
            w_new <- monowarp_ref(x_new, object$x_grid, w_grid_t, grid_index_new)
            if (settings$lite) {
              w_approx <- NULL
            } else {
             w_approx <- update_obs_in_approx(object$w_approx, w_t)
              w_approx <- add_pred_to_approx(w_approx, w_new, m, settings$ordering_new) 
            }
          } else {
            w_t <- object$w[[t]]
            w_new <- matrix(nrow = n_new, ncol = D)
            for (i in 1:D) {
              if (layers == 2) {
                if (object$settings$pmx) { # Optional prior mean of x
                  prior_mean_new <- x_new[, i]
                  prior_mean <- object$x[, i]
                  prior_tau2 <- object$settings$inner_tau2
                }
              } 
              k <- krig_vec(w_t[, i], object$theta_w[t, i], g = eps, tau2 = prior_tau2,
                            v = object$v, m = m, x = ifel(layers == 2, object$x, z_t), 
                            x_new = ifel(layers == 2, x_new, z_new),
                            NNarray_pred = ifel(layers == 2, NN_x_new, NULL),
                            prior_mean = prior_mean, 
                            prior_mean_new = prior_mean_new) # lite = TRUE version only
              w_new[, i] <- k$mean
            } # end of i for loop
            if (settings$store_latent) out$w_new[[j]] <- w_new
            
            if (settings$lite) {
              w_approx <- NULL
            } else {
              w_approx <- update_obs_in_approx(object$w_approx, w_t)
              w_approx <- add_pred_to_approx(w_approx, w_new, m, settings$ordering_new)
            }
          } # end of monowarp else statement
        }
        
        # 1 layer: map x_new to y
        # 2 and 3 layers: map w_new to y
        if (sep) { # only occurs in one layer or monowarped two-layer
          k <- krig_vec(object$y, 
                        ifel(monowarp, object$theta_y[t, ], object$theta[t, ]), 
                        object$g[t], object$tau2[t],
                        s2 = settings$lite, sigma = !settings$lite, v = object$v, 
                        m = m, 
                        x = ifel(monowarp, w_t, object$x), 
                        x_new = ifel(monowarp, w_new, x_new), 
                        NNarray_pred = ifel(monowarp, NULL, NN_x_new),
                        approx = ifel(monowarp, w_approx, x_approx), 
                        sep = TRUE, f_min = f_min)
        } else {
          k <- krig_vec(object$y, ifel(layers == 1, object$theta[t], object$theta_y[t]), 
                        object$g[t], object$tau2[t],
                        s2 = settings$lite, sigma = !settings$lite, v = object$v, 
                        m = m,
                        x = ifel(layers == 1, object$x, w_t), 
                        x_new = ifel(layers == 1, x_new, w_new), 
                        NNarray_pred = ifel(layers == 1, NN_x_new, NULL),
                        approx = ifel(layers == 1, x_approx, w_approx),
                        f_min = f_min)
        }
        out$mu_t[, j] <- k$mean
        if (settings$lite) {
          out$s2_sum <- out$s2_sum + k$s2
          if (settings$return_all) out$s2_t[, j] <- k$s2
        } else out$sigma_sum <- out$sigma_sum + k$sigma
        if (settings$EI) {
          if (settings$lite) {
            sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
          } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
          out$ei_sum <- out$ei_sum + exp_improv(k$mean, sig2, ifel(f_min, k$f_min, y_min))
        }
        if (!is.null(settings$entropy_limit)) {
          if (settings$lite) {
            sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
          } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
          out$ent_sum <- out$ent_sum + calc_entropy(k$mean, sig2, settings$entropy_limit)
        }
        j <- j + 1
      } # end of t for loop
      return(out)
    } # end of foreach loop
    
    stopCluster(cl)
      
    # Group elements out of the list
    mu_t <- do.call(cbind, lapply(result, with, eval(parse(text = "mu_t"))))
    if (settings$lite) {
      s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum"))))
      if (settings$return_all) s2_t <- do.call(cbind, lapply(result, with, eval(parse(text = "s2_t"))))
    } else {
      sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum"))))
    }
    if (layers >= 2) {
      if (settings$store_latent) {
        w_new_list <- unlist(lapply(result, with, eval(parse(text = "w_new"))), recursive = FALSE)
        if (layers == 3) z_new_list <- unlist(lapply(result, with, eval(parse(text = "z_new"))), recursive = FALSE)
      }
    }
    if (settings$EI) ei_sum <- Reduce("+", lapply(result, with, eval(parse(text = "ei_sum"))))
    if (!is.null(settings$entropy_limit)) ent_sum <- Reduce("+", lapply(result, with, 
                                                               eval(parse(text = "ent_sum"))))
    
  } # end of else statement

  # Add variables to the output list
  object$mean <- rowMeans(mu_t)
  if (layers >= 2) {
    if (settings$store_latent) {
      object$w_new <- w_new_list
      if (layers == 3) object$z_new <- z_new_list
    }
  }
  if (settings$lite) { 
    object$s2 <- s2_sum / object$nmcmc + apply(mu_t, 1, var)
    if (settings$return_all) {
      object$mean_all <- mu_t
      object$s2_all <- s2_t
    }
  } else object$Sigma <- sigma_sum / object$nmcmc + cov(t(mu_t))
  if (settings$EI) object$EI <- ei_sum / object$nmcmc
  if (!is.null(settings$entropy_limit)) object$entropy <- drop(ent_sum / object$nmcmc)
  toc <- proc.time()[[3]]
  object$time <- object$time + unname(toc - tic)
  
  return(object)
}