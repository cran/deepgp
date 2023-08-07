
# Function Contents -----------------------------------------------------------
# External: (see documentation in predict.R file)
#   predict.gpvec
#   predict.dgp2vec
#   predict.dgp3vec

# Predict One Layer Vecchia ---------------------------------------------------
#' @rdname predict
#' @export

predict.gpvec <- function(object, x_new, m = object$m, lite = TRUE,
                          return_all = FALSE, entropy_limit = NULL,
                          cores = 1, ...) {
  
  tic <- proc.time()[3]
  object <- clean_prediction(object) # remove previous predictions if present
  sep <- (is.matrix(object$theta))
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  if (return_all & !lite) stop("return_all only offered when lite = TRUE")
  if (!is.null(entropy_limit) & !is.numeric(entropy_limit))
    stop("entropy_limit must be numeric")
  
  # Pre-calculate nearest neighbors for x 
  if (lite) {
    NN_x_new <- FNN::get.knnx(object$x, x_new, m)$nn.index
    x_approx <-  NULL
  } else { 
    NN_x_new <- NULL
    x_approx <- add_pred_to_approx(object$x_approx, x_new, m)
  }
  
  if (cores == 1) { # run serial for loop
    
    mu_t <- matrix(nrow = n_new, ncol = object$nmcmc)
    if (lite) {
      s2_sum <- rep(0, times = n_new)
      if (return_all) s2_t <- matrix(nrow = n_new, ncol = object$nmcmc)
    } else sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    if (!is.null(entropy_limit)) ent_sum <- rep(0, times = n_new)
      
    for (t in 1:object$nmcmc) {
      if (sep) theta <- object$theta[t, ] else theta <- object$theta[t]
      k <- krig_vec(object$y, theta, object$g[t], object$tau2[t],
                    s2 = lite, sigma = !lite, v = object$v, m = m,
                    x = object$x, x_new = x_new, NNarray_pred = NN_x_new,
                    approx = x_approx, sep = sep)
      mu_t[, t] <- k$mean
      if (lite) {
        s2_sum <- s2_sum + k$s2
        if (return_all) s2_t[, t] <- k$s2
      } else sigma_sum <- sigma_sum + k$sigma
      if (!is.null(entropy_limit)) {
        if (lite) {
          sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
        ent_sum <- ent_sum + calc_entropy(k$mean, sig2, entropy_limit)
      }
    } # end of t for loop
    
  } else { # run in parallel using foreach
    
    iters <- 1:object$nmcmc
    chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
    if (cores > detectCores()) warning("cores is greater than available nodes")
    
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    thread <- NULL
    result <- foreach(thread = 1:cores) %dopar% {
      out <- list()
      out$mu_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      if (lite) {
        out$s2_sum <- rep(0, times = n_new)
        if (return_all) out$s2_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      } else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
      if (!is.null(entropy_limit)) out$ent_sum <- rep(0, times = n_new)
    
      # calculate predictions for each candidate MCMC iterations
      j <- 1
      for (t in chunks[[thread]]) {
        if (sep) theta <- object$theta[t, ] else theta <- object$theta[t]
        k <- krig_vec(object$y, theta, object$g[t], object$tau2[t],
                      s2 = lite, sigma = !lite, v = object$v, m = m,
                      x = object$x, x_new = x_new, NNarray_pred = NN_x_new,
                      approx = x_approx, sep = sep)
        out$mu_t[, j] <- k$mean
        if (lite) {
          out$s2_sum <- out$s2_sum + k$s2
          if (return_all) out$s2_t[, j] <- k$s2
        } else out$sigma_sum <- out$sigma_sum + k$sigma
        if (!is.null(entropy_limit)) {
          if (lite) {
            sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
          } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
          out$ent_sum <- out$ent_sum + calc_entropy(k$mean, sig2, entropy_limit)
        }
        j <- j + 1
      } # end of t for loop
      return(out)
    } # end of foreach loop
    
    stopCluster(cl)
      
    # Group elements out of the list
    mu_t <- do.call(cbind, lapply(result, with, eval(parse(text = "mu_t"))))
    if (lite) {
      s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum"))))
      if (return_all) s2_t <- do.call(cbind, lapply(result, with, eval(parse(text = "s2_t"))))
    } else {
      sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum"))))
    }
    if (!is.null(entropy_limit)) ent_sum <- Reduce("+", lapply(result, with, 
                                                               eval(parse(text = "ent_sum"))))
    
  } # end of else statement

  # Add variables to the output list
  mu_cov <- cov(t(mu_t))
  object$mean <- rowMeans(mu_t)
  if (lite) { 
    object$s2 <- s2_sum / object$nmcmc + diag(mu_cov)
    if (return_all) {
      object$mean_all <- mu_t
      object$s2_all <- s2_t
    }
  } else object$Sigma <- sigma_sum / object$nmcmc + mu_cov
  if (!is.null(entropy_limit)) object$entropy <- ent_sum / object$nmcmc
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  
  return(object)
}

# Predict Two Layer Vecchia ---------------------------------------------------
#' @rdname predict
#' @export
#' 
predict.dgp2vec <- function(object, x_new, m = object$m, lite = TRUE, 
                            store_latent = FALSE, mean_map = TRUE, 
                            return_all = FALSE, entropy_limit = NULL,
                            cores = 1, ...) {
  
  tic <- proc.time()[3]
  object <- clean_prediction(object) # remove previous predictions if present
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  D <- ncol(object$w[[1]])
  if (!mean_map) stop("mean_map = FALSE has not yet been implemented")
  if (return_all & !lite) stop("return_all only offered when lite = TRUE")
  if (!is.null(entropy_limit) & !is.numeric(entropy_limit))
    stop("entropy_limit must be numeric")

  # Pre-calculate nearest neighbors for x 
  if (lite) {
    NN_x_new <- FNN::get.knnx(object$x, x_new, m)$nn.index
    x_approx <-  NULL
  } else { 
    NN_x_new <- NULL
    x_approx <- add_pred_to_approx(object$x_approx, x_new, m)
  }
  
  if (cores == 1) { # run serial for loop

    if (store_latent) w_new_list <- list()
    mu_t <- matrix(nrow = n_new, ncol = object$nmcmc)
    if (lite) {
      s2_sum <- rep(0, times = n_new)
      if (return_all) s2_t <- matrix(nrow = n_new, ncol = object$nmcmc)
    } else sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    if (!is.null(entropy_limit)) ent_sum <- rep(0, times = n_new)

    for (t in 1:object$nmcmc) {
      w_t <- object$w[[t]]
        
      # Map x_new to w_new (separately for each dimension)
      w_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        if (object$settings$pmx) {
          prior_mean_new <- x_new[, i]
          prior_mean <- object$x[, i]
        } else {
          prior_mean_new <- 0
          prior_mean <- rep(0, length(object$y))
        }
        k <- krig_vec(w_t[, i], object$theta_w[t, i], g = eps, 
                      v = object$v, m = m, x = object$x, x_new = x_new,
                      NNarray_pred = NN_x_new,
                      prior_mean = prior_mean, prior_mean_new = prior_mean_new,
                      tau2 = object$settings$tau2)
        w_new[, i] <- k$mean
      } # end of i for loop
      if (store_latent) w_new_list[[t]] <- w_new
        
      if (lite) {
        w_approx <- NULL
      } else {
        w_approx <- update_obs_in_approx(object$w_approx, w_t)
        w_approx <- add_pred_to_approx(w_approx, w_new, m)
      }
        
      # Map w_new to y
      k <- krig_vec(object$y, object$theta_y[t], object$g[t], object$tau2[t],
                    s2 = lite, sigma = !lite, v = object$v, m = m,
                    x = w_t, x_new = w_new, approx = w_approx)
      mu_t[, t] <- k$mean
      if (lite) {
        s2_sum <- s2_sum + k$s2
        if (return_all) s2_t[, t] <- k$s2
      } else sigma_sum <- sigma_sum + k$sigma
      if (!is.null(entropy_limit)) {
        if (lite) {
          sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
        ent_sum <- ent_sum + calc_entropy(k$mean, sig2, entropy_limit)
      }
    } # end of t for loop
    
  } else { # run in parallel using foreach
    
    iters <- 1:object$nmcmc
    chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
    if (cores > detectCores()) warning('cores is greater than available nodes')
    
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    thread <- NULL
    result <- foreach(thread = 1:cores) %dopar% {
      out <- list()
      if (store_latent) out$w_new <- list()
      out$mu_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      if (lite) {
        out$s2_sum <- rep(0, times = n_new)
        if (return_all) out$s2_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      } else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
      if (!is.null(entropy_limit)) out$ent_sum <- rep(0, times = n_new)
  
      j <- 1
      for (t in chunks[[thread]]) {
        w_t <- object$w[[t]]
        
        # Map x_new to w_new (separately for each dimension)
        w_new <- matrix(nrow = n_new, ncol = D)
        for (i in 1:D) {
          if (object$settings$pmx) {
            prior_mean_new <- x_new[, i]
            prior_mean <- object$x[, i]
          } else {
            prior_mean_new <- 0
            prior_mean <- rep(0, length(object$y))
          }
          k <- krig_vec(w_t[, i], object$theta_w[t, i], g = eps, 
                        v = object$v, m = m, x = object$x, x_new = x_new,
                        NNarray_pred = NN_x_new, prior_mean = prior_mean, 
                        prior_mean_new = prior_mean_new,
                        tau2 = object$settings$tau2)
          w_new[, i] <- k$mean
        } # end of i for loop
        if (store_latent) out$w_new[[j]] <- w_new
        
        if (lite) {
          w_approx <- NULL
        } else {
          w_approx <- update_obs_in_approx(object$w_approx, w_t)
          w_approx <- add_pred_to_approx(w_approx, w_new, m)
        }
  
        # Map w_new to y
        k <- krig_vec(object$y, object$theta_y[t], object$g[t], object$tau2[t],
                      s2 = lite, sigma = !lite, v = object$v, m = m,
                      x = w_t, x_new = w_new, approx = w_approx)
        out$mu_t[, j] <- k$mean
        if (lite) {
          out$s2_sum <- out$s2_sum + k$s2
          if (return_all) out$s2_t[, j] <- k$s2
        } else out$sigma_sum <- out$sigma_sum + k$sigma
        if (!is.null(entropy_limit)) {
          if (lite) {
            sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
          } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
          out$ent_sum <- out$ent_sum + calc_entropy(k$mean, sig2, entropy_limit)
        }
        j <- j + 1
      } # end of t for loop
      return(out)
    } # end of foreach statement
      
    stopCluster(cl)
    
    # Group elements out of the list
    mu_t <- do.call(cbind, lapply(result, with, eval(parse(text = "mu_t"))))
    if (lite) {
      s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum"))))
      if (return_all) s2_t <- do.call(cbind, lapply(result, with, eval(parse(text = "s2_t"))))
    } else {
      sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum"))))
    }
    if (store_latent) w_new_list <- unlist(lapply(result, with, eval(parse(text = "w_new"))), 
                                                  recursive = FALSE) 
    if (!is.null(entropy_limit)) ent_sum <- Reduce("+", lapply(result, with, 
                                                               eval(parse(text = "ent_sum"))))
    
  } # end of else statement
  
  # Add variables to the output list
  mu_cov <- cov(t(mu_t))
  object$mean <- rowMeans(mu_t)
  if (store_latent) object$w_new <- w_new_list
  if (lite) {
    object$s2 <- s2_sum / object$nmcmc + diag(mu_cov)
    if (return_all) {
      object$mean_all <- mu_t
      object$s2_all <- s2_t
    }
  } else object$Sigma <- sigma_sum / object$nmcmc + mu_cov
  if (!is.null(entropy_limit)) object$entropy <- ent_sum / object$nmcmc
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  
  return(object)
}

# Predict Three Layer Vecchia -------------------------------------------------
#' @rdname predict
#' @export

predict.dgp3vec <- function(object, x_new, m = object$m, lite = TRUE, 
                            store_latent = FALSE, mean_map = TRUE,
                            return_all = FALSE, entropy_limit = NULL,
                            cores = 1, ...) {
  
  tic <- proc.time()[3]
  object <- clean_prediction(object) # remove previous predictions if present
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  D <- ncol(object$z[[1]])
  if (!mean_map) stop("mean_map = FALSE has not yet been implemented")
  if (return_all & !lite) stop("return_all only offered when lite = TRUE")
  if (!is.null(entropy_limit) & !is.numeric(entropy_limit))
    stop("entropy_limit must be numeric")
  
  # Pre-calculate nearest neighbors for x 
  if (lite) {
    NN_x_new <- FNN::get.knnx(object$x, x_new, m)$nn.index
    x_approx <-  NULL
  } else { 
    NN_x_new <- NULL
    x_approx <- add_pred_to_approx(object$x_approx, x_new, m)
  }
  
  if (cores == 1) { # run serial for loop 
    
    if (store_latent) {
      z_new_list <- list()
      w_new_list <- list()
    }
    mu_t <- matrix(nrow = n_new, ncol = object$nmcmc)
    if (lite) {
      s2_sum <- vector(length = n_new)
      if (return_all) s2_t <- matrix(nrow = n_new, ncol = object$nmcmc)
    } else sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    if (!is.null(entropy_limit)) ent_sum <- rep(0, times = n_new)
      
    for (t in 1:object$nmcmc) {
      z_t <- object$z[[t]]
      w_t <- object$w[[t]]
        
      # Map x_new to z_new (separately for each dimension)
      z_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) { # mean_map = TRUE only
        k <- krig_vec(z_t[, i], object$theta_z[t, i], g = eps, tau2 = 1,
                      v = object$v, m = m, x = object$x, x_new = x_new,
                      NNarray_pred = NN_x_new)
        z_new[, i] <- k$mean
      } # end of i for loop
      if (store_latent) z_new_list[[t]] <- z_new
        
      # Map z_new to w_new (separately for each dimension)
      w_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        k <- krig_vec(w_t[, i], object$theta_w[t, i], g = eps, tau2 = 1,
                      v = object$v, m = m, x = z_t, x_new = z_new)
        w_new[, i] <- k$mean
      } # end of i for loop
      if (store_latent) w_new_list[[t]] <- w_new
        
      if (lite) {
        w_approx <- NULL
      } else {
        w_approx <- update_obs_in_approx(object$w_approx, w_t)
        w_approx <- add_pred_to_approx(w_approx, w_new, m)
      }
        
      # Map w_new to y
      k <- krig_vec(object$y, object$theta_y[t], object$g[t], object$tau2[t],
                    s2 = lite, sigma = !lite, v = object$v, m = m,
                    x = w_t, x_new = w_new, approx = w_approx)
      mu_t[, t] <- k$mean
      if (lite) {
        s2_sum <- s2_sum + k$s2
        if (return_all) s2_t[, t] <- k$s2
      } else sigma_sum <- sigma_sum + k$sigma
      if (!is.null(entropy_limit)) {
        if (lite) {
          sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
        ent_sum <- ent_sum + calc_entropy(k$mean, sig2, entropy_limit)
      }
    } # end of t for loop

  } else { # run in parallel using foreach
  
    iters <- 1:object$nmcmc
    chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
    if (cores > detectCores()) warning('cores is greater than available nodes')
    
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    thread <- NULL
    result <- foreach(thread = 1:cores) %dopar% {
      out <- list()
      if (store_latent) {
        out$z_new <- list()
        out$w_new <- list()
      }
      out$mu_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      if (lite) {
        out$s2_sum <- vector(length = n_new)
        if (return_all) out$s2_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      } else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
      if (!is.null(entropy_limit)) ent_sum <- rep(0, times = n_new)
    
      j <- 1
      for (t in chunks[[thread]]) {
        z_t <- object$z[[t]]
        w_t <- object$w[[t]]
      
        # Map x_new to z_new (separately for each dimension)
        z_new <- matrix(nrow = n_new, ncol = D)
        for (i in 1:D) { # mean_map = TRUE only
          k <- krig_vec(z_t[, i], object$theta_z[t, i], g = eps, tau2 = 1,
                        v = object$v, m = m, x = object$x, x_new = x_new,
                        NNarray_pred = NN_x_new)
          z_new[, i] <- k$mean
        } # end of i for loop
        if (store_latent) out$z_new[[j]] <- z_new
      
        # Map z_new to w_new (separately for each dimension)
        w_new <- matrix(nrow = n_new, ncol = D)
        for (i in 1:D) {
          k <- krig_vec(w_t[, i], object$theta_w[t, i], g = eps, tau2 = 1,
                        v = object$v, m = m, x = z_t, x_new = z_new)
          w_new[, i] <- k$mean
        } # end of i for loop
        if (store_latent) out$w_new[[j]] <- w_new
      
        if (lite) {
          w_approx <- NULL
        } else {
          w_approx <- update_obs_in_approx(object$w_approx, w_t)
          w_approx <- add_pred_to_approx(w_approx, w_new, m)
        }
  
        # Map w_new to y
        k <- krig_vec(object$y, object$theta_y[t], object$g[t], object$tau2[t],
                      s2 = lite, sigma = !lite, v = object$v, m = m,
                      x = w_t, x_new = w_new, approx = w_approx)
        out$mu_t[, j] <- k$mean
        if (lite) {
          out$s2_sum <- out$s2_sum + k$s2
          if (return_all) out$s2_t[, j] <- k$s2
        } else out$sigma_sum <- out$sigma_sum + k$sigma
        if (!is.null(entropy_limit)) {
          if (lite) {
            sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
          } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
          out$ent_sum <- out$ent_sum + calc_entropy(k$mean, sig2, entropy_limit)
        }
        j <- j + 1
      } # end of t for loop
      return(out)
    } # end of foreach loop
    
    stopCluster(cl)
    
    # Group elements out of the list
    mu_t <- do.call(cbind, lapply(result, with, eval(parse(text = "mu_t"))))
    if (lite) {
      s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum"))))
      if (return_all) s2_t <- do.call(cbind, lapply(result, with, eval(parse(text = "s2_t"))))
    } else {
      sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum"))))
    }
    if (store_latent) {
      z_new_list <- unlist(lapply(result, with, eval(parse(text = "z_new"))), recursive = FALSE)
      w_new_list <- unlist(lapply(result, with, eval(parse(text = "w_new"))), recursive = FALSE)
    }
    if (!is.null(entropy_limit)) ent_sum <- Reduce("+", lapply(result, with, 
                                                               eval(parse(text = "ent_sum"))))
    
  } # end of else statement
  
  # Add variables to the output list
  mu_cov <- cov(t(mu_t))
  object$mean <- rowMeans(mu_t)
  if (store_latent) {
    object$z_new <- z_new_list
    object$w_new <- w_new_list
  }
  if (lite) {
    object$s2 <- s2_sum / object$nmcmc + diag(mu_cov)
    if (return_all) {
      object$mean_all <- mu_t
      object$s2_all <- s2_t
    }
  } else object$Sigma <- sigma_sum / object$nmcmc + mu_cov
  if (!is.null(entropy_limit)) object$entropy <- ent_sum / object$nmcmc
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  
  return(object)
}