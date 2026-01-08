
# Function Contents -----------------------------------------------------------
# External: (see documentation in predict.R file)
#   predict.gpvec
#   predict.dgp2vec
#   predict.dgp3vec
# Internal:
#   predict_shallow_vec: predictions for one layer Vecchia models
#   predict_deep_vec: predictions for two and three layer Vecchia models

# predict.gpvec ---------------------------------------------------------------
#' @rdname predict
#' @export

predict.gpvec <- function(object, x_new, m = NULL,
                          ord_new = NULL, lite = TRUE, grad = FALSE, 
                          return_all = FALSE, EI = FALSE, entropy_limit = NULL, 
                          cores = 1, ...) {
  
  if (return_all & !lite) 
    stop("return_all only offered when lite = TRUE")
  if (!is.null(entropy_limit) & !is.numeric(entropy_limit))
    stop("entropy_limit must be numeric")
  cores <- check_cores(cores, object$nmcmc)
  if (grad) { 
    if (object$v != 999) stop("grad only offered with cov = 'exp2'")
    if (!lite) stop("grad = TRUE requires lite = TRUE")
  }
  
  settings <- list(ord_new = ord_new, lite = lite, grad = grad,
                   return_all = return_all, EI = EI, 
                   entropy_limit = entropy_limit, cores = cores)
  object <- predict_shallow_vec(object, x_new, m, settings)
  return(object)
}

# predict.dgp2vec -------------------------------------------------------------
#' @rdname predict
#' @export
#' 
predict.dgp2vec <- function(object, x_new, m = NULL,
                            ord_new = NULL, lite = TRUE, grad = FALSE,
                            store_latent = FALSE, mean_map = TRUE, 
                            return_all = FALSE, EI = FALSE, 
                            entropy_limit = NULL, cores = 1, ...) {
  if (return_all & !lite) 
    stop("return_all only offered when lite = TRUE")
  if (!is.null(entropy_limit) & !is.numeric(entropy_limit))
    stop("entropy_limit must be numeric")
  cores <- check_cores(cores, object$nmcmc)
  if (grad) { 
    if (object$v != 999) stop("grad only offered with cov = 'exp2'")
    if (!lite) stop("grad = TRUE requires lite = TRUE")
    if (object$settings$monowarp) stop("grad not offered for monowarp = TRUE")
  }
  
  settings <- list(ord_new = ord_new, lite = lite, grad = grad,
                   store_latent = store_latent, mean_map = mean_map,
                   return_all = return_all, EI = EI, 
                   entropy_limit = entropy_limit, cores = cores)
  object <- predict_deep_vec(object, x_new, m, settings, layers = 2)
  return(object)
}
  
# predict.dgp3vec -------------------------------------------------------------
#' @rdname predict
#' @export
  
predict.dgp3vec <- function(object, x_new, m = NULL, 
                            ord_new = NULL, lite = TRUE,
                            store_latent = FALSE, 
                            mean_map = TRUE, return_all = FALSE, EI = FALSE, 
                            entropy_limit = NULL, cores = 1, ...) {
  
  if (return_all & !lite) 
    stop("return_all only offered when lite = TRUE")
  if (!is.null(entropy_limit) & !is.numeric(entropy_limit))
    stop("entropy_limit must be numeric")
  cores <- check_cores(cores, object$nmcmc)
  
  settings <- list(ord_new = ord_new, lite = lite, grad = FALSE,
                   store_latent = store_latent, mean_map = mean_map,
                   return_all = return_all, EI = EI, 
                   entropy_limit = entropy_limit, cores = cores)
  object <- predict_deep_vec(object, x_new, m, settings, layers = 3)
  return(object)
}

# predict_shallow_vec ---------------------------------------------------------

predict_shallow_vec <- function(object, x_new, m, settings, samples_only = FALSE) { 
  # One-layer only
  # If samples_only = TRUE, posterior sample draws are returned instead of
  # summarized posterior moments
  
  tic <- proc.time()[[3]]
  if (is.vector(x_new)) x_new <- as.matrix(x_new)
  if (ncol(x_new) != ncol(object$x)) stop("dimension of x_new does not match dimension of x")
  object$x_new <- x_new
  n_new <- nrow(x_new)
  d <- ncol(object$x) 
  sep <- is.matrix(object$theta)

  grad_enhance <- !is.null(object$dydx)
  if (grad_enhance) { 
    if (settings$EI) stop("EI not implemented for grad_enhance")
    if (!is.null(settings$entropy_limit)) stop("entropy_limit not implemented for grad_enhance")
    y <- c(object$y, as.vector(object$dydx))
  } else y <- object$y
  if (is.null(m)) {
    if (settings$lite) {
      m <- min(nrow(object$x_approx$x_ord), 2*object$x_approx$m)
    } else m <- min(nrow(object$x_approx$x_ord) + n_new - 1, 2*object$x_approx$m)
  }

  if (settings$EI) y_min <- min(object$y) # use smallest observed value
  
  if (!is.null(settings$ord_new)) {
    if (settings$lite) message("ord_new is only relevant when lite = FALSE")
    if (length(settings$ord_new) != n_new | 
        min(settings$ord_new) != 1 | 
        max(settings$ord_new) != n_new | 
        sum(duplicated(settings$ord_new)) > 0) {
      stop("ord_new must contain 1:nrow(x_new) without any duplicates")
    }
  }
  
  object$x_approx <- clean_pred_from_approx(object$x_approx)
  object$x_approx <- add_pred_to_approx(object$x_approx, x_new, m, lite = settings$lite,
                                        grad = settings$grad, ord_new = settings$ord_new)
  
  if (settings$cores == 1) { # run serial for loop
    
    if (samples_only) {
      samples <- matrix(nrow = settings$nper*object$nmcmc, ncol = n_new)
      if (settings$grad) grad_samples <- array(dim = c(settings$nper*object$nmcmc, n_new, d))
    } else {
      mu_t <- matrix(nrow = object$nmcmc, ncol = n_new)
      if (settings$lite) {
        s2_sum <- rep(0, times = n_new)
        if (settings$return_all) s2_t <- matrix(nrow = object$nmcmc, ncol = n_new)
      } else sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
      if (settings$grad) { # no return_all or lite = FALSE options
        grad_mu_t <- array(dim = c(object$nmcmc, n_new, d))
        grad_s2_sum <- matrix(0, nrow = n_new, ncol = d)
      }
      if (settings$EI) ei_sum <- rep(0, times = n_new)
      if (!is.null(settings$entropy_limit)) ent_sum <- rep(0, times = n_new)
    }
    
    for (t in 1:object$nmcmc) {
      g <- ifel(length(object$g) == 1, object$g, object$g[t])
      
      k <- krig_vec(y, approx = object$x_approx, 
                    tau2 = object$tau2[t], 
                    theta = ifel(sep, object$theta[t, ], object$theta[t]),
                    g = g, 
                    v = object$v,
                    sep = sep,
                    s2 = ifel(samples_only, FALSE, settings$lite),
                    sigma = ifel(samples_only, FALSE, !settings$lite),
                    grad = settings$grad,
                    nsamples = ifel(samples_only, settings$nper, 0))
      
      if (samples_only) { 
        row_indx <- (t*settings$nper - settings$nper + 1):(t*settings$nper)
        samples[row_indx, ] <- k$samples
        if (settings$grad) grad_samples[row_indx, , ] <- k$grad_samples
      } else {
        mu_t[t, ] <- k$mean
        if (settings$lite) {
          s2_sum <- s2_sum + k$s2
          if (settings$return_all) s2_t[t, ] <- k$s2
        } else sigma_sum <- sigma_sum + k$sigma
        if (settings$grad) {
          grad_mu_t[t, , ] <- k$grad_mean
          grad_s2_sum <- grad_s2_sum + k$grad_s2
        }
        if (settings$EI) {
          if (settings$lite) {
            sig2 <- k$s2 - (object$tau2[t] * g) 
          } else sig2 <- diag(k$sigma) - (object$tau2[t] * g)
          ei_sum <- ei_sum + exp_improv(k$mean, sig2, y_min)
        }
        if (!is.null(settings$entropy_limit)) {
          if (settings$lite) {
            sig2 <- k$s2 - (object$tau2[t] * g) 
          } else sig2 <- diag(k$sigma) - (object$tau2[t] * g)
          ent_sum <- ent_sum + calc_entropy(k$mean, sig2, settings$entropy_limit)
        }
      }
    } # end of t for loop
    
  } else { # run in parallel using foreach
    
    iters <- 1:object$nmcmc
    chunks <- split(iters, sort(cut(iters, settings$cores, labels = FALSE)))
    cl <- makeCluster(settings$cores)
    registerDoParallel(cl)
    
    thread <- NULL
    result <- foreach(thread = 1:settings$cores) %dopar% {
      out <- list()
      if (samples_only) {
        out$samples <- matrix(nrow = settings$nper*length(chunks[[thread]]), ncol = n_new)
        if (settings$grad) out$grad_samples <- array(dim = c(settings$nper*length(chunks[[thread]]), n_new, d))
      } else {
        out$mu_t <- matrix(nrow = length(chunks[[thread]]), ncol = n_new)
        if (settings$lite) {
          out$s2_sum <- rep(0, times = n_new)
          if (settings$return_all) 
            out$s2_t <- matrix(nrow = length(chunks[[thread]]), ncol = n_new)
        } else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
        if (settings$grad) {
          out$grad_mu_t <- array(dim = c( length(chunks[[thread]]), n_new, d))
          out$grad_s2_sum <- matrix(0, nrow = n_new, ncol = d)
        }
        if (settings$EI) out$ei_sum <- rep(0, times = n_new)
        if (!is.null(settings$entropy_limit)) out$ent_sum <- rep(0, times = n_new)
      }
      
      j <- 1
      for (t in chunks[[thread]]) {
        g <- ifel(length(object$g) == 1, object$g, object$g[t])
        
        k <- krig_vec(y, approx = object$x_approx, 
                      tau2 = object$tau2[t], 
                      theta = ifel(sep, object$theta[t, ], object$theta[t]),
                      g = g, 
                      v = object$v,
                      sep = sep,
                      s2 = ifel(samples_only, FALSE, settings$lite),
                      sigma = ifel(samples_only, FALSE, !settings$lite),
                      grad = settings$grad,
                      nsamples = ifel(samples_only, settings$nper, 0))
        
        if (samples_only) {
          row_indx <- (j*settings$nper - settings$nper + 1):(j*settings$nper)
          out$samples[row_indx, ] <- k$samples
          if (settings$grad) out$grad_samples[row_indx, , ] <- k$grad_samples
        } else {
          out$mu_t[j, ] <- k$mean
          if (settings$lite) {
            out$s2_sum <- out$s2_sum + k$s2
            if (settings$return_all) out$s2_t[j, ] <- k$s2
          } else out$sigma_sum <- out$sigma_sum + k$sigma
          if (settings$grad) {
            out$grad_mu_t[j, , ] <- k$grad_mean
            out$grad_s2_sum <- out$grad_s2_sum + k$grad_s2
          }
          if (settings$EI) {
            if (settings$lite) {
              sig2 <- k$s2 - (object$tau2[t] * g) 
            } else sig2 <- diag(k$sigma) - (object$tau2[t] * g)
            out$ei_sum <- out$ei_sum + exp_improv(k$mean, sig2, y_min)
          }
          if (!is.null(settings$entropy_limit)) {
            if (settings$lite) {
              sig2 <- k$s2 - (object$tau2[t] * g) 
            } else sig2 <- diag(k$sigma) - (object$tau2[t] * g)
            out$ent_sum <- out$ent_sum + calc_entropy(k$mean, sig2, settings$entropy_limit)
          }
        }
        j <- j + 1
      } # end of t for loop
      return(out)
    } # end of foreach loop
    
    stopCluster(cl)
    
    # Group elements out of the list
    if (samples_only) {
      samples <- do.call(rbind, lapply(result, with, eval(parse(text = "samples"))))
      if (settings$grad) grad_samples <- abind::abind(lapply(result, with, 
                                                             eval(parse(text = "grad_samples"))), along = 1)
    } else {
      mu_t <- do.call(rbind, lapply(result, with, eval(parse(text = "mu_t"))))
      if (settings$lite) {
        s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum"))))
        if (settings$return_all) s2_t <- do.call(rbind, lapply(result, with, eval(parse(text = "s2_t"))))
      } else {
        sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum"))))
      }
      if (settings$grad) {
        grad_mu_t <- abind::abind(lapply(result, with, eval(parse(text = "grad_mu_t"))), along = 1)
        grad_s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "grad_s2_sum"))))
      }
      if (settings$EI) ei_sum <- Reduce("+", lapply(result, with, eval(parse(text = "ei_sum"))))
      if (!is.null(settings$entropy_limit)) ent_sum <- Reduce("+", lapply(result, with, 
                                                                          eval(parse(text = "ent_sum"))))
    }
    
  } # end of else statement
  
  # Add variables to the output list
  if (samples_only) {
    if (settings$grad) {
      return(list(y = samples, dy = ifel(d == 1, grad_samples[, , 1], grad_samples)))
    } else return(samples)
  } else {
    object$mean <- colMeans(mu_t)
    if (settings$lite) { 
      object$s2 <- s2_sum / object$nmcmc + ifel(object$nmcmc > 1, apply(mu_t, 2, var), 0)
      if (settings$return_all) {
        object$mean_all <- mu_t
        object$s2_all <- s2_t
      }
    } else object$Sigma <- sigma_sum / object$nmcmc + ifel(object$nmcmc > 1, cov(mu_t), 0)
    if (settings$grad) { # TODO: speed up by removing this for loop?
      object$grad_mean <- matrix(nrow = n_new, ncol = d)
      object$grad_s2 <- grad_s2_sum / object$nmcmc
      for (i in 1:d) {
        object$grad_mean[, i] <- colMeans(grad_mu_t[, , i])
        object$grad_s2[, i] <- object$grad_s2[, i] + ifel(object$nmcmc > 1, apply(grad_mu_t[, , i], 2, var), 0)
      }
    }
    if (settings$EI) object$EI <- ei_sum / object$nmcmc
    if (!is.null(settings$entropy_limit)) object$entropy <- drop(ent_sum / object$nmcmc)
    toc <- proc.time()[[3]]
    object$time <- object$time + unname(toc - tic)
    return(object)
  }
}

# predict_deep_vec ------------------------------------------------------------

predict_deep_vec <- function(object, x_new, m, settings, layers, samples_only = FALSE) {
  # Two-layer or three-layer
  # If samples_only = TRUE, posterior sample draws are returned instead of
  # summarized posterior moments

  tic <- proc.time()[[3]]
  if (is.vector(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  n <- length(object$y)
  d <- ncol(object$x)
  D <- dim(object$w)[3]
  
  grad_enhance <- !is.null(object$dydx)
  if (grad_enhance) {
    if (settings$EI) stop("EI not implemented for grad_enhance")
    if (!is.null(settings$entropy_limit)) stop("entropy_limit not implemented for grad_enhance")
    # y will be recalculated for every w_t
  } else y <- object$y
  if (is.null(m)) {
    if (settings$lite) {
      m <- min(nrow(object$w_approx$x_ord), 2*object$w_approx$m)
    } else m <- min(nrow(object$w_approx$x_ord) + n_new - 1, 2*object$w_approx$m)
  }

  # Prespecify prior mean (if not zero)
  if (object$settings$pmx) {
    if (grad_enhance) {
      w_prior_mean <- object$settings$w_prior_mean
    } else w_prior_mean <- object$x
    if (settings$grad) {
      w_prior_mean_new <- get_prior_mean(x_new)
    } else w_prior_mean_new <- x_new
  }
  
  if (settings$EI) y_min <- min(object$y) # use smallest observed value
  
  if (!is.null(settings$ord_new)) {
    if (settings$lite) message("ord_new is only relevant when lite = FALSE")
    if (length(settings$ord_new) != n_new | 
        min(settings$ord_new) != 1 | 
        max(settings$ord_new) != n_new | 
        sum(duplicated(settings$ord_new)) > 0) {
      stop("ord_new must contain 1:nrow(x_new) without any duplicates")
    }
  }
  
  # Add x_new to x_approx (w and potentially z will be added later after mapping)  
  if (!object$settings$monowarp) {
    object$x_approx <- clean_pred_from_approx(object$x_approx)
    object$x_approx <- add_pred_to_approx(object$x_approx, x_new, 
                                          m = ifel(settings$mean_map, min(n, m), m), 
                                          lite = settings$mean_map,
                                          grad = settings$grad, 
                                          ord_new = settings$ord_new)
  }
  
  if (settings$cores == 1) { # run serial for loop
    
    if (samples_only) {
      samples <- matrix(nrow = settings$nper*object$nmcmc, ncol = n_new)
      if (settings$grad) grad_samples <- array(dim = c(settings$nper*object$nmcmc, n_new, d))
    } else {
      mu_t <- matrix(nrow = object$nmcmc, ncol = n_new)
      if (settings$lite) {
        s2_sum <- rep(0, times = n_new)
        if (settings$return_all) s2_t <- matrix(nrow = object$nmcmc, ncol = n_new)
      } else sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
      if (settings$grad) {
        grad_mu_t <- array(dim = c(object$nmcmc, n_new, d))
        grad_s2_sum <- matrix(0, nrow = n_new, ncol = d)
      }
      if (settings$EI) ei_sum <- rep(0, times = n_new)
      if (!is.null(settings$entropy_limit)) ent_sum <- rep(0, times = n_new)
      if (settings$store_latent) {
        w_new_store <- array(dim = c(object$nmcmc, n_new, D))
        if (layers == 3) z_new_store <- array(dim = c(object$nmcmc, n_new, D))
      }
    }
      
    for (t in 1:object$nmcmc) {
      g <- ifel(length(object$g) == 1, object$g, object$g[t])
      
      if (layers == 3) { # no pmx, monowarp, grad, or grad_enhance option

        # First map x_new to z_new
        z_t <- as.matrix(object$z[t, , ])
        z_new <- matrix(nrow = n_new, ncol = D)
        for (i in 1:D) {
          k <- krig_vec(z_t[, i], 
                        approx = object$x_approx,
                        tau2 = object$settings$tau2_z,
                        theta = object$theta_z[t, i], 
                        g = eps, 
                        v = object$v,
                        nsamples = ifel(settings$mean_map, 0, 1))
          if (settings$mean_map) {
            z_new[, i] <- k$mean
          } else z_new[, i] <- k$samples
        } # end of i for loop
        object$z_approx <- clean_pred_from_approx(object$z_approx)
        object$z_approx$x_ord <- z_t[object$z_approx$ord, , drop = FALSE]
        object$z_approx <- add_pred_to_approx(object$z_approx, z_new, m, 
                                              lite = settings$mean_map, grad = FALSE,
                                              ord_new = settings$ord_new)
        if (settings$store_latent) z_new_store[t, , ] <- z_new

        # Then map z_new to w_new
        w_t <- as.matrix(object$w[t, , ])
        w_new <- matrix(nrow = n_new, ncol = D)
        for (i in 1:D) {
          k <- krig_vec(w_t[, i], 
                        approx = object$z_approx,
                        tau2 = object$settings$tau2_w,
                        theta = object$theta_w[t, i], 
                        g = eps,
                        v = object$v,
                        nsamples = ifel(settings$mean_map, 0, 1))
          if (settings$mean_map) {
            w_new[, i] <- k$mean
          } else w_new[, i] <- k$samples
        } # end of i for loop
        object$w_approx <- clean_pred_from_approx(object$w_approx)
        object$w_approx$x_ord <- w_t[object$w_approx$ord, , drop = FALSE]
        object$w_approx <- add_pred_to_approx(object$w_approx, w_new, m, 
                                              lite = settings$lite, grad = FALSE,
                                              ord_new = settings$ord_new)
        if (settings$store_latent) w_new_store[t, , ] <- w_new

      } else if (layers == 2) {

        w_t <- as.matrix(object$w[t, , ]) # includes gradients if grad_enhance = TRUE
        if (object$settings$monowarp) { 
          w_new <- monotransform(x_new, object$x_grid, object$w_grid[t, , ])
        } else { # use kriging
          w_new <- matrix(nrow = n_new, ncol = D)
          if (settings$grad) dwdx <- array(dim = c(n_new, D, d))
          for (i in 1:D) { # no grad_enhance option here
            k <- krig_vec(w_t[, i], 
                          approx = object$x_approx,
                          tau2 = object$settings$tau2_w,
                          theta = object$theta_w[t, i],
                          g = eps, 
                          v = object$v, 
                          grad = settings$grad,
                          nsamples = ifel(settings$mean_map, 0, 1),
                          prior_mean = ifel(object$settings$pmx, w_prior_mean[, i], 0),
                          prior_mean_new = ifel(object$settings$pmx, w_prior_mean_new[, i], 0))
            if (settings$mean_map) {
              w_new[, i] <- k$mean
              if (settings$grad) dwdx[, i, ] <- k$grad_mean
            } else {
              w_new[, i] <- k$samples
              if (settings$grad) dwdx[, i, ] <- k$grad_samples[1, , ] # only one sample returned
            }
          } # end of i for loop
        } # end of else statement
        object$w_approx <- clean_pred_from_approx(object$w_approx)
        object$w_approx$x_ord <- w_t[object$w_approx$ord, , drop = FALSE]
        object$w_approx <- add_pred_to_approx(object$w_approx, w_new, m, 
                                              lite = settings$lite, grad = settings$grad,
                                              ord_new = settings$ord_new)
        if (settings$store_latent) w_new_store[t, , ] <- w_new
      
      } # end of layers == 2 else statement
      
      # Finally: map w_new to y
      if (grad_enhance) {
        dydw <- get_dydw(w_t, object$dydx)
        y <- c(object$y, as.vector(dydw))
      }
      k <- krig_vec(y,
                    approx = object$w_approx,
                    tau2 = object$tau2_y[t],
                    theta = object$theta_y[t], 
                    g = g, 
                    v = object$v,
                    s2 = ifel(samples_only, FALSE, settings$lite), 
                    sigma = ifel(samples_only, FALSE, !settings$lite),
                    grad = settings$grad,
                    nsamples = ifel(samples_only, settings$nper, 0))
      
      if (samples_only) {
        row_indx <- (t*settings$nper - settings$nper + 1):(t*settings$nper)
        samples[row_indx, ] <- k$samples
        if (settings$grad) {
          for (i in 1:d) {
            grad_samples[row_indx, , i] <- 0
            for (j in 1:D) { # grad_samples are dy/dw, must apply multivariate chain rule to get dy/dx
              grad_samples[row_indx, , i] <- grad_samples[row_indx, , i] + k$grad_samples[, , i] * dwdx[, j, i]
            }
          }
        }
      } else {
        mu_t[t, ] <- k$mean
        if (settings$lite) {
          s2_sum <- s2_sum + k$s2
          if (settings$return_all) s2_t[t, ] <- k$s2
        } else sigma_sum <- sigma_sum + k$sigma
        if (settings$grad) { # MUST USE MULTIVARIATE CHAIN RULE
          for (i in 1:d) {
            grad_mu_t[t, , i] <- rowSums(k$grad_mean * dwdx[, , i])
            grad_s2_sum[, i] <- grad_s2_sum[, i] + rowSums(dwdx[, , i]^2 * k$grad_s2)
          }
        }
        if (settings$EI) {
          if (settings$lite) {
            sig2 <- k$s2 - (object$tau2[t] * g) 
          } else sig2 <- diag(k$sigma) - (object$tau2[t] * g)
          ei_sum <- ei_sum + exp_improv(k$mean, sig2, y_min)
        }
        if (!is.null(settings$entropy_limit)) {
          if (settings$lite) {
            sig2 <- k$s2 - (object$tau2[t] * g) 
          } else sig2 <- diag(k$sigma) - (object$tau2[t] * g)
          ent_sum <- ent_sum + calc_entropy(k$mean, sig2, settings$entropy_limit)
        }
      }
    } # end of t for loop
    
  } else { # run in parallel using foreach
    
    iters <- 1:object$nmcmc
    chunks <- split(iters, sort(cut(iters, settings$cores, labels = FALSE)))
    cl <- makeCluster(settings$cores)
    registerDoParallel(cl)
    
    thread <- NULL
    result <- foreach(thread = 1:settings$cores) %dopar% {
      out <- list()
      if (samples_only) {
        out$samples <- matrix(nrow = settings$nper*length(chunks[[thread]]), ncol = n_new)
        if (settings$grad) out$grad_samples <- array(dim = c(settings$nper*length(chunks[[thread]]), n_new, d))
      } else {
        out$mu_t <- matrix(nrow = length(chunks[[thread]]), ncol = n_new)
        if (settings$lite) {
          out$s2_sum <- rep(0, times = n_new)
          if (settings$return_all) 
            out$s2_t <- matrix(nrow = length(chunks[[thread]]), ncol = n_new)
        } else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
        if (settings$grad) {
          out$grad_mu_t <- array(dim = c(length(chunks[[thread]]), n_new, d))
          out$grad_s2_sum <- matrix(0, nrow = n_new, ncol = d)
        }
        if (settings$EI) out$ei_sum <- rep(0, times = n_new)
        if (!is.null(settings$entropy_limit)) out$ent_sum <- rep(0, times = n_new)
        if (settings$store_latent) {
          out$w_new <- array(dim = c(length(chunks[[thread]]), n_new, D))
          if (layers == 3) out$z_new <- array(dim = c(length(chunks[[thread]]), n_new, D))
        }
      }
      
      j <- 1
      for (t in chunks[[thread]]) {
        g <- ifel(length(object$g) == 1, object$g, object$g[t])
        
        if (layers == 3) { # no pmx, monowarp, grad, or grad_enhance option

          # First map x_new to z_new
          z_t <- as.matrix(object$z[t, , ])
          z_new <- matrix(nrow = n_new, ncol = D)
          for (i in 1:D) {
            k <- krig_vec(z_t[, i], 
                          approx = object$x_approx,
                          tau2 = object$settings$tau2_z,
                          theta = object$theta_z[t, i], 
                          g = eps, 
                          v = object$v,
                          nsamples = ifel(settings$mean_map, 0, 1))
            if (settings$mean_map) {
              z_new[, i] <- k$mean
            } else z_new[, i] <- k$samples
          } # end of i for loop
          object$z_approx <- clean_pred_from_approx(object$z_approx)
          object$z_approx$x_ord <- z_t[object$z_approx$ord, , drop = FALSE]
          object$z_approx <- add_pred_to_approx(object$z_approx, z_new, m, 
                                                lite = settings$mean_map, grad = FALSE,
                                                ord_new = settings$ord_new)
          if (settings$store_latent) out$z_new[j, , ] <- z_new

          # Then map z_new to w_new
          w_t <- as.matrix(object$w[t, , ])
          w_new <- matrix(nrow = n_new, ncol = D)
          for (i in 1:D) {
            k <- krig_vec(w_t[, i], 
                          approx = object$z_approx,
                          tau2 = object$settings$tau2_w,
                          theta = object$theta_w[t, i], 
                          g = eps,
                          v = object$v,
                          nsamples = ifel(settings$mean_map, 0, 1))
            if (settings$mean_map) {
              w_new[, i] <- k$mean
            } else w_new[, i] <- k$samples
          } # end of i for loop
          object$w_approx <- clean_pred_from_approx(object$w_approx)
          object$w_approx$x_ord <- w_t[object$w_approx$ord, , drop = FALSE]
          object$w_approx <- add_pred_to_approx(object$w_approx, w_new, m, 
                                                lite = settings$lite, grad = FALSE,
                                                ord_new = settings$ord_new)
          if (settings$store_latent) out$w_new[j, , ] <- w_new

        } else if (layers == 2) {

          w_t <- as.matrix(object$w[t, , ]) # includes gradients if grad_enhance = TRUE
          if (object$settings$monowarp) { 
            w_new <- monotransform(x_new, object$x_grid, object$w_grid[t, , ])
          } else { # use kriging
            w_new <- matrix(nrow = n_new, ncol = D)
            if (settings$grad) dwdx <- array(dim = c(n_new, D, d))
            for (i in 1:D) { # no grad_enhance option here
              k <- krig_vec(w_t[, i], 
                            approx = object$x_approx,
                            tau2 = object$settings$tau2_w,
                            theta = object$theta_w[t, i],
                            g = eps, 
                            v = object$v, 
                            grad = settings$grad,
                            nsamples = ifel(settings$mean_map, 0, 1),
                            prior_mean = ifel(object$settings$pmx, w_prior_mean[, i], 0),
                            prior_mean_new = ifel(object$settings$pmx, w_prior_mean_new[, i], 0))
              if (settings$mean_map) {
                w_new[, i] <- k$mean
                if (settings$grad) dwdx[, i, ] <- k$grad_mean
              } else {
                w_new[, i] <- k$samples
                if (settings$grad) dwdx[, i, ] <- k$grad_samples[1, , ] # only one sample returned
              }
            } # end of i for loop
          } # end of else statement
          object$w_approx <- clean_pred_from_approx(object$w_approx)
          object$w_approx$x_ord <- w_t[object$w_approx$ord, , drop = FALSE]
          object$w_approx <- add_pred_to_approx(object$w_approx, w_new, m, 
                                                lite = settings$lite, grad = settings$grad,
                                                ord_new = settings$ord_new)
          if (settings$store_latent) out$w_new[j, , ] <- w_new
      
      } # end of layers == 2 else statement
        
      # Finally: map w_new to y
      if (grad_enhance) {
        dydw <- get_dydw(w_t, object$dydx)
        y <- c(object$y, as.vector(dydw))
      }
      k <- krig_vec(y,
                    approx = object$w_approx,
                    tau2 = object$tau2_y[t],
                    theta = object$theta_y[t], 
                    g = g, 
                    v = object$v,
                    s2 = ifel(samples_only, FALSE, settings$lite), 
                    sigma = ifel(samples_only, FALSE, !settings$lite),
                    grad = settings$grad,
                    nsamples = ifel(samples_only, settings$nper, 0))
        
        if (samples_only) {
          row_indx <- (j*settings$nper - settings$nper + 1):(j*settings$nper)
          out$samples[row_indx, ] <- k$samples
          if (settings$grad) {
            for (i in 1:d) {
              out$grad_samples[row_indx, , i] <- 0
              for (l in 1:D) { # grad_samples are dy/dw, must apply multivariate chain rule to get dy/dx
                out$grad_samples[row_indx, , i] <- out$grad_samples[row_indx, , i] + k$grad_samples[, , i] * dwdx[, l, i]
              }
            }
          }
        } else {
          out$mu_t[j, ] <- k$mean
          if (settings$lite) {
            out$s2_sum <- out$s2_sum + k$s2
            if (settings$return_all) out$s2_t[j, ] <- k$s2
          } else out$sigma_sum <- out$sigma_sum + k$sigma
          if (settings$grad) { # MUST USE MULTIVARIATE CHAIN RULE
            for (i in 1:d) {
              out$grad_mu_t[j, , i] <- rowSums(k$grad_mean * dwdx[, , i])
              out$grad_s2_sum[, i] <- out$grad_s2_sum[, i] + rowSums(dwdx[, , i]^2 * k$grad_s2)
            }
          }
          if (settings$EI) {
            if (settings$lite) {
              sig2 <- k$s2 - (object$tau2[t] * g) 
            } else sig2 <- diag(k$sigma) - (object$tau2[t] * g)
            out$ei_sum <- out$ei_sum + exp_improv(k$mean, sig2, y_min)
          }
          if (!is.null(settings$entropy_limit)) {
            if (settings$lite) {
              sig2 <- k$s2 - (object$tau2[t] * g) 
            } else sig2 <- diag(k$sigma) - (object$tau2[t] * g)
            out$ent_sum <- out$ent_sum + calc_entropy(k$mean, sig2, settings$entropy_limit)
          }
        }
        j <- j + 1
      } # end of t for loop
      return(out)
    } # end of foreach loop
    
    stopCluster(cl)
      
    # Group elements out of the list
    if (samples_only) {
      samples <- do.call(rbind, lapply(result, with, eval(parse(text = "samples"))))
      if (settings$grad) grad_samples <- abind::abind(lapply(result, with, 
                                                             eval(parse(text = "grad_samples"))), along = 1)
    } else {
      mu_t <- do.call(rbind, lapply(result, with, eval(parse(text = "mu_t"))))
      if (settings$lite) {
        s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum"))))
        if (settings$return_all) s2_t <- do.call(rbind, lapply(result, with, eval(parse(text = "s2_t"))))
      } else {
        sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum"))))
      }
      if (settings$grad) {
        grad_mu_t <- abind::abind(lapply(result, with, eval(parse(text = "grad_mu_t"))), along = 1)
        grad_s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "grad_s2_sum"))))
      }
      if (settings$store_latent) {
        w_new_store <- abind::abind(lapply(result, with, eval(parse(text = "w_new"))), along = 1)
        if (layers == 3) z_new_store <- abind::abind(lapply(result, with, eval(parse(text = "z_new"))), along = 1)
      }
      if (settings$EI) ei_sum <- Reduce("+", lapply(result, with, eval(parse(text = "ei_sum"))))
      if (!is.null(settings$entropy_limit)) ent_sum <- Reduce("+", lapply(result, with, 
                                                                 eval(parse(text = "ent_sum"))))
    }
    
  } # end of else statement

  # Add variables to the output list
  if (samples_only) {
    if (settings$grad) {
      return(list(y = samples, dy = ifel(d == 1, grad_samples[, , 1], grad_samples)))
    } else return(samples)
  } else {
    object$mean <- colMeans(mu_t)
    if (settings$store_latent) {
      object$w_new <- w_new_store
      if (layers == 3) object$z_new <- z_new_store
    }
    if (settings$lite) { 
      object$s2 <- s2_sum / object$nmcmc + ifel(object$nmcmc > 1, apply(mu_t, 2, var), 0)
      if (settings$return_all) {
        object$mean_all <- mu_t
        object$s2_all <- s2_t
      }
    } else object$Sigma <- sigma_sum / object$nmcmc + ifel(object$nmcmc > 1, cov(mu_t), 0)
    if (settings$grad) {
      object$grad_mean <- matrix(nrow = n_new, ncol = d)
      object$grad_s2 <- grad_s2_sum / object$nmcmc
      for (i in 1:d) {
        object$grad_mean[, i] <- colMeans(grad_mu_t[, , i])
        object$grad_s2[, i] <- object$grad_s2[, i] + ifel(object$nmcmc > 1, apply(grad_mu_t[, , i], 2, var), 0)
      }
    }
    if (settings$EI) object$EI <- ei_sum / object$nmcmc
    if (!is.null(settings$entropy_limit)) object$entropy <- drop(ent_sum / object$nmcmc)
    toc <- proc.time()[[3]]
    object$time <- object$time + unname(toc - tic)
    return(object)
  }
}