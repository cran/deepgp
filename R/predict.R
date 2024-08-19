
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   predict.gp
#   predict.dgp2
#   predict.dgp3
# Internal:
#   predict_nonvec (used in all three of the above)

# Define Predict for S3 Objects -----------------------------------------------
#' @name predict
#' @title Predict posterior mean and variance/covariance
#' @description Acts on a \code{gp}, \code{dgp2}, or \code{dgp3} object.
#'     Calculates posterior mean and variance/covariance over specified input 
#'     locations.  Optionally calculates expected improvement (EI) or entropy 
#'     over candidate inputs.  Optionally utilizes SNOW parallelization.
#' 
#' @details All iterations in the object are used for prediction, so samples 
#'     should be burned-in.  Thinning the samples using \code{trim} will speed 
#'     up computation.  Posterior moments are calculated using conditional 
#'     expectation and variance.  As a default, only point-wise variance is 
#'     calculated.  Full covariance may be calculated using \code{lite = FALSE}. 
#'     
#'     Expected improvement is calculated with the goal of minimizing the 
#'     response.  See Chapter 7 of Gramacy (2020) for details.  Entropy is 
#'     calculated based on two classes separated by the specified limit.  
#'     See Sauer (2023, Chapter 3) for details.
#'     
#'     SNOW parallelization reduces computation time but requires 
#'     more memory storage.
#' 
#' @param object object from \code{fit_one_layer}, \code{fit_two_layer}, or 
#'        \code{fit_three_layer} with burn-in already removed
#' @param x_new matrix of predictive input locations
#' @param lite logical indicating whether to calculate only point-wise 
#'        variances (\code{lite = TRUE}) or full covariance 
#'        (\code{lite = FALSE})
#' @param store_latent logical indicating whether to store and return mapped 
#'        values of latent layers (two or three layer models only)
#' @param mean_map logical indicating whether to map hidden layers using 
#'        conditional mean (\code{mean_map = TRUE}) or using a random sample
#'        from the full MVN distribution (two or three layer models only),
#'        \code{mean_map = FALSE} is not yet implemented for fits with 
#'        \code{vecchia = TRUE} 
#' @param return_all logical indicating whether to return mean and point-wise
#'        variance prediction for ALL samples (only available for \code{lite = TRUE})
#' @param EI logical indicating whether to calculate expected improvement 
#'        (for minimizing the response)
#' @param entropy_limit optional limit state for entropy calculations (separating
#'        passes and failures), default value of \code{NULL} bypasses entropy
#'        calculations
#' @param cores number of cores to utilize in parallel
#' @param m size of Vecchia conditioning sets (only for fits with 
#'        \code{vecchia = TRUE}), defaults to the \code{m} used for MCMC
#' @param ordering_new optional ordering for Vecchia approximation, must correspond
#'        to rows of \code{x_new}, defaults to random, is applied to all layers
#'        in deeper models
#' @param ... N/A
#' @return object of the same class with the following additional elements:
#' \itemize{
#'   \item \code{x_new}: copy of predictive input locations
#'   \item \code{mean}: predicted posterior mean, indices correspond to 
#'         \code{x_new} locations
#'   \item \code{s2}: predicted point-wise variances, indices correspond to 
#'         \code{x_new} locations (only returned when \code{lite = TRUE})
#'   \item \code{mean_all}: predicted posterior mean for each sample (column
#'         indices), only returned when \code{return_all = TRUE}
#'   \item \code{s2_all} predicted point-wise variances for each sample (column
#'         indices), only returned when \code{return-all = TRUE}
#'   \item \code{Sigma}: predicted posterior covariance, indices correspond to 
#'         \code{x_new} locations (only returned when \code{lite = FALSE})
#'   \item \code{EI}: vector of expected improvement values, indices correspond 
#'         to \code{x_new} locations (only returned when \code{EI = TRUE})
#'   \item \code{entropy}: vector of entropy values, indices correspond to 
#'         \code{x_new} locations (only returned when \code{entropy_limit} is
#'         numeric)
#'   \item \code{w_new}: list of hidden layer mappings (only returned when 
#'         \code{store_latent = TRUE}), list index corresponds to iteration and 
#'         row index corresponds to \code{x_new} location (two or three layer 
#'         models only)
#'   \item \code{z_new}: list of hidden layer mappings (only returned when 
#'         \code{store_latent = TRUE}), list index corresponds to iteration and 
#'         row index corresponds to \code{x_new} location (three layer models only) 
#' }
#' Computation time is added to the computation time of the existing object.
#' 
#' @references 
#' Sauer, A. (2023). Deep Gaussian process surrogates for computer experiments. 
#'      *Ph.D. Dissertation, Department of Statistics, Virginia Polytechnic Institute and State University.*
#'      \cr\cr
#' Sauer, A., Gramacy, R.B., & Higdon, D. (2023). Active learning for deep 
#'      Gaussian process surrogates. *Technometrics, 65,* 4-18.  arXiv:2012.08015
#'      \cr\cr
#' Sauer, A., Cooper, A., & Gramacy, R. B. (2023). Vecchia-approximated deep Gaussian 
#'      processes for computer experiments. 
#'      *Journal of Computational and Graphical Statistics, 32*(3), 824-837.  arXiv:2204.02904
#'      \cr\cr
#' Barnett, S., Beesley, L. J., Booth, A. S., Gramacy, R. B., & Osthus D. (2024). Monotonic 
#'      warpings for additive and deep Gaussian processes. *In Review.* arXiv:2408.01540
#'
#' @examples 
#' # See ?fit_one_layer, ?fit_two_layer, or ?fit_three_layer
#' # for examples
#'
#' @rdname predict
NULL

# Predict One Layer -----------------------------------------------------------
#' @rdname predict
#' @export

predict.gp <- function(object, x_new, lite = TRUE, return_all = FALSE, 
                       EI = FALSE, entropy_limit = NULL, cores = 1, ...) {
  
  settings <- list(lite = lite, return_all = return_all, EI = EI, 
                   entropy_limit = entropy_limit, cores = cores)
  object <- predict_nonvec(object, x_new, settings, layers = 1)
  return(object)
}

# Predict Two Layer -----------------------------------------------------------
#' @rdname predict
#' @export

predict.dgp2 <- function(object, x_new, lite = TRUE, store_latent = FALSE, 
                         mean_map = TRUE, return_all = FALSE, EI = FALSE, 
                         entropy_limit = NULL, cores = 1, ...) {
  
  settings <- list(lite = lite, store_latent = store_latent, mean_map = mean_map,
                   return_all = return_all, EI = EI, 
                   entropy_limit = entropy_limit, cores = cores)
  object <- predict_nonvec(object, x_new, settings, layers = 2)
  return(object)
}
  
# Predict Three Layer ---------------------------------------------------------
#' @rdname predict
#' @export

predict.dgp3 <- function(object, x_new, lite = TRUE, store_latent = FALSE, 
                         mean_map = TRUE, return_all = FALSE, EI = FALSE, 
                         entropy_limit = NULL, cores = 1, ...) {
  
  settings <- list(lite = lite, store_latent = store_latent, mean_map = mean_map,
                   return_all = return_all, EI = EI, 
                   entropy_limit = entropy_limit, cores = cores)
  object <- predict_nonvec(object, x_new, settings, layers = 3)
  return(object)
}

# Predict Non-Vecchia ---------------------------------------------------------
  
predict_nonvec <- function(object, x_new, settings, layers) {
  
  tic <- proc.time()[[3]]
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  monowarp <- (!is.null(object$x_grid))
  object$x_new <- x_new
  n_new <- nrow(object$x_new)

  if (monowarp) {
     grid_index <- fo_approx_init(object$x_grid, object$x)
     grid_index_new <- fo_approx_init(object$x_grid, x_new)
  }

  if (layers >= 2) {
    if (monowarp) D <- ncol(object$w_grid[[1]]) else D <- ncol(object$w[[1]]) 
    if (settings$lite & !settings$mean_map)
      stop("mean_map = FALSE requires lite = FALSE")
    if (!settings$mean_map) 
      message("mean_map = FALSE may cause numerical instability in latent layer mapping")
  }

  if (layers == 1) {
    sep <- is.matrix(object$theta)
  } else sep <- monowarp # separable lengthscales only used in monotonic warpings

  if (!sep) {
    dx <- sq_dist(object$x)
    dx_cross <- sq_dist(object$x_new, object$x)
    if (!settings$lite) dx_new <- sq_dist(object$x_new) else dx_new <- NULL
  }

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
        for (i in 1:D) {
          k <- krig(z_t[, i], dx, dx_new, dx_cross, object$theta_z[t, i], 
                    g = eps, sigma = !settings$mean_map, v = object$v)
          if (settings$mean_map) {
            z_new[, i] <- k$mean
          } else z_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        } # end of i for loop
        if (settings$store_latent) z_new_list[[t]] <- z_new
        dz <- sq_dist(z_t)
        dz_cross <- sq_dist(z_new, z_t)
        if (!settings$mean_map) dz_new <- sq_dist(z_new) else dz_new <- NULL
      }
      
      if (layers >= 2) {
        # 2 layers: map x_new to w_new (separately for each dimension)
        # 3 layers: map z_new to w_new (separately for each dimension)
        if (monowarp) { # simply get the monowarped values at x_new locations
          w_grid_t <- object$w_grid[[t]]
          w_t <- monowarp_ref(object$x, object$x_grid, w_grid_t, grid_index)
          w_new <- monowarp_ref(x_new, object$x_grid, w_grid_t, grid_index_new)
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
            k <- krig(w_t[, i], 
                    ifel(layers == 2, dx, dz), 
                    ifel(layers == 2, dx_new, dz_new),
                    ifel(layers == 2, dx_cross, dz_cross),
                    object$theta_w[t, i], g = eps, tau2 = prior_tau2,
                    sigma = !settings$mean_map,
                    v = object$v, prior_mean = prior_mean,
                    prior_mean_new = prior_mean_new)
            if (settings$mean_map) {
              w_new[, i] <- k$mean
            } else w_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
          } # end of i for loop
          dw <- sq_dist(w_t)
          dw_cross <- sq_dist(w_new, w_t)
          if (!settings$lite) dw_new <- sq_dist(w_new) else dw_new <- NULL
        } # end of monowarp else statement
        if (settings$store_latent) w_new_list[[t]] <- w_new
      }
      
      # 1 layer: map x_new to y
      # 2 and 3 layers: map w_new to y
      if (sep) { # only occurs in one layer or monowarped two-layer
        k <- krig_sep(object$y, 
                      ifel(monowarp, w_t, object$x), 
                      ifel(monowarp, w_new, x_new), 
                      ifel(monowarp, object$theta_y[t, ], object$theta[t, ]), 
                      object$g[t], 
                      object$tau2[t], s2 = settings$lite, 
                      sigma = !settings$lite, f_min = f_min, 
                      v = object$v)
      } else { 
        k <- krig(object$y, ifel(layers == 1, dx, dw), 
                  ifel(layers == 1, dx_new, dw_new), 
                  ifel(layers == 1, dx_cross, dw_cross),
                  ifel(layers == 1, object$theta[t], object$theta_y[t]), 
                  object$g[t], object$tau2[t], s2 = settings$lite, 
                  sigma = !settings$lite, f_min = f_min, v = object$v)
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
        if (settings$return_all) out$s2_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
      } else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
      if (settings$EI) out$ei_sum <- rep(0, times = n_new)
      if (!is.null(settings$entropy_limit)) out$ent_sum <- rep(0, times = n_new)
      if (layers >= 2) {
        if (settings$store_latent) {
          out$w_new <- list()
          if (layers == 3) out$z_new <- list()
        }
      }
      
      j <- 1
      for (t in chunks[[thread]]) {

        if (layers == 3) { 
          # 3 layers: map x_new to z_new (separately for each dimension)
          z_t <- object$z[[t]]
          z_new <- matrix(nrow = n_new, ncol = D)
          for (i in 1:D) {
            k <- krig(z_t[, i], dx, dx_new, dx_cross, object$theta_z[t, i], 
                      g = eps, sigma = !settings$mean_map, v = object$v)
            if (settings$mean_map) {
              z_new[, i] <- k$mean
            } else z_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
          } # end of i for loop
          if (settings$store_latent) out$z_new[[j]] <- z_new
          dz <- sq_dist(z_t)
          dz_cross <- sq_dist(z_new, z_t)
          if (!settings$mean_map) dz_new <- sq_dist(z_new) else dz_new <- NULL
        }
        
        if (layers >= 2) {
          # 2 layers: map x_new to w_new (separately for each dimension)
          # 3 layers: map z_new to w_new (separately for each dimension)
          if (monowarp) { # simply get the monowarped values at x_new locations
            w_grid_t <- object$w_grid[[t]]
            w_t <- monowarp_ref(object$x, object$x_grid, w_grid_t, grid_index)
            w_new <- monowarp_ref(x_new, object$x_grid, w_grid_t, grid_index_new)
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
              k <- krig(w_t[, i], ifel(layers == 2, dx, dz), 
                        ifel(layers == 2, dx_new, dz_new),
                        ifel(layers == 2, dx_cross, dz_cross),
                        object$theta_w[t, i], g = eps, tau2 = prior_tau2,
                        sigma = !settings$mean_map,
                        v = object$v, prior_mean = prior_mean,
                        prior_mean_new = prior_mean_new)
              if (settings$mean_map) {
                w_new[, i] <- k$mean
              } else w_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
            } # end of i for loop
            dw <- sq_dist(w_t)
            dw_cross <- sq_dist(w_new, w_t)
            if (!settings$lite) dw_new <- sq_dist(w_new) else dw_new <- NULL
          } # end of monowarp if statement
          if (settings$store_latent) out$w_new[[t]] <- w_new
        }
        
        # 1 layer: map x_new to y
        # 2 and 3 layers: map w_new to y
        if (sep) { # only occurs in one layer or monowarped two-layer
          k <- krig_sep(object$y, 
                        ifel(monowarp, w_t, object$x), 
                        ifel(monowarp, w_new, x_new), 
                        ifel(monowarp, object$theta_y[t, ], object$theta[t, ]), 
                        object$g[t], 
                        object$tau2[t], s2 = settings$lite, 
                        sigma = !settings$lite, f_min = f_min, 
                        v = object$v)
        } else { 
          k <- krig(object$y, ifel(layers == 1, dx, dw), 
                    ifel(layers == 1, dx_new, dw_new), 
                    ifel(layers == 1, dx_cross, dw_cross),
                    ifel(layers == 1, object$theta[t], object$theta_y[t]), 
                    object$g[t], object$tau2[t], s2 = settings$lite, 
                    sigma = !settings$lite, f_min = f_min, v = object$v)
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
    } # end of foreach statement
    
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
        if (layers == 3)  z_new_list <- unlist(lapply(result, with, eval(parse(text = "z_new"))), recursive = FALSE)
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
