
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   predict.gp
#   predict.dgp2
#   predict.dgp3
# Internal:
#   predict_shallow: predictions for one layer models
#   predict_deep: predictions for two and three layer models

# predict S3 class -----------------------------------------------
#' @name predict
#' @title Predict posterior mean and variance/covariance
#' @description Acts on a \code{gp}, \code{gpvec}, \code{dgp2}, \code{dgp2vec},
#'     \code{dgp3}, or \code{dgp3vec} object.
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
#' @param x_new vector or matrix of predictive input locations
#' @param lite logical indicating whether to calculate only point-wise 
#'        variances (\code{lite = TRUE}) or full covariance 
#'        (\code{lite = FALSE})
#' @param grad logical indicating whether to additionally calculate/return 
#'        predictions of the gradient (one and two layer models only)
#' @param store_latent logical indicating whether to store and return mapped 
#'        values of latent layers (two or three layer models only)
#' @param mean_map logical indicating whether to map hidden layers using 
#'        conditional mean (\code{mean_map = TRUE}) or using a random sample
#'        from the full MVN distribution (two or three layer models only)
#' @param return_all logical indicating whether to return mean and point-wise
#'        variance prediction for ALL samples (only available for \code{lite = TRUE})
#' @param EI logical indicating whether to calculate expected improvement 
#'        (for minimizing the response)
#' @param entropy_limit optional limit state for entropy calculations (separating
#'        passes and failures), default value of \code{NULL} bypasses entropy
#'        calculations
#' @param cores number of cores to utilize for SNOW parallelization
#' @param m size of Vecchia conditioning sets, defaults to the lower of twice the 
#'        \code{m} used for MCMC or the maximum available (only for fits with 
#'        \code{vecchia = TRUE}), 
#' @param ord_new optional ordering for Vecchia approximation with \code{lite = FALSE}, 
#'        must correspond to rows of \code{x_new}, defaults to random, is 
#'        applied to all layers in deeper models
#' @param ... N/A
#' @return object of the same class with the following additional elements:
#' \itemize{
#'   \item \code{x_new}: copy of predictive input locations
#'   \item \code{mean}: predicted posterior mean, indices correspond to 
#'         \code{x_new} locations
#'   \item \code{s2}: predicted point-wise variances, indices correspond to 
#'         \code{x_new} locations (only returned when \code{lite = TRUE})
#'   \item \code{mean_all}: predicted posterior mean for each sample (rows correspond
#'         to iterations), only returned when \code{return_all = TRUE}
#'   \item \code{s2_all}: predicted point-wise variances for each sample (rows correspond
#'         to iterations), only returned when \code{return_all = TRUE}
#'   \item \code{Sigma}: predicted posterior covariance, indices correspond to 
#'         \code{x_new} locations (only returned when \code{lite = FALSE})
#'   \item \code{grad_mean}: predicted posterior mean of the gradient (rows correspond
#'         to \code{x_new}, columns correspond to dimension, only returned when
#'         \code{grad = TRUE})
#'   \item \code{grad_s2}: predicted point-wise variances of the gradient (rows correspond
#'         to \code{x_new}, columns correspond to dimension, only returned when
#'         \code{grad = TRUE})
#'   \item \code{EI}: vector of expected improvement values, indices correspond 
#'         to \code{x_new} locations (only returned when \code{EI = TRUE})
#'   \item \code{entropy}: vector of entropy values, indices correspond to 
#'         \code{x_new} locations (only returned when \code{entropy_limit} is
#'         numeric)
#'   \item \code{w_new}: array of hidden layer mappings, with dimensions corresponding
#'         to iteration, then \code{x_new} location, then dimension (only returned when 
#'         \code{store_latent = TRUE})
#'   \item \code{z_new}: array of hidden layer mappings, with dimensions corresponding
#'         to iteration, then \code{x_new} location, then dimension (only returned when 
#'         \code{store_latent = TRUE})
#' }
#' Computation time is added to the computation time of the existing object.
#' 
#' @references 
#' Sauer, A. (2023). Deep Gaussian process surrogates for computer experiments. 
#'      *Ph.D. Dissertation, Department of Statistics, Virginia Polytechnic Institute and State University.*
#'      \url{http://hdl.handle.net/10919/114845}
#'      \cr\cr
#' Booth, A. S. (2025). Deep Gaussian processes with gradients. arXiv:2512.18066
#'      \cr\cr
#' Sauer, A., Gramacy, R.B., & Higdon, D. (2023). Active learning for deep 
#'      Gaussian process surrogates. *Technometrics, 65,* 4-18.  arXiv:2012.08015
#'      \cr\cr
#' Sauer, A., Cooper, A., & Gramacy, R. B. (2023). Vecchia-approximated deep Gaussian 
#'      processes for computer experiments. 
#'      *Journal of Computational and Graphical Statistics, 32*(3), 824-837.  arXiv:2204.02904
#'      \cr\cr
#' Gramacy, R. B., Sauer, A. & Wycoff, N. (2022). Triangulation candidates for Bayesian 
#'     optimization.  *Advances in Neural Information Processing Systems (NeurIPS), 35,* 
#'     35933-35945.  arXiv:2112.07457
#'     \cr\cr
#' Booth, A., Renganathan, S. A. & Gramacy, R. B. (2025). Contour location for 
#'     reliability in airfoil simulation experiments using deep Gaussian 
#'     processes. *Annals of Applied Statistics, 19*(1), 191-211. arXiv:2308.04420
#'     \cr\cr
#' Barnett, S., Beesley, L. J., Booth, A. S., Gramacy, R. B., & Osthus D. (2025). 
#'     Monotonic warpings for additive and deep Gaussian processes. 
#'     *Statistics and Computing, 35*(3), 65. arXiv:2408.01540
#'
#' @examples 
#' # See ?fit_one_layer, ?fit_two_layer, or ?fit_three_layer
#' # for examples
#'
#' @rdname predict
NULL

# predict.gp ------------------------------------------------------------------
#' @rdname predict
#' @export

predict.gp <- function(object, x_new, lite = TRUE, grad = FALSE, return_all = FALSE, 
                       EI = FALSE, entropy_limit = NULL, cores = 1, ...) {
  
  if (return_all & !lite) 
    stop("return_all only offered when lite = TRUE")
  if (!is.null(entropy_limit) & !is.numeric(entropy_limit))
    stop("entropy_limit must be numeric")
  cores <- check_cores(cores, object$nmcmc)
  if (grad) { 
    if (object$v != 999) stop("grad only offered with cov = 'exp2'")
    if (!lite) stop("grad = TRUE requires lite = TRUE")
  }
  
  settings <- list(lite = lite, grad = grad, return_all = return_all, EI = EI, 
                   entropy_limit = entropy_limit, cores = cores)
  object <- predict_shallow(object, x_new, settings)
  return(object)
}

# predict.dgp2 ----------------------------------------------------------------
#' @rdname predict
#' @export

predict.dgp2 <- function(object, x_new, lite = TRUE, grad = FALSE, 
                         store_latent = FALSE, 
                         mean_map = TRUE, return_all = FALSE, EI = FALSE, 
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
  
  settings <- list(lite = lite, grad = grad, store_latent = store_latent, 
                   mean_map = mean_map, return_all = return_all, EI = EI, 
                   entropy_limit = entropy_limit, cores = cores)
  object <- predict_deep(object, x_new, settings, layers = 2)
  return(object)
}
  
# predict.dgp3 ----------------------------------------------------------------
#' @rdname predict
#' @export

predict.dgp3 <- function(object, x_new, lite = TRUE,
                         store_latent = FALSE, 
                         mean_map = TRUE, return_all = FALSE, EI = FALSE, 
                         entropy_limit = NULL, cores = 1, ...) {
  
  if (return_all & !lite) 
    stop("return_all only offered when lite = TRUE")
  if (!is.null(entropy_limit) & !is.numeric(entropy_limit))
    stop("entropy_limit must be numeric")
  cores <- check_cores(cores, object$nmcmc)
  
  settings <- list(lite = lite, grad = FALSE, store_latent = store_latent, 
                   mean_map = mean_map, return_all = return_all, EI = EI, 
                   entropy_limit = entropy_limit, cores = cores)
  object <- predict_deep(object, x_new, settings, layers = 3)
  return(object)
}

# predict_shallow -------------------------------------------------------------

predict_shallow <- function(object, x_new, settings, samples_only = FALSE) { 
  # One-layer only
  # If samples_only = TRUE, posterior sample draws are returned instead of
  # summarized posterior moments

  tic <- proc.time()[[3]]
  if (is.vector(x_new)) x_new <- as.matrix(x_new)
  if (ncol(x_new) != ncol(object$x)) stop("dimension of x_new does not match dimension of x")
  object$x_new <- x_new
  n_new <- nrow(x_new)
  x <- object$x
  d <- ncol(x)
  sep <- is.matrix(object$theta)

  grad_enhance <- !is.null(object$dydx)
  if (grad_enhance) { 
    if (settings$EI) stop("EI not implemented for grad_enhance")
    if (!is.null(settings$entropy_limit)) stop("entropy_limit not implemented for grad_enhance")
    y <- c(object$y, as.vector(object$dydx))
  } else y <- object$y
  
  # Pre-calculate distance matrices?
  if (sep | grad_enhance | settings$grad) {
    xdmat <- NULL
    xdmat_cross <- NULL
    xdmat_new <- NULL
  } else {
    xdmat <- sq_dist(x)
    xdmat_cross <- sq_dist(x_new, x)
    if (settings$lite) xdmat_new <- NULL else xdmat_new <- sq_dist(x_new)
  }

  if (settings$EI) y_min <- min(object$y) # use smallest observed value
  
  if (settings$cores == 1) { # run serial for loop
    
    if (samples_only) {
      samples <- matrix(nrow = settings$nper*object$nmcmc, ncol = n_new)
      if (settings$grad) grad_samples <- array(dim = c(settings$nper*object$nmcmc, n_new, d))
    } else {
      mu_t <- matrix(nrow = object$nmcmc, ncol = n_new)
      if (settings$lite) {
        s2_sum <- rep(0, times = n_new)
        if (settings$return_all) s2_t <- matrix(nrow = n_new, ncol = object$nmcmc)
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

      k <- krig(y, 
                xdmat = xdmat, 
                xdmat_new = xdmat_new, 
                xdmat_cross = xdmat_cross,
                x = x,
                x_new = x_new,
                tau2 = object$tau2[t],
                theta = ifel(sep, object$theta[t, ], object$theta[t]),
                g = g, 
                v = object$v,
                sep = sep,
                s2 = settings$lite, 
                sigma = !settings$lite, 
                grad = settings$grad,
                grad_enhance = grad_enhance,
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
          if (settings$return_all) out$s2_t <- matrix(nrow = length(chunks[[thread]]), ncol = n_new)
        } else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
        if (settings$grad) {
          out$grad_mu_t <- array(dim = c(length(chunks[[thread]]), n_new, d))
          out$grad_s2_sum <- matrix(0, nrow = n_new, ncol = d)
        }
        if (settings$EI) out$ei_sum <- rep(0, times = n_new)
        if (!is.null(settings$entropy_limit)) out$ent_sum <- rep(0, times = n_new)
      }
      
      j <- 1
      for (t in chunks[[thread]]) {
        g <- ifel(length(object$g) == 1, object$g, object$g[t])

        k <- krig(y, 
                  xdmat = xdmat, 
                  xdmat_new = xdmat_new, 
                  xdmat_cross = xdmat_cross,
                  x = x,
                  x_new = x_new,
                  tau2 = object$tau2[t],
                  theta = ifel(sep, object$theta[t, ], object$theta[t]),
                  g = g, 
                  v = object$v,
                  sep = sep,
                  s2 = settings$lite, 
                  sigma = !settings$lite, 
                  grad = settings$grad,
                  grad_enhance = grad_enhance,
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
    } # end of foreach statement
    
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
  
  # Add variables to the output list (post process if grad = TRUE)  
  if (samples_only) {
    if (settings$grad) {
      return(list(y = samples, dydx = ifel(d == 1, grad_samples[, , 1], grad_samples)))
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
    if (settings$grad) {
      object$grad_mean <- matrix(nrow = n_new, ncol = d)
      object$grad_s2 <- grad_s2_sum / object$nmcmc
      for (i in 1:d) {
        object$grad_mean[, i] <- colMeans(grad_mu_t[, , i]) # average over MCMC iterations
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


# predict_deep ----------------------------------------------------------------
  
predict_deep <- function(object, x_new, settings, layers, samples_only = FALSE) { 
  # Two-layer or three-layer
  # If samples_only = TRUE, posterior sample draws are returned instead of
  # summarized posterior moments
  
  tic <- proc.time()[[3]]
  if (is.vector(x_new)) x_new <- as.matrix(x_new)
  if (ncol(x_new) != ncol(object$x)) stop("dimension of x_new does not match dimension of x")
  object$x_new <- x_new
  n_new <- nrow(x_new)
  n <- length(object$y)
  x <- object$x
  d <- ncol(x)
  D <- dim(object$w)[3]
  
  grad_enhance <- !is.null(object$dydx)
  if (grad_enhance) {
    if (settings$EI) stop("EI not implemented for grad_enhance")
    if (!is.null(settings$entropy_limit)) stop("entropy_limit not implemented for grad_enhance")
    # y will be recalculated for every w_t
  } else y <- object$y

  # Pre-calculate distance matrices?
  if (grad_enhance | settings$grad | object$settings$monowarp) {
    xdmat <- NULL
    xdmat_cross <- NULL
    xdmat_new <- NULL
    wdmat <- NULL
    wdmat_cross <- NULL
    wdmat_new <- NULL
  } else {
    xdmat <- sq_dist(x)
    xdmat_cross <- sq_dist(x_new, x)
    if (settings$mean_map) xdmat_new <- NULL else xdmat_new <- sq_dist(x_new)
    # wdmat will be recalculated for every w_t
  }

  # Prespecify prior mean (if not zero)
  if (object$settings$pmx) {
    if (grad_enhance) {
      w_prior_mean <- object$settings$w_prior_mean
    } else w_prior_mean <- x
    if (settings$grad) {
      w_prior_mean_new <- get_prior_mean(x_new)
    } else w_prior_mean_new <- x_new
  }

  if (settings$EI) y_min <- min(object$y) # use smallest observed value
  
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
      if (settings$store_latent) {
        w_new_store <- array(dim = c(object$nmcmc, nrow(x_new), D))
        if (layers == 3) z_new_store <- array(dim = c(object$nmcmc, nrow(x_new), D))
      }
    }
    
    for (t in 1:object$nmcmc) {
      g <- ifel(length(object$g) == 1, object$g, object$g[t])
      
      if (layers == 3) { # no pmx, monowarp, grad, or grad_enhance option

        # First map x_new to z_new
        z_t <- as.matrix(object$z[t, , ])
        z_new <- matrix(nrow = n_new, ncol = D)
        for (i in 1:D) {
          k <- krig(z_t[, i], 
                    xdmat = xdmat, 
                    xdmat_new = xdmat_new, 
                    xdmat_cross = xdmat_cross, 
                    tau2 = object$settings$tau2_z,
                    theta = object$theta_z[t, i], 
                    g = eps, 
                    v = object$v, 
                    nsamples = ifel(settings$mean_map, 0, 1))
          if (settings$mean_map) {
            z_new[, i] <- k$mean
          } else z_new[, i] <- k$samples
        } # end of i for loop
        zdmat <- sq_dist(z_t)
        zdmat_cross <- sq_dist(z_new, z_t)
        if (!settings$mean_map) zdmat_new <- sq_dist(z_new) else zdmat_new <- NULL
        if (settings$store_latent) z_new_store[t, , ] <- z_new

        # Then map z_new to w_new
        w_t <- as.matrix(object$w[t, , ])
        w_new <- matrix(nrow = n_new, ncol = D)
        for (i in 1:D) {
          k <- krig(w_t[, i], 
                    xdmat = zdmat, 
                    xdmat_new = zdmat_new,
                    xdmat_cross = zdmat_cross,
                    tau2 = object$settings$tau2_w,
                    theta = object$theta_w[t, i], 
                    g = eps,
                    v = object$v,
                    nsamples = ifel(settings$mean_map, 0, 1))
          if (settings$mean_map) {
            w_new[, i] <- k$mean
          } else w_new[, i] <- k$samples
        }
        wdmat <- sq_dist(w_t)
        wdmat_cross <- sq_dist(w_new, w_t)
        if (!settings$lite) wdmat_new <- sq_dist(w_new) else wdmat_new <- NULL
        if (settings$store_latent) w_new_store[t, , ] <- w_new
      
      } else if (layers == 2) {

        w_t <- as.matrix(object$w[t, , ]) # includes gradients if grad_enhance = TRUE
        if (object$settings$monowarp) {
          w_new <- monotransform(x_new, object$x_grid, object$w_grid[t, , ])
        } else { # use kriging
          w_new <- matrix(nrow = n_new, ncol = D)
          if (settings$grad) dwdx <- array(dim = c(n_new, D, d))
          for (i in 1:D) {
            k <- krig(w_t[, i], 
                      xdmat = xdmat, 
                      xdmat_new = xdmat_new,
                      xdmat_cross = xdmat_cross,
                      x = x,
                      x_new = x_new,
                      tau2 = object$settings$tau2_w,
                      theta = object$theta_w[t, i], 
                      g = eps,
                      v = object$v,
                      grad = settings$grad,
                      grad_enhance = grad_enhance,
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
        } # end of monowarp else statement
        if (!grad_enhance & !settings$grad) {
          wdmat <- sq_dist(w_t) 
          wdmat_cross <- sq_dist(w_new, w_t)
          if (settings$lite) wdmat_new <- NULL else wdmat_new <- sq_dist(w_new)
        }
        if (settings$store_latent) w_new_store[t, , ] <- w_new

      } # end of layers == 2 else statement
      
      # Finally: map w_new to y
      if (grad_enhance) {
        dydw <- get_dydw(w_t, object$dydx)
        y <- c(object$y, as.vector(dydw))
        w_t <- w_t[1:n, , drop = FALSE] # remove gradients for prediction
      } 
      k <- krig(y, 
                xdmat = wdmat, 
                xdmat_new = wdmat_new, 
                xdmat_cross = wdmat_cross,
                x = w_t,
                x_new = w_new,
                tau2 = object$tau2_y[t],
                theta = object$theta_y[t], 
                g = g, 
                v = object$v,
                s2 = settings$lite, 
                sigma = !settings$lite, 
                grad = settings$grad,
                grad_enhance = grad_enhance,
                nsamples = ifel(samples_only, settings$nper, 0))
     
      if (samples_only) {
        row_indx <- (t*settings$nper - settings$nper + 1):(t*settings$nper)
        samples[row_indx, ] <- k$samples
        if (settings$grad) {
          grad_samples[row_indx, , i] <- 0
          for (j in 1:D) {
            grad_samples[row_indx, , i] <- grad_samples[row_indx, , i] + k$grad_samples[, , i] * dwdx[, j, i]
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
        if (settings$grad) out$grad_samples <- array(dim = c(settings$nper*object$nmcmc, n_new, d))
      } else {
        out$mu_t <- matrix(nrow = length(chunks[[thread]]), ncol = n_new)
        if (settings$lite) {
          out$s2_sum <- rep(0, times = n_new)
          if (settings$return_all) out$s2_t <- matrix(nrow = length(chunks[[thread]]), ncol = n_new)
        } else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
        if (settings$grad) { # no return_all or lite = FALSE options
          out$grad_mu_t <- array(dim = c(length(chunks[[thread]]), n_new, d))
          out$grad_s2_sum <- matrix(0, nrow = n_new, ncol = d)
        }
        if (settings$EI) out$ei_sum <- rep(0, times = n_new)
        if (!is.null(settings$entropy_limit)) out$ent_sum <- rep(0, times = n_new)
        if (settings$store_latent) {
          out$w_new <- array(dim = c(length(chunks[[thread]]), nrow(x_new), D))
          if (layers == 3) out$z_new <- array(dim = c(length(chunks[[thread]]), nrow(x_new), D))
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
            k <- krig(z_t[, i], 
                      xdmat = xdmat, 
                      xdmat_new = xdmat_new, 
                      xdmat_cross = xdmat_cross, 
                      tau2 = object$settings$tau2_z,
                      theta = object$theta_z[t, i], 
                      g = eps, 
                      v = object$v, 
                      nsamples = ifel(settings$mean_map, 0, 1))
            if (settings$mean_map) {
              z_new[, i] <- k$mean
            } else z_new[, i] <- k$samples
          } # end of i for loop
          zdmat <- sq_dist(z_t)
          zdmat_cross <- sq_dist(z_new, z_t)
          if (!settings$mean_map) zdmat_new <- sq_dist(z_new) else zdmat_new <- NULL
          if (settings$store_latent) out$z_new[j, , ] <- z_new

          # Then map z_new to w_new
          w_t <- as.matrix(object$w[t, , ])
          w_new <- matrix(nrow = n_new, ncol = D)
          for (i in 1:D) {
            k <- krig(w_t[, i], 
                      xdmat = zdmat, 
                      xdmat_new = zdmat_new,
                      xdmat_cross = zdmat_cross,
                      tau2 = object$settings$tau2_w,
                      theta = object$theta_w[t, i], 
                      g = eps,
                      v = object$v,
                      nsamples = ifel(settings$mean_map, 0, 1))
            if (settings$mean_map) {
              w_new[, i] <- k$mean
            } else w_new[, i] <- k$samples
          }
          wdmat <- sq_dist(w_t)
          wdmat_cross <- sq_dist(w_new, w_t)
          if (!settings$lite) wdmat_new <- sq_dist(w_new) else wdmat_new <- NULL
          if (settings$store_latent) out$w_new[j, , ] <- w_new
        
        } else if (layers == 2) {
        
          w_t <- as.matrix(object$w[t, , ]) # includes gradients if grad_enhance = TRUE
          if (object$settings$monowarp) {
            w_new <- monotransform(x_new, object$x_grid, object$w_grid[t, , ])
          } else { # use kriging
            w_new <- matrix(nrow = n_new, ncol = D)
            if (settings$grad) dwdx <- array(dim = c(n_new, D, d))
            for (i in 1:D) {
              k <- krig(w_t[, i], 
                        xdmat = xdmat, 
                        xdmat_new = xdmat_new,
                        xdmat_cross = xdmat_cross,
                        x = x,
                        x_new = x_new,
                        tau2 = object$settings$tau2_w,
                        theta = object$theta_w[t, i], 
                        g = eps, 
                        v = object$v, 
                        grad = settings$grad,
                        grad_enhance = grad_enhance,
                        nsamples = ifel(settings$mean_map, 0, 1),
                        prior_mean = ifel(object$settings$pmx, w_prior_mean[, i], 0),
                        prior_mean_new = ifel(object$settings$pmx, w_prior_mean_new[, i], 0))
              if (settings$mean_map) {
                w_new[, i] <- k$mean
                if (settings$grad) dwdx[, i, ] <- k$grad_mean
              } else {
                w_new[, i] <- k$samples
                if (settings$grad) dwdx[, i, ] <- k$grad_samples[1, , ]
              }
            } # end of i for loop
          } # end of monowarp else statement 
          if (!grad_enhance & !settings$grad) {
            wdmat <- sq_dist(w_t)
            wdmat_cross <- sq_dist(w_new, w_t)
            if (!settings$lite) wdmat_new <- sq_dist(w_new) else wdmat_new <- NULL
          }
          if (settings$store_latent) out$w_new[j, , ] <- w_new
        } # end of layers == 2 else statement
        
        # Finally: map w_new to y
        if (grad_enhance) {
          dydw <- get_dydw(w_t, object$dydx)
          y <- c(object$y, as.vector(dydw))
          w_t <- w_t[1:n, , drop = FALSE] # remove gradients for prediction
        }
        k <- krig(y, 
                  xdmat = wdmat, 
                  xdmat_new = wdmat_new, 
                  xdmat_cross = wdmat_cross,
                  x = w_t,
                  x_new = w_new,
                  tau2 = object$tau2_y[t],
                  theta = object$theta_y[t], 
                  g = g, 
                  v = object$v,
                  s2 = settings$lite, 
                  sigma = !settings$lite, 
                  grad = settings$grad,
                  grad_enhance = grad_enhance,
                  nsamples = ifel(samples_only, settings$nper, 0))
        
        if (samples_only) {
          row_indx <- (j*settings$nper - settings$nper + 1):(j*settings$nper)
          out$samples[row_indx, ] <- k$samples
          if (settings$grad) {
            out$grad_samples[row_indx, , i] <- 0
            for (j in 1:D) {
              out$grad_samples[row_indx, , i] <- out$grad_samples[row_indx, , i] + k$grad_samples[, , i] * dwdx[, j, i]
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
    } # end of foreach statement
    
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
        if (layers == 3)  z_new_store <- abind::abind(lapply(result, with, eval(parse(text = "z_new"))), along = 1)
      }
      if (settings$EI) ei_sum <- Reduce("+", lapply(result, with, eval(parse(text = "ei_sum"))))
      if (!is.null(settings$entropy_limit)) ent_sum <- Reduce("+", lapply(result, with, 
                                                                 eval(parse(text = "ent_sum"))))
    }
  } # end of else statement
  
  # Add variables to the output list
  if (samples_only) {
    if (settings$grad) {
      return(list(y = samples, dydx = ifel(d == 1, grad_samples[, , 1], grad_samples)))
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
    if (settings$grad) {
      object$grad_mean <- matrix(nrow = n_new, ncol = d)
      object$grad_s2 <- grad_s2_sum / object$nmcmc
      for (i in 1:d) {
        object$grad_mean[, i] <- colMeans(as.matrix(grad_mu_t[, , i]))
        object$grad_s2[, i] <- object$grad_s2[, i] + ifel(object$nmcmc > 1, apply(grad_mu_t[, , i], 2, var), 0)
      }
    }
    if (settings$store_latent) {
      object$w_new <- w_new_store
      if (layers == 3) object$z_new <- z_new_store
    }
    if (settings$EI) object$EI <- ei_sum / object$nmcmc
    if (!is.null(settings$entropy_limit)) object$entropy <- drop(ent_sum / object$nmcmc)
    toc <- proc.time()[[3]]
    object$time <- object$time + unname(toc - tic)
    return(object)
  }
}
