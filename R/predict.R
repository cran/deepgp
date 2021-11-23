
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   predict.gp
#   predict.dgp2
#   predict.dgp3

# Define Predict for S3 Objects --------------------------------------------------
#' @name predict
#' @title Predict posterior mean and variance/covariance
#' @description Acts on a \code{gp}, \code{dgp2}, or \code{dgp3} object.
#'     Calculates posterior mean and variance/covariance over specified input 
#'     locations.  Optionally calculates expected improvement (EI) over candidate
#'     inputs.  Optionally utilizes SNOW parallelization.
#' 
#' @details All iterations in the object are used for prediction, so samples 
#'     should be burned-in.  Thinning the samples using \code{trim} will speed 
#'     up computation.  Posterior moments are calculated using conditional 
#'     expectation and variance.  As a default, only point-wise variance is 
#'     calculated.  Full covariance may be calculated using \code{lite = FALSE}. 
#'     
#'     Expected improvement is calculated with the goal of minimizing the 
#'     response.  See Chapter 7 of Gramacy (2020) for details.
#'     
#'     SNOW parallelization reduces computation time but requires significantly 
#'     more memory storage.  Use \code{cores = 1} if memory is limited.
#' 
#' @param object object from \code{fit_one_layer}, \code{fit_two_layer}, or 
#'        \code{fit_three_layer} with burn-in already removed
#' @param x_new matrix of predictive input locations
#' @param lite logical indicating whether to calculate only point-wise 
#'        variances (\code{lite = TRUE}) or full covariance 
#'        (\code{lite = FALSE})
#' @param EI logical indicating whether to calculate expected improvement 
#'        (for minimizing the response)
#' @param store_latent logical indicating whether to store and return mapped 
#'        values of latent layers (\code{dgp2} or \code{dgp3} only)
#' @param mean_map logical indicating whether to map hidden layers using 
#'        conditional mean (\code{mean_map = TRUE}) or using a random sample
#'        from the full MVN distribution (\code{dgp2} or \code{dgp3} only) 
#' @param cores number of cores to utilize in parallel, by default no 
#'        parallelization is used
#' @param ... N/A
#' @return object of the same class with the following additional elements:
#' \itemize{
#'   \item \code{x_new}: copy of predictive input locations
#'   \item \code{tau2}: vector of tau2 estimates (governing the magnitude of 
#'         the covariance)
#'   \item \code{mean}: predicted posterior mean, indices correspond to 
#'         \code{x_new} location
#'   \item \code{s2}: predicted point-wise variances, indices correspond to 
#'         \code{x_new} location (only returned when \code{lite = TRUE})
#'   \item \code{s2_smooth}: predicted point-wise variances with \code{g} 
#'         removed, indices correspond to \code{x_new} location (only returned 
#'         when \code{lite = TRUE})
#'   \item \code{Sigma}: predicted posterior covariance, indices correspond to 
#'         \code{x_new} location (only returned when \code{lite = FALSE})
#'   \item \code{Sigma_smooth}: predicted posterior covariance with \code{g} 
#'         removed from the diagonal (only returned when \code{lite = FALSE})
#'   \item \code{EI}: vector of expected improvement values, indices correspond 
#'         to \code{x_new} location (only returned when \code{EI = TRUE})
#'   \item \code{w_new}: list of hidden layer mappings (only returned when 
#'         \code{store_latent = TRUE}), list index corresponds to iteration and 
#'         row index corresponds to \code{x_new} location (\code{dgp2} and 
#'         \code{dgp3} only)
#'   \item \code{z_new}: list of hidden layer mappings (only returned when 
#'         \code{store_latent = TRUE}), list index corresponds to iteration and 
#'         row index corresponds to \code{x_new} location (\code{dgp3} only) 
#' }
#' Computation time is added to the computation time of the existing object.
#' 
#' @references 
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." \emph{Technometrics, to appear;} arXiv:2012.08015. 
#'     \cr\cr
#' Gramacy, RB. \emph{Surrogates: Gaussian Process Modeling, Design, and 
#'     Optimization for the Applied Sciences}. Chapman Hall, 2020.
#' 
#' @examples 
#' # See "deepgp-package", "fit_one_layer", "fit_two_layer", or "fit_three_layer"
#' # for an example
#' 
#' @rdname predict
NULL

# Predict One Layer Function --------------------------------------------------
#' @rdname predict
#' @export

predict.gp <- function(object, x_new, lite = TRUE, EI = FALSE, cores = 1, ...) {
  
  tic <- proc.time()[3]
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  dx <- sq_dist(object$x)
  d_new <- sq_dist(object$x_new)
  d_cross <- sq_dist(object$x_new, object$x)
  
  if (cores == 1) { # running on a single core, no parallelization
    
    tau2 <- vector(length = object$nmcmc)
    mu_t <- matrix(nrow = n_new, ncol = object$nmcmc)
    if (lite) {
      s2_sum <- vector(length = n_new)
    } else sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    if (EI) ei_sum <- vector(length = n_new)
    
    # calculate predictions for each candidate MCMC iteration
    for(t in 1:object$nmcmc) {
      k <- krig(object$y, dx, d_new, d_cross, object$theta[t], object$g[t], 
                s2 = lite, sigma = !lite, tau2 = TRUE, f_min = EI, 
                cov = object$cov, v = object$v)
      tau2[t] <- k$tau2
      mu_t[, t] <- k$mean
      if (lite) {
        s2_sum <- s2_sum + k$s2
      } else sigma_sum <- sigma_sum + k$sigma
      if (EI) {
        if (lite) {
          sig2 <- k$s2 - (k$tau2 * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (k$tau2 * object$g[t])
        ei_sum <- ei_sum + exp_improv(k$mean, sig2, k$f_min)
      }
    }
    
  } else { # use foreach to run in parallel
    
    # prepare parallel clusters
    if (cores > detectCores()) warning('cores is greater than available nodes')
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    # calculate predictions for each candidate MCMC iteration
    result <- foreach(t = 1:object$nmcmc) %dopar% {
      k <- krig(object$y, dx, d_new, d_cross, object$theta[t], object$g[t], 
                s2 = lite, sigma = !lite, tau2 = TRUE, f_min = EI,
                cov = object$cov, v = object$v)
      if (EI) {
        if (lite) {
          sig2 <- k$s2 - (k$tau2 * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (k$tau2 * object$g[t])
        k$ei <- exp_improv(k$mean, sig2, k$f_min)
      }
      return(k)
    }
    
    stopCluster(cl)
    
    # group elements out of the list
    mu_t <- sapply(result, with, eval(parse(text = "mean")))
    tau2 <- sapply(result, with, eval(parse(text = "tau2")))
    if (lite) {
      s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2"))))
    } else {
      sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma"))))
    }
    if (EI) ei_sum <- Reduce("+", lapply(result, with, eval(parse(text = "ei"))))
  } # end of else statement
  
  # add variables to the output list
  mu_cov <- cov(t(mu_t))
  object$mean <- rowMeans(mu_t)
  object$tau2 <- tau2
  if (lite) { 
    object$s2 <- s2_sum / object$nmcmc + diag(mu_cov)
    object$s2_smooth <- object$s2 - mean(object$g * object$tau2)
  } else {
    object$Sigma <- sigma_sum / object$nmcmc + mu_cov
    object$Sigma_smooth <- object$Sigma - diag(mean(object$g * object$tau2), n_new)
  }
  if (EI) object$EI <- ei_sum / object$nmcmc
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  
  return(object)
}

# Predict Two Layer Function --------------------------------------------------
#' @rdname predict
#' @export

predict.dgp2 <- function(object, x_new, lite = TRUE, store_latent = FALSE, 
                         mean_map = TRUE, EI = FALSE, cores = 1, ...) {
  
  tic <- proc.time()[3]
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  D <- ncol(object$w[[1]])
  dx <- sq_dist(object$x)
  d_new <- sq_dist(object$x_new)
  d_cross <- sq_dist(object$x_new, object$x)
  
  if (cores == 1) { # running on a single core, no parallelization
    
    if (store_latent) w_new_list <- list()
    tau2 <- vector(length = object$nmcmc)
    mu_t <- matrix(nrow = n_new, ncol = object$nmcmc)
    if (lite) {
      s2_sum <- vector(length = n_new)
    } else sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    if (EI) ei_sum <- vector(length = n_new)
    
    # calculate predictions for each candidate MCMC iteration
    for(t in 1:object$nmcmc) {
      w_t <- object$w[[t]]
      
      # map x_new to w_new (separately for each dimension)
      w_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        if (mean_map) {
          k <- krig(w_t[, i], dx, NULL, d_cross, object$theta_w[t, i], 
                    g = eps, cov = object$cov, v = object$v)
          w_new[, i] <- k$mean
        } else {
          k <- krig(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], 
                    g = eps, sigma = TRUE, cov = object$cov, v = object$v)
          w_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      if (store_latent) w_new_list[[t]] <- w_new
      
      # map w_new to y
      k <- krig(object$y, sq_dist(w_t), sq_dist(w_new), 
                sq_dist(w_new, w_t), object$theta_y[t], object$g[t], 
                s2 = lite, sigma = !lite, tau2 = TRUE, f_min = EI,
                cov = object$cov, v = object$v)
      tau2[t] <- k$tau2
      mu_t[, t] <- k$mean
      if (lite) {
        s2_sum <- s2_sum + k$s2
      } else sigma_sum <- sigma_sum + k$sigma
      if (EI) {
        if (lite) {
          sig2 <- k$s2 - (k$tau2 * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (k$tau2 * object$g[t])
        ei_sum <- ei_sum + exp_improv(k$mean, sig2, k$f_min)
      }
    } # end of t for loop
    
  } else { # use foreach to run in parallel
    
    # prepare parallel clusters
    if (cores > detectCores()) warning('cores is greater than available nodes')
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    result <- foreach(t = 1:object$nmcmc) %dopar% {
      w_t <- object$w[[t]]
      
      # map x_new to w_new (separately for each dimension)
      w_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        if (mean_map) {
          k <- krig(w_t[, i], dx, NULL, d_cross, object$theta_w[t, i], 
                    g = eps, cov = object$cov, v = object$v)
          w_new[, i] <- k$mean
        } else {
          k <- krig(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], 
                    g = eps, sigma = TRUE, cov = object$cov, v = object$v)
          w_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      
      # map w_new to y
      k <- krig(object$y, sq_dist(w_t), sq_dist(w_new), 
                sq_dist(w_new, w_t), object$theta_y[t], object$g[t], 
                s2 = lite, sigma = !lite, tau2 = TRUE, f_min = EI,
                cov = object$cov, v = object$v)
      if (store_latent) k$w_new <- w_new
      if (EI) {
        if (lite) {
          sig2 <- k$s2 - (k$tau2 * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (k$tau2 * object$g[t])
        k$ei <- exp_improv(k$mean, sig2, k$f_min)
      }
      return(k)
    } # end of foreach statement
    
    stopCluster(cl)
    
    # group elements out of the list
    mu_t <- sapply(result, with, eval(parse(text = "mean")))
    tau2 <- sapply(result, with, eval(parse(text = "tau2")))
    if (lite) {
      s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2"))))
    } else {
      sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma"))))
    }
    if (store_latent) w_new_list <- lapply(result, with, eval(parse(text = "w_new")))
    if (EI) ei_sum <- Reduce("+", lapply(result, with, eval(parse(text = "ei"))))
  } # end of else statement
  
  # add variables to the output list
  mu_cov <- cov(t(mu_t))
  object$mean <- rowMeans(mu_t)
  object$tau2 <- tau2
  if (store_latent) object$w_new <- w_new_list
  if (lite) {
    object$s2 <- s2_sum / object$nmcmc + diag(mu_cov)
    object$s2_smooth <- object$s2 - mean(object$g * object$tau2)
  } else {
    object$Sigma <- sigma_sum / object$nmcmc + mu_cov
    object$Sigma_smooth <- object$Sigma - diag(mean(object$g * object$tau2), n_new)
  }
  if (EI) object$EI <- ei_sum / object$nmcmc
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  
  return(object)
}

# Predict Three Layer Function ------------------------------------------------
#' @rdname predict
#' @export

predict.dgp3 <- function(object, x_new, lite = TRUE, store_latent = FALSE, 
                         mean_map = TRUE, EI = FALSE, cores = 1, ...) {
  
  tic <- proc.time()[3]
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  D <- ncol(object$z[[1]])
  dx <- sq_dist(object$x)
  d_new <- sq_dist(object$x_new)
  d_cross <- sq_dist(object$x_new, object$x)
  
  if (cores == 1) { # running on a single core, no parallelization
    if (store_latent) {
      z_new_list <- list()
      w_new_list <- list()
    }
    tau2 <- vector(length = object$nmcmc)
    mu_t <- matrix(nrow = n_new, ncol = object$nmcmc)
    if (lite) {
      s2_sum <- vector(length = n_new)
    } else sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    if (EI) ei_sum <- vector(length = n_new)
    
    # calculate predictions for each candidate MCMC iteration
    for (t in 1:object$nmcmc) {
      z_t <- object$z[[t]]
      w_t <- object$w[[t]]
      
      # map x_new to z_new (separately for each dimension)
      z_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        if (mean_map) {
          k <- krig(z_t[, i], dx, NULL, d_cross, object$theta_z[t, i], 
                    g = eps, cov = object$cov, v = object$v)
          z_new[, i] <- k$mean
        } else {
          k <- krig(z_t[, i], dx, d_new, d_cross, object$theta_z[t, i], 
                    g = eps, sigma = TRUE, cov = object$cov, v = object$v)
          z_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      if (store_latent) z_new_list[[t]] <- z_new
      
      # map z_new to w_new (separately for each dimension)
      w_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        if (mean_map) { 
          k <- krig(w_t[, i], sq_dist(z_t), NULL, sq_dist(z_new, z_t),
                    object$theta_w[t, i], g = eps, cov = object$cov, v = object$v)
          w_new[, i] <- k$mean
        } else {
          k <- krig(w_t[, i], sq_dist(z_t), sq_dist(z_new), sq_dist(z_new, z_t),
                    object$theta_w[t, i], g = eps, sigma = TRUE,
                    cov = object$cov, v = object$v)
          w_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      if (store_latent) w_new_list[[t]] <- w_new
      
      # map w_new to y
      k <- krig(object$y, sq_dist(w_t), sq_dist(w_new), sq_dist(w_new, w_t),
                object$theta_y[t], object$g[t], s2 = lite, sigma = !lite,
                tau2 = TRUE, f_min = EI, cov = object$cov, v = object$v)
      tau2[t] <- k$tau2
      mu_t[, t] <- k$mean
      if (lite) {
        s2_sum <- s2_sum + k$s2
      } else sigma_sum <- sigma_sum + k$sigma
      if (EI) {
        if (lite) {
          sig2 <- k$s2 - (k$tau2 * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (k$tau2 * object$g[t])
        ei_sum <- ei_sum + exp_improv(k$mean, sig2, k$f_min)
      }
    } # end of t for loop
    
  } else { # use foreach to run in parallel
    
    # prepare parallel clusters
    if (cores > detectCores()) warning('cores is greater than available nodes')
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    result <- foreach(t = 1:object$nmcmc) %dopar% {
      
      z_t <- object$z[[t]]
      w_t <- object$w[[t]]
      
      # map x_new to z_new (separately for each dimension)
      z_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        if (mean_map) {
          k <- krig(z_t[, i], dx, NULL, d_cross, object$theta_z[t, i], 
                    g = eps, cov = object$cov, v = object$v)
          z_new[, i] <- k$mean
        } else {
          k <- krig(z_t[, i], dx, d_new, d_cross, object$theta_z[t, i], 
                    g = eps, sigma = TRUE, cov = object$cov, v = object$v)
          z_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      
      # map z_new to w_new (separately for each dimension)
      w_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        if (mean_map) { 
          k <- krig(w_t[, i], sq_dist(z_t), NULL, sq_dist(z_new, z_t),
                    object$theta_w[t, i], g = eps, cov = object$cov, v = object$v)
          w_new[, i] <- k$mean
        } else {
          k <- krig(w_t[, i], sq_dist(z_t), sq_dist(z_new), sq_dist(z_new, z_t),
                    object$theta_w[t, i], g = eps, sigma = TRUE,
                    cov = object$cov, v = object$v)
          w_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      
      # map w_new to y
      k <- krig(object$y, sq_dist(w_t), sq_dist(w_new), sq_dist(w_new, w_t),
                object$theta_y[t], object$g[t], s2 = lite, sigma = !lite,
                tau2 = TRUE, f_min = EI, cov = object$cov, v = object$v)
      if (store_latent) {
        k$z_new <- z_new
        k$w_new <- w_new
      }
      if (EI) {
        if (lite) {
          sig2 <- k$s2 - (k$tau2 * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (k$tau2 * object$g[t])
        k$ei <- exp_improv(k$mean, sig2, k$f_min)
      }
      return(k)
    } # end of foreach statement
    
    stopCluster(cl)
    
    # group elements out of the list
    mu_t <- sapply(result, with, eval(parse(text = "mean")))
    tau2 <- sapply(result, with, eval(parse(text = "tau2")))
    if (lite) {
      s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2"))))
    } else {
      sigma_t_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma"))))
    }
    if (store_latent) {
      z_new_list <- lapply(result, with, eval(parse(text = "z_new")))
      w_new_list <- lapply(result, with, eval(parse(text = "w_new")))
    }
    if (EI) ei_sum <- Reduce("+", lapply(result, with, eval(parse(text = "ei"))))
  } # end of else statement
  
  # add variables to the output list
  mu_cov <- cov(t(mu_t))
  object$mean <- rowMeans(mu_t)
  object$tau2 <- tau2
  if (store_latent) {
    object$z_new <- z_new_list
    object$w_new <- w_new_list
  }
  if (lite) {
    object$s2 <- s2_sum / object$nmcmc + diag(mu_cov)
    object$s2_smooth <- object$s2 - mean(object$g * object$tau2)
  } else {
    object$Sigma <- sigma_sum / object$nmcmc + mu_cov
    object$Sigma_smooth <- object$Sigma - diag(mean(object$g * object$tau2), n_new)
  }
  if (EI) object$EI <- ei_sum / object$nmcmc
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  
  return(object)
}
