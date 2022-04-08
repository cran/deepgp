
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   predict.gp
#   predict.dgp2
#   predict.dgp3

# Define Predict for S3 Objects -----------------------------------------------
#' @name predict
#' @title Predict posterior mean and variance/covariance
#' @description Acts on a \code{gp}, \code{dgp2}, or \code{dgp3} object.
#'     Calculates posterior mean and variance/covariance over specified input 
#'     locations.  Optionally calculates expected improvement (EI) over 
#'     candidate inputs.  Optionally utilizes SNOW parallelization.
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
#' @param EI logical indicating whether to calculate expected improvement 
#'        (for minimizing the response)
#' @param cores number of cores to utilize in parallel, defaults to available
#'        cores minus one
#' @param m size of Vecchia conditioning sets (only for fits with 
#'        \code{vecchia = TRUE}), defaults to the \code{m} used for MCMC
#' @param ... N/A
#' @return object of the same class with the following additional elements:
#' \itemize{
#'   \item \code{x_new}: copy of predictive input locations
#'   \item \code{mean}: predicted posterior mean, indices correspond to 
#'         \code{x_new} locations
#'   \item \code{s2}: predicted point-wise variances, indices correspond to 
#'         \code{x_new} locations (only returned when \code{lite = TRUE})
#'   \item \code{s2_smooth}: predicted point-wise variances with \code{g} 
#'         removed, indices correspond to \code{x_new} locations (only returned 
#'         when \code{lite = TRUE})
#'   \item \code{Sigma}: predicted posterior covariance, indices correspond to 
#'         \code{x_new} locations (only returned when \code{lite = FALSE})
#'   \item \code{Sigma_smooth}: predicted posterior covariance with \code{g} 
#'         removed from the diagonal (only returned when \code{lite = FALSE})
#'   \item \code{EI}: vector of expected improvement values, indices correspond 
#'         to \code{x_new} locations (only returned when \code{EI = TRUE})
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
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." \emph{Technometrics, to appear;} arXiv:2012.08015. 
#'     \cr\cr
#' Sauer, A, A Cooper, and RB Gramacy. 2022. "Vecchia-approximated Deep Gaussian
#'     Processes for Computer Experiments." \emph{pre-print on arXiv:2204.02904} 
#'     \cr\cr   
#' Gramacy, RB. \emph{Surrogates: Gaussian Process Modeling, Design, and 
#'     Optimization for the Applied Sciences}. Chapman Hall, 2020.
#' 
#' @examples 
#' # See "fit_one_layer", "fit_two_layer", or "fit_three_layer"
#' # for an example
#' 
#' @rdname predict
NULL

# Predict One Layer -----------------------------------------------------------
#' @rdname predict
#' @export

predict.gp <- function(object, x_new, lite = TRUE, EI = FALSE, 
                       cores = detectCores() - 1, ...) {
  
  tic <- proc.time()[3]
  object <- clean_prediction(object) # remove previous predictions if present
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  dx <- sq_dist(object$x)
  d_new <- sq_dist(object$x_new)
  d_cross <- sq_dist(object$x_new, object$x)
  
  # Prepare clusters
  iters <- 1:object$nmcmc
  if (cores == 1) {
    chunks <- list(iters)
  } else chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
  if (cores > detectCores()) warning('cores is greater than available nodes')
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  thread <- NULL
  result <- foreach(thread = 1:cores) %dopar% {
    out <- list()
    out$mu_t <- matrix(nrow = n_new, ncol = length(chunks[[thread]]))
    if (lite) {
      out$s2_sum <- rep(0, times = n_new)
    } else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    if (EI) out$ei_sum <- rep(0, times = n_new)
    
    j <- 1
    for(t in chunks[[thread]]) {
      k <- krig(object$y, dx, d_new, d_cross, object$theta[t], object$g[t], 
                object$tau2[t], s2 = lite, sigma = !lite, f_min = EI, 
                v = object$v)
      out$mu_t[, j] <- k$mean
      if (lite) {
        out$s2_sum <- out$s2_sum + k$s2
      } else out$sigma_sum <- out$sigma_sum + k$sigma
      if (EI) {
        if (lite) {
          sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
        out$ei_sum <- out$ei_sum + exp_improv(k$mean, sig2, k$f_min)
      }
      j <- j + 1
    } # end of t for loop
    return(out)
  } # end of foreach loop
  
  stopCluster(cl)
    
  # Group elements out of list
  mu_t <- do.call(cbind, lapply(result, with, eval(parse(text = "mu_t"))))
  if (lite) {
    s2_sum <- Reduce("+", lapply(result, with, eval(parse(text = "s2_sum"))))
  } else {
    sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum"))))
  }
  if (EI) ei_sum <- Reduce("+", lapply(result, with, eval(parse(text = "ei_sum"))))
  
  # Add variables to the output list
  mu_cov <- cov(t(mu_t))
  object$mean <- rowMeans(mu_t)
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

# Predict Two Layer -----------------------------------------------------------
#' @rdname predict
#' @export

predict.dgp2 <- function(object, x_new, lite = TRUE, store_latent = FALSE, 
                         mean_map = TRUE, EI = FALSE, 
                         cores = detectCores() - 1, ...) {
  
  tic <- proc.time()[3]
  object <- clean_prediction(object) # remove previous predictions if present
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  D <- ncol(object$w[[1]])
  dx <- sq_dist(object$x)
  d_new <- sq_dist(object$x_new)
  d_cross <- sq_dist(object$x_new, object$x)
  
  # Prepare clusters
  iters <- 1:object$nmcmc
  if (cores == 1) {
    chunks <- list(iters)
  } else chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
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
    } else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    if (EI) out$ei_sum <- rep(0, times = n_new)
  
    j <- 1
    for(t in chunks[[thread]]) {
      w_t <- object$w[[t]]
      
      # Map x_new to w_new (separately for each dimension)
      w_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        if (mean_map) {
          k <- krig(w_t[, i], dx, NULL, d_cross, object$theta_w[t, i], 
                    g = eps, v = object$v)
          w_new[, i] <- k$mean
        } else {
          k <- krig(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], 
                    g = eps, sigma = TRUE, v = object$v)
          w_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      if (store_latent) out$w_new[[j]] <- w_new
      
      # Map w_new to y
      k <- krig(object$y, sq_dist(w_t), sq_dist(w_new), 
                sq_dist(w_new, w_t), object$theta_y[t], object$g[t], 
                object$tau2[t], s2 = lite, sigma = !lite, f_min = EI,
                v = object$v)
      out$mu_t[, j] <- k$mean
      if (lite) {
        out$s2_sum <- out$s2_sum + k$s2
      } else out$sigma_sum <- out$sigma_sum + k$sigma
      if (EI) {
        if (lite) {
          sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
        out$ei_sum <- out$ei_sum + exp_improv(k$mean, sig2, k$f_min)
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
  } else {
    sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum"))))
  }
  if (store_latent) w_new <- unlist(lapply(result, with, eval(parse(text = "w_new"))), 
                                         recursive = FALSE)
  if (EI) ei_sum <- Reduce("+", lapply(result, with, eval(parse(text = "ei_sum"))))

  # Add variables to the output list
  mu_cov <- cov(t(mu_t))
  object$mean <- rowMeans(mu_t)
  if (store_latent) object$w_new <- w_new
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

# Predict Three Layer ---------------------------------------------------------
#' @rdname predict
#' @export

predict.dgp3 <- function(object, x_new, lite = TRUE, store_latent = FALSE, 
                         mean_map = TRUE, EI = FALSE, 
                         cores = detectCores() - 1, ...) {
  
  tic <- proc.time()[3]
  object <- clean_prediction(object) # remove previous predictions if present
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  object$x_new <- x_new
  n_new <- nrow(object$x_new)
  D <- ncol(object$z[[1]])
  dx <- sq_dist(object$x)
  d_new <- sq_dist(object$x_new)
  d_cross <- sq_dist(object$x_new, object$x)
  
  # Prepare clusters
  iters <- 1:object$nmcmc
  if (cores == 1) {
    chunks <- list(iters)
  } else chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
  if (cores > detectCores()) warning("cores is greater than available nodes")
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
      out$s2_sum <- rep(0, times = n_new)
    } else out$sigma_sum <- matrix(0, nrow = n_new, ncol = n_new)
    if (EI) out$ei_sum <- rep(0, times = n_new)
    
    j <- 1
    for (t in chunks[[thread]]) {
      z_t <- object$z[[t]]
      w_t <- object$w[[t]]
      
      # Map x_new to z_new (separately for each dimension)
      z_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        if (mean_map) {
          k <- krig(z_t[, i], dx, NULL, d_cross, object$theta_z[t, i], 
                    g = eps, v = object$v)
          z_new[, i] <- k$mean
        } else {
          k <- krig(z_t[, i], dx, d_new, d_cross, object$theta_z[t, i], 
                    g = eps, sigma = TRUE, v = object$v)
          z_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      if (store_latent) out$z_new[[j]] <- z_new
      
      # Map z_new to w_new (separately for each dimension)
      w_new <- matrix(nrow = n_new, ncol = D)
      for (i in 1:D) {
        if (mean_map) { 
          k <- krig(w_t[, i], sq_dist(z_t), NULL, sq_dist(z_new, z_t),
                    object$theta_w[t, i], g = eps, v = object$v)
          w_new[, i] <- k$mean
        } else {
          k <- krig(w_t[, i], sq_dist(z_t), sq_dist(z_new), sq_dist(z_new, z_t),
                    object$theta_w[t, i], g = eps, sigma = TRUE,
                    v = object$v)
          w_new[, i] <- mvtnorm::rmvnorm(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      if (store_latent) out$w_new[[j]] <- w_new
      
      # Map w_new to y
      k <- krig(object$y, sq_dist(w_t), sq_dist(w_new), sq_dist(w_new, w_t),
                object$theta_y[t], object$g[t], object$tau2[t], s2 = lite, sigma = !lite,
                f_min = EI, v = object$v)
      out$mu_t[, j] <- k$mean
      if (lite) {
        out$s2_sum <- out$s2_sum + k$s2
      } else out$sigma_sum <- out$sigma_sum + k$sigma
      if (EI) {
        if (lite) {
          sig2 <- k$s2 - (object$tau2[t] * object$g[t]) 
        } else sig2 <- diag(k$sigma) - (object$tau2[t] * object$g[t])
        out$ei_sum <- out$ei_sum + exp_improv(k$mean, sig2, k$f_min)
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
  } else {
    sigma_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma_sum"))))
  }
  if (store_latent) {
    z_new <- unlist(lapply(result, with, eval(parse(text = "z_new"))), recursive = FALSE)
    w_new <- unlist(lapply(result, with, eval(parse(text = "w_new"))), recursive = FALSE)
  }
  if (EI) ei_sum <- Reduce("+", lapply(result, with, eval(parse(text = "ei_sum"))))
  
  # Add variables to the output list
  mu_cov <- cov(t(mu_t))
  object$mean <- rowMeans(mu_t)
  if (store_latent) {
    object$z_new <- z_new
    object$w_new <- w_new
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