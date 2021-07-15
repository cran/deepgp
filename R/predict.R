
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   predict.gp
#   predict.dgp2
#   predict.dgp3

# Define Predict for S3 Objects -----------------------------------------------
#' @name predict
#' @title Predict posterior mean and variance/covariance
#' @description Acts on a "\code{gp}", "\code{dgp2}", or "\code{dgp3}" object.
#'     Calculates posterior mean and variance/covariance over specified input 
#'     locations.  Optionally utilizes SNOW parallelization.
#' 
#' @details All iterations in the object are used for prediction, so samples 
#'     should be burned-in.  Thinning the samples using \code{trim} will speed 
#'     up computation.  Posterior moments are calculated using conditional 
#'     expectation and variance.  As a default, only point-wise variance is 
#'     calculated.  Full covariance may be calculated using \code{lite = FALSE}.  
#'     The storage of means and point-wise variances for each individual 
#'     iteration (specified using \code{store_all = TRUE}) is required in order 
#'     to use \code{EI}.
#'     
#'     SNOW parallelization reduces computation time but requires significantly 
#'     more memory storage.  Use \code{cores = 1} if memory is limited.
#' 
#' @param object object from \code{fit_one_layer}, \code{fit_two_layer}, or 
#'        \code{fit_three_layer} with burn-in already removed
#' @param x_new matrix of predictive input locations
#' @param lite logical indicating whether to calculate only point-wise 
#'        variances (\code{lite = TRUE}), or full covariance 
#'        (\code{lite = FALSE})
#' @param store_all logical indicating whether to store mean and variance for 
#'        each iteration
#' @param mean_map denotes whether to map hidden layers using conditional mean 
#'        or a random sample from the full MVN distribution 
#'        ("\code{dgp2}" or "\code{dgp3}" only) 
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
#'   \item \code{mu_t}: matrix of posterior mean for each iteration, column 
#'         index corresponds to iteration and row index corresponds to 
#'         \code{x_new} location (only returned when \code{store_all = TRUE})
#'   \item \code{s2_t}: matrix of posterior point-wise variance for each 
#'         iteration, column index corresponds to iteration and row index 
#'         corresponds to \code{x_new} location (only returned when 
#'         \code{store_all = TRUE})
#'   \item \code{w_new}: list of hidden layer mappings, list index corresponds 
#'         to iteration and row index corresponds to \code{x_new} location 
#'         ("\code{dgp2}" and "\code{dgp3}" only)
#'   \item \code{z_new}: list of hidden layer mappings, list index corresponds 
#'         to iteration and row index corresponds to \code{x_new} location 
#'         ("\code{dgp3}" only) 
#' }
#' Computation time is added to the computation time of the existing object.
#' 
#' @references 
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." arXiv:2012.08015. \cr\cr
#' Gramacy, RB. \emph{Surrogates: Gaussian Process Modeling, Design, and 
#'     Optimization for the Applied Sciences}. Chapman Hall, 2020.
#' 
#' @examples 
#' # See "deepgp-package", "fit_one_layer", "fit_two_layer", or 
#' # "fit_three_layer" for an example
#' 
#' @rdname predict
NULL

# Predict One Layer Function --------------------------------------------------
#' @rdname predict
#' @export

predict.gp <- function(object, x_new, lite = TRUE, store_all = FALSE, 
                       cores = 1, ...) {
  
  tic <- proc.time()[3]
  
  # check that x_new is a matrix
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  
  object$x_new <- x_new
  
  m <- nrow(object$x_new)
  dx <- sq_dist(object$x)
  d_new <- sq_dist(object$x_new)
  d_cross <- sq_dist(object$x_new, object$x)
  
  if (cores == 1) { # running on a single core, no parallelization
    
    tau2 <- vector(length = object$nmcmc)
    mu_t <- matrix(nrow = m, ncol = object$nmcmc)
    if (lite) {
      s2_t <- matrix(nrow = m, ncol = object$nmcmc)
    } else sigma_t_sum <- matrix(0, nrow = m, ncol = m) # full covariance
    
    for(t in 1:object$nmcmc) {
      # map x_new to mu_t and sigma_t
      k <- krig(object$y, dx, d_new, d_cross, object$theta[t], object$g[t], 
                s2 = lite, sigma = !lite)
      tau2[t] <- k$tau2
      mu_t[, t] <- k$mean
      if (lite) {
        s2_t[, t] <- k$s2
      } else sigma_t_sum <- sigma_t_sum + k$sigma
    }
    
  } else { # use foreach to run in parallel
    
    # prepare parallel clusters
    if (cores > detectCores()) warning('cores is greater than available nodes')
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    result <- foreach(t = 1:object$nmcmc) %dopar% {
      k <- krig(object$y, dx, d_new, d_cross, object$theta[t], object$g[t], 
                s2 = lite, sigma = !lite)
      return(k)
    }
    
    stopCluster(cl)
    
    # group elements out of the list
    mu_t <- sapply(result, with, eval(parse(text = "mean")))
    tau2 <- sapply(result, with, eval(parse(text = "tau2")))
    if (lite) {
      s2_t <- sapply(result, with, eval(parse(text = "s2")))
    } else sigma_t_sum <- Reduce("+", lapply(result, with, 
                                             eval(parse(text = "sigma"))))
  } # end of else statement
  
  # calculate expected value and covariance of means
  mu_y <- rowMeans(mu_t)
  mu_cov <- cov(t(mu_t))
  
  # store mean and s2 across iterations if requested
  if (store_all) {
    object$mu_t <- mu_t
    if (lite) object$s2_t <- s2_t
    # no option available to store all sigma matrices at this time
  }
  
  # add variables to the output list
  object$tau2 <- tau2
  object$mean <- mu_y
  if (lite) { 
    object$s2 <- rowSums(s2_t) / object$nmcmc + diag(mu_cov)
    object$s2_smooth <- object$s2 - mean(object$g * object$tau2)
  } else {
    object$Sigma <- sigma_t_sum / object$nmcmc + mu_cov
    object$Sigma_smooth <- object$Sigma - diag(mean(object$g * object$tau2), m)
  }
  
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  
  return(object)
}

# Predict Two Layer Function --------------------------------------------------
#' @rdname predict
#' @export

predict.dgp2 <- function(object, x_new, lite = TRUE, store_all = FALSE,
                         mean_map = TRUE, cores = 1, ...) {
  
  tic <- proc.time()[3]
  
  # check that x_new is a matrix
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  
  object$x_new <- x_new
  
  m <- nrow(object$x_new)
  D <- ncol(object$w[[1]])
  dx <- sq_dist(object$x)
  d_new <- sq_dist(object$x_new)
  d_cross <- sq_dist(object$x_new, object$x)
  
  if (cores == 1) { # running on a single core, no parallelization
    w_new <- list()
    tau2 <- vector(length = object$nmcmc)
    mu_t <- matrix(nrow = m, ncol = object$nmcmc)
    if (lite) {
      s2_t <- matrix(nrow = m, ncol = object$nmcmc)
    } else sigma_t_sum <- matrix(0, nrow = m, ncol = m)
    
    for(t in 1:object$nmcmc) {
      
      w_t <- object$w[[t]]
      
      # map x_new to w_new (separately for each dimension)
      w_new[[t]] <- matrix(nrow = m, ncol = D)
      for (i in 1:D){
        if (mean_map) {
          k <- krig(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], 
                    g = NULL, s2 = FALSE, sigma = FALSE,
                    tau2 = FALSE)
          w_new[[t]][, i] <- k$mean
        } else {
          k <- krig(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], 
                    g = NULL, tau2 = FALSE)
          w_new[[t]][, i] <- rand_mvn(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      
      # map w_new to mu_t and sigma_t
      k <- krig(object$y, sq_dist(w_t), sq_dist(w_new[[t]]), 
                sq_dist(w_new[[t]], w_t), object$theta_y[t], object$g[t], 
                s2 = lite, sigma = !lite)
      tau2[t] <- k$tau2
      mu_t[, t] <- k$mean
      if (lite) {
        s2_t[, t] <- k$s2
      } else sigma_t_sum <- sigma_t_sum + k$sigma
    } # end of t for loop
    
  } else { # use foreach to run in parallel
    
    # prepare parallel clusters
    if (cores > detectCores()) warning('cores is greater than available nodes')
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    result <- foreach(t = 1:object$nmcmc) %dopar% {
      
      w_t <- object$w[[t]]
      
      # map x_new to w_new (separately for each dimension)
      w_new <- matrix(nrow = m, ncol = D)
      for (i in 1:D) {
        if (mean_map) {
          k <- krig(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], 
                    g = NULL, sigma = FALSE, tau2 = FALSE)
          w_new[, i] <- k$mean
        } else {
          k <- krig(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], 
                    g = NULL, tau2 = FALSE)
          w_new[, i] <- rand_mvn(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      
      # map w_new to mu_t and sigma_t
      k <- krig(object$y, sq_dist(w_t), sq_dist(w_new), sq_dist(w_new, w_t),
                object$theta_y[t], object$g[t], s2 = lite, sigma = !lite)
      k$w_new <- w_new
      return(k)
    } # end of foreach statement
    
    stopCluster(cl)
    
    # group elements out of the list
    mu_t <- sapply(result, with, eval(parse(text = "mean")))
    tau2 <- sapply(result, with, eval(parse(text = "tau2")))
    if (lite) {
      s2_t <- sapply(result, with, eval(parse(text = "s2")))
    } else sigma_t_sum <- Reduce("+", lapply(result, with, 
                                             eval(parse(text = "sigma"))))
    w_new <- lapply(result, with, eval(parse(text = "w_new")))
  } # end of else statement
  
  # calculate expected value and covariance of means
  mu_y <- rowMeans(mu_t)
  mu_cov <- cov(t(mu_t))
  
  # store mean and s2 across iterations if requested
  if (store_all) {
    object$mu_t <- mu_t
    if (lite) object$s2_t <- s2_t
    # no option available to store all sigma matrices at this time
  }
  
  # add variables to the output list
  object$w_new <- w_new
  object$tau2 <- tau2
  object$mean <- mu_y
  if (lite) {
    object$s2 <- rowSums(s2_t) / object$nmcmc + diag(mu_cov)
    object$s2_smooth <- object$s2 - mean(object$g * object$tau2)
  } else {
    object$Sigma <- sigma_t_sum
    object$Sigma_smooth <- object$Sigma - diag(mean(object$g * object$tau2), m)
  }
  
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  
  return(object)
}

# Predict Three Layer Function ------------------------------------------------
#' @rdname predict
#' @export

predict.dgp3 <- function(object, x_new, lite = TRUE, store_all = FALSE, 
                         mean_map = TRUE, cores = 1, ...) {
  
  tic <- proc.time()[3]
  
  # check that x_new is a matrix
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  
  object$x_new <- x_new
  
  m <- nrow(object$x_new)
  D <- ncol(object$z[[1]])
  dx <- sq_dist(object$x)
  d_new <- sq_dist(object$x_new)
  d_cross <- sq_dist(object$x_new, object$x)
  
  if (cores == 1) { # running on a single core, no parallelization
    z_new <- list()
    w_new <- list()
    tau2 <- vector(length = object$nmcmc)
    mu_t <- matrix(nrow = m, ncol = object$nmcmc)
    if (lite) {
      s2_t <- matrix(nrow = m, ncol = object$nmcmc)
    } else sigma_t_sum <- matrix(0, nrow = m, ncol = m)
    
    for(t in 1:object$nmcmc) {
      
      z_t <- object$z[[t]]
      w_t <- object$w[[t]]
      
      # map x_new to z_new (separately for each dimension)
      z_new[[t]] <- matrix(nrow = m, ncol = D)
      for (i in 1:D) {
        if (mean_map) {
          k <- krig(z_t[, i], dx, d_new, d_cross, object$theta_z[t, i], 
                    g = NULL, sigma = FALSE, tau2 = FALSE)
          z_new[[t]][, i] <- k$mean
        } else {
          k <- krig(z_t[, i], dx, d_new, d_cross, object$theta_z[t, i], 
                    g = NULL, tau2 = FALSE)
          z_new[[t]][, i] <- rand_mvn(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      
      # map z_new to w_new (separately for each dimension)
      w_new[[t]] <- matrix(nrow = m, ncol = D)
      for (i in 1:D) {
        if (mean_map) { 
          k <- krig(w_t[, i], sq_dist(z_t), sq_dist(z_new[[t]]), 
                    sq_dist(z_new[[t]], z_t), object$theta_w[t, i], g = NULL, 
                    sigma = FALSE, tau2 = FALSE)
          w_new[[t]][, i] <- k$mean
        } else {
          k <- krig(w_t[, i], sq_dist(z_t), sq_dist(z_new[[t]]), 
                    sq_dist(z_new[[t]], z_t), object$theta_w[t, i], g = NULL, 
                    tau2 = FALSE)
          w_new[[t]][, i] <- rand_mvn(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      
      # map w_new to mu_t and sigma_t
      k <- krig(object$y, sq_dist(w_t), sq_dist(w_new[[t]]), 
                sq_dist(w_new[[t]], w_t), object$theta_y[t], object$g[t], 
                s2 = lite, sigma = !lite)
      tau2[t] <- k$tau2
      mu_t[, t] <- k$mean
      if (lite) {
        s2_t[, t] <- k$s2
      } else sigma_t_sum <- sigma_t_sum + k$sigma
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
      z_new <- matrix(nrow = m, ncol = D)
      for (i in 1:D) {
        if (mean_map) {
          k <- krig(z_t[, i], dx, d_new, d_cross, object$theta_z[t, i], 
                    g = NULL, sigma = FALSE, tau2 = FALSE)
          z_new[, i] <- k$mean
        } else {
          k <- krig(z_t[, i], dx, d_new, d_cross, object$theta_z[t, i], 
                    g = NULL, tau2 = FALSE)
          z_new[, i] <- rand_mvn(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      
      # map z_new to w_new (separately for each dimension)
      w_new <- matrix(nrow = m, ncol = D)
      for (i in 1:D) {
        if (mean_map) {
          k <- krig(w_t[, i], sq_dist(z_t), sq_dist(z_new), sq_dist(z_new, z_t),
                    object$theta_w[t, i], g = NULL, sigma = FALSE, tau2 = FALSE)
          w_new[, i] <- k$mean
        } else {
          k <- krig(w_t[, i], sq_dist(z_t), sq_dist(z_new), sq_dist(z_new, z_t),
                    object$theta_w[t, i], g = NULL, tau2 = FALSE)
          w_new[, i] <- rand_mvn(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      
      # map w_new to mu_t and sigma_t
      k <- krig(object$y, sq_dist(w_t), sq_dist(w_new), sq_dist(w_new, w_t),
                object$theta_y[t], object$g[t], s2 = lite, sigma = !lite)
      k$z_new <- z_new
      k$w_new <- w_new
      return(k)
    } # end of foreach statement
    
    stopCluster(cl)
    
    # group elements out of the list
    mu_t <- sapply(result, with, eval(parse(text = "mean")))
    tau2 <- sapply(result, with, eval(parse(text = "tau2")))
    if (lite) {
      s2_t <- sapply(result, with, eval(parse(text = "s2")))
    } else sigma_t_sum <- Reduce("+", lapply(result, with, 
                                             eval(parse(text = "sigma"))))
    z_new <- lapply(result, with, eval(parse(text = "z_new")))
    w_new <- lapply(result, with, eval(parse(text = "w_new")))
    
  } # end of else statement
  
  # calculate expected value and covariance of means
  mu_y <- rowMeans(mu_t)
  mu_cov <- cov(t(mu_t))
  
  # store mean and s2 across iterations if requested
  if (store_all) {
    object$mu_t <- mu_t
    if (lite) object$s2_t <- s2_t
    # no option available to store all sigma matrices at this time
  }
  
  # add variables to the output list
  object$z_new <- z_new
  object$w_new <- w_new
  object$tau2 <- tau2
  object$mean <- mu_y
  if (lite) {
    object$s2 <- rowSums(s2_t) / object$nmcmc + diag(mu_cov)
    object$s2_smooth <- object$s2 - mean(object$g * object$tau2)
  } else {
    object$Sigma <- sigma_t_sum / object$nmcmc + mu_cov
    object$Sigma_smooth <- object$Sigma - diag(mean(object$g * object$tau2), m)
  }
  
  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  return(object)
}
