
# Function Contents -----------------------------------------------------------
# Internal:
#   logl: evaluates MVN log likelihood with zero mean
#   sample_g: conducts Metropolis Hastings sampling for nugget
#   sample_theta: conducts Metropolis Hastings sampling for theta
#   sample_w: conducts Elliptical Slice Sampling for w layer
#   sample_z: conducts Elliptical Slice Sampling for z layer
# External (see documentation below):
#   fit_one_layer
#   fit_two_layer
#   fit_three_layer
#   continue (S3 method for gp, dgp2, dgp3 classes)
#   trim (S3 method for gp, dgp2, dgp3 classes)
#   predict (S3 method for gp, dgp2, dgp3 classes)
#   plot (S3 method for gp, dgp2, dgp3 classes)

# Imported Functions ----------------------------------------------------------
#' @importFrom grDevices heat.colors
#' @importFrom graphics image lines matlines par plot points
#' @importFrom stats cov dgamma dnorm pnorm qnorm rnorm runif var
#' @importFrom parallel makeCluster detectCores stopCluster

# Log Likelihood Function -----------------------------------------------------
# Calculates log likelihood for multivariate normal distribution with zero mean

logl <- function(out_vec, in_dmat, g, theta, outer = TRUE) {
  
  K <- calc_K(in_dmat, theta, g)
  id <- invdet(K)

  if (outer) { # use profile log likelihood (with tau2 integrated out)
    n <- length(out_vec)
    logl <- (- n * 0.5) * log(t(out_vec) %*% id$Mi %*% (out_vec)) - 0.5 * id$ldet
  } else {
    logl <- (- 0.5) * id$ldet - 0.5 * (t(out_vec) %*% id$Mi %*% (out_vec))
  }
  return(c(logl))
}

# Sample Nugget Function ------------------------------------------------------
# Completes one iteration of Metropolis Hastings sampling for nugget

sample_g <- function(out_vec, in_dmat, g_t, theta, alpha, beta, l, u, 
                     ll_prev = NULL) {

  # Propose value
  g_star <- runif(1, min = l * g_t / u, max = u * g_t / l)

  # Compute acceptance threshold
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev)) 
    ll_prev <- logl(out_vec, in_dmat, g_t, theta, outer = TRUE)
  lpost_threshold <-  ll_prev + dgamma(g_t, alpha, beta, log = TRUE) + 
                      log(ru) - log(g_t) + log(g_star)
 
  ll_new <- logl(out_vec, in_dmat, g_star, theta, outer = TRUE)
  
  # Accept or reject
  if (ll_new + dgamma(g_star, alpha, beta, log = TRUE) > lpost_threshold) {
    return(list(g = g_star, ll = ll_new))
  } else{
    return(list(g = g_t, ll = ll_prev))
  }
}

# Sample Theta Function -------------------------------------------------------
# Completes one iteration of Metropolis Hastings sampling for length scale

sample_theta <- function(out_vec, in_dmat, g, theta_t, alpha, beta, l, u, 
                         outer, ll_prev = NULL) {

  # Propose value
  theta_star <- runif(1, min = l * theta_t / u, max = u * theta_t / l)

  # Compute acceptance threshold
  ru <- runif(1, min = 0, max = 1)
  if (is.null(ll_prev))
    ll_prev <- logl(out_vec, in_dmat, g, theta_t, outer)
  lpost_threshold <- ll_prev + dgamma(theta_t, alpha, beta, log = TRUE) + 
                      log(ru) - log(theta_t) + log(theta_star)
  
  ll_new <- logl(out_vec, in_dmat, g, theta_star, outer)

  # Accept or reject
  if (ll_new + dgamma(theta_star, alpha, beta, log = TRUE) > lpost_threshold) {
    return(list(theta = theta_star, ll = ll_new))
  } else{
    return(list(theta = theta_t, ll = ll_prev))
  }
}

# Elliptical Slice W Function -------------------------------------------------
# Completes one iteration of Elliptical Slice Sampling for a hidden layer

sample_w <- function(out_vec, w_t, w_t_dmat, in_dmat, g, theta_y, theta_w,
                     ll_prev = NULL) {

  D <- ncol(w_t) # dimension of hidden layer

  if (is.null(ll_prev)) 
    ll_prev <- logl(out_vec, w_t_dmat, g, theta_y, outer = TRUE)
  
  for (i in 1:D) { # separate sampling for each dimension of hidden layer

    # Draw from prior distribution
    w_prior <- rand_mvn(1, sigma = calc_K(in_dmat, theta_w[i]))

    # Initialize a and bounds on a
    a <- runif(1, min = 0, max = 2 * pi)
    amin <- a - 2 * pi
    amax <- a

    # Compute acceptance threshold - based on all dimensions of previous w
    ru <- runif(1, min = 0, max = 1)
    ll_threshold <- ll_prev + log(ru)

    # Calculate proposed values, accept or reject, repeat if necessary
    accept <- FALSE
    count <- 0
    w_prev <- w_t[, i] # store for re-proposal
    
    while (accept == FALSE) {
      count <- count + 1

      # Calculate proposed values and new likelihood
      w_t[, i] <- w_prev * cos(a) + w_prior * sin(a)
      dw <- sq_dist(w_t)
      new_logl <- logl(out_vec, dw, g, theta_y, outer = TRUE)
      
      # Accept or reject
      if (new_logl > ll_threshold) { # accept
        ll_prev <- new_logl
        accept <- TRUE
      } else { # reject
        # update the bounds on a and repeat
        if (a < 0) {
          amin <- a
        } else {
          amax <- a
        }
        a <- runif(1, amin, amax)
        if (count > 100) stop('reached maximum iterations of ESS')
      } # end of else statement
    } # end of while loop
  } # end of i for loop

  return(list(w = w_t, ll = ll_prev, dw = dw))
}

# Elliptical Slice Z Function -------------------------------------------------
# Completes one iteration of Elliptical Slice Sampling for a hidden layer

sample_z <- function(out_mat, z_t, z_t_dmat, in_dmat, g, theta_w, theta_z,
                     ll_prev = NULL) {

  D <- ncol(z_t) # dimension of hidden layer

  if (is.null(ll_prev)) {
    ll_prev <- 0
    for (j in 1:D)
      ll_prev <- ll_prev + logl(out_mat[, j], z_t_dmat, g, theta_w[j], outer = FALSE)
  }
  
  for (i in 1:D) { # separate sampling for each dimension of hidden layer

    # Draw from prior distribution
    z_prior <- rand_mvn(1, sigma = calc_K(in_dmat, theta_z[i]))

    # Initialize a and bounds on a
    a <- runif(1, min = 0, max = 2 * pi)
    amin <- a - 2 * pi
    amax <- a

    # Compute acceptance threshold - based on all dimensions of previous z
    ru <- runif(1, min = 0, max = 1)
    ll_threshold <- ll_prev + log(ru)

    # Calculate proposed values, accept or reject, repeat if necessary
    accept <- FALSE
    count <- 0
    z_prev <- z_t[, i] # store for re-proposal

    while (accept == FALSE) {
      count <- count + 1

      # Calculate proposed values and new likelihood
      z_t[, i] <- z_prev * cos(a) + z_prior * sin(a)
      dz <- sq_dist(z_t)
      new_logl <- 0
      for (j in 1:D) new_logl <- new_logl + logl(out_mat[, j], dz, g, 
                                                 theta_w[j], outer = FALSE)
      
      # Accept or reject
      if (new_logl > ll_threshold) { # accept
        ll_prev <- new_logl
        accept <- TRUE
      } else { # reject
        # update the bounds on a and repeat
        if (a < 0) {
          amin <- a
        } else {
          amax <- a
        }
        a <- runif(1, amin, amax)
        if (count > 100) stop('reached maximum iterations of ESS')
      } # end of else statement
    } # end of while loop
  } # end of i for loop

  return(list(z = z_t, ll = ll_prev, dz = dz))
}

# Fit One Layer Function ------------------------------------------------------
#' @title MCMC sampling for one layer GP
#' @description Conducts MCMC sampling of hyperparameters for a one layer 
#'     GP.  Covariance structure is based on inverse exponentiated squared 
#'     euclidean distance with length scale parameter "\code{theta}" governing 
#'     the strength of the correlation and nugget parameter "\code{g}" governing 
#'     noise.
#'
#' @details Utilizes Metropolis Hastings sampling of the length scale and
#'     nugget parameters with proposals and priors controlled by \code{settings}.
#'     Proposals for \code{g} and \code{theta} follow a uniform sliding window scheme,
#'     e.g. 
#'     
#'     \code{g_star <- runif(1, l * g_t / u, u * g_t / l)}, 
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{settings}.

#'     Priors on \code{g} and \code{theta} follow Gamma distributions with shape parameter
#'     (\code{alpha}) and rate parameter (\code{beta}) provided in \code{settings}.  
#'     These priors are designed for "\code{x}" scaled to [0,1] and "\code{y}" 
#'     scaled to have mean zero and variance 1.  
#'     
#'     The output object of class 
#'     "\code{gp}" is designed for use with \code{continue}, \code{trim}, and 
#'     \code{predict}.
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param nmcmc number of MCMC iterations
#' @param trace logical indicating whether to print iteration progress
#' @param g_0 initial value for \code{g}
#' @param theta_0 initial value for \code{theta}
#' @param true_g if true nugget is known it may be specified here (set to a small
#'        value to make fit deterministic)
#' @param settings hyperparameters for proposals and priors on \code{g} and \code{theta}
#' @return a list of the S3 class "\code{gp}" with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{settings}: copy of proposal/prior settings
#'   \item \code{g}: vector of MCMC samples for \code{g}
#'   \item \code{theta}: vector of MCMC samples for \code{theta}
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." arXiv:2012.08015. \cr\cr
#' Gramacy, RB. \emph{Surrogates: Gaussian Process Modeling, Design, and Optimization 
#'     for the Applied Sciences}. Chapman Hall, 2020.
#' 
#' @examples 
#' # Toy example (runs in less than 5 seconds) -----------------------------------
#' # This example uses a small number of MCMC iterations in order to run quickly
#' # More iterations are required to get appropriate fits
#' # Function defaults are recommended (see additional example below)
#' 
#' f <- function(x) {
#'   if (x <= 0.4) return(-1)
#'   if (x >= 0.6) return(1)
#'   if (x > 0.4 & x < 0.6) return(10*(x-0.5))
#' }
#' x <- seq(0.05, 0.95, length = 7)
#' y <- sapply(x, f)
#' x_new <- seq(0, 1, length = 100)
#' 
#' # Fit model and calculate EI
#' fit <- fit_one_layer(x, y, nmcmc = 500)
#' fit <- trim(fit, 400)
#' fit <- predict(fit, x_new, lite = TRUE, store_all = TRUE)
#' ei <- EI(fit)
#' 
#' \donttest{
#' # One Layer and EI ------------------------------------------------------------
#' 
#' f <- function(x) {
#'   sin(5 * pi * x) / (2 * x) + (x - 1) ^ 4
#' }
#'   
#' # Training data
#' x <- seq(0.5, 2, length = 30)
#' y <- f(x) + rnorm(30, 0, 0.01)
#'   
#' # Testing data
#' xx <- seq(0.5, 2, length = 100)
#' yy <- f(xx)
#'   
#' # Standardize inputs and outputs
#' xx <- (xx - min(x)) / (max(x) - min(x))
#' x <- (x - min(x)) / (max(x) - min(x))
#' yy <- (yy - mean(y)) / sd(y)
#' y <- (y - mean(y)) / sd(y)
#'   
#' # Conduct MCMC
#' fit <- fit_one_layer(x, y, nmcmc = 10000)
#' plot(fit) # investigate trace plots
#' fit <- trim(fit, 8000, 2)
#'   
#' # Predict and calculate EI
#' fit <- predict(fit, xx, lite = TRUE, store_all = TRUE)
#' ei <- EI(fit)
#'   
#' # Visualize Fit
#' plot(fit)
#' par(new = TRUE) # overlay EI
#' plot(xx, ei$value, type = 'l', lty = 2, axes = FALSE, xlab = '', ylab = '')
#' 
#' # Select next design point
#' x_new <- xx[which.max(ei$value)]
#' 
#' # Evaluate fit
#' rmse(yy, fit$mean) # lower is better
#' }
#' 
#' @export

fit_one_layer <- function(x, y, nmcmc = 10000, trace = TRUE, g_0 = 0.01, theta_0 = 0.5, true_g = NULL,
                          settings = list(l = 1, u = 2, alpha = list(g = 1.5, theta = 1.5), 
                                          beta = list(g = 3.9, theta = 3.9/1.5))) {

  tic <- proc.time()[3]

  # check that x is a matrix
  if (is.numeric(x)) x <- as.matrix(x)
  
  # check that y is a vector
  if (!is.numeric(y)) stop('y must be numeric')

  # check that x and y have matching dimension
  if (nrow(x) != length(y)) stop('dimensions of x and y do not match')
  
  # check that x is scaled properly
  if (min(x) < -0.5 | min(x) > 0.5 | max(x) < 0.5 | max(x) > 1.5) 
    warning('this function is designed for x over the range [0,1]')
  
  # check that y is scaled properly (only if nugget is not specified)
  if (is.null(true_g) & (mean(y) < -0.9 | mean(y) > 0.9 | var(y) < 0.1 | var(y) > 1.9))
    warning('this function is designed for y scaled to mean zero and variance 1')
  
  # check that all settings have been defined
  if (is.null(settings$l)) settings$l <- 1
  if (is.null(settings$u)) settings$u <- 2
  if (is.null(settings$alpha$g)) settings$alpha$g <- 1.5
  if (is.null(settings$alpha$theta)) settings$alpha$theta <- 1.5
  if (is.null(settings$beta$g)) settings$beta$g <- 3.9
  if (is.null(settings$beta$theta)) settings$beta$theta <- 3.9/1.5

  # create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings)
  class(out) <- 'gp'

  n <- length(y) # sample size
  dx <- sq_dist(x)

  # Create locations to store results, initialize mcmc values
  g <- vector(length = nmcmc)
  if (is.null(true_g)) g[1] <- g_0 else g[1] <- true_g
  theta <- vector(length = nmcmc)
  theta[1] <- theta_0
  ll <- NULL

  # Run Gibbs sampling iterations
  for (j in 2:nmcmc) {
    if(trace) if(j %% 500 == 0) cat(j, '\n')

    # sample g
    if (is.null(true_g)) {
      samp <- sample_g(y, dx, g[j-1], theta[j-1], alpha = settings$alpha$g, 
                       beta = settings$beta$g, l = settings$l, u = settings$u, ll_prev = ll)
      g[j] <- samp$g
      ll <- samp$ll
    } else g[j] <- true_g

    # sample theta
    samp <- sample_theta(y, dx, g[j], theta[j-1], alpha = settings$alpha$theta,
                         beta = settings$beta$theta, l = settings$l, u = settings$u,
                         outer = TRUE, ll_prev = ll)
    theta[j] <- samp$theta
    ll <- samp$ll
  } # end of j for loop

  out$g <- g
  out$theta <- theta
  toc <- proc.time()[3]
  out$time <- toc - tic

  return(out)
}

# Fit Two Layer Function ------------------------------------------------------
#' @title MCMC sampling for two layer deep GP
#' @description Conducts MCMC sampling of hyperparameters and hidden layer 
#'     "\code{w}" for a two layer deep GP.  Covariance structure is based on 
#'     inverse exponentiated squared euclidean distance.  Separate length scale 
#'     parameters "\code{theta_w}" and "\code{theta_y}" govern the correlation 
#'     strength of the hidden layer and outer layer respectively.  Nugget 
#'     parameter "\code{g}" governs noise on the outer layer.
#'
#' @details Maps inputs "\code{x}" through hidden layer "\code{w}" to outputs 
#'     "\code{y}".  Conducts sampling of the hidden layer using Elliptical Slice 
#'     sampling.  Utilizes Metropolis Hastings sampling of the length scale and
#'     nugget parameters with proposals and priors controlled by \code{settings}.
#'     Proposals for \code{g}, \code{theta_y}, and \code{theta_w} follow a uniform
#'     sliding window scheme, e.g.
#'     
#'     \code{g_star <- runif(1, l * g_t / u, u * g_t / l)}, 
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{settings}.   Priors 
#'     on \code{g} and \code{theta} follow Gamma distributions with shape parameter
#'     (\code{alpha}) and rate parameter (\code{beta}) provided in \code{settings}.  
#'     These priors are designed for "\code{x}" scaled to [0,1] and "\code{y}" 
#'     scaled to have mean zero and variance 1.  
#'     
#'     The output object of class 
#'     "\code{dgp2}" is designed for use with \code{continue}, \code{trim}, 
#'     and \code{predict}. If \code{w_0} is of dimension \code{nrow(x) - 1} by 
#'     \code{D}, the final row is predicted using kriging.  This is helpful in 
#'     sequential design when adding a new input location and starting the MCMC 
#'     at the place where the previous MCMC left off.
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param D integer designating dimension of hidden layer, defaults to dimension of \code{x}
#' @param nmcmc number of MCMC iterations
#' @param trace logical indicating whether to print iteration progress
#' @param w_0 initial value for hidden layer \code{w}, defaults to identity mapping (must be 
#'        matrix of dimension \code{nrow(x)} by \code{D} or dimension \code{nrow(x) - 1} by \code{D})
#' @param g_0 initial value for \code{g}
#' @param theta_y_0 initial value for \code{theta_y} (length scale of outer layer)
#' @param theta_w_0 initial value for \code{theta_w} (length scale of inner layer), 
#'        may be single value or vector of length \code{D}
#' @param true_g if true nugget is known it may be specified here (set to a small
#'        value to make fit deterministic)
#' @param settings hyperparameters for proposals and priors on \code{g}, \code{theta_y}, and \code{theta_w}
#' @return a list of the S3 class "\code{dgp2}" with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{settings}: copy of proposal/prior settings
#'   \item \code{g}: vector of MCMC samples for \code{g}
#'   \item \code{theta_y}: vector of MCMC samples for \code{theta_y} (length
#'         scale of outer layer)
#'   \item \code{theta_w}: matrix of MCMC samples for \code{theta_w} (length 
#'         scale of inner layer)
#'   \item \code{w}: list of MCMC samples for hidden layer \code{w}
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." arXiv:2012.08015. \cr\cr
#' Murray, I, RP Adams, and D MacKay. 2010. "Elliptical slice sampling." 
#'     \emph{Journal of Machine Learning Research 9}, 541-548.
#' 
#' @examples 
#' # Toy example (runs in less than 5 seconds) -----------------------------------
#' # This example uses a small number of MCMC iterations in order to run quickly
#' # More iterations are required to get appropriate fits
#' # Function defaults are recommended (see additional example below)
#' 
#' f <- function(x) {
#'   if (x <= 0.4) return(-1)
#'   if (x >= 0.6) return(1)
#'   if (x > 0.4 & x < 0.6) return(10*(x-0.5))
#' }
#' x <- seq(0.05, 0.95, length = 7)
#' y <- sapply(x, f)
#' x_new <- seq(0, 1, length = 100)
#' 
#' # Fit model and calculate ALC
#' fit <- fit_two_layer(x, y, nmcmc = 500)
#' fit <- trim(fit, 400)
#' fit <- predict(fit, x_new)
#' alc <- ALC(fit)
#' 
#' \donttest{
#' # Two Layer and ALC -----------------------------------------------------------
#' 
#' f <- function(x) {
#'   exp(-10 * x) * (cos(10 * pi * x - 1) + sin(10 * pi * x - 1)) * 5 - 0.2
#' }
#' 
#' # Training data
#' x <- seq(0, 1, length = 30)
#' y <- f(x) + rnorm(30, 0, 0.05)
#' 
#' # Testing data
#' xx <- seq(0, 1, length = 100)
#' yy <- f(xx)
#' 
#' # Conduct MCMC
#' fit <- fit_two_layer(x, y, D = 1, nmcmc = 9000)
#' fit <- continue(fit, 1000)
#' plot(fit) # investigate trace plots
#' fit <- trim(fit, 8000, 2)
#' 
#' # Option 1 - calculate ALC from MCMC iterations
#' alc <- ALC(fit, xx)
#' 
#' # Option 2 - calculate ALC after predictions
#' fit <- predict(fit, xx)
#' alc <- ALC(fit)
#' 
#' # Visualize fit
#' plot(fit)
#' par(new = TRUE) # overlay ALC
#' plot(xx, alc$value, type = 'l', lty = 2, axes = FALSE, xlab = '', ylab = '')
#' 
#' # Select next design point
#' x_new <- xx[which.max(alc$value)]
#' 
#' # Evaluate fit
#' rmse(yy, fit$mean) # lower is better
#' }
#' 
#' @export

fit_two_layer <- function(x, y, D = ifelse(is.matrix(x), ncol(x), 1), nmcmc = 10000, trace = TRUE,
                         w_0 = suppressWarnings(matrix(x, nrow = length(y), ncol = D)), 
                         g_0 = 0.01, theta_y_0 = 0.5, theta_w_0 = 1, true_g = NULL,
                         settings = list(l = 1, u = 2, alpha = list(g = 1.5, theta_w = 1.5, theta_y = 1.5), 
                                         beta = list(g = 3.9, theta_w = 3.9/4, theta_y = 3.9/6))) {
  
  tic <- proc.time()[3]

  # check that x is a matrix
  if (is.numeric(x)) x <- as.matrix(x)

  # check that y is a vector
  if (!is.numeric(y)) stop('y must be numeric')

  # check that x and y have matching dimension
  if (nrow(x) != length(y)) stop('dimensions of x and y do not match')
  
  # check that x is scaled properly
  if (min(x) < -0.5 | min(x) > 0.5 | max(x) < 0.5 | max(x) > 1.5) 
    warning('this function is designed for x over the range [0,1]')
  
  # check that y is scaled properly (only if nugget is not specified)
  if (is.null(true_g) & (mean(y) < -0.9 | mean(y) > 0.9 | var(y) < 0.1 | var(y) > 1.9))
    warning('this function is designed for y scaled to mean zero and variance 1')
  
  # check that all settings have been defined
  if (is.null(settings$l)) settings$l <- 1
  if (is.null(settings$u)) settings$u <- 2
  if (is.null(settings$alpha$g)) settings$alpha$g <- 1.5
  if (is.null(settings$alpha$theta_w)) settings$alpha$theta_w <- 1.5
  if (is.null(settings$alpha$theta_y)) settings$alpha$theta_y <- 1.5
  if (is.null(settings$beta$g)) settings$beta$g <- 3.9
  if (is.null(settings$beta$theta_w)) settings$beta$theta_w <- 3.9/4
  if (is.null(settings$beta$theta_y)) settings$beta$theta_y <- 3.9/6

  # check that w_0 is a matrix
  if (!is.matrix(w_0)) {
    warning('w_0 must be a matrix, trying to coerce')
    w_0 <- matrix(w_0)
  }

  # check that the dimension of w_0 matches the specified dimension
  if (ncol(w_0) != D) stop('dimension of w_0 does not match D')

  # if theta_0_w is a single value, expand to dimension of W
  if (length(theta_w_0) != D & length(theta_w_0) == 1) theta_w_0 <- rep(theta_w_0, D)

  # create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings)
  class(out) <- 'dgp2'

  n <- length(y) # sample size
  dx <- sq_dist(x)

  # if w_0 is from previous sequential design iteration, predict at new point
  if (nrow(w_0) == n - 1) {
    # predict w for final row of x (separately for each dimension)
    new_w <- vector(length = D)
    old_x <- x[1:(n-1), ]
    new_x <- matrix(x[n, ], nrow = 1)
    for (i in 1:D) {
      new_w[i] <- krig(w_0[, i], dx[1:(n-1), 1:(n-1)], d_cross = sq_dist(new_x, old_x),
                       theta = theta_w_0[i], g = NULL, mean = TRUE, s2 = FALSE, sigma = FALSE,
                       tau2 = FALSE)$mean
    }
    w_0 <- rbind(w_0, new_w)
  }

  # Create locations to store results, initialize mcmc values
  g <- vector(length = nmcmc)
  if (is.null(true_g)) g[1] <- g_0 else g[1] <- true_g
  theta_y <- vector(length = nmcmc)
  theta_y[1] <- theta_y_0
  theta_w <- matrix(nrow = nmcmc, ncol = D)
  theta_w[1,] <- theta_w_0
  w_out <- list()
  w_out[[1]] <- w_0
  dw <- sq_dist(w_out[[1]])
  ll_outer <- NULL

  # Run Gibbs sampling iterations
  for (j in 2:nmcmc) {
    
    if(trace) if(j %% 500 == 0) cat(j,'\n')
    
    dw <- sq_dist(w_out[[j-1]])

    # sample g
    if (is.null(true_g)) {
      samp <- sample_g(y, dw, g[j - 1], theta_y[j - 1], alpha = settings$alpha$g, 
                       beta = settings$beta$g, l = settings$l, u = settings$u, ll_prev = ll_outer)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else g[j] <- true_g

    # sample theta_y
    samp <- sample_theta(y, dw, g[j], theta_y[j - 1], alpha = settings$alpha$theta_y,
                         beta = settings$beta$theta_y, l = settings$l, u = settings$u,
                                outer = TRUE, ll_prev = ll_outer)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    
    # sample theta_w
    for (i in 1:D) {
      samp <- sample_theta(w_out[[j - 1]][, i], dx, NULL, theta_w[j-1, i],
                           alpha = settings$alpha$theta_w, beta = settings$beta$theta_w,
                           l = settings$l, u = settings$u, outer = FALSE)
      theta_w[j, i] <- samp$theta
    }

    # sample w
    samp <- sample_w(y, w_out[[j-1]], dw, dx, g[j], theta_y[j], theta_w[j, ], 
                           ll_prev = ll_outer)
    w_out[[j]] <- samp$w
    ll_outer <- samp$ll
    dw <- samp$dw
  } # end of j for loop

  out$g <- g
  out$theta_y <- theta_y
  out$theta_w <- theta_w
  out$w <- w_out
  toc <- proc.time()[3]
  out$time <- toc - tic

  return(out)
}

# Fit Three Layer Function ----------------------------------------------------
#' @title MCMC sampling for three layer deep GP
#' @description Conducts MCMC sampling of hyperparameters, hidden layer 
#'     "\code{z}", and hidden layer "\code{w}" for a three layer deep GP.  
#'     Covariance structure is based on inverse exponentiated squared euclidean 
#'     distance.  Separate length scale parameters "\code{theta_z}", 
#'     "\code{theta_w}", and "\code{theta_y}" govern the correlation 
#'     strength of the inner layer, middle layer, and outer layer respectively.  
#'     Nugget parameter "\code{g}" governs noise on the outer layer.
#'
#' @details Maps inputs "\code{x}" through hidden layer "\code{z}" then hidden
#'     layer "\code{w}" to outputs "\code{y}".  Conducts sampling of the hidden 
#'     layers using Elliptical Slice sampling.  Utilizes Metropolis Hastings 
#'     sampling of the length scale and nugget parameters with proposals and priors
#'     controlled by \code{settings}.  Proposals for \code{g}, \code{theta_y}, 
#'     \code{theta_w}, and \code{theta_z} follow a uniform sliding window scheme, e.g.
#'     
#'     \code{g_star <- runif(1, l * g_t / u, u * g_t / l)},
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{settings}.  Priors
#'     on \code{g}, \code{theta_y}, \code{theta_w}, and \code{theta_z} follow Gamma
#'     distributions with shape parameter (\code{alpha}) and rate parameter (\code{beta}) 
#'     provided in \code{settings}.  These priors are designed for "\code{x}" scaled to 
#'     [0,1] and "\code{y}" scaled to have mean zero and variance 1.  
#'     
#'     The output object of class "\code{dgp3}" is designed for use with \code{continue}, \code{trim}, 
#'     and \code{predict}. If \code{z_0} and \code{w_0} are of dimension 
#'     \code{nrow(x) - 1} by \code{D}, the final rows are predicted using 
#'     kriging.  This is helpful in sequential design when adding a new input 
#'     location and starting the MCMC at the place where the previous MCMC 
#'     left off.
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param D integer designating dimension of hidden layers, defaults to dimension of \code{x}
#' @param nmcmc number of MCMC iterations
#' @param trace logical indicating whether to print iteration progress
#' @param w_0 initial value for hidden layer \code{w}, defaults to identity mapping (must be 
#'        matrix of dimension \code{nrow(x)} by \code{D} or dimension \code{nrow(x) - 1} by \code{D})
#' @param z_0 initial value for hidden layer \code{z}, defaults to identity mapping (must be 
#'        matrix of dimension \code{nrow(x)} by \code{D} or dimension \code{nrow(x) - 1} by \code{D})
#' @param g_0 initial value for \code{g}
#' @param theta_y_0 initial value for \code{theta_y} (length scale of outer layer)
#' @param theta_w_0 initial value for \code{theta_w} (length scale of middle layer),
#'        may be single value or vector of length \code{D}
#' @param theta_z_0 initial value for \code{theta_z} (length scale of inner layer),
#'        may be single value or vector of length \code{D}
#' @param true_g if true nugget is known it may be specified here (set to a small value
#'        to make fit deterministic)
#' @param settings hyperparameters for proposals and priors on \code{g}, \code{theta_y},
#'        \code{theta_w}, and \code{theta_z}
#' @return a list of the S3 class "\code{dgp3}" with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{settings}: copy of proposal/prior settings
#'   \item \code{g}: vector of MCMC samples for \code{g}
#'   \item \code{theta_y}: vector of MCMC samples for \code{theta_y} (length scale of outer layer)
#'   \item \code{theta_w}: matrix of MCMC samples for \code{theta_w} (length scale of middle layer)
#'   \item \code{theta_z}: matrix of MCMC samples for \code{theta_z} (length scale of inner layer)
#'   \item \code{w}: list of MCMC samples for middle hidden layer \code{w}
#'   \item \code{z}: list of MCMC samples for inner hidden layer \code{z}
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." arXiv:2012.08015. \cr\cr
#' Murray, I, RP Adams, and D MacKay. 2010. "Elliptical slice sampling."
#'      \emph{Journal of Machine Learning Research 9}, 541-548.
#' 
#' @examples 
#' # Toy example (runs in less than 5 seconds) -----------------------------------
#' # This example uses a small number of MCMC iterations in order to run quickly
#' # More iterations are required to get appropriate fits
#' # Function defaults are recommended (see additional example below)
#' 
#' f <- function(x) {
#'   if (x <= 0.4) return(-1)
#'   if (x >= 0.6) return(1)
#'   if (x > 0.4 & x < 0.6) return(10*(x-0.5))
#' }
#' x <- seq(0.05, 0.95, length = 7)
#' y <- sapply(x, f)
#' x_new <- seq(0, 1, length = 100)
#' 
#' # Fit model and calculate IMSPE
#' fit <- fit_three_layer(x, y, nmcmc = 500)
#' fit <- trim(fit, 400)
#' fit <- predict(fit, x_new)
#' imse <- IMSE(fit)
#' 
#' \donttest{
#' # Three Layer and IMSE -------------------------------------------------------
#' 
#' f <- function(x) {
#'   i <- which(x <= 0.48)
#'   x[i] <- 2 * sin(pi * x[i] * 4) + 0.4 * cos(pi * x[i] * 16)
#'   x[-i] <- 2 * x[-i] - 1
#'   return(x)
#' }
#' 
#' # Training data
#' x <- seq(0, 1, length = 30)
#' y <- f(x) + rnorm(30, 0, 0.05)
#' 
#' # Testing data
#' xx <- seq(0, 1, length = 100)
#' yy <- f(xx)
#' 
#' # Conduct MCMC
#' fit <- fit_three_layer(x, y, D = 1, nmcmc = 10000)
#' plot(fit) # investigate trace plots
#' fit <- trim(fit, 8000, 2)
#' 
#' # Option 1 - calculate IMSE from only MCMC iterations
#' imse <- IMSE(fit, xx)
#' 
#' # Option 2 - calculate IMSE after predictions
#' fit <- predict(fit, xx)
#' imse <- IMSE(fit)
#' 
#' # Visualize fit
#' plot(fit)
#' par(new = TRUE) # overlay IMSPE
#' plot(xx, imse$value, type = 'l', lty = 2, axes = FALSE, xlab = '', ylab = '')
#' 
#' # Select next design point
#' x_new <- xx[which.min(imse$value)]
#' 
#' # Evaluate fit
#' rmse(yy, fit$mean) # lower is better
#' }
#' 
#' @export

fit_three_layer <- function(x, y, D = ifelse(is.matrix(x), ncol(x), 1), nmcmc = 10000, trace = TRUE,
                           w_0 = suppressWarnings(matrix(x, nrow = length(y), ncol = D)),
                           z_0 = suppressWarnings(matrix(x, nrow = length(y), ncol = D)), 
                           g_0 = 0.01, theta_y_0 = 0.5, theta_w_0 = 1, theta_z_0 = 1, true_g = NULL,
                           settings = list(l = 1, u = 2, alpha = list(g = 1.5, theta_z = 1.5, theta_w = 1.5, theta_y = 1.5), 
                                           beta = list(g = 3.9, theta_z = 3.9/4, theta_w = 3.9/12, theta_y = 3.9/6))) {

  tic <- proc.time()[3]

  # check that x is a matrix
  if (is.numeric(x)) x <- as.matrix(x)

  # check that y is a vector
  if (!is.numeric(y)) stop('y must be numeric')

  # check that x and y have matching dimension
  if (nrow(x) != length(y)) stop('dimensions of x and y do not match')
  
  # check that x is scaled properly
  if (min(x) < -0.5 | min(x) > 0.5 | max(x) < 0.5 | max(x) > 1.5) 
    warning('this function is designed for x over the range [0,1]')
  
  # check that y is scaled properly (only if nugget is not specified)
  if (is.null(true_g) & (mean(y) < -0.9 | mean(y) > 0.9 | var(y) < 0.1 | var(y) > 1.9))
    warning('this function is designed for y scaled to mean zero and variance 1')
  
  # check that all settings have been defined
  if (is.null(settings$l)) settings$l <- 1
  if (is.null(settings$u)) settings$u <- 2
  if (is.null(settings$alpha$g)) settings$alpha$g <- 1.5
  if (is.null(settings$alpha$theta_z)) settings$alpha$theta_z <- 1.5
  if (is.null(settings$alpha$theta_w)) settings$alpha$theta_w <- 1.5
  if (is.null(settings$alpha$theta_y)) settings$alpha$theta_y <- 1.5
  if (is.null(settings$beta$g)) settings$beta$g <- 3.9
  if (is.null(settings$beta$theta_z)) settings$beta$theta_z <- 3.9/4
  if (is.null(settings$beta$theta_w)) settings$beta$theta_w <- 3.9/12
  if (is.null(settings$beta$theta_y)) settings$beta$theta_y <- 3.9/6

  # check that w_0 is a matrix
  if(!is.matrix(w_0)) {
    warning('w_0 must be a matrix, trying to coerce')
    w_0 <- matrix(w_0)
  }

  # check that z_0 is a matrix
  if(!is.matrix(z_0)) {
    warning('z_0 must be a matrix, trying to coerce')
    z_0 <- matrix(z_0)
  }

  # check that the dimension of w_0 matches the specified dimension
  if (ncol(w_0) != D) stop('dimension of w_0 does not match D')

  # check that the dimension of z_0 matches the specified dimension
  if (ncol(z_0) != D) stop('dimension of z_0 does not match D')

  # if theta_0_w is a single value, expand to dimension of W (same for Z)
  if (length(theta_w_0) != D & length(theta_w_0) == 1) theta_w_0 <- rep(theta_w_0, D)
  if (length(theta_z_0) != D & length(theta_z_0) == 1) theta_z_0 <- rep(theta_z_0, D)

  # create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings)
  class(out) <- 'dgp3'

  n <- length(y) # sample size
  dx <- sq_dist(x)

  # if z_0 is from previous sequential design iteration, predict at new point
  if (nrow(z_0) == n - 1) {
    # predict z for final row of x (separately for each dimension)
    new_z <- vector(length = D)
    old_x <- x[1:(n-1), ]
    new_x = matrix(x[n,], nrow = 1)
    for (i in 1:D) {
      new_z[i] <- krig(z_0[, i], dx[1:(n-1), 1:(n-1)], d_cross = sq_dist(new_x, old_x),
                       theta = theta_z_0[i],g = NULL, mean = TRUE, s2 = FALSE, sigma = FALSE,
                       tau2 = FALSE)$mean
    }
    z_0 <- rbind(z_0, new_z)
  }

  # if w_0 is from previous sequential design iteration, predict at new point
  if (nrow(w_0) == n - 1) {
    # predict w for final row of z (separately for each dimension)
    new_w <- vector(length = D)
    old_z <- z_0[1:(n-1), ]
    new_z = matrix(z_0[n, ], nrow = 1)
    for (i in 1:D) {
      new_w[i] <- krig(w_0[, i], sq_dist(old_z), d_cross = sq_dist(new_z, old_z),
                       theta = theta_w_0[i], g = NULL, mean = TRUE, s2 = FALSE, sigma = FALSE,
                       tau2 = FALSE)$mean
    }
    w_0 <- rbind(w_0, new_w)
  }

  # Create locations to store results
  g <- vector(length = nmcmc)
  if (is.null(true_g)) g[1] <- g_0 else g[1] <- true_g
  theta_y <- vector(length = nmcmc)
  theta_y[1] <- theta_y_0
  theta_w <- matrix(nrow = nmcmc, ncol = D)
  theta_w[1, ] <- theta_w_0 
  theta_z <- matrix(nrow = nmcmc, ncol = D)
  theta_z[1, ] <- theta_z_0
  z_out <- list()
  z_out[[1]] <- z_0
  dz <- sq_dist(z_out[[1]])
  w_out <- list()
  w_out[[1]] <- w_0
  dw <- sq_dist(w_out[[1]])
  ll_outer <- NULL
  
  # Run Gibbs sampling iterations
  for (j in 2:nmcmc) {
    if(trace) if(j %% 500 == 0) cat(j,'\n')

    # sample g
    if (is.null(true_g)) {
      samp <- sample_g(y, dw, g[j-1], theta_y[j-1], alpha = settings$alpha$g, beta = settings$beta$g,
                       l = settings$l, u = settings$u, ll_prev = ll_outer)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else g[j] <- true_g

    # sample theta_y
    samp <- sample_theta(y, dw, g[j], theta_y[j-1], alpha = settings$alpha$theta_y,
                         beta = settings$beta$theta_y, l = settings$l, u = settings$u, outer = TRUE,
                         ll_prev = ll_outer)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    
    # sample theta_w
    ll_mid <- 0 # re-calculated each time since we have a new z
    for (i in 1:D) {
      samp <- sample_theta(w_out[[j-1]][, i], dz, g = NULL, theta_w[j-1, i],
                            alpha = settings$alpha$theta_w, beta = settings$beta$theta_w,
                           l = settings$l, u = settings$u, outer = FALSE)
      theta_w[j, i] <- samp$theta
      ll_mid <- ll_mid + samp$ll
    }

    # sample theta_z
    for (i in 1:D) {
      samp <- sample_theta(z_out[[j-1]][, i], dx, g = NULL, theta_z[j-1, i], 
                           alpha = settings$alpha$theta_z, beta = settings$beta$theta_z,
                           l = settings$l, u = settings$u, outer = FALSE)
      theta_z[j, i] <- samp$theta
    }

    # sample z
    samp <- sample_z(w_out[[j-1]], z_out[[j-1]], dz, dx, g = NULL, theta_w[j, ], theta_z[j, ],
                     ll_prev = ll_mid)
    z_out[[j]] <- samp$z
    dz <- samp$dz

    # sample w
    samp <- sample_w(y, w_out[[j-1]], dw, dz, g = g[j], theta_y[j], theta_w[j, ], ll_prev = ll_outer)
    w_out[[j]] <- samp$w
    ll_outer <- samp$ll
    dw <- samp$dw
  } # end of j for loop

  out$g <- g
  out$theta_y <- theta_y
  out$theta_w <- theta_w
  out$theta_z <- theta_z
  out$w <- w_out
  out$z <- z_out
  toc <- proc.time()[3]
  out$time <- toc - tic
  
  return(out)
}

# Define Continue for S3 Objects ----------------------------------------------
#' @title Continues MCMC sampling
#' @description Acts on a "\code{gp}", "\code{dgp2}", or "\code{dgp3}" object.  
#'     Continues MCMC sampling of hyperparameters and hidden layers and appends
#'     results to existing object.
#'     
#' @details See "\code{fit_one_layer}", "\code{fit_two_layer}", or 
#'     "\code{fit_three_layer}" for details on MCMC.  The resulting object 
#'     will have \code{nmcmc} equal to the previous \code{nmcmc} plus 
#'     \code{new_mcmc}.  It is recommended to start an MCMC fit then 
#'     investigate trace plots to assess burn-in.  The primary use of this 
#'     function is to gather more MCMC iterations in order to obtain burned-in 
#'     samples.
#' 
#' @param object object from \code{fit_one_layer}, \code{fit_two_layer}, or 
#'        \code{fit_three_layer}
#' @param new_mcmc number of MCMC iterations to conduct and append
#' @param trace logical indicating whether to print iteration progress
#' @return object of the same class with the new iterations appended
#' 
#' @examples 
#' # See "deepgp-package" or "fit_two_layer" for an example
#' 
#' @rdname continue
#' @export

continue <- function(object, new_mcmc, trace)
  UseMethod("continue", object)

# Continue One Layer Function -------------------------------------------------
#' @rdname continue
#' @export

continue.gp <- function(object, new_mcmc = 1000, trace = TRUE) {
  
  tic <- proc.time()[3]
  
  # Use true nugget if it was specified
  if (length(unique(object$g)) == 1) {
    true_g = object$g[1] 
  } else true_g = NULL

  # Run continuing MCMC iterations
  new_fit <- fit_one_layer(object$x, object$y, nmcmc = new_mcmc, trace = trace,
                    g_0 = object$g[object$nmcmc], theta_0 = object$theta[object$nmcmc],
                    true_g = true_g, settings = object$settings)

  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$g <- c(object$g, new_fit$g)
  object$theta <- c(object$theta, new_fit$theta)

  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  return(object)
}

# Continue Two Layer Function -------------------------------------------------
#' @rdname continue
#' @export

continue.dgp2 <- function(object, new_mcmc = 1000, trace = TRUE) {

  tic <- proc.time()[3]
  
  # Use true nugget if it was specified
  if (length(unique(object$g)) == 1) {
    true_g = object$g[1] 
  } else true_g = NULL

  # Run continuing MCMC iterations
  new_fit <- fit_two_layer(object$x, object$y, D = ncol(object$w[[1]]),
                           nmcmc = new_mcmc, trace = trace,
                           w_0 = object$w[[object$nmcmc]], g_0 = object$g[object$nmcmc],
                           theta_y_0 = object$theta_y[object$nmcmc],
                           theta_w_0 = object$theta_w[object$nmcmc, ], true_g = true_g,
                           settings = object$settings)

  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$g <- c(object$g, new_fit$g)
  object$theta_y <- c(object$theta_y, new_fit$theta_y)
  object$theta_w <- rbind(object$theta_w, new_fit$theta_w)
  object$w <- c(object$w, new_fit$w)

  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  return(object)
}

# Continue Three Layer Function -----------------------------------------------
#' @rdname continue
#' @export

continue.dgp3 <- function(object, new_mcmc = 1000, trace = TRUE) {

  tic <- proc.time()[3]
  
  # Use true nugget if it was specified
  if (length(unique(object$g)) == 1) {
    true_g = object$g[1] 
  } else true_g = NULL

  # Run continuing MCMC iterations
  new_fit <- fit_three_layer(object$x, object$y, D = ncol(object$w[[1]]),
                             nmcmc = new_mcmc, trace = trace,
                             w_0 = object$w[[object$nmcmc]], z_0 = object$z[[object$nmcmc]],
                             g_0 = object$g[object$nmcmc], theta_y_0 = object$theta_y[object$nmcmc],
                             theta_w_0 = object$theta_w[object$nmcmc, ],
                             theta_z_0 = object$theta_z[object$nmcmc, ], true_g = true_g,
                             settings = object$settings)

  # Append new information to original fit
  object$nmcmc <- object$nmcmc + new_mcmc
  object$g <- c(object$g, new_fit$g)
  object$theta_y <- c(object$theta_y, new_fit$theta_y)
  object$theta_w <- rbind(object$theta_w, new_fit$theta_w)
  object$theta_z <- rbind(object$theta_z, new_fit$theta_z)
  object$w <- c(object$w, new_fit$w)
  object$z <- c(object$z, new_fit$z)

  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)
  return(object)
}

# Define Trim for S3 Objects --------------------------------------------------
#' @title Trim/Thin MCMC iterations
#' @description Acts on a "\code{gp}", "\code{dgp2}", or "\code{dgp3}" object.
#'    Removes the specified number of MCMC iterations (starting at the first 
#'    iteration).  After these samples are removed, the remaining samples may 
#'    be thinned.
#' 
#' @details The resulting object will have \code{nmcmc} equal to the previous 
#'     \code{nmcmc} minus \code{burn} divided by \code{thin}.  It is recommended 
#'     to start an MCMC fit then investigate trace plots to assess burn-in.  Once
#'     burn-in has been achieved, use this function to remove the starting 
#'     iterations.  Thinning reduces the size of the resulting object while 
#'     accounting for the high correlation between consecutive iterations.
#'     
#' @param object object from \code{fit_one_layer}, \code{fit_two_layer}, or 
#'        \code{fit_three_layer}
#' @param burn integer specifying number of iterations to cut off as burn-in
#' @param thin integer specifying amount of thinning (\code{thin = 1} keeps all 
#'        iterations, \code{thin = 2} keeps every other iteration, 
#'        \code{thin = 10} keeps every tenth iteration, etc.)
#' @return object of the same class with the selected iterations removed
#' 
#' @examples 
#' # See "deepgp-package", "fit_one_layer", "fit_two_layer", or "fit_three_layer"
#' # for an example
#' 
#' @rdname trim
#' @export

trim <- function(object, burn, thin)
  UseMethod("trim", object)


# Trim One Layer Function -----------------------------------------------------
#' @rdname trim
#' @export

trim.gp <- function(object, burn, thin = 1) {

  tic <- proc.time()[3]

  if (burn >= object$nmcmc) stop('burn must be less than nmcmc')

  nmcmc <- object$nmcmc
  indx <- (burn + 1):nmcmc
  indx <- indx[which(indx %% thin == 0)]

  object$nmcmc <- length(indx)
  object$g <- object$g[indx]
  object$theta <- object$theta[indx]

  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)

  return(object)
}

# Trim Two Layer Function -----------------------------------------------------
#' @rdname trim
#' @export

trim.dgp2 <- function(object, burn, thin = 1) {

  tic <- proc.time()[3]

  if (burn >= object$nmcmc) stop('burn must be less than nmcmc')

  nmcmc <- object$nmcmc
  indx <- (burn + 1):nmcmc
  indx <- indx[which(indx %% thin == 0)]

  object$nmcmc <- length(indx)
  object$g <- object$g[indx]
  object$theta_y <- object$theta_y[indx]
  object$theta_w <- as.matrix(object$theta_w[indx,])
  object$w <- object$w[indx]

  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)

  return(object)
}

# Trim Three Layer Function ---------------------------------------------------
#' @rdname trim
#' @export

trim.dgp3 <- function(object, burn, thin = 1) {

  tic <- proc.time()[3]

  if (burn >= object$nmcmc) stop('burn must be less than nmcmc')

  nmcmc <- object$nmcmc
  indx <- (burn + 1):nmcmc
  indx <- indx[which(indx %% thin == 0)]

  object$nmcmc <- length(indx)
  object$g <- object$g[indx]
  object$theta_y <- object$theta_y[indx]
  object$theta_w <- as.matrix(object$theta_w[indx,])
  object$theta_z <- as.matrix(object$theta_z[indx,])
  object$w <- object$w[indx]
  object$z <- object$z[indx]

  toc <- proc.time()[3]
  object$time <- object$time + (toc - tic)

  return(object)
}

# Define Predict for S3 Objects --------------------------------------------------
#' @name predict
#' @title Predict posterior mean and variance/covariance
#' @description Acts on a "\code{gp}", "\code{dgp2}", or "\code{dgp3}" object.
#'     Calculates posterior mean and variance/covariance over specified input 
#'     locations.  Optionally utilizes SNOW parallelization.
#' 
#' @details All iterations in the object are used for prediction, so samples 
#'     should be burned-in.  Thinning the samples using \code{trim} will speed up 
#'     computation.  Posterior moments are calculated using conditional expectation
#'     and variance.  As a default, only point-wise variance is calculated.  Full covariance
#'     may be calculated using \code{lite = FALSE}.  The storage of means and point-wise
#'     variances for each individual iteration (specified using \code{store_all = TRUE})
#'     is required in order to use \code{EI}.
#'     
#'     SNOW parallelization reduces computation time but requires significantly more 
#'     memory storage.  Use \code{cores = 1} if memory is limited.
#' 
#' @param object object from \code{fit_one_layer}, \code{fit_two_layer}, or 
#'        \code{fit_three_layer} with burn-in already removed
#' @param x_new matrix of predictive input locations
#' @param lite logical indicating whether to calculate only point-wise variances (\code{lite = TRUE}),
#'        or full covariance (\code{lite = FALSE})
#' @param store_all logical indicating whether to store mean and variance for each iteration
#' @param mean_map denotes whether to map hidden layers using conditional mean or a random sample
#'        from the full MVN distribution ("\code{dgp2}" or "\code{dgp3}" only) 
#' @param cores number of cores to utilize in parallel, by default no parallelization
#'        is used
#' @param ... N/A
#' @return object of the same class with the following additional elements:
#' \itemize{
#'   \item \code{x_new}: copy of predictive input locations
#'   \item \code{tau2}: vector of tau2 estimates (governing the magnitude of 
#'         the covariance)
#'   \item \code{mean}: predicted posterior mean, indices correspond to 
#'         \code{x_new} location
#'   \item \code{s2}: predicted point-wise variances, indices correspond to \code{x_new} location
#'         (only returned when \code{lite = TRUE})
#'   \item \code{s2_smooth}: predicted point-wise variances with \code{g} removed, indices 
#'         correspond to \code{x_new} location (only returned when \code{lite = TRUE})
#'   \item \code{Sigma}: predicted posterior covariance, indices correspond to 
#'         \code{x_new} location (only returned when \code{lite = FALSE})
#'   \item \code{Sigma_smooth}: predicted posterior covariance with \code{g} removed 
#'         from the diagonal (only returned when \code{lite = FALSE})
#'   \item \code{mu_t}: matrix of posterior mean for each iteration, column index corresponds 
#'         to iteration and row index corresponds to \code{x_new} location (only returned when 
#'         \code{store_all = TRUE})
#'   \item \code{s2_t}: matrix of posterior point-wise variance for each iteration, column index
#'         corresponds to iteration and row index corresponds to \code{x_new} location (only
#'         returned when \code{store_all = TRUE})
#'   \item \code{w_new}: list of hidden layer mappings, list index corresponds to iteration and 
#'         row index corresponds to \code{x_new} location ("\code{dgp2}" and "\code{dgp3}" only)
#'   \item \code{z_new}: list of hidden layer mappings, list index corresponds to iteration and 
#'         row index corresponds to \code{x_new} location ("\code{dgp3}" only) 
#' }
#' Computation time is added to the computation time of the existing object.
#' 
#' @references 
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." arXiv:2012.08015. \cr\cr
#' Gramacy, RB. \emph{Surrogates: Gaussian Process Modeling, Design, and Optimization 
#'     for the Applied Sciences}. Chapman Hall, 2020.
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

predict.gp <- function(object, x_new, lite = TRUE, store_all = FALSE, cores = 1, ...) {
  
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
      k <- krig(object$y, dx, d_new, d_cross, object$theta[t], object$g[t], s2 = lite, sigma = !lite)
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
      k <- krig(object$y, dx, d_new, d_cross, object$theta[t], object$g[t], s2 = lite, sigma = !lite)
      return(k)
    }
  
    stopCluster(cl)
  
    # group elements out of the list
    mu_t <- sapply(result, with, eval(parse(text = "mean")))
    tau2 <- sapply(result, with, eval(parse(text = "tau2")))
    if (lite) {
        s2_t <- sapply(result, with, eval(parse(text = "s2")))
    } else sigma_t_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma"))))
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
          k <- krig(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], g = NULL, s2 = FALSE, sigma = FALSE,
                    tau2 = FALSE)
          w_new[[t]][, i] <- k$mean
        } else {
          k <- krig(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], g = NULL, tau2 = FALSE)
          w_new[[t]][, i] <- rand_mvn(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      
      # map w_new to mu_t and sigma_t
      k <- krig(object$y, sq_dist(w_t), sq_dist(w_new[[t]]), sq_dist(w_new[[t]], w_t),
                object$theta_y[t], object$g[t], s2 = lite, sigma = !lite)
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
          k <- krig(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], g = NULL, sigma = FALSE,
                    tau2 = FALSE)
          w_new[, i] <- k$mean
        } else {
          k <- krig(w_t[, i], dx, d_new, d_cross, object$theta_w[t, i], g = NULL, tau2 = FALSE)
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
    } else sigma_t_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma"))))
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
          k <- krig(z_t[, i], dx, d_new, d_cross, object$theta_z[t, i], g = NULL, sigma = FALSE,
                    tau2 = FALSE)
          z_new[[t]][, i] <- k$mean
        } else {
          k <- krig(z_t[, i], dx, d_new, d_cross, object$theta_z[t, i], g = NULL, tau2 = FALSE)
          z_new[[t]][, i] <- rand_mvn(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      
      # map z_new to w_new (separately for each dimension)
      w_new[[t]] <- matrix(nrow = m, ncol = D)
      for (i in 1:D) {
        if (mean_map) { 
          k <- krig(w_t[, i], sq_dist(z_t), sq_dist(z_new[[t]]), sq_dist(z_new[[t]], z_t),
                    object$theta_w[t, i], g = NULL, sigma = FALSE, tau2 = FALSE)
          w_new[[t]][, i] <- k$mean
        } else {
          k <- krig(w_t[, i], sq_dist(z_t), sq_dist(z_new[[t]]), sq_dist(z_new[[t]], z_t),
                    object$theta_w[t, i], g = NULL, tau2 = FALSE)
          w_new[[t]][, i] <- rand_mvn(1, k$mean, k$sigma)
        } 
      } # end of i for loop
      
      # map w_new to mu_t and sigma_t
      k <- krig(object$y, sq_dist(w_t), sq_dist(w_new[[t]]), sq_dist(w_new[[t]], w_t),
                object$theta_y[t], object$g[t], s2 = lite, sigma = !lite)
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
          k <- krig(z_t[, i], dx, d_new, d_cross, object$theta_z[t, i], g = NULL, sigma = FALSE,
                    tau2 = FALSE)
          z_new[, i] <- k$mean
        } else {
          k <- krig(z_t[, i], dx, d_new, d_cross, object$theta_z[t, i], g = NULL, tau2 = FALSE)
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
    } else sigma_t_sum <- Reduce("+", lapply(result, with, eval(parse(text = "sigma"))))
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

# Define Plot for S3 Objects --------------------------------------------------
#' @name plot
#' @title Plots object from "\code{deepgp}" package
#' 
#' @description Acts on a "\code{gp}", "\code{dgp2}", or "\code{dgp3}" object.  
#'     Generates trace plots for length scale and nugget hyperparameters.
#'     Generates plots of hidden layers for one-dimensional inputs.  Generates
#'     plots of the posterior mean and estimated 95\% prediction intervals for 
#'     one-dimensional inputs; generates heat maps of the posterior mean and 
#'     point-wise variance for two-dimensional inputs.
#'     
#' @details Trace plots are useful in assessing burn-in.  Hidden layer plots 
#'     are colored on a gradient - red lines represent earlier iterations and 
#'     yellow lines represent later iterations - to help assess burn-in of the 
#'     hidden layers.  These plots are meant to help in model fitting and 
#'     visualization.
#' 
#' @param x object of class \code{gp}, \code{dgp2}, or \code{dgp3}
#' @param trace logical indicating whether to generate trace plots
#' @param hidden logical indicating whether to generate plots of hidden layers
#'        ("\code{dgp2}" or "\code{dgp3}" only)
#' @param predict logical indicating whether to generate posterior predictive plot
#' @param ... N/A
#' 
#' @examples 
#' # See "deepgp-package", "fit_one_layer", "fit_two_layer", or "fit_three_layer"
#' # for an example
#' 
#' @rdname plot
NULL

# Plot One Layer Function -----------------------------------------------------
#' @rdname plot
#' @export

plot.gp <- function(x, trace = TRUE, predict = TRUE, ...) {
  
  # save and restore par settings
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # extract dimensions of x
  Dx <- ncol(x$x)

  # if mcmc only, change predict to FALSE
  if (is.null(x$mean)) predict <- FALSE

  if(trace) {
    par(mfrow = c(1, 2), mar = c(5, 4, 2, 2))
    plot(x$g, type = 'l', ylab = 'g', xlab = 'Iteration',
         main = 'Trace Plot of g')
    plot(x$theta, type = 'l', ylab = 'theta_y', xlab = 'Iteration',
         main = 'Trace Plot of theta')
  }

  if(predict) {
    if (Dx == 1) {
      par(mfrow = c(1, 1))
      if (is.null(x$Sigma)) {
        q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
      } else {
        y_samples <- rand_mvn(50, x$mean, x$Sigma_smooth)
        q1 <- x$mean + qnorm(0.05, 0, sqrt(diag(x$Sigma)))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(diag(x$Sigma)))
      }
      o <- order(x$x_new)
      plot(x$x_new[o], x$mean[o], type = 'l', xlab = 'X', ylab = 'Y', ylim = c(min(q1), max(q3)),
           col = 'blue', main = 'Posterior Mean and 95% PI')
      if (!is.null(x$Sigma)) {
        matlines(x$x_new[o], y_samples[o,], col = 'grey', lty = 1)
        lines(x$x_new[o], x$mean[o], col = 'blue')
      }
      lines(x$x_new[o], q1[o], col = 'blue', lty = 2)
      lines(x$x_new[o], q3[o], col = 'blue', lty = 2)
      points(x$x, x$y, pch = 20)
    } else if (Dx == 2) {
      if (!requireNamespace("akima", quietly = TRUE)) {
        stop("Package \"akima\" needed for this plot. Please install it.",
             call. = FALSE)
      }
      cols <- heat.colors(128)
      i1 <- akima::interp(x$x_new[, 1], x$x_new[, 2], x$mean)
      if (is.null(x$Sigma)) {
        i2 <- akima::interp(x$x_new[, 1], x$x_new[, 2], sqrt(x$s2))
      } else i2 <- akima::interp(x$x_new[, 1], x$x_new[, 2], sqrt(diag(x$Sigma)))
      par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
      image(i1, col = cols, main = 'Posterior Mean', xlab = 'X1', ylab = 'X2')
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
      image(i2, col = cols, main = 'Posterior Variance', xlab = 'X1', 
            ylab = 'X2')
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
    } else cat('Dimension of X too large for default plotting')
  }
}

# Plot Two Layer Function -----------------------------------------------------
#' @rdname plot
#' @export

plot.dgp2 <- function(x, trace = TRUE, hidden = FALSE, predict = TRUE, ...) {
  
  # save and restore par settings
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # extract dimensions of x and w
  Dx <- ncol(x$x)
  D <- ncol(x$w[[1]])

  # if mcmc only, change predict to FALSE
  if (is.null(x$mean)) predict <- FALSE

  if(trace) {
    par(mfrow = c(1, D + 2), mar = c(5, 4, 2, 2))
    plot(x$g, type = 'l', ylab = 'g', xlab = 'Iteration',
         main = 'Trace Plot of g')
    plot(x$theta_y, type = 'l', ylab = 'theta_y', xlab = 'Iteration',
         main = 'Trace Plot of theta_y')
    for (i in 1:D)
      plot(x$theta_w[, i], type = 'l', ylab = 'theta_w', xlab = 'Iteration',
           main = paste0('Trace Plot of theta_w [', i, ']'))
  }

  if (hidden) {
    if (Dx == 1 & D == 1) {
      # specify the hidden layers to plot
      indx <- floor(seq(from = 1, to = x$nmcmc, length = 50))
      if (indx[1] == 0) indx[1] = 1
      col <- heat.colors(50 + 10) # add ten to avoid using colors that are too light
      par(mfrow = c(2, D), mar = c(4, 4, 2, 2))
      
      # plot w to y
      o <- order(x$w[[indx[1]]])
      plot(x$w[[indx[1]]][o] - mean(x$w[[indx[1]]]), x$y, type = 'l', xlab = 'W', 
           ylab = 'Y', col = col[1], main = paste0('MCMC samples of W to Y'), 
           xlim = c(min(unlist(x$w)), max(unlist(x$w))))
      for (j in 2:length(indx)){
        o <- order(x$w[[indx[j]]])
        lines(x$w[[indx[j]]][o] - mean(x$w[[indx[j]]]), x$y, col = col[j])
      }
      
      # plot x to w
      o <- order(x$x)
      plot(x$x[o], x$w[[indx[1]]][o] - mean(x$w[[indx[1]]]), type = 'l', xlab = 'X', 
           ylab = 'W', col = col[1], main = paste0('MCMC samples of X to W'), 
           ylim = c(min(unlist(x$w)), max(unlist(x$w))))
      for (j in 2:length(indx)) {
        lines(x$x[o], x$w[[indx[j]]][o] - mean(x$w[[indx[j]]]), col = col[j])
      }
    } else cat('Default plotting not prepared for these dimensions')
  }

  if(predict) {
    if (Dx == 1){
      par(mfrow = c(1, 1), mar = c(5, 4, 2, 2))
      if (is.null(x$Sigma)) {
        q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
      } else {
        y_samples <- rand_mvn(50, x$mean, x$Sigma_smooth)
        q1 <- x$mean + qnorm(0.05, 0, sqrt(diag(x$Sigma)))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(diag(x$Sigma)))
      }
      o <- order(x$x_new)
      plot(x$x_new[o], x$mean[o], type = 'l', xlab = 'X', ylab = 'Y', ylim = c(min(q1), max(q3)),
           col = 'blue', main = 'Posterior Mean and 95% PI')
      if (!is.null(x$Sigma)) {
        matlines(x$x_new[o], y_samples[o,], col = 'grey', lty = 1)
        lines(x$x_new[o], x$mean[o], col = 'blue')
      }
      lines(x$x_new[o], q1[o], col = 'blue', lty = 2)
      lines(x$x_new[o], q3[o], col = 'blue', lty = 2)
      points(x$x, x$y, pch = 20)

    } else if (Dx == 2) {
      if (!requireNamespace("akima", quietly = TRUE)) {
        stop("Package \"akima\" needed for this plot. Please install it.",
             call. = FALSE)
      }
      cols <- heat.colors(128)
      i1 <- akima::interp(x$x_new[, 1], x$x_new[, 2], x$mean)
      if (is.null(x$Sigma)) {
        i2 <- akima::interp(x$x_new[, 1], x$x_new[, 2], sqrt(x$s2))
      } else i2 <- akima::interp(x$x_new[, 1], x$x_new[, 2], sqrt(diag(x$Sigma)))
      par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
      image(i1, col = cols, main = 'Posterior Mean', xlab = 'X1', ylab = 'X2')
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
      image(i2, col = cols, main = 'Posterior Variance', xlab = 'X1', ylab = 'X2')
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
    } else cat('Dimension of X too large for default plotting')
  }
}

# Plot Three Layer Function ---------------------------------------------------
#' @rdname plot
#' @export

plot.dgp3 <- function(x, trace = TRUE, hidden = FALSE, predict = TRUE, ...) {
  
  # save and restore par settings
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # extract dimensions of x, w, and z
  Dx <- ncol(x$x)
  D <- ncol(x$w[[1]])

  # if mcmc only, change predict to FALSE
  if (is.null(x$mean)) predict <- FALSE

  if(trace) {
    par(mfrow = c(2, D + 1), mar = c(5, 4, 2, 2))
    plot(x$g, type = 'l', ylab = 'g', xlab = 'Iteration',
         main = 'Trace Plot of g')
    plot(x$theta_y, type = 'l', ylab = 'theta_y', xlab = 'Iteration',
         main = 'Trace Plot of theta_y')
    for (i in 1:D)
      plot(x$theta_w[,i], type = 'l', ylab = 'theta_w', xlab = 'Iteration',
           main = paste0('Trace Plot of theta_w [', i, ']'))
    for (i in 1:D)
      plot(x$theta_z[,i], type = 'l', ylab = 'theta_z', xlab = 'Iteration',
           main = paste0('Trace Plot of theta_z [', i, ']'))
  }

  if(hidden) {
    if (Dx == 1 & D == 1) {
      # specify the hidden layers to plot
      indx <- floor(seq(from = 1, to = x$nmcmc, length = 50))
      if (indx[1] == 0) indx[1] <- 1
      col <- heat.colors(50 + 10)
      par(mfrow = c(3, 1), mar = c(4, 4, 2, 2))
      
      # plot w to y
      o <- order(x$w[[indx[1]]])
      plot(x$w[[indx[1]]][o] - mean(x$w[[indx[1]]]), x$y, type = 'l', xlab = 'W', 
           ylab = 'Y', col = col[1], main = paste0('MCMC samples of W to Y'), 
           xlim = c(min(unlist(x$w)), max(unlist(x$w))))
      for (j in 2:length(indx)) {
        o <- order(x$w[[indx[j]]])
        lines(x$w[[indx[j]]][o] - mean(x$w[[indx[j]]]), x$y, col = col[j])
      }

      # plot z to w
      o <- order(x$z[[indx[1]]])
      plot(x$z[[indx[1]]][o] - mean(x$z[[indx[1]]]), 
           x$w[[indx[1]]][o] - mean(x$w[[indx[1]]]), type = 'l', xlab = 'Z', 
           ylab = 'W', col = col[1], main = paste0('MCMC samples of Z to W'), 
           xlim = c(min(unlist(x$z)), max(unlist(x$z))),
           ylim = c(min(unlist(x$w)), max(unlist(x$w))))
      for (j in 2:length(indx)) {
        o <- order(x$z[[indx[j]]])
        lines(x$z[[indx[j]]][o] - mean(x$z[[indx[j]]]), 
              x$w[[indx[j]]][o] - mean(x$w[[indx[j]]]), col = col[j])
      }
      
      # plot x to z
      o <- order(x$x)
      plot(x$x[o], x$z[[indx[1]]][o] - mean(x$z[[indx[1]]]), type = 'l', 
           xlab = 'X', ylab = 'Z', col = col[1], 
           main = paste0('MCMC samples of X to Z'), 
           ylim = c(min(unlist(x$z)), max(unlist(x$z))))
      for (j in 2:length(indx)) {
        lines(x$x[o], x$z[[indx[j]]][o] - mean(x$z[[indx[j]]]), col = col[j])
      }
    } else cat('Default plotting not prepared for these dimensions')
  }

  if(predict) {
    if (Dx == 1) {
      par(mfrow = c(1, 1), mar = c(5, 4, 2, 2))
      if (is.null(x$Sigma)) {
        q1 <- x$mean + qnorm(0.05, 0, sqrt(x$s2))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(x$s2))
      } else {
        y_samples <- rand_mvn(50, x$mean, x$Sigma_smooth)
        q1 <- x$mean + qnorm(0.05, 0, sqrt(diag(x$Sigma)))
        q3 <- x$mean + qnorm(0.95, 0, sqrt(diag(x$Sigma)))
      }
      o <- order(x$x_new)
      plot(x$x_new[o], x$mean[o], type = 'l', xlab = 'X', ylab = 'Y', ylim = c(min(q1), max(q3)),
           col = 'blue', main = 'Posterior Mean and 95% PI')
      if (!is.null(x$Sigma)) {
        matlines(x$x_new[o], y_samples[o,], col = 'grey', lty = 1)
        lines(x$x_new[o], x$mean[o], col = 'blue')
      }
      lines(x$x_new[o], q1[o], col = 'blue', lty = 2)
      lines(x$x_new[o], q3[o], col = 'blue', lty = 2)
      points(x$x, x$y, pch = 20)
    } else if (Dx == 2) {
      if (!requireNamespace("akima", quietly = TRUE)) {
        stop("Package \"akima\" needed for this function to work. Please install it.",
             call. = FALSE)
      }
      cols <- heat.colors(128)
      i1 <- akima::interp(x$x_new[, 1], x$x_new[, 2], x$mean)
      if (is.null(x$Sigma)) {
        i2 <- akima::interp(x$x_new[, 1], x$x_new[, 2], sqrt(x$s2))
      } else i2 <- akima::interp(x$x_new[, 1], x$x_new[, 2], sqrt(diag(x$Sigma)))
      par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
      image(i1, col = cols, main = 'Posterior Mean', xlab = 'X1', ylab = 'X2')
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
      image(i2, col = cols, main = 'Posterior Variance', xlab = 'X1', ylab = 'X2')
      points(x$x[, 1], x$x[, 2], pch = 20, cex = 0.5)
    } else cat('Dimension of X too large for default plotting')
  }
}

