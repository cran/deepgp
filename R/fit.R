
# Function Contents -----------------------------------------------------------
# External (see documentation below):
#   fit_one_layer
#   fit_two_layer
#   fit_three_layer

# Imported Functions ----------------------------------------------------------
#' @importFrom grDevices heat.colors
#' @importFrom graphics image lines matlines par plot points
#' @importFrom stats cov dgamma dnorm pnorm qnorm rnorm runif var
#' @importFrom parallel makeCluster detectCores stopCluster

# Fit One Layer Function ------------------------------------------------------
#' @title MCMC sampling for one layer GP
#' @description Conducts MCMC sampling of hyperparameters for a one layer 
#'     GP.  Covariance structure is based on inverse exponentiated squared 
#'     euclidean distance with length scale parameter "\code{theta}" governing 
#'     the strength of the correlation and nugget parameter "\code{g}" 
#'     governing noise.
#'
#' @details Utilizes Metropolis Hastings sampling of the length scale and
#'     nugget parameters with proposals and priors controlled by \code{settings}.
#'     Proposals for \code{g} and \code{theta} follow a uniform sliding window 
#'     scheme, e.g. 
#'     
#'     \code{g_star <- runif(1, l * g_t / u, u * g_t / l)}, 
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{settings}.

#'     Priors on \code{g} and \code{theta} follow Gamma distributions with 
#'     shape parameter (\code{alpha}) and rate parameter (\code{beta}) provided 
#'     in \code{settings}.  These priors are designed for "\code{x}" scaled 
#'     to [0,1] and "\code{y}" scaled to have mean 0 and variance 1.  
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
#' @param true_g if true nugget is known it may be specified here (set to a 
#'        small value to make fit deterministic).  Note - values that are too 
#'        small may cause numerical issues in matrix inversions.
#' @param settings hyperparameters for proposals and priors on \code{g} and 
#'        \code{theta}
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
#' Gramacy, RB. \emph{Surrogates: Gaussian Process Modeling, Design, and 
#'     Optimization for the Applied Sciences}. Chapman Hall, 2020.
#' 
#' @examples 
#' # Toy example (runs in less than 5 seconds) --------------------------------
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
#' # One Layer and EI ---------------------------------------------------------
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

fit_one_layer <- function(x, y, nmcmc = 10000, trace = TRUE, g_0 = 0.01, 
                          theta_0 = 0.1, true_g = NULL, 
                          settings = list(l = 1, u = 2, 
                                          alpha = list(g = 1.5, theta = 1.5), 
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
  if (is.null(true_g) & (mean(y) < -0.9 | mean(y) > 0.9 | 
                         var(y) < 0.1 | var(y) > 1.9))
    warning('this function is designed for y scaled to mean 0 and variance 1')
  
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
                       beta = settings$beta$g, l = settings$l, u = settings$u, 
                       ll_prev = ll)
      g[j] <- samp$g
      ll <- samp$ll
    } else g[j] <- true_g

    # sample theta
    samp <- sample_theta(y, dx, g[j], theta[j-1], alpha = settings$alpha$theta,
                         beta = settings$beta$theta, l = settings$l, 
                         u = settings$u, outer = TRUE, ll_prev = ll)
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
#'     "\code{y}".  Conducts sampling of the hidden layer using Elliptical 
#'     Slice sampling.  Utilizes Metropolis Hastings sampling of the length 
#'     scale and nugget parameters with proposals and priors controlled by 
#'     \code{settings}.  Proposals for \code{g}, \code{theta_y}, and 
#'     \code{theta_w} follow a uniform sliding window scheme, e.g.
#'     
#'     \code{g_star <- runif(1, l * g_t / u, u * g_t / l)}, 
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{settings}.   
#'     Priors on \code{g} and \code{theta} follow Gamma distributions with 
#'     shape parameter (\code{alpha}) and rate parameter (\code{beta}) provided 
#'     in \code{settings}.  These priors are designed for "\code{x}" scaled to 
#'     [0,1] and "\code{y}" scaled to have mean 0 and variance 1.  
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
#' @param D integer designating dimension of hidden layer, defaults to 
#'        dimension of \code{x}
#' @param nmcmc number of MCMC iterations
#' @param trace logical indicating whether to print iteration progress
#' @param w_0 initial value for hidden layer \code{w}, defaults to identity 
#'        mapping (must be matrix of dimension \code{nrow(x)} by \code{D} or 
#'        dimension \code{nrow(x) - 1} by \code{D})
#' @param g_0 initial value for \code{g}
#' @param theta_y_0 initial value for \code{theta_y} (length scale of outer 
#'        layer)
#' @param theta_w_0 initial value for \code{theta_w} (length scale of inner 
#'        layer), may be single value or vector of length \code{D}
#' @param true_g if true nugget is known it may be specified here (set to a 
#'        small value to make fit deterministic).  Note - values that are too 
#'        small may cause numerical issues in matrix inversions.
#' @param settings hyperparameters for proposals and priors on \code{g}, 
#'        \code{theta_y}, and \code{theta_w}
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
#' # Toy example (runs in less than 5 seconds) --------------------------------
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
#' # Two Layer and ALC --------------------------------------------------------
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

fit_two_layer <- function(x, y, D = ifelse(is.matrix(x), ncol(x), 1), 
                          nmcmc = 10000, trace = TRUE,
                          w_0 = suppressWarnings(matrix(x, nrow = length(y), ncol = D)), 
                          g_0 = 0.01, theta_y_0 = 0.1, theta_w_0 = 0.1, 
                          true_g = NULL,
                          settings = list(l = 1, u = 2, 
                                          alpha = list(g = 1.5, theta_w = 1.5, theta_y = 1.5), 
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
  if (is.null(true_g) & (mean(y) < -0.9 | mean(y) > 0.9 | 
                         var(y) < 0.1 | var(y) > 1.9))
    warning('this function is designed for y scaled to mean 0 and variance 1')
  
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
  if (length(theta_w_0) != D & length(theta_w_0) == 1) 
    theta_w_0 <- rep(theta_w_0, D)

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
      new_w[i] <- krig(w_0[, i], dx[1:(n-1), 1:(n-1)], 
                       d_cross = sq_dist(new_x, old_x),
                       theta = theta_w_0[i], g = NULL, mean = TRUE, s2 = FALSE, 
                       sigma = FALSE, tau2 = FALSE)$mean
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
      samp <- sample_g(y, dw, g[j - 1], theta_y[j - 1], 
                       alpha = settings$alpha$g, beta = settings$beta$g, 
                       l = settings$l, u = settings$u, ll_prev = ll_outer)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else g[j] <- true_g

    # sample theta_y
    samp <- sample_theta(y, dw, g[j], theta_y[j - 1], 
                         alpha = settings$alpha$theta_y, 
                         beta = settings$beta$theta_y, l = settings$l, 
                         u = settings$u, outer = TRUE, ll_prev = ll_outer)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    
    # sample theta_w
    for (i in 1:D) {
      samp <- sample_theta(w_out[[j - 1]][, i], dx, NULL, theta_w[j-1, i],
                           alpha = settings$alpha$theta_w, 
                           beta = settings$beta$theta_w,
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
#'     sampling of the length scale and nugget parameters with proposals and 
#'     priors controlled by \code{settings}.  Proposals for \code{g}, 
#'     \code{theta_y}, \code{theta_w}, and \code{theta_z} follow a uniform 
#'     sliding window scheme, e.g.
#'     
#'     \code{g_star <- runif(1, l * g_t / u, u * g_t / l)},
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{settings}.  
#'     Priors on \code{g}, \code{theta_y}, \code{theta_w}, and \code{theta_z} 
#'     follow Gamma distributions with shape parameter (\code{alpha}) and rate 
#'     parameter (\code{beta}) provided in \code{settings}.  These priors are 
#'     designed for "\code{x}" scaled to [0,1] and "\code{y}" scaled to have 
#'     mean 0 and variance 1.  
#'     
#'     The output object of class "\code{dgp3}" is designed for use with 
#'     \code{continue}, \code{trim}, and \code{predict}. If \code{z_0} and 
#'     \code{w_0} are of dimension \code{nrow(x) - 1} by \code{D}, the final 
#'     rows are predicted using kriging.  This is helpful in sequential design 
#'     when adding a new input location and starting the MCMC at the place 
#'     where the previous MCMC left off.
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param D integer designating dimension of hidden layers, defaults to 
#'        dimension of \code{x}
#' @param nmcmc number of MCMC iterations
#' @param trace logical indicating whether to print iteration progress
#' @param w_0 initial value for hidden layer \code{w}, defaults to identity 
#'        mapping (must be matrix of dimension \code{nrow(x)} by \code{D} or 
#'        dimension \code{nrow(x) - 1} by \code{D})
#' @param z_0 initial value for hidden layer \code{z}, defaults to identity 
#'        mapping (must be matrix of dimension \code{nrow(x)} by \code{D} or 
#'        dimension \code{nrow(x) - 1} by \code{D})
#' @param g_0 initial value for \code{g}
#' @param theta_y_0 initial value for \code{theta_y} (length scale of outer 
#'        layer)
#' @param theta_w_0 initial value for \code{theta_w} (length scale of middle 
#'        layer), may be single value or vector of length \code{D}
#' @param theta_z_0 initial value for \code{theta_z} (length scale of inner 
#'        layer), may be single value or vector of length \code{D}
#' @param true_g if true nugget is known it may be specified here (set to a 
#'        small value to make fit deterministic).  Note - values that are too 
#'        small may cause numerical issues in matrix inversions.
#' @param settings hyperparameters for proposals and priors on \code{g}, 
#'        \code{theta_y}, \code{theta_w}, and \code{theta_z}
#' @return a list of the S3 class "\code{dgp3}" with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{settings}: copy of proposal/prior settings
#'   \item \code{g}: vector of MCMC samples for \code{g}
#'   \item \code{theta_y}: vector of MCMC samples for \code{theta_y} (length 
#'         scale of outer layer)
#'   \item \code{theta_w}: matrix of MCMC samples for \code{theta_w} (length 
#'         scale of middle layer)
#'   \item \code{theta_z}: matrix of MCMC samples for \code{theta_z} (length 
#'         scale of inner layer)
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
#' # Toy example (runs in less than 5 seconds) --------------------------------
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
#' # Three Layer and IMSE -----------------------------------------------------
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

fit_three_layer <- function(x, y, D = ifelse(is.matrix(x), ncol(x), 1), 
                            nmcmc = 10000, trace = TRUE,
                            w_0 = suppressWarnings(matrix(x, nrow = length(y), ncol = D)),
                            z_0 = suppressWarnings(matrix(x, nrow = length(y), ncol = D)), 
                            g_0 = 0.01, theta_y_0 = 0.1, theta_w_0 = 0.1, 
                            theta_z_0 = 0.1, true_g = NULL,
                            settings = list(l = 1, u = 2, 
                                            alpha = list(g = 1.5, theta_z = 1.5, theta_w = 1.5, theta_y = 1.5), 
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
  if (is.null(true_g) & (mean(y) < -0.9 | mean(y) > 0.9 | 
                         var(y) < 0.1 | var(y) > 1.9))
    warning('this function is designed for y scaled to mean 0 and variance 1')
  
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
  if (length(theta_w_0) != D & length(theta_w_0) == 1) 
    theta_w_0 <- rep(theta_w_0, D)
  if (length(theta_z_0) != D & length(theta_z_0) == 1) 
    theta_z_0 <- rep(theta_z_0, D)

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
      new_z[i] <- krig(z_0[, i], dx[1:(n-1), 1:(n-1)], 
                       d_cross = sq_dist(new_x, old_x), theta = theta_z_0[i],
                       g = NULL, mean = TRUE, s2 = FALSE, sigma = FALSE,
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
      new_w[i] <- krig(w_0[, i], sq_dist(old_z), 
                       d_cross = sq_dist(new_z, old_z), theta = theta_w_0[i], 
                       g = NULL, mean = TRUE, s2 = FALSE, sigma = FALSE,
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
      samp <- sample_g(y, dw, g[j-1], theta_y[j-1], alpha = settings$alpha$g, 
                       beta = settings$beta$g, l = settings$l, u = settings$u, 
                       ll_prev = ll_outer)
      g[j] <- samp$g
      ll_outer <- samp$ll
    } else g[j] <- true_g

    # sample theta_y
    samp <- sample_theta(y, dw, g[j], theta_y[j-1], 
                         alpha = settings$alpha$theta_y,
                         beta = settings$beta$theta_y, l = settings$l, 
                         u = settings$u, outer = TRUE, ll_prev = ll_outer)
    theta_y[j] <- samp$theta
    ll_outer <- samp$ll
    
    # sample theta_w
    ll_mid <- 0 # re-calculated each time since we have a new z
    for (i in 1:D) {
      samp <- sample_theta(w_out[[j-1]][, i], dz, g = NULL, theta_w[j-1, i],
                           alpha = settings$alpha$theta_w, 
                           beta = settings$beta$theta_w, l = settings$l, 
                           u = settings$u, outer = FALSE)
      theta_w[j, i] <- samp$theta
      ll_mid <- ll_mid + samp$ll
    }

    # sample theta_z
    for (i in 1:D) {
      samp <- sample_theta(z_out[[j-1]][, i], dx, g = NULL, theta_z[j-1, i], 
                           alpha = settings$alpha$theta_z, 
                           beta = settings$beta$theta_z, l = settings$l, 
                           u = settings$u, outer = FALSE)
      theta_z[j, i] <- samp$theta
    }

    # sample z
    samp <- sample_z(w_out[[j-1]], z_out[[j-1]], dz, dx, g = NULL, theta_w[j, ], 
                     theta_z[j, ], ll_prev = ll_mid)
    z_out[[j]] <- samp$z
    dz <- samp$dz

    # sample w
    samp <- sample_w(y, w_out[[j-1]], dw, dz, g = g[j], theta_y[j], 
                     theta_w[j, ], ll_prev = ll_outer)
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
