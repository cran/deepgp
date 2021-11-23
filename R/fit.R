
# Function Contents -----------------------------------------------------------
# External (see documentation below):
#   fit_one_layer
#   fit_two_layer
#   fit_three_layer

# Fit One Layer Function ------------------------------------------------------
#' @title MCMC sampling for one layer GP
#' @description Conducts MCMC sampling of hyperparameters for a one layer 
#'     GP.  Length scale parameter \code{theta} governs 
#'     the strength of the correlation and nugget parameter \code{g} 
#'     governs noise.  In Matern covariance, \code{v} governs smoothness.
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
#'     in \code{settings}.  These priors are designed for \code{x} scaled 
#'     to [0, 1] and \code{y} scaled to have mean 0 and variance 1.  
#'     
#'     The output object of class \code{gp} is designed for use with 
#'     \code{continue}, \code{trim}, and \code{predict}.
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param nmcmc number of MCMC iterations
#' @param verb logical indicating whether to print iteration progress
#' @param g_0 initial value for \code{g}
#' @param theta_0 initial value for \code{theta}
#' @param true_g if true nugget is known it may be specified here (set to a 
#'        small value to make fit deterministic).  Note - values that are too 
#'        small may cause numerical issues in matrix inversions.
#' @param settings hyperparameters for proposals and priors on \code{g} and 
#'        \code{theta}
#' @param cov covariance kernel, either Matern or squared exponential (\code{exp2})
#' @param v Matern smoothness parameter (only used if \code{cov = "matern"})
#' @return a list of the S3 class \code{gp} with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{settings}: copy of proposal/prior settings
#'   \item \code{cov}: copy of covariance kernel setting
#'   \item \code{v}: copy of Matern smoothness parameter (if \code{cov = "matern"})
#'   \item \code{g}: vector of MCMC samples for \code{g}
#'   \item \code{theta}: vector of MCMC samples for \code{theta}
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." \emph{Technometrics, to appear;} arXiv:2012.08015. 
#'     \cr\cr
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
#'   if (x > 0.4 & x < 0.6) return(10 * (x - 0.5))
#' }
#' x <- seq(0.05, 0.95, length = 7)
#' y <- sapply(x, f)
#' x_new <- seq(0, 1, length = 100)
#' 
#' # Fit model and calculate EI
#' fit <- fit_one_layer(x, y, nmcmc = 500)
#' fit <- trim(fit, 400)
#' fit <- predict(fit, x_new, EI = TRUE)
#' 
#' \donttest{
#' # One Layer and EI -----------------------------------------------------------
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
#' fit <- predict(fit, xx, EI = TRUE)
#' 
#' # Visualize Fit
#' plot(fit)
#' par(new = TRUE) # overlay EI
#' plot(xx, fit$EI, type = 'l', lty = 2, axes = FALSE, xlab = '', ylab = '')
#' 
#' # Select next design point
#' x_new <- xx[which.max(fit$EI)]
#' 
#' # Evaluate fit
#' rmse(yy, fit$mean) # lower is better
#' }
#' 
#' @export

fit_one_layer <- function(x, y, nmcmc = 10000, verb = TRUE, g_0 = 0.01, 
                          theta_0 = 0.1, true_g = NULL, 
                          settings = list(l = 1, u = 2, 
                                          alpha = list(g = 1.5, theta = 1.5), 
                                          beta = list(g = 3.9, theta = 3.9/1.5)),
                          cov = c("matern", "exp2"), v = 2.5) {

  tic <- proc.time()[3]
  cov <- match.arg(cov)

  # Check inputs
  if (is.numeric(x)) x <- as.matrix(x)
  settings <- check_settings(settings, layers = 1)
  test <- check_inputs(x, y, true_g) # returns NULL if all checks pass
  if (cov == "matern")
    if(!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")

  # Create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, cov = cov)
  if (cov == "matern") out$v <- v

  # Conduct MCMC
  samples <- gibbs_one_layer(x, y, nmcmc, verb, g_0, theta_0, true_g,
                              settings, cov, v)
  
  out <- c(out, samples)
  toc <- proc.time()[3]
  out$time <- toc - tic
  class(out) <- "gp"
  return(out)
}

# Fit Two Layer Function ------------------------------------------------------
#' @title MCMC sampling for two layer deep GP
#' @description Conducts MCMC sampling of hyperparameters and hidden layer 
#'     \code{w} for a two layer deep GP.  Separate length scale 
#'     parameters \code{theta_w} and \code{theta_y} govern the correlation 
#'     strength of the hidden layer and outer layer respectively.  Nugget 
#'     parameter \code{g} governs noise on the outer layer.  In Matern 
#'     covariance, \code{v} governs smoothness.
#'
#' @details Maps inputs \code{x} through hidden layer \code{w} to outputs 
#'     \code{y}.  Conducts sampling of the hidden layer using Elliptical 
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
#'     in \code{settings}.  These priors are designed for \code{x} scaled to 
#'     [0, 1] and \code{y} scaled to have mean 0 and variance 1.  
#'     
#'     The output object of class \code{dgp2} is designed for use with 
#'     \code{continue}, \code{trim}, and \code{predict}. If \code{w_0} is 
#'     of dimension \code{nrow(x) - 1} by 
#'     \code{D}, the final row is predicted using kriging.  This is helpful in 
#'     sequential design when adding a new input location and starting the MCMC 
#'     at the place where the previous MCMC left off.
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param D integer designating dimension of hidden layer, defaults to 
#'        dimension of \code{x}
#' @param nmcmc number of MCMC iterations
#' @param verb logical indicating whether to print iteration progress
#' @param w_0 initial value for hidden layer \code{w} (must be matrix 
#'        of dimension \code{nrow(x)} by \code{D} or  dimension 
#'        \code{nrow(x) - 1} by \code{D}).  Defaults to the identity mapping.
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
#' @param cov covariance kernel, either Matern or squared exponential (\code{exp2})
#' @param v Matern smoothness parameter (only used if \code{cov = "matern"})
#' @return a list of the S3 class \code{dgp2} with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{settings}: copy of proposal/prior settings
#'   \item \code{cov}: copy of covariance kernel setting
#'   \item \code{v}: copy of Matern smoothness parameter (if \code{cov = "matern"}) 
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
#'     Process Surrogates." \emph{Technometrics, to appear;} arXiv:2012.08015. 
#'     \cr\cr
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
#'   if (x > 0.4 & x < 0.6) return(10 * (x - 0.5))
#' }
#' x <- seq(0.05, 0.95, length = 7)
#' y <- sapply(x, f)
#' x_new <- seq(0, 1, length = 100)
#' 
#' # Fit model and calculate ALC
#' fit <- fit_two_layer(x, y, nmcmc = 500, cov = "exp2")
#' fit <- trim(fit, 400)
#' fit <- predict(fit, x_new, store_latent = TRUE)
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
#' fit <- fit_two_layer(x, y, D = 1, nmcmc = 9000, cov = "exp2")
#' fit <- continue(fit, 1000)
#' plot(fit) # investigate trace plots
#' fit <- trim(fit, 8000, 2)
#' 
#' # Option 1 - calculate ALC from MCMC iterations
#' alc <- ALC(fit, xx)
#' 
#' # Option 2 - calculate ALC after predictions
#' fit <- predict(fit, xx, store_latent = TRUE)
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
                          nmcmc = 10000, verb = TRUE, 
                          w_0 = suppressWarnings(matrix(x, nrow = length(y), ncol = D)), 
                          g_0 = 0.01, theta_y_0 = 0.1, theta_w_0 = 0.1, true_g = NULL,
                          settings = list(l = 1, u = 2, 
                                          alpha = list(g = 1.5, theta_w = 1.5, theta_y = 1.5), 
                                          beta = list(g = 3.9, theta_w = 3.9/4, theta_y = 3.9/6)),
                          cov = c("matern", "exp2"), v = 2.5) {
  
  tic <- proc.time()[3]
  cov <- match.arg(cov)
  
  # Check inputs
  if (is.numeric(x)) x <- as.matrix(x)
  if (!is.matrix(w_0)) w_0 <- as.matrix(w_0)
  settings <- check_settings(settings, layers = 2)
  test <- check_inputs(x, y, true_g, w_0, D) # returns NULL if all checks pass
  if (cov == "matern")
    if(!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  if (length(theta_w_0) == 1) theta_w_0 <- rep(theta_w_0, D)
  
  # Create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, cov = cov)
  if (cov == "matern") out$v <- v
  
  # If w_0 is from previous sequential design iteration, predict at new point
  if (nrow(w_0) == nrow(x) - 1) 
    w_0 <- fill_final_row(x, w_0, D, theta_w_0, cov, v)
  
  # Conduct MCMC
  samples <- gibbs_two_layer(x, y, nmcmc, D, verb, w_0, g_0, theta_y_0,
                             theta_w_0, true_g, settings, cov, v)
  
  out <- c(out, samples)
  toc <- proc.time()[3]
  out$time <- toc - tic
  class(out) <- "dgp2"
  return(out)
}

# Fit Three Layer Function ----------------------------------------------------
#' @title MCMC sampling for three layer deep GP
#' @description Conducts MCMC sampling of hyperparameters, hidden layer 
#'     \code{z}, and hidden layer \code{w} for a three layer deep GP.  
#'     Separate length scale parameters \code{theta_z}, 
#'     \code{theta_w}, and \code{theta_y} govern the correlation 
#'     strength of the inner layer, middle layer, and outer layer respectively.  
#'     Nugget parameter \code{g} governs noise on the outer layer.  In Matern 
#'     covariance, \code{v} governs smoothness.
#'
#' @details Maps inputs \code{x} through hidden layer \code{z} then hidden
#'     layer \code{w} to outputs \code{y}.  Conducts sampling of the hidden 
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
#'     designed for \code{x} scaled to [0, 1] and \code{y} scaled to have 
#'     mean 0 and variance 1.  
#'     
#'     The output object of class \code{dgp3} is designed for use with 
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
#' @param verb logical indicating whether to print iteration progress
#' @param w_0 initial value for hidden layer \code{w} (must be matrix 
#'        of dimension \code{nrow(x)} by \code{D} or  dimension 
#'        \code{nrow(x) - 1} by \code{D}).  Defaults to the identity mapping.
#' @param z_0 initial value for hidden layer \code{z} (must be matrix 
#'        of dimension \code{nrow(x)} by \code{D} or  dimension 
#'        \code{nrow(x) - 1} by \code{D}).  Defaults to the identity mapping.
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
#' @param cov covariance kernel, either Matern or squared exponential (\code{exp2})
#' @param v Matern smoothness parameter (only used if \code{cov = "matern"})
#' @return a list of the S3 class \code{dgp3} with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{settings}: copy of proposal/prior settings
#'   \item \code{cov}: copy of covariance kernel setting
#'   \item \code{v}: copy of Matern smoothness parameter (if \code{cov = "matern"}) 
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
#'     Process Surrogates." \emph{Technometrics, to appear;} arXiv:2012.08015. 
#'     \cr\cr
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
#'   if (x > 0.4 & x < 0.6) return(10 * (x - 0.5))
#' }
#' x <- seq(0.05, 0.95, length = 7)
#' y <- sapply(x, f)
#' x_new <- seq(0, 1, length = 100)
#' 
#' # Fit model and calculate IMSPE
#' fit <- fit_three_layer(x, y, nmcmc = 500, cov = "exp2")
#' fit <- trim(fit, 400)
#' fit <- predict(fit, x_new, store_latent = TRUE)
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
#' fit <- fit_three_layer(x, y, D = 1, nmcmc = 10000, cov = "exp2")
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
                            nmcmc = 10000, verb = TRUE, 
                            w_0 = suppressWarnings(matrix(x, nrow = length(y), ncol = D)), 
                            z_0 = suppressWarnings(matrix(x, nrow = length(y), ncol = D)), 
                            g_0 = 0.01, theta_y_0 = 0.1, theta_w_0 = 0.1, 
                            theta_z_0 = 0.1, true_g = NULL,
                            settings = list(l = 1, u = 2, 
                                            alpha = list(g = 1.5, theta_z = 1.5, theta_w = 1.5, theta_y = 1.5), 
                                            beta = list(g = 3.9, theta_z = 3.9/4, theta_w = 3.9/12, theta_y = 3.9/6)),
                            cov = c("matern", "exp2"), v = 2.5) {

  tic <- proc.time()[3]
  cov <- match.arg(cov)
  
  # Check inputs
  if (is.numeric(x)) x <- as.matrix(x)
  if (!is.matrix(w_0)) w_0 <- as.matrix(w_0)
  if (!is.matrix(z_0)) z_0 <- as.matrix(z_0)
  settings <- check_settings(settings, layers = 3)
  test <- check_inputs(x, y, true_g, w_0, D, z_0)
  if (cov == "matern")
    if(!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  if (length(theta_w_0) == 1) theta_w_0 <- rep(theta_w_0, D)
  if (length(theta_z_0) == 1) theta_z_0 <- rep(theta_z_0, D)
  
  # Create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, cov = cov)
  if (cov == "matern") out$v <- v
  
  # If z_0/w_0 are from previous sequential design iteration, predict at new point
  if (nrow(z_0) == nrow(x) - 1) 
    z_0 <- fill_final_row(x, z_0, D, theta_z_0, cov, v)
  if (nrow(w_0) == nrow(x) - 1) 
    w_0 <- fill_final_row(z_0, w_0, D, theta_w_0, cov, v)
  
  # Conduct MCMC
  samples <- gibbs_three_layer(x, y, nmcmc, D, verb, w_0, z_0, g_0, theta_y_0,
                               theta_w_0, theta_z_0, true_g, settings, cov, v)
  
  out <- c(out, samples)
  toc <- proc.time()[3]
  out$time <- toc - tic
  class(out) <- "dgp3"
  return(out)
}
