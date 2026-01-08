
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   post_sample.gp
#   post_sample.gpvec
#   post_sample.dgp2
#   post_sample.dgp2vec
#   post_sample.dgp3
#   post_sample.dgp3vec

# post_sample S3 class --------------------------------------------------------
#' @name post_sample
#' @title Generates joint posterior samples from a trained GP/DGP
#' @description Acts on a \code{gp}, \code{gpvec}, \code{dgp2}, \code{dgp2vec},
#'     \code{dgp3}, or \code{dgp3vec} object.  Generates joint samples from the 
#'     posterior distribution at the provided locations.
#' 
#' @details By default, one sample is generated per each MCMC iteration.  This
#'     may be increased with the \code{nper} argument.
#'     
#'     SNOW parallelization reduces computation time but requires 
#'     more memory storage.
#' 
#' @param object object from \code{fit_one_layer}, \code{fit_two_layer}, or 
#'        \code{fit_three_layer} with burn-in already removed
#' @param x_new vector or matrix of predictive input locations
#' @param nper the number of samples to generate from each MCMC iteration.  The
#'        total number of samples will equal \code{nper*object$nmcmc}
#' @param m size of Vecchia conditioning sets (only for fits with 
#'        \code{vecchia = TRUE}), defaults to the twice the \code{m} used for MCMC
#' @param ord_new optional ordering for Vecchia approximation with \code{lite = FALSE}, 
#'        must correspond to rows of \code{x_new}, defaults to random, is 
#'        applied to all layers in deeper models
#' @param grad logical indicating whether to additionally calculate/return 
#'        samples of the gradient (one and two layer models only)
#' @param mean_map logical indicating whether to map hidden layers using 
#'        conditional mean (\code{mean_map = TRUE}) or using a random sample
#'        from the full MVN distribution (two or three layer models only)
#' @param cores number of cores to use for SNOW parallelization
#' @param ... N/A
#'
#' @return If \code{grad = FALSE}, returns matrix of samples.  Rows correspond to 
#'         \code{x_new} locations.  If \code{grad = TRUE}, returns a list with 
#'         \code{y} and \code{dydx} containing the respective samples.
#' @references 
#' Sauer, A. (2023). Deep Gaussian process surrogates for computer experiments. 
#'      *Ph.D. Dissertation, Department of Statistics, Virginia Polytechnic Institute and State University.*
#'      \cr\cr

#' @examples 
#' \donttest{
#' # Simple step function
#' f <- function(x) {
#'   return(pnorm((x - 0.5) / 0.065))
#' }
#'
#' # Training data
#' x <- seq(0, 1, length = 5)
#' y <- f(x)
#'
#' # Testing data
#' xx <- seq(0, 1, length = 100)
#' yy <- f(xx)
#'
#' plot(xx, yy, type = "l")
#' points(x, y, col = 2)
#'
#' # Conduct MCMC
#' fit <- fit_two_layer(x, y, nmcmc = 2000, true_g = 1e-6, cov = "exp2")
#' plot(fit, hidden = TRUE)
#' fit <- trim(fit, 1000, 2)
#' 
#' # Generate posterior samples, including gradients
#' samples <- post_sample(fit, xx, grad = TRUE, cores = 1)
#'
#' # Plot samples
#' par(mfrow = c(1, 2))
#' matplot(xx, t(samples$y), type = "l")
#' points(x, y, pch = 20)
#' matplot(xx, t(samples$dy), type = "l")
#' }
#'
#' @rdname post_sample
#' @export

post_sample <- function(object, x_new, nper = 1, ...)
  UseMethod("post_sample", object)

# post_sample.gp --------------------------------------------------------------
#' @rdname post_sample
#' @export

post_sample.gp <- function(object, x_new, nper = 1, grad = FALSE, cores = 1, ...) {
  
  cores <- check_cores(cores, object$nmcmc)
  if (grad) {
    if (object$v != 999) stop("grad only offered with cov = 'exp2'")
    if (min(object$g) < 1e-6) message("Warning: small g may cause numerical issues")
    if (cores > 1 & !requireNamespace("abind", quietly = TRUE))
      stop("Package \"abind\" needed for cores > 1 with grad = TRUE. Please install it.",
         call. = FALSE)
  }
  
  settings <- list(lite = FALSE, grad = grad, return_all = FALSE, EI = FALSE,
                   entropy_limit = NULL, cores = cores, nper = nper)
  samples <- predict_shallow(object, x_new, settings, samples_only = TRUE)
  return(samples)
}

# post_sample.gpvec -----------------------------------------------------------
#' @rdname post_sample
#' @export

post_sample.gpvec <- function(object, x_new, nper = 1, m = NULL,
                              ord_new = NULL, grad = FALSE, cores = 1, ...) {
  
  cores <- check_cores(cores, object$nmcmc)
  if (grad) {
    if (object$v != 999) stop("grad only offered with cov = 'exp2'")
    if (min(object$g) < 1e-6) message("Warning: small g may cause numerical issues")
    if (cores > 1 & !requireNamespace("abind", quietly = TRUE))
      stop("Package \"abind\" needed for cores > 1 with grad = TRUE. Please install it.",
           call. = FALSE)
  }
  
  settings <- list(ord_new = ord_new, lite = FALSE, grad = grad, 
                   return_all = FALSE, EI = FALSE, entropy_limit = NULL, 
                   cores = cores, nper = nper)
  samples <- predict_shallow_vec(object, x_new, m, settings, samples_only = TRUE)
  return(samples)
}

# post_sample.dgp2 ------------------------------------------------------------
#' @rdname post_sample
#' @export

post_sample.dgp2 <- function(object, x_new, nper = 1, grad = FALSE, 
                             mean_map = TRUE, cores = 1, ...) {
  
  cores <- check_cores(cores, object$nmcmc)
  if (grad) {
    if (object$v != 999) stop("grad only offered with cov = 'exp2'")
    if (min(object$g) < 1e-6) message("Warning: small g may cause numerical issues")
    if (object$settings$monowarp) stop("grad not offered for monowarp = TRUE")
    if (cores > 1 & !requireNamespace("abind", quietly = TRUE))
      stop("Package \"abind\" needed for cores > 1 with grad = TRUE. Please install it.",
           call. = FALSE)
  }

  settings <- list(lite = FALSE, grad = grad, store_latent = FALSE, 
                   mean_map = mean_map, return_all = FALSE, EI = FALSE,
                   entropy_limit = NULL, cores = cores, nper = nper)
  samples <- predict_deep(object, x_new, settings, layers = 2, samples_only = TRUE)
  return(samples) 
}
  
# post_sample.dgp2vec ---------------------------------------------------------
#' @rdname post_sample
#' @export

post_sample.dgp2vec <- function(object, x_new, nper = 1, m = NULL,
                                ord_new = NULL, grad = FALSE, 
                                mean_map = TRUE, cores = 1, ...) {
  
  cores <- check_cores(cores, object$nmcmc)
  if (grad) {
    if (object$v != 999) stop("grad only offered with cov = 'exp2'")
    if (min(object$g) < 1e-6) message("Warning: small g may cause numerical issues")
    if (object$settings$monowarp) stop("grad not offered for monowarp = TRUE")
    if (cores > 1 & !requireNamespace("abind", quietly = TRUE))
      stop("Package \"abind\" needed for cores > 1 with grad = TRUE. Please install it.",
           call. = FALSE)
  }

  settings <- list(ord_new = ord_new, lite = FALSE, grad = grad, 
                   store_latent = FALSE, mean_map = mean_map,
                   return_all = FALSE, EI = FALSE, entropy_limit = NULL, 
                   cores = cores, nper = nper)
  samples <- predict_deep_vec(object, x_new, m, settings, layers = 2, samples_only = TRUE)
  return(samples)
}

# post_sample.dgp3 ----------------------------------------------------------
#' @rdname post_sample
#' @export

post_sample.dgp3 <- function(object, x_new, nper = 1, mean_map = TRUE, 
                             cores = 1, ...) {
  
  cores <- check_cores(cores, object$nmcmc)

  settings <- list(lite = FALSE, grad = FALSE, store_latent = FALSE, 
                   mean_map = mean_map, return_all = FALSE, EI = FALSE,
                   entropy_limit = NULL, cores = cores, nper = nper)
  samples <- predict_deep(object, x_new, settings, layers = 3, samples_only = TRUE)
  return(samples) 
}

# post_sample.dgp3vec ---------------------------------------------------------
#' @rdname post_sample
#' @export

post_sample.dgp3vec <- function(object, x_new, nper = 1, m = NULL, ord_new = NULL, 
                                mean_map = TRUE, cores = 1, ...) {
  
  cores <- check_cores(cores, object$nmcmc)
  
  settings <- list(ord_new = ord_new, lite = FALSE, grad = FALSE, 
                   store_latent = FALSE, mean_map = mean_map,
                   return_all = FALSE, EI = FALSE, entropy_limit = NULL, 
                   cores = cores, nper = nper)
  samples <- predict_deep_vec(object, x_new, m, settings, layers = 3, samples_only = TRUE)
  return(samples)
}
