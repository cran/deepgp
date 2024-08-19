
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   fit_one_layer
#   fit_two_layer
#   fit_three_layer

# Fit One Layer ---------------------------------------------------------------
#' @title MCMC sampling for one layer GP
#' @description Conducts MCMC sampling of hyperparameters for a one layer 
#'     GP.  Length scale parameter \code{theta} governs 
#'     the strength of the correlation and nugget parameter \code{g} 
#'     governs noise.  In Matern covariance, \code{v} governs smoothness.
#'
#' @details Utilizes Metropolis Hastings sampling of the length scale and
#'     nugget parameters with proposals and priors controlled by 
#'     \code{settings}.  When \code{true_g} is set to a specific value, the 
#'     nugget is not estimated.  When \code{vecchia = TRUE}, all calculations 
#'     leverage the Vecchia approximation with specified conditioning set size 
#'     \code{m}.  Vecchia approximation is only implemented for 
#'     \code{cov = "matern"}.
#'     
#'     NOTE on OpenMP: The Vecchia implementation relies on OpenMP parallelization
#'     for efficient computation.  This function will produce a warning message 
#'     if the package was installed without OpenMP (this is the default for 
#'     CRAN packages installed on Apple machines).  To set up OpenMP 
#'     parallelization, download the package source code and install 
#'     using the gcc/g++ compiler.  
#'     
#'     Proposals for \code{g} and \code{theta} follow a uniform sliding window 
#'     scheme, e.g. 
#'     
#'     \code{g_star <- runif(1, l * g_t / u, u * g_t / l)}, 
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{settings}.
#'     To adjust these, set \code{settings = list(l = new_l, u = new_u)}.

#'     Priors on \code{g} and \code{theta} follow Gamma distributions with 
#'     shape parameters (\code{alpha}) and rate parameters (\code{beta}) 
#'     controlled within the \code{settings} list object.  
#'     Defaults have been updated with package version 1.1.3.  Default priors differ
#'     for noisy/deterministic settings and depend on whether \code{monowarp = TRUE}.  
#'     All default values are visible in the internal
#'     \code{deepgp:::check_settings} function.
#'     These priors are designed for \code{x} scaled 
#'     to [0, 1] and \code{y} scaled to have mean 0 and variance 1.  These may
#'     be adjusted using the \code{settings} input.
#'     
#'     The output object of class \code{gp} is designed for use with 
#'     \code{continue}, \code{trim}, and \code{predict}.
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param nmcmc number of MCMC iterations
#' @param sep logical indicating whether to use separable (\code{sep = TRUE})
#'        or isotropic (\code{sep = FALSE}) lengthscales
#' @param verb logical indicating whether to print iteration progress
#' @param g_0 initial value for \code{g}
#' @param theta_0 initial value for \code{theta}
#' @param true_g if true nugget is known it may be specified here (set to a 
#'        small value to make fit deterministic).  Note - values that are too 
#'        small may cause numerical issues in matrix inversions.
#' @param settings hyperparameters for proposals and priors (see details)
#' @param cov covariance kernel, either Matern or squared exponential 
#'        (\code{"exp2"})
#' @param v Matern smoothness parameter (only used if \code{cov = "matern"})
#' @param vecchia logical indicating whether to use Vecchia approximation
#' @param m size of Vecchia conditioning sets (only used if 
#'        \code{vecchia = TRUE})
#' @param ordering optional ordering for Vecchia approximation, must correspond
#'        to rows of \code{x}, defaults to random
#' 
#' @return a list of the S3 class \code{gp} or \code{gpvec} with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{settings}: copy of proposal/prior settings
#'   \item \code{v}: copy of Matern smoothness parameter (\code{v = 999} 
#'         indicates \code{cov = "exp2"})
#'   \item \code{g}: vector of MCMC samples for \code{g}
#'   \item \code{theta}: vector of MCMC samples for \code{theta}
#'   \item \code{tau2}: vector of MLE estimates for \code{tau2} 
#'         (scale parameter)
#'   \item \code{ll}: vector of MVN log likelihood for each Gibbs iteration
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Sauer, A. (2023). Deep Gaussian process surrogates for computer experiments. 
#'     *Ph.D. Dissertation, Department of Statistics, Virginia Polytechnic Institute and State University.*
#'      \cr\cr
#' Sauer, A., Gramacy, R.B., & Higdon, D. (2023). Active learning for deep 
#'     Gaussian process surrogates. *Technometrics, 65,* 4-18.  arXiv:2012.08015
#'     \cr\cr
#' Sauer, A., Cooper, A., & Gramacy, R. B. (2023). Vecchia-approximated deep Gaussian 
#'    processes for computer experiments. 
#'    *Journal of Computational and Graphical Statistics, 32*(3), 824-837.  arXiv:2204.02904
#' 
#' @examples 
#' # Additional examples including real-world computer experiments are available at: 
#' # https://bitbucket.org/gramacylab/deepgp-ex/
#' \donttest{
#' # Booth function (inspired by the Higdon function)
#' f <- function(x) {
#'   i <- which(x <= 0.58)
#'   x[i] <- sin(pi * x[i] * 6) + cos(pi * x[i] * 12)
#'   x[-i] <- 5 * x[-i] - 4.9
#'   return(x)
#' }
#' 
#' # Training data
#' x <- seq(0, 1, length = 25)
#' y <- f(x)
#' 
#' # Testing data
#' xx <- seq(0, 1, length = 200)
#' yy <- f(xx)
#' 
#' plot(xx, yy, type = "l")
#' points(x, y, col = 2)
#' 
#' # Example 1: full model (nugget fixed)
#' fit <- fit_one_layer(x, y, nmcmc = 2000, true_g = 1e-6)
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#' 
#' # Example 2: full model (nugget estimated, EI calculated)
#' fit <- fit_one_layer(x, y, nmcmc = 2000)
#' plot(fit) 
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1, EI = TRUE)
#' plot(fit)
#' par(new = TRUE) # overlay EI
#' plot(xx[order(xx)], fit$EI[order(xx)], type = 'l', lty = 2, 
#' axes = FALSE, xlab = '', ylab = '')
#' 
#' # Example 3: Vecchia approximated model (nugget estimated)
#' fit <- fit_one_layer(x, y, nmcmc = 2000, vecchia = TRUE, m = 10)
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#' }
#' 
#' @export

fit_one_layer <- function(x, y, nmcmc = 10000, sep = FALSE, verb = TRUE, g_0 = 0.001, 
                          theta_0 = 0.1, true_g = NULL, settings = NULL,
                          cov = c("matern", "exp2"), v = 2.5, 
                          vecchia = FALSE, m = min(25, length(y) - 1),
                          ordering = NULL) {
  
  tic <- proc.time()[[3]]
  cov <- match.arg(cov)
  if (vecchia) check_omp()
  if (cov == "exp2") v <- 999 # solely used as an indicator
  if (!vecchia & length(y) > 300) 
    message("'vecchia = TRUE' is recommended for faster computation.")
  if (nmcmc <= 1) stop("nmcmc must be greater than 1")

  # Check inputs
  if (is.numeric(x)) x <- as.matrix(x)
  if (sep & ncol(x) == 1) sep <- FALSE # no need for separable theta in one dimension
  test <- check_inputs(x, y, true_g) # returns NULL if all checks pass
  settings <- check_settings(settings, layers = 1, noisy = is.null(true_g))
  if (vecchia & m >= length(y)) stop("m must be less than the length of y")
  if (cov == "matern")
    if(!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  if (!is.null(ordering)) {
    if (!vecchia) message("ordering only used when vecchia = TRUE")
    test <- check_ordering(ordering, nrow(x)) # returns NULL if all checks pass
  }
  initial <- list(theta = theta_0, g = g_0, tau2 = 1)
  if (sep & (length(initial$theta) == 1)) initial$theta <- rep(initial$theta, ncol(x))
  
  # Create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, v = v)
  if (vecchia) {
    out$m <- m
    out$ordering <- ordering
  }

  # Conduct MCMC
  if (sep) {
    if (vecchia) {
      samples <- gibbs_one_layer_vec_sep(x, y, nmcmc, verb, initial, true_g,
                                         settings, v, m, ordering)
    } else{
      samples <- gibbs_one_layer_sep(x, y, nmcmc, verb, initial, true_g,
                                     settings, v)
    }
  } else {
    if (vecchia) {
      samples <- gibbs_one_layer_vec(x, y, nmcmc, verb, initial, true_g, 
                                     settings, v, m, ordering)
    } else { 
      samples <- gibbs_one_layer(x, y, nmcmc, verb, initial, true_g,
                                 settings, v)
    }
  }
  
  out <- c(out, samples)
  toc <- proc.time()[[3]]
  out$time <- unname(toc - tic)
  if (vecchia) class(out) <- "gpvec" else class(out) <- "gp"
  return(out)
}

# Fit Two Layer ---------------------------------------------------------------
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
#'     \code{settings}.  When \code{true_g} is set to a specific value, the 
#'     nugget is not estimated.  When \code{vecchia = TRUE}, all calculations
#'     leverage the Vecchia approximation with specified conditioning set size
#'     \code{m}.  Vecchia approximation is only implemented for 
#'     \code{cov = "matern"}.
#'   
#'     When \code{monowarp = TRUE}, each input dimension is warped separately and
#'     monotonically.  This requires \code{D = ncol(x)}.  Monotonic warpings trigger
#'     separable lengthscales on the outer layer (\code{theta_y}).  As a default, monotonic 
#'     warpings use the reference grid: \code{seq(0, 1, length = 50)}.  The grid size 
#'     may be controlled by passing a numeric integer to \code{monowarp}
#'     (i.e., \code{monowarp = 100} uses the grid \code{seq(0, 1, length = 100)}).
#'     Alternatively, any user-specified grid may be passed as the argument to 
#'     \code{monowarp}.
#'     
#'     When \code{pmx = TRUE}, the prior on the latent layer is set at \code{x} 
#'     (rather than the default of zero).  This requires \code{D = ncol(x)}.  If
#'     \code{pmx} is set to a numeric value, then that value is used as the scale
#'     parameter on the latent layer.  Specifying a small value here encourages
#'     an identity mapping.
#'     
#'     NOTE on OpenMP: The Vecchia implementation relies on OpenMP parallelization
#'     for efficient computation.  This function will produce a warning message 
#'     if the package was installed without OpenMP (this is the default for 
#'     CRAN packages installed on Apple machines).  To set up OpenMP 
#'     parallelization, download the package source code and install 
#'     using the gcc/g++ compiler.  
#'     
#'     Proposals for \code{g}, \code{theta_y}, and 
#'     \code{theta_w} follow a uniform sliding window scheme, e.g.
#'     
#'     \code{g_star <- runif(1, l * g_t / u, u * g_t / l)}, 
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{settings}.
#'     To adjust these, set \code{settings = list(l = new_l, u = new_u)}.    
#'     Priors on \code{g}, \code{theta_y}, and \code{theta_w} follow Gamma 
#'     distributions with shape parameters (\code{alpha}) and rate parameters 
#'     (\code{beta}) controlled within the \code{settings} list object.  
#'     Defaults have been updated with package version 1.1.3.  Default priors differ
#'     for noisy/deterministic settings and depend on whether \code{monowarp = TRUE}.  
#'     All default values are visible in the internal
#'     \code{deepgp:::check_settings} function.
#'     These priors are designed for \code{x} scaled to 
#'     [0, 1] and \code{y} scaled to have mean 0 and variance 1.  These may be 
#'     adjusted using the \code{settings} input.
#'     
#'     When \code{w_0 = NULL}, the hidden layer is initialized at \code{x} 
#'     (i.e. the identity mapping).  If \code{w_0} is of dimension 
#'     \code{nrow(x) - 1} by \code{D}, the final row is predicted using kriging. 
#'     This is helpful in sequential design when adding a new input location 
#'     and starting the MCMC at the place where the previous MCMC left off.
#'     
#'     The output object of class \code{dgp2} or \code{dgp2vec} is designed for 
#'     use with \code{continue}, \code{trim}, and \code{predict}.   
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param nmcmc number of MCMC iterations
#' @param D integer designating dimension of hidden layer, defaults to 
#'        dimension of \code{x}
#' @param monowarp indicates whether warpings should be forced to be 
#'        monotonic.  Input may be a matrix of
#'        grid points (or a vector which will be applied to every dimension)
#'        for interpolation of the cumulative sum, an integer
#'        specifying the number of grid points to use over the range [0, 1],
#'        or simply the boolean \code{TRUE} which triggers 50 grid points
#'        over the range [0, 1].
#' @param pmx "prior mean x", logical indicating whether \code{w} should have 
#'        prior mean of \code{x} (\code{TRUE}, requires \code{D = ncol(x)}) or prior 
#'        mean zero (\code{FALSE}).  \code{pmx = TRUE} is recommended for
#'        higher dimensions.  May be numeric, in which case the specified
#'        argument is used as the scale (\code{tau2}) in the latent \code{w}
#'        layer (default is 1).  Small values encourage identity mappings.
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
#' @param settings hyperparameters for proposals and priors (see details)
#' @param cov covariance kernel, either Matern or squared exponential 
#'        (\code{"exp2"})
#' @param v Matern smoothness parameter (only used if \code{cov = "matern"})
#' @param vecchia logical indicating whether to use Vecchia approximation
#' @param m size of Vecchia conditioning sets (only used if 
#'        \code{vecchia = TRUE})
#' @param ordering optional ordering for Vecchia approximation, must correspond
#'        to rows of \code{x}, defaults to random, is applied to \code{x}
#'        and \code{w}
#' 
#' @return a list of the S3 class \code{dgp2} or \code{dgp2vec} with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{settings}: copy of proposal/prior settings
#'   \item \code{v}: copy of Matern smoothness parameter (\code{v = 999} 
#'         indicates \code{cov = "exp2"}) 
#'   \item \code{g}: vector of MCMC samples for \code{g}
#'   \item \code{theta_y}: vector of MCMC samples for \code{theta_y} (length
#'         scale of outer layer)
#'   \item \code{theta_w}: matrix of MCMC samples for \code{theta_w} (length 
#'         scale of inner layer)
#'   \item \code{tau2}: vector of MLE estimates for \code{tau2} (scale 
#'         parameter of outer layer)
#'   \item \code{w}: list of MCMC samples for hidden layer \code{w}
#'   \item \code{ll}: vector of MVN log likelihood of the outer layer 
#'         for reach Gibbs iteration
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Sauer, A. (2023). Deep Gaussian process surrogates for computer experiments. 
#'     *Ph.D. Dissertation, Department of Statistics, Virginia Polytechnic Institute and State University.*
#'      \cr\cr
#' Sauer, A., Gramacy, R.B., & Higdon, D. (2023). Active learning for deep 
#'     Gaussian process surrogates. *Technometrics, 65,* 4-18.  arXiv:2012.08015
#'     \cr\cr
#' Sauer, A., Cooper, A., & Gramacy, R. B. (2023). Vecchia-approximated deep Gaussian 
#'      processes for computer experiments. 
#'      *Journal of Computational and Graphical Statistics, 32*(3), 824-837.  arXiv:2204.02904
#'      \cr\cr
#' Barnett, S., Beesley, L. J., Booth, A. S., Gramacy, R. B., & Osthus D. (2024). Monotonic 
#'      warpings for additive and deep Gaussian processes. *In Review.* arXiv:2408.01540
#' 
#' @examples 
#' # Additional examples including real-world computer experiments are available at: 
#' # https://bitbucket.org/gramacylab/deepgp-ex/
#' \donttest{
#' # Booth function (inspired by the Higdon function)
#' f <- function(x) {
#'   i <- which(x <= 0.58)
#'   x[i] <- sin(pi * x[i] * 6) + cos(pi * x[i] * 12)
#'   x[-i] <- 5 * x[-i] - 4.9
#'   return(x)
#' }
#' 
#' # Training data
#' x <- seq(0, 1, length = 25)
#' y <- f(x)
#' 
#' # Testing data
#' xx <- seq(0, 1, length = 200)
#' yy <- f(xx)
#' 
#' plot(xx, yy, type = "l")
#' points(x, y, col = 2)
#' 
#' # Example 1: full model (nugget estimated, using continue)
#' fit <- fit_two_layer(x, y, nmcmc = 1000)
#' plot(fit)
#' fit <- continue(fit, 1000) 
#' plot(fit, hidden = TRUE) # trace plots and ESS samples 
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#' 
#' # Example 2: Vecchia approximated model (nugget estimated)
#' # (Vecchia approximation is faster for larger data sizes)
#' fit <- fit_two_layer(x, y, nmcmc = 2000, vecchia = TRUE, m = 10)
#' plot(fit, hidden = TRUE) # trace plots and ESS samples
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#' 
#' # Example 3: Vecchia approximated model, re-approximated after burn-in 
#' fit <- fit_two_layer(x, y, nmcmc = 1000, vecchia = TRUE, m = 10)
#' fit <- continue(fit, 1000, re_approx = TRUE)
#' plot(fit, hidden = TRUE) # trace plots and ESS samples
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#' 
#' # Example 4: full model with monotonic warpings (nugget estimated)
#' fit <- fit_two_layer(x, y, nmcmc = 2000, monowarp = TRUE)
#' plot(fit, hidden = TRUE) # trace plots and ESS samples
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#' }
#' 
#' @export

fit_two_layer <- function(x, y, nmcmc = 10000, D = ifelse(is.matrix(x), ncol(x), 1), 
                          monowarp = FALSE, pmx = FALSE, 
                          verb = TRUE, w_0 = NULL, g_0 = 0.001,
                          theta_y_0 = 0.1, theta_w_0 = 0.1, true_g = NULL,
                          settings = NULL, cov = c("matern", "exp2"), v = 2.5,
                          vecchia = FALSE, m = min(25, length(y) - 1),
                          ordering = NULL) {

  tic <- proc.time()[[3]]
  cov <- match.arg(cov)
  if (vecchia) check_omp()
  if (cov == "exp2") v <- 999 # solely used as an indicator
  if (!vecchia & length(y) > 300) 
    message("'vecchia = TRUE' is recommended for faster computation.")
  if (nmcmc <= 1) stop("nmcmc must be greater than 1")

  # Check inputs
  if (is.numeric(x)) x <- as.matrix(x)
  test <- check_inputs(x, y, true_g) # returns NULL if all checks pass
  if (vecchia & m >= length(y)) stop("m must be less than the length of y")
  if (cov == "matern")
    if(!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  if (!is.null(ordering)) {
    if (!vecchia) message("ordering only used when vecchia = TRUE")
    test <- check_ordering(ordering, nrow(x)) # returns NULL if all checks pass
  }
  
  # Check monotonic settings
  if (is.numeric(monowarp)) { # anything except FALSE triggers monotonic warpings
    if (length(monowarp) == 1) { # monowarp is integer specifying length of grid
      grid <- seq(0, 1, length = monowarp)
      x_grid <- matrix(grid, nrow = monowarp, ncol = ncol(x), byrow = FALSE)
    } else { # monowarp may be a vector or a matrix
      if (!is.matrix(monowarp)) monowarp <- matrix(monowarp, ncol = 1)
      if (ncol(monowarp) == 1) { # monowarp is the grid for all dimensions
        x_grid <- matrix(sort(monowarp), nrow = nrow(monowarp), ncol = ncol(x), byrow = FALSE)
      } else { # monowarp is a matrix, must check dimensions, make sure ordered
        if (ncol(monowarp) != ncol(x)) stop("dimension of monowarp does not match x")
        x_grid <- monowarp
        for (i in 1:ncol(x_grid)) x_grid[, i] <- sort(x_grid[, i])
      }
    }
    monowarp <- TRUE
  } else if (monowarp) {
      x_grid <- matrix(seq(0, 1, length = 50), nrow = 50, ncol = ncol(x), byrow = FALSE)
  } else {
    x_grid <- NULL
    monowarp <- FALSE
  }
  if (monowarp) {
    for (i in 1:D) {
      if (min(x[, i]) < min(x_grid[, i]) | max(x[, i]) > max(x_grid[, i]))
        message("monowarp grid is narrower than range of x, consider scaling x to [0, 1]
                  or manually specifying the monowarp grid")
      }
  }
  if (monowarp & (ncol(x) != D)) stop("monowarp requires D = ncol(x)")

  # Check settings and initial values
  settings <- check_settings(settings, layers = 2, monowarp, is.null(true_g))
  initial <- list(w = w_0, theta_y = theta_y_0, theta_w = theta_w_0, 
                  g = g_0, tau2 = 1)
  initial <- check_initialization(initial, layers = 2, x = x, D = D, 
                                  vecchia = vecchia, v = v, m = m, 
                                  monowarp = monowarp)
  
  # Check prior mean setting
  if (pmx == 0) pmx <- FALSE
  if (is.numeric(pmx)) {
    settings$inner_tau2 <- pmx
    pmx <- TRUE
  } else {
    settings$inner_tau2 <- 1 # default value for pmx and prior mean zero
  } 
  if (pmx & (ncol(x) != D)) stop("pmx = TRUE requires D = ncol(x)")
  settings$pmx <- pmx
  
  # Create output object
  out <- list(x = x, y = y, x_grid = x_grid, nmcmc = nmcmc, settings = settings, v = v)
  if (vecchia) {
    out$m <- m
    out$ordering <- ordering
  }

  # Conduct MCMC
  if (vecchia) {
    if (monowarp) {
      samples <- gibbs_two_layer_vec_mono(x, y, x_grid, nmcmc, D, verb, initial,
                                          true_g, settings, v, m, ordering)
    } else {
      samples <- gibbs_two_layer_vec(x, y, nmcmc, D, verb, initial,
                                     true_g, settings, v, m, ordering)
    }
  } else { 
    if (monowarp) {
      samples <- gibbs_two_layer_mono(x, y, x_grid, nmcmc, D, verb, initial,
                                      true_g, settings, v)
    } else {
      samples <- gibbs_two_layer(x, y, nmcmc, D, verb, initial,
                                 true_g, settings, v)
    }
  } 
  
  out <- c(out, samples)
  toc <- proc.time()[[3]]
  out$time <- unname(toc - tic)
  if (vecchia) class(out) <- "dgp2vec" else class(out) <- "dgp2"
  return(out)
}

# Fit Three Layer -------------------------------------------------------------
#' @title MCMC sampling for three layer deep GP
#' @description Conducts MCMC sampling of hyperparameters, hidden layer 
#'     \code{z}, and hidden layer \code{w} for a three layer deep GP.  
#'     Separate length scale parameters \code{theta_z}, 
#'     \code{theta_w}, and \code{theta_y} govern the correlation 
#'     strength of the inner layer, middle layer, and outer layer respectively.  
#'     Nugget parameter \code{g} governs noise on the outer layer.  In Matern 
#'     covariance, \code{v} governs smoothness.
#'
#' @details \code{pmx = TRUE} option not yet implemented for three-layer DGP.
#' 
#'     Maps inputs \code{x} through hidden layer \code{z} then hidden
#'     layer \code{w} to outputs \code{y}.  Conducts sampling of the hidden 
#'     layers using Elliptical Slice sampling.  Utilizes Metropolis Hastings 
#'     sampling of the length scale and nugget parameters with proposals and 
#'     priors controlled by \code{settings}.  When \code{true_g} is set to a 
#'     specific value, the nugget is not estimated.  When 
#'     \code{vecchia = TRUE}, all calculations leverage the Vecchia 
#'     approximation with specified conditioning set size \code{m}.  Vecchia 
#'     approximation is only implemented for \code{cov = "matern"}.
#'     
#'     NOTE on OpenMP: The Vecchia implementation relies on OpenMP parallelization
#'     for efficient computation.  This function will produce a warning message 
#'     if the package was installed without OpenMP (this is the default for 
#'     CRAN packages installed on Apple machines).  To set up OpenMP 
#'     parallelization, download the package source code and install 
#'     using the gcc/g++ compiler.  
#'     
#'     Proposals for \code{g}, 
#'     \code{theta_y}, \code{theta_w}, and \code{theta_z} follow a uniform 
#'     sliding window scheme, e.g.
#'     
#'     \code{g_star <- runif(1, l * g_t / u, u * g_t / l)},
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{settings}.
#'     To adjust these, set \code{settings = list(l = new_l, u = new_u)}.  
#'     Priors on \code{g}, \code{theta_y}, \code{theta_w}, and \code{theta_z} 
#'     follow Gamma distributions with shape parameters (\code{alpha}) and rate 
#'     parameters (\code{beta}) controlled within the \code{settings} list 
#'     object.  Defaults have been updated with package version 1.1.3.  
#'     Default priors differ for noisy/deterministic settings and 
#'     depend on whether \code{monowarp = TRUE}.  All default values are 
#'     visible in the internal \code{deepgp:::check_settings} function.
#'     These priors are designed for \code{x} scaled to [0, 1] and \code{y} 
#'     scaled to have mean 0 and variance 1.  These may be adjusted using the 
#'     \code{settings} input.
#'     
#'     In the current version, the three-layer does not have any equivalent
#'     setting for \code{monowarp = TRUE} or \code{pmx = TRUE} as in 
#'     \code{fit_two_layer}.
#'     
#'     When \code{w_0 = NULL} and/or \code{z_0 = NULL}, the hidden layers are 
#'     initialized at \code{x} (i.e. the identity mapping).  The default prior 
#'     mean of the inner hidden layer \code{z} is zero, but may be adjusted to \code{x} 
#'     using \code{settings = list(z_prior_mean = x)}.  The prior mean of the
#'     middle hidden layer \code{w} is set at zero is is not user adjustable.
#'     If \code{w_0} and/or \code{z_0} is of dimension \code{nrow(x) - 1} by 
#'     \code{D}, the final row is predicted using kriging. This is helpful in 
#'     sequential design when adding a new input location and starting the MCMC 
#'     at the place where the previous MCMC left off.
#'     
#'     The output object of class \code{dgp3} or \code{dgp3vec} is designed for 
#'     use with \code{continue}, \code{trim}, and \code{predict}. 
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param nmcmc number of MCMC iterations
#' @param D integer designating dimension of hidden layers, defaults to 
#'        dimension of \code{x}
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
#' @param settings hyperparameters for proposals and priors (see details)
#' @param cov covariance kernel, either Matern or squared exponential 
#'        (\code{"exp2"})
#' @param v Matern smoothness parameter (only used if \code{cov = "matern"})
#' @param vecchia logical indicating whether to use Vecchia approximation
#' @param m size of Vecchia conditioning sets (only used if 
#'        \code{vecchia = TRUE})
#' @param ordering optional ordering for Vecchia approximation, must correspond
#'        to rows of \code{x}, defaults to random, is applied to \code{x},
#'        \code{w}, and \code{z}
#'        
#' @return a list of the S3 class \code{dgp3} or \code{dgp3vec} with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{settings}: copy of proposal/prior settings
#'   \item \code{v}: copy of Matern smoothness parameter (\code{v = 999} 
#'         indicates \code{cov = "exp2"}) 
#'   \item \code{g}: vector of MCMC samples for \code{g}
#'   \item \code{theta_y}: vector of MCMC samples for \code{theta_y} (length 
#'         scale of outer layer)
#'   \item \code{theta_w}: matrix of MCMC samples for \code{theta_w} (length 
#'         scale of middle layer)
#'   \item \code{theta_z}: matrix of MCMC samples for \code{theta_z} (length 
#'         scale of inner layer)
#'   \item \code{tau2}: vector of MLE estimates for \code{tau2} (scale 
#'         parameter of outer layer)
#'   \item \code{w}: list of MCMC samples for middle hidden layer \code{w}
#'   \item \code{z}: list of MCMC samples for inner hidden layer \code{z}
#'   \item \code{ll}: vector of MVN log likelihood of the outer layer 
#'         for reach Gibbs iteration
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Sauer, A. (2023). Deep Gaussian process surrogates for computer experiments. 
#'     *Ph.D. Dissertation, Department of Statistics, Virginia Polytechnic Institute and State University.*
#'      \cr\cr
#' Sauer, A., Gramacy, R.B., & Higdon, D. (2023). Active learning for deep 
#'     Gaussian process surrogates. *Technometrics, 65,* 4-18.  arXiv:2012.08015
#'     \cr\cr
#' Sauer, A., Cooper, A., & Gramacy, R. B. (2023). Vecchia-approximated deep Gaussian 
#'      processes for computer experiments. 
#'      *Journal of Computational and Graphical Statistics, 32*(3), 824-837.  arXiv:2204.02904
#' 
#' @examples 
#' # Additional examples including real-world computer experiments are available at: 
#' # https://bitbucket.org/gramacylab/deepgp-ex/
#' \donttest{
#' # G function in 2 dimensions (https://www.sfu.ca/~ssurjano/gfunc.html)
#' f <- function(xx, a = (c(1:length(xx)) - 1) / 2) { 
#'   new1 <- abs(4 * xx - 2) + a
#'   new2 <- 1 + a
#'   prod <- prod(new1 / new2)
#'   return((prod - 1) / 0.86)
#' }
#' 
#' # Training data
#' d <- 2
#' n <- 30
#' x <- matrix(runif(n * d), ncol = d)
#' y <- apply(x, 1, f)
#' 
#' # Testing data
#' n_test <- 500
#' xx <- matrix(runif(n_test * d), ncol = d)
#' yy <- apply(xx, 1, f)
#' 
#' i <- interp::interp(xx[, 1], xx[, 2], yy)
#' image(i, col = heat.colors(128))
#' contour(i, add = TRUE)
#' contour(i, level = -0.5, col = 4, add = TRUE) # potential failure limit
#' points(x)
#' 
#' # Example 1: full model (nugget estimated, entropy calculated)
#' fit <- fit_three_layer(x, y, nmcmc = 2000)
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, entropy_limit = -0.5, cores = 1)
#' plot(fit)
#' i <- interp::interp(xx[, 1], xx[, 2], fit$entropy)
#' image(i, col = heat.colors(128), main = "Entropy")
#' 
#' # Example 2: Vecchia approximated model (nugget fixed)
#' # (Vecchia approximation is faster for larger data sizes)
#' fit <- fit_three_layer(x, y, nmcmc = 2000, vecchia = TRUE, 
#'                        m = 10, true_g = 1e-6)
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#' }
#' 
#' @export

fit_three_layer <- function(x, y, nmcmc = 10000, D = ifelse(is.matrix(x), ncol(x), 1), 
                            verb = TRUE, w_0 = NULL, z_0 = NULL,
                            g_0 = 0.001, theta_y_0 = 0.1, theta_w_0 = 0.1, 
                            theta_z_0 = 0.1, true_g = NULL, settings = NULL,
                            cov = c("matern", "exp2"), v = 2.5, vecchia = FALSE, 
                            m = min(25, length(y) - 1), ordering = NULL) {
  
  tic <- proc.time()[[3]]
  cov <- match.arg(cov)
  if (vecchia) check_omp()
  if (cov == "exp2") v <- 999 # solely used as an indicator
  if (!vecchia & length(y) > 300) 
    message("'vecchia = TRUE' is recommended for faster computation.")
  if (nmcmc <= 1) stop("nmcmc must be greater than 1")

  # Check inputs
  if (is.numeric(x)) x <- as.matrix(x)
  test <- check_inputs(x, y, true_g)
  settings <- check_settings(settings, layers = 3, noisy = is.null(true_g))
  if (vecchia & m >= length(y)) stop("m must be less than the length of y")
  if (cov == "matern")
    if(!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  if (!is.null(ordering)) {
    if (!vecchia) message("ordering only used when vecchia = TRUE")
    test <- check_ordering(ordering, nrow(x)) # returns NULL if all checks pass
  }
  initial <- list(w = w_0, z = z_0, theta_y = theta_y_0, theta_w = theta_w_0, 
                  theta_z = theta_z_0, g = g_0, tau2 = 1)
  initial <- check_initialization(initial, layers = 3, x = x, D = D, 
                                  vecchia = vecchia, v = v, m = m)
  
  # Create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, v = v)
  if (vecchia) {
    out$m <- m
    out$ordering <- ordering
  }
  
  # Conduct MCMC
  if (vecchia) {
    samples <- gibbs_three_layer_vec(x, y, nmcmc, D, verb, initial, true_g,
                                     settings, v, m, ordering)
  } else {
    samples <- gibbs_three_layer(x, y, nmcmc, D, verb, initial, true_g, 
                                 settings, v)
  }
  
  out <- c(out, samples)
  toc <- proc.time()[[3]]
  out$time <- unname(toc - tic)
  if (vecchia) class(out) <- "dgp3vec" else class(out) <- "dgp3"
  return(out)
}
