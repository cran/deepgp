
# Function Contents -----------------------------------------------------------
# External: (see documentation below)
#   fit_one_layer
#   fit_two_layer
#   fit_three_layer

# fit_one_layer ---------------------------------------------------------------
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
#'     \code{m}.
#'     
#'     NOTE on OpenMP: The Vecchia implementation relies on OpenMP parallelization
#'     for efficient computation.  This function will produce a warning message 
#'     if the package was installed without OpenMP (this is the default for 
#'     CRAN packages installed on Apple machines).  To set up OpenMP 
#'     parallelization, download the package source code and install 
#'     using the gcc/g++ compiler.  
#'     
#'     Proposals for \code{g} and \code{theta} follow a uniform sliding window 
#'     scheme, e.g.,
#'     
#'     \code{g_star <- runif(1, l * g_t / u, u * g_t / l)}, 
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{settings}.
#'     To adjust these, set \code{settings = list(l = new_l, u = new_u)}.
#'
#'     Priors on \code{g} and \code{theta} follow Gamma distributions with 
#'     shape parameters (\code{alpha}) and rate parameters (\code{beta}) 
#'     controlled within the \code{settings} list object.  
#'     Default priors differ for noisy/deterministic settings.  
#'     All default values are visible in the internal
#'     \code{deepgp:::check_settings} function.
#'     These priors are designed for \code{x} scaled 
#'     to [0, 1] and \code{y} scaled to have mean 0 and variance 1.  These may
#'     be adjusted using the \code{settings} input.
#'     
#'     The output object of class \code{gp} is designed for use with
#'     \code{continue}, \code{trim}, \code{plot}, and \code{predict}.
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param dydx optional matrix of observed gradients, rows correspond to
#'        \code{x} locations, columns contain partial derivatives with 
#'        respect to that input dimension (\code{dim(dy)} must match \code{dim(x)})
#' @param nmcmc number of MCMC iterations
#' @param sep logical indicating whether to use separable (\code{sep = TRUE})
#'        or isotropic (\code{sep = FALSE}) lengthscales
#' @param verb logical indicating whether to print iteration progress
#' @param theta_0 initial value for \code{theta}
#' @param g_0 initial value for \code{g} (only used if \code{true_g = NULL})
#' @param true_g if true nugget is known it may be specified here (set to a 
#'        small value to make fit deterministic).  Note - values that are too 
#'        small may cause numerical issues in matrix inversions.
#' @param v Matern smoothness parameter (only used if \code{cov = "matern"})
#' @param settings hyperparameters for proposals and priors (see details)
#' @param cov covariance kernel, either Matern (\code{"matern"}) or squared 
#'        exponential (\code{"exp2"})
#' @param vecchia logical indicating whether to use Vecchia approximation
#' @param m size of Vecchia conditioning sets, defaults to the lower of 25 or
#'        the maximum available (only used if \code{vecchia = TRUE})
#' @param ord optional ordering for Vecchia approximation, must correspond
#'        to rows of \code{x}, defaults to random
#' @param cores number of cores to use for OpenMP parallelization 
#'        (\code{vecchia = TRUE} only).  Defaults to \code{min(4, maxcores - 1)} 
#'        where \code{maxcores} is the number of detectable available cores.
#' 
#' @return a list of the S3 class \code{gp} or \code{gpvec} with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{settings}: copy of proposal/prior settings
#'   \item \code{v}: copy of Matern smoothness parameter (\code{v = 999} 
#'         indicates \code{cov = "exp2"})
#'   \item \code{dydx}: copy of dydx (if not NULL)
#'   \item \code{grad_indx}: stacked partial derivative indices (only if \code{dydx} is provided)
#'   \item \code{g}: vector of MCMC samples for \code{g}
#'   \item \code{theta}: vector of MCMC samples for \code{theta}
#'   \item \code{tau2}: vector of MLE estimates for \code{tau2} (scale parameter)
#'   \item \code{x_approx}: Vecchia approximation object (\code{vecchia = TRUE} only)
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
#' Booth, A. S. (2025). Deep Gaussian processes with gradients. arXiv:2512.18066
#'      \cr\cr
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
#' xx <- seq(0, 1, length = 100)
#' yy <- f(xx)
#' 
#' plot(xx, yy, type = "l")
#' points(x, y, col = 2)
#'
#' # Example 1: nugget fixed, calculating EI
#' fit <- fit_one_layer(x, y, nmcmc = 2000, true_g = 1e-6)
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1, EI = TRUE)
#' plot(fit)
#' par(new = TRUE) # overlay EI
#' plot(xx[order(xx)], fit$EI[order(xx)], type = 'l', lty = 2, 
#'    axes = FALSE, xlab = '', ylab = '')
#'
#' # Example 2: convert fit to Vecchia object before predicting
#' # (this is faster if the training data set is large)
#' fit <- to_vec(fit)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#' 
#' # Example 3: using Vecchia for training and testing
#' fit <- fit_one_layer(x, y, nmcmc = 2000, true_g = 1e-6, vecchia = TRUE, m = 10)
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#' }
#' 
#' @export

fit_one_layer <- function(x, y, dydx = NULL, nmcmc = 10000, sep = FALSE, verb = TRUE, 
                          theta_0 = 0.01, g_0 = 0.001, true_g = NULL, v = 2.5,
                          settings = NULL, cov = c("matern", "exp2"),
                          vecchia = FALSE, m = NULL, ord = NULL, cores = NULL) {
  
  tic <- proc.time()[[3]]
  cov <- match.arg(cov)
  if (is.vector(x)) x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  if (sep & d == 1) sep <- FALSE # no need for separable theta in one dimension
  
  # Check inputs and settings
  test <- check_inputs(x, y, true_g, nmcmc) # returns NULL if all checks pass
  settings <- check_settings(settings, layers = 1, noisy = is.null(true_g))
  settings$sep <- sep
  
  # Check gradients
  if (!is.null(dydx)) {
    grad_enhance <- TRUE
    if (is.vector(dydx)) dydx <- as.matrix(dydx) # one dimension only
    test <- check_gradients(n, d, dydx, cov, true_g, vecchia = vecchia) # returns NULL
  } else grad_enhance <- FALSE
  
  # Check covariance
  if (cov == "matern") {
    if(!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  } else if (cov == "exp2") {
    v <- 999 # indicator for "exp2" kernel
  } else stop("cov must be 'matern' or 'exp2'")
  
  # Check vecchia 
  if (vecchia) {
    cores <- check_cores(cores)
    if (is.null(m)) m <- min(25, ifel(grad_enhance, n*(d + 1) - 1, n - 1))
    test <- check_vecchia(n, d, m, ord, grad_enhance) 
  } else {
    if (n > 200) message("'vecchia = TRUE' is recommended for faster computation.")
    if (!is.null(cores)) message("cores is only used when 'vecchia = TRUE'")
  }
  
  # Create initial list
  if (sep & (length(theta_0) == 1)) theta_0 <- rep(theta_0, d)
  initial <- list(theta = theta_0, g = g_0)
  
  # Create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, v = v)
  if (grad_enhance) {
    out$dydx <- dydx
    out$grad_indx <- rep(0:d, each = n)
  }

  # Conduct MCMC
  if (vecchia) {
    samples <- gibbs_one_layer_vec(x, y, dydx, nmcmc, verb, initial, true_g,
                                   settings, v, m, ord, cores)
  } else {
    samples <- gibbs_one_layer(x, y, dydx, nmcmc, verb, initial, true_g,
                               settings, v)
  }
  
  out <- c(out, samples)
  toc <- proc.time()[[3]]
  out$time <- unname(toc - tic)
  if (vecchia) class(out) <- "gpvec" else class(out) <- "gp"
  return(out)
}

# fit_two_layer ---------------------------------------------------------------
#' @title MCMC sampling for two layer deep GP
#' @description Conducts MCMC sampling of hyperparameters and hidden layer 
#'     \code{w} for a two layer deep GP.  Separate length scale 
#'     parameters \code{theta_w} and \code{theta_y} govern the correlation 
#'     strength of the hidden layer and outer layer respectively.  Nugget 
#'     parameter \code{g} governs noise on the outer layer.  In Matern 
#'     covariance, \code{v} governs smoothness.
#'
#' @details Maps inputs \code{x} through hidden layer \code{w} to outputs 
#'     \code{y}.  Conducts sampling of the hidden layer using elliptical 
#'     slice sampling.  Utilizes Metropolis Hastings sampling of the length 
#'     scale and nugget parameters with proposals and priors controlled by 
#'     \code{settings}.  When \code{true_g} is set to a specific value, the 
#'     nugget is not estimated.  When \code{vecchia = TRUE}, all calculations
#'     leverage the Vecchia approximation with specified conditioning set size
#'     \code{m}.
#'   
#'     When \code{monowarp = TRUE}, each input dimension is warped separately and
#'     monotonically.  This requires \code{D = ncol(x)} with \code{x} scaled to the 
#'     unit cube.  New in version 1.2.0 - monotonic warpings estimate separate
#'     scale parameters (\code{tau2_w}) on each latent node and use an isotropic 
#'     lengthscale on the outer layer.  As a default, monotonic 
#'     warpings use the reference grid: \code{seq(0, 1, length = 50)}.  The grid size 
#'     may be controlled by passing a numeric integer to \code{monowarp}
#'     (i.e., \code{monowarp = 100} uses the grid \code{seq(0, 1, length = 100)}).
#'     
#'     When \code{pmx = TRUE}, the prior on the latent layer is set at \code{x} 
#'     (rather than the default of zero).  This requires \code{D = ncol(x)}.
#'     
#'     NOTE on OpenMP: The Vecchia implementation relies on OpenMP parallelization
#'     for efficient computation.  This function will produce a warning message 
#'     if the package was installed without OpenMP (this is the default for 
#'     CRAN packages installed on Apple machines).  To set up OpenMP 
#'     parallelization, download the package source code and install 
#'     using the gcc/g++ compiler.  
#'     
#'     Proposals for \code{g}, \code{theta_y}, and 
#'     \code{theta_w} follow a uniform sliding window scheme, e.g.,
#'     
#'     \code{g_star <- runif(1, l * g_t / u, u * g_t / l)}, 
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{settings}.
#'     To adjust these, set \code{settings = list(l = new_l, u = new_u)}.    
#'     Priors on \code{g}, \code{theta_y}, and \code{theta_w} follow Gamma 
#'     distributions with shape parameters (\code{alpha}) and rate parameters 
#'     (\code{beta}) controlled within the \code{settings} list object.  
#'     Default priors differ for noisy/deterministic settings and depend on 
#'     whether \code{monowarp = TRUE}.  
#'     All default values are visible in the internal
#'     \code{deepgp:::check_settings} function.
#'     These priors are designed for \code{x} scaled to 
#'     [0, 1] and \code{y} scaled to have mean 0 and variance 1.  These may be 
#'     adjusted using the \code{settings} input.
#'
#'     The scale on the latent layer (\code{tau2_w}) may also be specified in 
#'     \code{settings}.  Defaults to 1. 
#'     
#'     When \code{w_0 = NULL}, the hidden layer is initialized at \code{x} 
#'     (i.e., the identity mapping).  If \code{w_0} is of dimension 
#'     \code{nrow(x) - 1} by \code{D}, the final row is filled-in using the GP
#'     posterior mean. 
#'     This is helpful in sequential design when adding a new input location 
#'     and starting the MCMC at the place where the previous MCMC left off.
#'     
#'     The output object of class \code{dgp2} or \code{dgp2vec} is designed for 
#'     use with \code{continue}, \code{trim}, and \code{predict}.   
#'
#' @param x vector or matrix of input locations
#' @param y vector of response values
#' @param dydx optional matrix of observed gradients, rows correspond to
#'        \code{x} locations, columns contain partial derivatives with 
#'        respect to that input dimension (\code{dim(dy)} must match \code{dim(x)})
#' @param nmcmc number of MCMC iterations
#' @param D integer designating dimension of hidden layer, defaults to 
#'        dimension of \code{x}
#' @param monowarp logical or numeric.  If \code{FALSE}, warpings are not forced to
#'        be monotonic.  If \code{TRUE}, each input dimension is individually monotonically
#'        warped with a default grid size of 50.  If numeric, triggers monotonic 
#'        warpings with the provided grid size.
#' @param pmx "prior mean x", logical indicating whether \code{w} should have 
#'        prior mean of \code{x} (\code{TRUE}, requires \code{D = ncol(x)}) or prior 
#'        mean zero (\code{FALSE}).  \code{pmx = TRUE} is recommended for
#'        higher dimensions.
#' @param verb logical indicating whether to print iteration progress
#' @param w_0 initial value for hidden layer \code{w} (rows must correspond to
#'        rows of \code{x}, requires \code{ncol(w_0) = D}.  Defaults to the 
#'        identity mapping.  If \code{nrow(w_0) < nrow(x)}, missing initial values
#'        are filled-in with the GP posterior mean.
#' @param theta_y_0 initial value for \code{theta_y} (length scale of outer 
#'        layer)
#' @param theta_w_0 initial value for \code{theta_w} (length scale of inner 
#'        layer), may be single value or vector of length \code{D}
#' @param g_0 initial value for \code{g} (only used if \code{true_g = NULL})
#' @param true_g if true nugget is known it may be specified here (set to a 
#'        small value to make fit deterministic).  Note - values that are too 
#'        small may cause numerical issues in matrix inversions.
#' @param v Matern smoothness parameter (only used if \code{cov = "matern"})
#' @param settings hyperparameters for proposals and priors (see details)
#' @param cov covariance kernel, either Matern (\code{"matern"}) or squared 
#'        exponential (\code{"exp2"})
#' @param vecchia logical indicating whether to use Vecchia approximation
#' @param m size of Vecchia conditioning sets, defaults to the lower of 25 or
#'        the maximum available (only used if \code{vecchia = TRUE})
#' @param ord optional ordering for Vecchia approximation, must correspond
#'        to rows of \code{x}, defaults to random, is applied to both \code{x}
#'        and \code{w}
#' @param cores number of cores to use for OpenMP parallelization 
#'        (\code{vecchia = TRUE} only).  Defaults to \code{min(4, maxcores - 1)} 
#'        where \code{maxcores} is the number of detectable available cores.
#' 
#' @return a list of the S3 class \code{dgp2} or \code{dgp2vec} with elements:
#' \itemize{
#'   \item \code{x}: copy of input matrix
#'   \item \code{y}: copy of response vector
#'   \item \code{nmcmc}: number of MCMC iterations
#'   \item \code{settings}: copy of proposal/prior settings
#'   \item \code{v}: copy of Matern smoothness parameter (\code{v = 999} 
#'         indicates \code{cov = "exp2"}) 
#'   \item \code{x_grid}: grid used for monotonic warpings (\code{monowarp = TRUE} only)
#'   \item \code{dydx}: copy of dydx (if not NULL)
#'   \item \code{grad_indx}: stacked partial derivative indices (only if \code{dydx} is provided)
#'   \item \code{g}: vector of MCMC samples for \code{g}
#'   \item \code{tau2_y}: vector of MLE estimates for \code{tau2} on the outer layer
#'   \item \code{theta_y}: vector of MCMC samples for \code{theta_y} (length
#'         scale of outer layer)
#'   \item \code{tau2_w}: matrix of MLE estimates for \code{tau2} on inner layer
#'         (only returned if \code{monowarp = TRUE}, otherwise this is fixed in \code{settings})
#'   \item \code{theta_w}: matrix of MCMC samples for \code{theta_w} (length 
#'         scale of inner layer)
#'   \item \code{w}: list of MCMC samples for hidden layer \code{w}
#'   \item \code{w_grid}: \code{w} values at \code{x_grid} locations (\code{monowarp = TRUE} only)
#'   \item \code{w_approx}: Vecchia approximation object for outer layer (\code{vecchia = TRUE} only)
#'   \item \code{x_approx}: Vecchia approximation object for inner layer (\code{vecchia = TRUE} only)
#'   \item \code{ll}: vector of MVN log likelihood of the outer layer 
#'         for reach Gibbs iteration
#'   \item \code{time}: computation time in seconds
#' }
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
#' Barnett, S., Beesley, L. J., Booth, A. S., Gramacy, R. B., & Osthus D. (2025). 
#'     Monotonic warpings for additive and deep Gaussian processes. 
#'     *Statistics and Computing, 35*(3), 65. arXiv:2408.01540
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
#' xx <- seq(0, 1, length = 100)
#' yy <- f(xx)
#' 
#' plot(xx, yy, type = "l")
#' points(x, y, col = 2)
#' 
#' # Example 1: nugget fixed, using continue
#' fit <- fit_two_layer(x, y, nmcmc = 1000, true_g = 1e-6)
#' plot(fit)
#' fit <- continue(fit, 1000) 
#' plot(fit, hidden = TRUE) # trace plots and ESS samples 
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#'
#' # Example 2: using Vecchia, re-approximated after burn-in 
#' fit <- fit_two_layer(x, y, nmcmc = 1000, true_g = 1e-6, vecchia = TRUE, m = 10)
#' fit <- continue(fit, 1000, re_approx = TRUE)
#' plot(fit, hidden = TRUE) # trace plots and ESS samples
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#'
#' # Example 3: using monotonic warpings
#' fit <- fit_two_layer(x, y, nmcmc = 2000, true_g = 1e-6, monowarp = TRUE)
#' plot(fit, hidden = TRUE) # trace plots and ESS samples
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#' }
#'
#' @export

fit_two_layer <- function(x, y, dydx = NULL, nmcmc = 10000, 
                          D = ifelse(is.matrix(x), ncol(x), 1), 
                          monowarp = FALSE, pmx = FALSE, 
                          verb = TRUE, w_0 = NULL,
                          theta_y_0 = 0.01, theta_w_0 = 0.1, g_0 = 0.001, 
                          true_g = NULL, v = 2.5,
                          settings = NULL, cov = c("matern", "exp2"), 
                          vecchia = FALSE, m = NULL, ord = NULL, cores = NULL) {

  tic <- proc.time()[[3]]
  cov <- match.arg(cov)
  if (is.vector(x)) x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)

  # Check gradients
  if (!is.null(dydx)) {
    grad_enhance <- TRUE
    if (is.vector(dydx)) dydx <- as.matrix(dydx) # one dimension only
    test <- check_gradients(n, d, dydx, cov, true_g, D, vecchia, monowarp) # returns NULL
  } else grad_enhance <- FALSE
  
  # Check inputs and settings
  test <- check_inputs(x, y, true_g, nmcmc) # returns NULL
  settings <- check_settings(settings, layers = 2, noisy = is.null(true_g),
                             monowarp = (monowarp | is.numeric(monowarp)))
  
  # Check covariance
  if (cov == "matern") {
    if (!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  } else if (cov == "exp2") {
    v <- 999 # indicator for "exp2" kernel
  } else stop("cov must be 'matern' or 'exp2'")
  
  # Check vecchia
  if (vecchia) {
    cores <- check_cores(cores)
    if (is.null(m)) m <- min(25, ifel(grad_enhance, n*(d + 1) - 1, n - 1))
    test <- check_vecchia(n, d, m, ord, grad_enhance)
  } else {
    if (n > 200) message("'vecchia = TRUE' is recommended for faster computation.")
    if (!is.null(cores)) message("cores is only used when 'vecchia = TRUE'")
  }

  # Check monowarp 
  if (is.numeric(monowarp)) {
    ng <- monowarp
    if (ng > 200) message(paste0("Warning: Vecchia approximation is not implemented for monowarp",
                                  " grid, we recommend a smaller value"))
    monowarp <- TRUE
  } else ng <- 50
  if (monowarp) {
    if (min(x) < 0 | max(x) > 1) stop("monowarp requires x be scaled to [0, 1]^d")
    if (d != D) stop("monowarp = TRUE requires D = ncol(x)")
    if (!is.null(w_0)) message("monowarp overwrites w_0 with x")
    x_grid <- seq(0, 1, length = ng)
  } 
  settings$monowarp <- monowarp

  # Check pmx
  if (pmx & (d != D)) stop("pmx = TRUE requires D = ncol(x)")
  settings$pmx <- pmx
  if (pmx & grad_enhance) settings$w_prior_mean <- get_prior_mean(x)

  # Create and check initial list
  initial <- list(w = w_0, theta_y = theta_y_0, theta_w = theta_w_0, g = g_0)
  initial <- check_initialization(initial, layers = 2, D = D, grad_enhance = grad_enhance, 
                                  x = x, v = v, pmx = pmx, vecchia = vecchia, m = m)
  
  # Create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, v = v)
  if (monowarp) out$x_grid <- x_grid
  if (grad_enhance) {
    out$dydx <- dydx
    out$grad_indx <- rep(0:d, each = n) 
  }

  # Conduct MCMC
  if (vecchia) {
    if (monowarp) {
      samples <- gibbs_two_layer_vec_mono(x, y, x_grid, nmcmc, verb, initial, true_g,
                                          settings, v, m, ord, cores)
    } else if (grad_enhance) {
      samples <- gibbs_two_layer_vec_grad(x, y, dydx, nmcmc, verb, initial,
                                          true_g, settings, v, m, ord, cores)
    } else {
      samples <- gibbs_two_layer_vec(x, y, nmcmc, D, verb, initial,
                                     true_g, settings, v, m, ord, cores)
    }
  } else { 
    if (monowarp) {
      samples <- gibbs_two_layer_mono(x, y, x_grid, nmcmc, verb, initial, true_g, 
                                      settings, v)
    } else if (grad_enhance) {
      samples <- gibbs_two_layer_grad(x, y, dydx, nmcmc, verb, initial, true_g, 
                                      settings, v)
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

# fit_three_layer -------------------------------------------------------------
#' @title MCMC sampling for three layer deep GP
#' @description Conducts MCMC sampling of hyperparameters, hidden layer 
#'     \code{z}, and hidden layer \code{w} for a three layer deep GP.  
#'     Separate length scale parameters \code{theta_z}, 
#'     \code{theta_w}, and \code{theta_y} govern the correlation 
#'     strength of the inner layer, middle layer, and outer layer respectively.  
#'     Nugget parameter \code{g} governs noise on the outer layer.  In Matern 
#'     covariance, \code{v} governs smoothness.
#'
#'     Currently, there are no \code{pmx}, \code{monowarp}, or \code{dydx} 
#'     options.  
#'
#' @details Maps inputs \code{x} through hidden layer \code{z} then hidden
#'     layer \code{w} to outputs \code{y}.  Conducts sampling of the hidden 
#'     layers using elliptical slice sampling.  Utilizes Metropolis Hastings 
#'     sampling of the length scale and nugget parameters with proposals and 
#'     priors controlled by \code{settings}.  When \code{true_g} is set to a 
#'     specific value, the nugget is not estimated.  When 
#'     \code{vecchia = TRUE}, all calculations leverage the Vecchia 
#'     approximation with specified conditioning set size \code{m}.
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
#'     sliding window scheme, e.g.,
#'     
#'     \code{g_star <- runif(1, l * g_t / u, u * g_t / l)},
#'     
#'     with defaults \code{l = 1} and \code{u = 2} provided in \code{settings}.
#'     To adjust these, set \code{settings = list(l = new_l, u = new_u)}.  
#'     Priors on \code{g}, \code{theta_y}, \code{theta_w}, and \code{theta_z} 
#'     follow Gamma distributions with shape parameters (\code{alpha}) and rate 
#'     parameters (\code{beta}) controlled within the \code{settings} list 
#'     object.  Default priors differ for noisy/deterministic settings.  
#'     All default values are 
#'     visible in the internal \code{deepgp:::check_settings} function.
#'     These priors are designed for \code{x} scaled to [0, 1] and \code{y} 
#'     scaled to have mean 0 and variance 1.  These may be adjusted using the 
#'     \code{settings} input.
#'
#'     The scale on the latent layers (\code{tau2_z} and \code{tau2_w}) may also 
#'     be specified in \code{settings}.  Defaults to 1. 
#'     
#'     When \code{w_0 = NULL} and/or \code{z_0 = NULL}, the hidden layers are 
#'     initialized at \code{x} (i.e., the identity mapping).
#'     If \code{w_0} and/or \code{z_0} is of dimension \code{nrow(x) - 1} by 
#'     \code{D}, the final row is filled-in using the GP posterior mean. This is helpful in 
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
#' @param w_0 initial value for hidden layer \code{w} (rows must correspond to
#'        rows of \code{x}, requires \code{ncol(w_0) = D}.  Defaults to the 
#'        identity mapping.  If \code{nrow(w_0) < nrow(x)}, missing initial values
#'        are filled-in with the GP posterior mean.
#' @param z_0 initial value for hidden layer \code{z} (rows must correspond to
#'        rows of \code{x}, requires \code{ncol(z_0) = D}.  Defaults to the 
#'        identity mapping.  If \code{nrow(z_0) < nrow(x)}, missing initial values
#'        are filled-in with the GP posterior mean.
#' @param theta_y_0 initial value for \code{theta_y} (length scale of outer 
#'        layer)
#' @param theta_w_0 initial value for \code{theta_w} (length scale of middle 
#'        layer), may be single value or vector of length \code{D}
#' @param theta_z_0 initial value for \code{theta_z} (length scale of inner 
#'        layer), may be single value or vector of length \code{D}
#' @param g_0 initial value for \code{g}
#' @param true_g if true nugget is known it may be specified here (set to a 
#'        small value to make fit deterministic).  Note - values that are too 
#'        small may cause numerical issues in matrix inversions.
#' @param v Matern smoothness parameter (only used if \code{cov = "matern"})
#' @param settings hyperparameters for proposals and priors (see details)
#' @param cov covariance kernel, either Matern (\code{"matern"}) or squared exponential 
#'        (\code{"exp2"})
#' @param vecchia logical indicating whether to use Vecchia approximation
#' @param m size of Vecchia conditioning sets, defaults to the lower of 25 or
#'        the maximum available (only used if \code{vecchia = TRUE})
#' @param ord optional ordering for Vecchia approximation, must correspond
#'        to rows of \code{x}, defaults to random, is applied to \code{x},
#'        \code{w}, and \code{z}
#' @param cores number of cores to use for OpenMP parallelization 
#'        (\code{vecchia = TRUE} only).  Defaults to \code{min(4, maxcores - 1)} 
#'        where \code{maxcores} is the number of detectable available cores.
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
#'   \item \code{tau2_y}: vector of MLE estimates for \code{tau2} on outer layer
#'   \item \code{theta_y}: vector of MCMC samples for \code{theta_y} (length 
#'         scale of outer layer)
#'   \item \code{theta_w}: matrix of MCMC samples for \code{theta_w} (length 
#'         scale of middle layer)
#'   \item \code{theta_z}: matrix of MCMC samples for \code{theta_z} (length 
#'         scale of inner layer)
#'   \item \code{w}: list of MCMC samples for middle hidden layer \code{w}
#'   \item \code{z}: list of MCMC samples for inner hidden layer \code{z}
#'   \item \code{w_approx}: Vecchia approximation object for outer layer (\code{vecchia = TRUE} only)
#'   \item \code{z_approx}: Vecchia approximation object for middle layer (\code{vecchia = TRUE} only)
#'   \item \code{x_approx}: Vecchia approximation object for inner layer (\code{vecchia = TRUE} only)
#'   \item \code{ll}: vector of MVN log likelihood of the outer layer 
#'         for reach Gibbs iteration
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Sauer, A. (2023). Deep Gaussian process surrogates for computer experiments. 
#'      *Ph.D. Dissertation, Department of Statistics, Virginia Polytechnic Institute and State University.*
#'      \url{http://hdl.handle.net/10919/114845}
#'      \cr\cr
#' Sauer, A., Gramacy, R.B., & Higdon, D. (2023). Active learning for deep 
#'      Gaussian process surrogates. *Technometrics, 65,* 4-18.  arXiv:2012.08015
#'      \cr\cr
#' Sauer, A., Cooper, A., & Gramacy, R. B. (2023). Vecchia-approximated deep Gaussian 
#'      processes for computer experiments. 
#'      *Journal of Computational and Graphical Statistics, 32*(3), 824-837.  arXiv:2204.02904
#'      \cr\cr
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
#' # Example 1: nugget fixed, calculating entropy
#' fit <- fit_three_layer(x, y, nmcmc = 2000, true_g = 1e-6)
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, entropy_limit = -0.5, cores = 1)
#' plot(fit)
#' i <- interp::interp(xx[, 1], xx[, 2], fit$entropy)
#' image(i, col = heat.colors(128), main = "Entropy")
#'
#' # Example 2: using Vecchia
#' fit <- fit_three_layer(x, y, nmcmc = 2000, true_g = 1e-6, vecchia = TRUE, m = 10)
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#' }
#' 
#' @export

fit_three_layer <- function(x, y, nmcmc = 10000, 
                            D = ifelse(is.matrix(x), ncol(x), 1), 
                            verb = TRUE, w_0 = NULL, z_0 = NULL, theta_y_0 = 0.01, 
                            theta_w_0 = 0.1, theta_z_0 = 0.1, g_0 = 0.001, 
                            true_g = NULL, v = 2.5, settings = NULL,
                            cov = c("matern", "exp2"), vecchia = FALSE, 
                            m = NULL, ord = NULL, cores = NULL) {
  
  tic <- proc.time()[[3]]
  cov <- match.arg(cov)
  if (is.vector(x)) x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  
  # Check inputs and settings
  test <- check_inputs(x, y, true_g, nmcmc) # returns NULL if all checks pass
  settings <- check_settings(settings, layers = 3, noisy = is.null(true_g))
  settings$pmx <- FALSE
  settings$monowarp <- FALSE
  
  # Check covariance
  if (cov == "matern") {
    if(!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  } else if (cov == "exp2") {
    v <- 999 # indicator for "exp2" kernel
  } else stop("cov must be 'matern' or 'exp2'")
  
  # Check vecchia
  if (vecchia) {
    cores <- check_cores(cores)
    if (is.null(m)) m <- min(25, n - 1)
    test <- check_vecchia(n, d, m, ord)
  } else { 
    if (n > 200) message("'vecchia = TRUE' is recommended for faster computation.")
    if (!is.null(cores)) message("cores is only used when 'vecchia = TRUE'")
  }  

  # Create and check initial list
  initial <- list(w = w_0, z = z_0, theta_y = theta_y_0, 
                  theta_w = theta_w_0, theta_z = theta_z_0, g = g_0)
  initial <- check_initialization(initial, layers = 3, D = D, x = x, v = v,
                                  vecchia = vecchia, m = m)
  
  # Create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, v = v)
  
  # Conduct MCMC
  if (vecchia) {
    samples <- gibbs_three_layer_vec(x, y, nmcmc, D, verb, initial, true_g,
                                     settings, v, m, ord, cores)
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
