
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
#'     controlled within the \code{settings} list object.  Defaults are
#'     \itemize{
#'         \item \code{settings$alpha$g <- 1.5}
#'         \item \code{settings$beta$g <- 3.9}
#'         \item \code{settings$alpha$theta <- 1.5}
#'         \item \code{settings$beta$theta <- 3.9 / 1.5}
#'     }
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
#' Sauer, A., Cooper, A., & Gramacy, R. B. (2022). Vecchia-approximated deep 
#'     Gaussian processes for computer experiments. 
#'     *Journal of Computational and Graphical Statistics,* 1-14.  arXiv:2204.02904
#' 
#' @examples 
#' # Examples of real-world implementations are available at: 
#' # https://bitbucket.org/gramacylab/deepgp-ex/
#' \donttest{
#' # G function (https://www.sfu.ca/~ssurjano/gfunc.html)
#' f <- function(xx, a = (c(1:length(xx)) - 1) / 2) { 
#'     new1 <- abs(4 * xx - 2) + a
#'     new2 <- 1 + a
#'     prod <- prod(new1 / new2)
#'     return((prod - 1) / 0.86)
#' }
#' 
#' # Training data
#' d <- 1 
#' n <- 20
#' x <- matrix(runif(n * d), ncol = d)
#' y <- apply(x, 1, f)
#' 
#' # Testing data
#' n_test <- 100
#' xx <- matrix(runif(n_test * d), ncol = d)
#' yy <- apply(xx, 1, f)
#' 
#' plot(xx[order(xx)], yy[order(xx)], type = "l")
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
#'       axes = FALSE, xlab = '', ylab = '')
#'       
#' # Example 3: Vecchia approximated model
#' fit <- fit_one_layer(x, y, nmcmc = 2000, vecchia = TRUE, m = 10) 
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
#' }
#' 
#' @export

fit_one_layer <- function(x, y, nmcmc = 10000, sep = FALSE, verb = TRUE, g_0 = 0.01, 
                          theta_0 = 0.1, true_g = NULL, settings = NULL,
                          cov = c("matern", "exp2"), v = 2.5, 
                          vecchia = FALSE, m = min(25, length(y) - 1)) {
  
  tic <- proc.time()[3]
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
  settings <- check_settings(settings, layers = 1)
  initial <- list(theta = theta_0, g = g_0, tau2 = 1)
  if (m >= length(y)) stop("m must be less than the length of y")
  if (cov == "matern")
    if(!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  
  # Create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, v = v)
  if (vecchia) out$m <- m

  # Conduct MCMC
  if (sep) {
    if (vecchia) {
      samples <- gibbs_one_layer_vec_sep(x, y, nmcmc, verb, initial, true_g,
                                         settings, v, m)
    } else{
      samples <- gibbs_one_layer_sep(x, y, nmcmc, verb, initial, true_g,
                                     settings, v)
    }
  } else {
    if (vecchia) {
      samples <- gibbs_one_layer_vec(x, y, nmcmc, verb, initial, true_g, 
                                     settings, v, m)
    } else { 
      samples <- gibbs_one_layer(x, y, nmcmc, verb, initial, true_g,
                                 settings, v)
    }
  }
  
  out <- c(out, samples)
  toc <- proc.time()[3]
  out$time <- toc - tic
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
#'     Defaults are
#'     \itemize{
#'         \item \code{settings$alpha$g <- 1.5}
#'         \item \code{settings$beta$g <- 3.9}
#'         \item \code{settings$alpha$theta_w <- 1.5}
#'         \item \code{settings$beta$theta_w <- 3.9 / 4}
#'         \item \code{settings$alpha$theta_y <- 1.5}
#'         \item \code{settings$beta$theta_y <- 3.9 / 6}
#'     }
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
#' @param pmx "prior mean X", logical indicating whether W should have prior
#'        mean of X (\code{TRUE}, requires \code{D = ncol(X)}) or prior 
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
#' Sauer, A., Cooper, A., & Gramacy, R. B. (2022). Vecchia-approximated deep 
#'     Gaussian processes for computer experiments. 
#'     *Journal of Computational and Graphical Statistics,* 1-14.  arXiv:2204.02904
#' 
#' @examples 
#' # Examples of real-world implementations are available at: 
#' # https://bitbucket.org/gramacylab/deepgp-ex/
#' \donttest{
#' # G function (https://www.sfu.ca/~ssurjano/gfunc.html)
#' f <- function(xx, a = (c(1:length(xx)) - 1) / 2) { 
#'     new1 <- abs(4 * xx - 2) + a
#'     new2 <- 1 + a
#'     prod <- prod(new1 / new2)
#'     return((prod - 1) / 0.86)
#' }
#' 
#' # Training data
#' d <- 1 
#' n <- 20
#' x <- matrix(runif(n * d), ncol = d)
#' y <- apply(x, 1, f)
#' 
#' # Testing data
#' n_test <- 100
#' xx <- matrix(runif(n_test * d), ncol = d)
#' yy <- apply(xx, 1, f)
#' 
#' plot(xx[order(xx)], yy[order(xx)], type = "l")
#' points(x, y, col = 2)
#' 
#' # Example 1: full model (nugget estimated, using continue)
#' fit <- fit_two_layer(x, y, nmcmc = 1000)
#' plot(fit)
#' fit <- continue(fit, 1000) 
#' plot(fit) 
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit, hidden = TRUE)
#' 
#' # Example 2: Vecchia approximated model
#' # (Vecchia approximation is faster for larger data sizes)
#' fit <- fit_two_layer(x, y, nmcmc = 2000, vecchia = TRUE, m = 10)
#' plot(fit) 
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit, hidden = TRUE)
#' 
#' # Example 3: Vecchia approximated model (re-approximated after burn-in)
#' fit <- fit_two_layer(x, y, nmcmc = 1000, vecchia = TRUE, m = 10)
#' fit <- continue(fit, 1000, re_approx = TRUE)
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit, hidden = TRUE)
#' }
#' 
#' @export

fit_two_layer <- function(x, y, nmcmc = 10000, D = ifelse(is.matrix(x), ncol(x), 1), 
                          pmx = FALSE, verb = TRUE, w_0 = NULL, g_0 = 0.01,
                          theta_y_0 = 0.1, theta_w_0 = 0.1, true_g = NULL,
                          settings = NULL, cov = c("matern", "exp2"), v = 2.5,
                          vecchia = FALSE, m = min(25, length(y) - 1)) {

  tic <- proc.time()[3]
  cov <- match.arg(cov)
  if (vecchia) check_omp()
  if (cov == "exp2") v <- 999 # solely used as an indicator
  if (!vecchia & length(y) > 300) 
    message("'vecchia = TRUE' is recommended for faster computation.")
  if (nmcmc <= 1) stop("nmcmc must be greater than 1")

  # Check inputs
  if (is.numeric(x)) x <- as.matrix(x)
  test <- check_inputs(x, y, true_g) # returns NULL if all checks pass
  settings <- check_settings(settings, layers = 2, D, length(y))
  initial <- list(w = w_0, theta_y = theta_y_0, theta_w = theta_w_0, 
                  g = g_0, tau2 = 1)
  initial <- check_initialization(initial, layers = 2, x = x, D = D, 
                                  vecchia = vecchia, v = v, m = m)
  if (m >= length(y)) stop("m must be less than the length of y")
  if (cov == "matern")
    if(!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  
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
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, v = v)
  if (vecchia) out$m <- m

  # Conduct MCMC
  if (vecchia) {
    samples <- gibbs_two_layer_vec(x, y, nmcmc, D, verb, initial,
                                   true_g, settings, v, m)
  } else { 
    samples <- gibbs_two_layer(x, y, nmcmc, D, verb, initial,
                               true_g, settings, v)
  } 
  
  out <- c(out, samples)
  toc <- proc.time()[3]
  out$time <- toc - tic
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
#'     object.  Defaults are 
#'     \itemize{
#'         \item \code{settings$alpha$g <- 1.5}
#'         \item \code{settings$beta$g <- 3.9}
#'         \item \code{settings$alpha$theta_z <- 1.5}
#'         \item \code{settings$beta$theta_z <- 3.9 / 4}
#'         \item \code{settings$alpha$theta_w <- 1.5}
#'         \item \code{settings$beta$theta_w <- 3.9 / 12}
#'         \item \code{settings$alpha$theta_y <- 1.5}
#'         \item \code{settings$beta$theta_y <- 3.9 / 6}
#'     }
#'     These priors are designed for \code{x} scaled to [0, 1] and \code{y} 
#'     scaled to have mean 0 and variance 1.  These may be adjusted using the 
#'     \code{settings} input.
#'     
#'     In the current version, the three-layer does not have any equivalent
#'     setting for \code{pmx = TRUE} as in \code{fit_two_layer}.
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
#' Sauer, A., Cooper, A., & Gramacy, R. B. (2022). Vecchia-approximated deep 
#'     Gaussian processes for computer experiments. 
#'     *Journal of Computational and Graphical Statistics,* 1-14.  arXiv:2204.02904
#' 
#' @examples 
#' # Examples of real-world implementations are available at: 
#' # https://bitbucket.org/gramacylab/deepgp-ex/
#' \donttest{
#' # G function (https://www.sfu.ca/~ssurjano/gfunc.html)
#' f <- function(xx, a = (c(1:length(xx)) - 1) / 2) { 
#'     new1 <- abs(4 * xx - 2) + a
#'     new2 <- 1 + a
#'     prod <- prod(new1 / new2)
#'     return((prod - 1) / 0.86)
#' }
#' 
#' # Training data
#' d <- 2
#' n <- 30
#' x <- matrix(runif(n * d), ncol = d)
#' y <- apply(x, 1, f)
#' 
#' # Testing data
#' n_test <- 100
#' xx <- matrix(runif(n_test * d), ncol = d)
#' yy <- apply(xx, 1, f)
#' 
#' i <- interp::interp(xx[, 1], xx[, 2], yy)
#' image(i, col = heat.colors(128))
#' contour(i, add = TRUE)
#' points(x)
#' 
#' # Example 1: full model (nugget estimated)
#' fit <- fit_three_layer(x, y, nmcmc = 2000)
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' fit <- predict(fit, xx, cores = 1)
#' plot(fit)
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
                            g_0 = 0.01, theta_y_0 = 0.1, theta_w_0 = 0.1, 
                            theta_z_0 = 0.1, true_g = NULL, settings = NULL,
                            cov = c("matern", "exp2"), v = 2.5, vecchia = FALSE, 
                            m = min(25, length(y) - 1)) {
  
  tic <- proc.time()[3]
  cov <- match.arg(cov)
  if (vecchia) check_omp()
  if (cov == "exp2") v <- 999 # solely used as an indicator
  if (!vecchia & length(y) > 300) 
    message("'vecchia = TRUE' is recommended for faster computation.")
  if (nmcmc <= 1) stop("nmcmc must be greater than 1")

  # Check inputs
  if (is.numeric(x)) x <- as.matrix(x)
  test <- check_inputs(x, y, true_g)
  settings <- check_settings(settings, layers = 3, D, length(y))
  initial <- list(w = w_0, z = z_0, theta_y = theta_y_0, theta_w = theta_w_0, 
                  theta_z = theta_z_0, g = g_0, tau2 = 1)
  initial <- check_initialization(initial, layers = 3, x = x, D = D, 
                                  vecchia = vecchia, v = v, m = m)
  if (m >= length(y)) stop("m must be less than the length of y")
  if (cov == "matern")
    if(!(v %in% c(0.5, 1.5, 2.5))) 
      stop("v must be one of 0.5, 1.5, or 2.5")
  
  # Create output object
  out <- list(x = x, y = y, nmcmc = nmcmc, settings = settings, v = v)
  if (vecchia) out$m <- m
  
  # Conduct MCMC
  if (vecchia) {
    samples <- gibbs_three_layer_vec(x, y, nmcmc, D, verb, initial, true_g,
                                     settings, v, m)
  } else {
    samples <- gibbs_three_layer(x, y, nmcmc, D, verb, initial, true_g, 
                                 settings, v)
  }
  
  out <- c(out, samples)
  toc <- proc.time()[3]
  out$time <- toc - tic
  if (vecchia) class(out) <- "dgp3vec" else class(out) <- "dgp3"
  return(out)
}
