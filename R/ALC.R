
# Function Contents -----------------------------------------------------------
# Internal:
#   alc.C: calculates ALC using C
# External (see documentation below):
#   ALC (S3 method for gp, dgp2, dgp3 classes)

# ALC C -----------------------------------------------------------------------

alc.C <- function(X, Ki, theta, g, Xcand, Xref, tau2, verb = 0) {
  result <- .C("alcGP_R",
               X = as.double(t(X)),
               n = as.integer(nrow(X)),
               col = as.integer(ncol(X)),
               Ki = as.double(t(Ki)),
               d = as.double(theta),
               g = as.double(g),
               ncand = as.integer(nrow(Xcand)),
               Xcand = as.double(t(Xcand)),
               nref = as.integer(nrow(Xref)),
               Xref = as.double(t(Xref)),
               phi = as.double(tau2),
               verb = as.integer(verb),
               alc = double(nrow(Xcand)),
               PACKAGE = "deepgp")
  return(result$alc)
}

# Define ALC for S3 Objects ---------------------------------------------------
#' @title Active Learning Cohn for Sequential Design
#' @description Acts on a \code{gp}, \code{dgp2}, or \code{dgp3} object. 
#'    Current version requires squared exponential covariance 
#'    (\code{cov = "exp2"}).  Calculates ALC over the input locations 
#'    \code{x_new} using specified reference grid.  If no reference grid is 
#'    specified, \code{x_new} is used as the reference.  Optionally utilizes 
#'    SNOW parallelization.  User should 
#'    select the point with the highest ALC to add to the design.   
#'    
#' @details Not yet implemented for Vecchia-approximated fits.
#' 
#'     All iterations in the object are used in the calculation, so samples 
#'     should be burned-in.  Thinning the samples using \code{trim} will
#'     speed up computation.  This function may be used in two ways:
#'     \itemize{
#'         \item Option 1: called on an object with only MCMC iterations, in 
#'               which case \code{x_new} must be specified
#'         \item Option 2: called on an object that has been predicted over, in 
#'               which case the \code{x_new} from \code{predict} is used
#'     }
#'     In Option 2,  it is recommended to set \code{store_latent = TRUE} for 
#'     \code{dgp2} and \code{dgp3} objects so
#'     latent mappings do not have to be re-calculated.  Through \code{predict}, 
#'     the user may specify a mean mapping (\code{mean_map = TRUE}) or a full 
#'     sample from the MVN distribution over \code{w_new} 
#'     (\code{mean_map = FALSE}).  When the object has not yet been predicted
#'     over (Option 1), the mean mapping is used.
#'     
#'     SNOW parallelization reduces computation time but requires more memory
#'     storage.  C code derived from the "laGP" package (Robert B Gramacy and 
#'     Furong Sun).
#' 
#' @param object object of class \code{gp}, \code{dgp2}, or \code{dgp3}
#' @param x_new matrix of possible input locations, if object has been run 
#'        through \code{predict} the previously stored \code{x_new} is used
#' @param ref optional reference grid for ALC approximation, if \code{ref = NULL} 
#'        then \code{x_new} is used
#' @param cores number of cores to utilize in parallel, by default no 
#'        parallelization is used
#' @return list with elements:
#' \itemize{
#'   \item \code{value}: vector of ALC values, indices correspond to \code{x_new}
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Sauer, A., Gramacy, R.B., & Higdon, D. (2023). Active learning for deep 
#'     Gaussian process surrogates. *Technometrics, 65,* 4-18.  arXiv:2012.08015
#'     \cr\cr
#' Seo, S, M Wallat, T Graepel, and K Obermayer. 2000. Gaussian Process Regression:
#'     Active Data Selection and Test Point Rejection. In Mustererkennung 2000, 
#'     2734. New York, NY: SpringerVerlag.\cr\cr
#' Gramacy, RB and F Sun. (2016). laGP: Large-Scale Spatial Modeling via Local 
#'     Approximate Gaussian Processes in R. \emph{Journal of Statistical Software 
#'     72} (1), 1-46. doi:10.18637/jss.v072.i01
#' 
#' @examples
#' # --------------------------------------------------------
#' # Example 1: toy step function, runs in less than 5 seconds
#' # --------------------------------------------------------
#' 
#' f <- function(x) {
#'     if (x <= 0.4) return(-1)
#'     if (x >= 0.6) return(1)
#'     if (x > 0.4 & x < 0.6) return(10*(x-0.5))
#' }
#' 
#' x <- seq(0.05, 0.95, length = 7)
#' y <- sapply(x, f)
#' x_new <- seq(0, 1, length = 100)
#' 
#' # Fit model and calculate ALC
#' fit <- fit_two_layer(x, y, nmcmc = 100, cov = "exp2")
#' fit <- trim(fit, 50)
#' fit <- predict(fit, x_new, cores = 1, store_latent = TRUE)
#' alc <- ALC(fit)
#' 
#' \donttest{
#' # --------------------------------------------------------
#' # Example 2: damped sine wave
#' # --------------------------------------------------------
#' 
#' f <- function(x) {
#'     exp(-10*x) * (cos(10*pi*x - 1) + sin(10*pi*x - 1)) * 5 - 0.2
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
#' plot(xx, yy, type = "l")
#' points(x, y, col = 2)
#' 
#' # Conduct MCMC (can replace fit_two_layer with fit_one_layer/fit_three_layer)
#' fit <- fit_two_layer(x, y, D = 1, nmcmc = 2000, cov = "exp2")
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' 
#' # Option 1 - calculate ALC from MCMC iterations
#' alc <- ALC(fit, xx)
#' 
#' # Option 2 - calculate ALC after predictions
#' fit <- predict(fit, xx, cores = 1, store_latent = TRUE)
#' alc <- ALC(fit)
#' 
#' # Visualize fit
#' plot(fit)
#' par(new = TRUE) # overlay ALC
#' plot(xx, alc$value, type = 'l', lty = 2, 
#'      axes = FALSE, xlab = '', ylab = '')
#' 
#' # Select next design point
#' x_new <- xx[which.max(alc$value)]
#' }
#' 
#' @rdname ALC
#' @export

ALC <- function(object, x_new, ref, cores)
  UseMethod("ALC", object)

# ALC One Layer ---------------------------------------------------------------
#' @rdname ALC
#' @export

ALC.gp <- function(object, x_new = NULL, ref = NULL, cores = 1) {
  
  tic <- proc.time()[3]
  
  if (object$v != 999) stop("Currently, ALC is only implemented for the
                             un-approximated squared exponential kernel.  
                             Re-fit model with 'vecchia = FALSE' and  
                             cov = 'exp2' in order to use ALC.")
  
  if (is.null(x_new)) {
    if (is.null(object$x_new)) {
      stop("x_new has not been specified")
    } else x_new <- object$x_new
  }
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  if (!is.null(ref) & !is.matrix(ref)) stop("ref must be a matrix")
  if (is.null(ref)) ref <- x_new
  
  dx <- sq_dist(object$x)
  n_new <- nrow(x_new)
  
  if (cores == 1) { # run serial for loop 
    
    alc <- rep(0, times = n_new)
    for (t in 1:object$nmcmc) {
      K <- Exp2(dx, 1, object$theta[t], object$g[t])
      Ki <- invdet(K)$Mi
      alc <- alc + alc.C(object$x, Ki, object$theta[t], object$g[t], 
                         x_new, ref, object$tau2[t])
    } # end of t for loop

  } else { # run in parallel using foreach

    iters <- 1:object$nmcmc
    chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
    if (cores > detectCores()) warning("cores is greater than available nodes")
    
    cl <- makeCluster(cores)
    registerDoParallel(cl)
      
    thread <- NULL
    alc <- foreach(thread = 1:cores, .combine = "+") %dopar% {
      alc_sum <- rep(0, times = n_new)
        
      for (t in chunks[[thread]]) {
        K <- Exp2(dx, 1, object$theta[t], object$g[t])
        Ki <- invdet(K)$Mi
        alc_sum <- alc_sum + alc.C(object$x, Ki, object$theta[t], object$g[t], 
                                   x_new, ref, object$tau2[t])
      } # end of t for loop
      return(alc_sum)
    } # end of foreach statement
    
    stopCluster(cl)
    
  } # end of else statement
  
  toc <- proc.time()[3]
  return(list(value = alc / object$nmcmc, time = toc - tic))
}

# ALC Two Layer Function ------------------------------------------------------
#' @rdname ALC
#' @export

ALC.dgp2 <- function(object, x_new = NULL, ref = NULL, cores = 1) {
  
  tic <- proc.time()[3]
  
  if (object$v != 999) stop("Currently, ALC is only implemented for the
                             un-approximated squared exponential kernel.  
                             Re-fit model with 'vecchia = FALSE' and  
                             cov = 'exp2' in order to use ALC.")
  
  if (is.null(x_new)) {
    if (is.null(object$x_new)) {
      stop("x_new has not been specified")
    } else {
      x_new <- object$x_new
      if (is.null(object$w_new)) {
        predicted <- FALSE 
        message("next time, use store_latent = TRUE inside prediction to 
                speed up computation")
      } else predicted <- TRUE
    }
  } else predicted <- FALSE
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  if (!is.null(ref) & !is.matrix(ref)) stop('ref must be a matrix')
  
  # Specify pre-calculations if predicted is FALSE
  n_new <- nrow(x_new)
  if (!predicted) {
    D <- ncol(object$w[[1]])
    dx <- sq_dist(object$x)
    d_cross <- sq_dist(x_new, object$x)
  }
  
  if (cores == 1) { # run serial for loop 
    
    alc <- rep(0, times = n_new)
    for (t in 1:object$nmcmc) {
      w <- object$w[[t]]
        
      if (predicted) {
        w_new <- object$w_new[[t]]
      } else {
        w_new <- matrix(nrow = n_new, ncol = D)
        for (i in 1:D)
          w_new[, i] <- krig(w[, i], dx, NULL, d_cross, object$theta_w[t, i], 
                             g = eps, v = 999)$mean
      } 
        
      if (is.null(ref)) ref <- w_new
        
      K <- Exp2(sq_dist(w), 1, object$theta_y[t], object$g[t])
      Ki <- invdet(K)$Mi
      alc <- alc + alc.C(w, Ki, object$theta_y[t], object$g[t], w_new, 
                         ref, object$tau2[t])
    } # end of t for loop
    
  } else { # run in parallel using foreach 
    
    message("WARNING - recommend cores = 1.  Odd behavior noticed when cores > 1")
    
    iters <- 1:object$nmcmc
    chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
    if (cores > detectCores()) warning("cores is greater than available nodes")
    
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    
    thread <- NULL
    alc <- foreach(thread = 1:cores, .combine = "+") %dopar% {
      alc_sum <- rep(0, times = n_new)
      
      for (t in chunks[[thread]]) {
        w <- object$w[[t]]
        
        if (predicted) {
          w_new <- object$w_new[[t]]
        } else {
          w_new <- matrix(nrow = n_new, ncol = D)
          for (i in 1:D)
            w_new[, i] <- krig(w[, i], dx, NULL, d_cross, object$theta_w[t, i], 
                               g = eps, v = 999)$mean
        } 
        
        if (is.null(ref)) ref <- w_new
        
        K <- Exp2(sq_dist(w), 1, object$theta_y[t], object$g[t])
        Ki <- invdet(K)$Mi
        alc_sum <- alc_sum + alc.C(w, Ki, object$theta_y[t], object$g[t], w_new, 
                                   ref, object$tau2[t])
      } # end of t for loop
      return(alc_sum)
    } # end of foreach statement
    
    stopCluster(cl)
    
  } # end of else statement
  
  toc <- proc.time()[3]
  return(list(value = alc / object$nmcmc, time = toc - tic))
}

# ALC Three Layer Function ----------------------------------------------------
#' @rdname ALC
#' @export

ALC.dgp3 <- function(object, x_new = NULL, ref = NULL, cores = 1) {
  
  tic <- proc.time()[3]
  
  if (object$v != 999) stop("Currently, ALC is only implemented for the
                             un-approximated squared exponential kernel.  
                             Re-fit model with 'vecchia = FALSE' and  
                             cov = 'exp2' in order to use ALC.")
  
  if (is.null(x_new)) {
    if (is.null(object$x_new)) {
      stop("x_new has not been specified")
    } else {
      x_new <- object$x_new
      if (is.null(object$w_new)) {
        predicted <- FALSE 
        message("next time, use store_latent = TRUE inside prediction to 
                speed up computation")
      } else predicted <- TRUE
    } 
  } else predicted <- FALSE
  if (is.numeric(x_new)) x_new <- as.matrix(x_new)
  if (!is.null(ref) & !is.matrix(ref)) stop('ref must be a matrix')
  
  # Specify pre-calculations if predicted is FALSE
  n_new <- nrow(x_new)
  if (!predicted) {
    D <- ncol(object$w[[1]])
    dx <- sq_dist(object$x)
    d_cross <- sq_dist(x_new, object$x)
  }
  
  if (cores == 1) { # run serial for loop
    
    alc <- rep(0, times = n_new)
    for (t in 1:object$nmcmc) {
      w <- object$w[[t]]
        
      if (predicted) {
        w_new <- object$w_new[[t]]
      } else {
        z <- object$z[[t]]
        z_new <- matrix(nrow = n_new, ncol = D)
        for (i in 1:D)
          z_new[, i] <- krig(z[, i], dx, NULL, d_cross, object$theta_z[t, i], 
                             g = eps, v = 999)$mean
        w_new <- matrix(nrow = n_new, ncol = D)
        for (i in 1:D) 
          w_new[, i] <- krig(w[, i], sq_dist(z), NULL, sq_dist(z_new, z), 
                             object$theta_w[t, i], g = eps, v = 999)$mean
      } 
        
      if (is.null(ref)) ref <- w_new
        
      K <- Exp2(sq_dist(w), 1, object$theta_y[t], object$g[t])
      Ki <- invdet(K)$Mi
      alc <- alc + alc.C(w, Ki, object$theta_y[t], object$g[t], w_new, 
                         ref, object$tau2[t])
    } # end of t for loop
     
  } else { # run in parallel using foreach
    
    message("WARNING - recommend cores = 1.  Odd behavior noticed when cores > 1")
    
    iters <- 1:object$nmcmc
    chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
    if (cores > detectCores()) warning("cores is greater than available nodes")
    
    cl <- makeCluster(cores)
    registerDoParallel(cl)
     
    thread <- NULL
    alc <- foreach(thread = 1:cores, .combine = '+') %dopar% {
      alc_sum <- rep(0, times = n_new)
      
      for (t in chunks[[thread]]) {
        w <- object$w[[t]]
        
        if (predicted) {
          w_new <- object$w_new[[t]]
        } else {
          z <- object$z[[t]]
          z_new <- matrix(nrow = n_new, ncol = D)
          for (i in 1:D)
            z_new[, i] <- krig(z[, i], dx, NULL, d_cross, object$theta_z[t, i], 
                               g = eps, v = 999)$mean
          w_new <- matrix(nrow = n_new, ncol = D)
          for (i in 1:D) 
            w_new[, i] <- krig(w[, i], sq_dist(z), NULL, sq_dist(z_new, z), 
                               object$theta_w[t, i], g = eps, v = 999)$mean
        } 
        
        if (is.null(ref)) ref <- w_new
        
        K <- Exp2(sq_dist(w), 1, object$theta_y[t], object$g[t])
        Ki <- invdet(K)$Mi
        alc_sum <- alc_sum + alc.C(w, Ki, object$theta_y[t], object$g[t], w_new, 
                                   ref, object$tau2[t])
      } # end of t for loop
      return(alc_sum)
    } # end of foreach statement
    
    stopCluster(cl)
  
  } # end of else statement
  
  toc <- proc.time()[3]
  return(list(value = alc / object$nmcmc, time = toc - tic))
}
