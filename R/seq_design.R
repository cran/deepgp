
# Function Contents -----------------------------------------------------------
# Internal:
#   alc.C: calculates ALC using C
#   Wij.C: calculates Wij matrix using C
#   exp_improv: calculates expected improvement
# External (see documentation below):
#   ALC (S3 method for gp, dgp2, dgp3 classes)
#   IMSE (S3 method for gp, dgp2, dgp3 classes)

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
               alc = double(nrow(Xcand)))
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
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." \emph{Technometrics, to appear;} arXiv:2012.08015. 
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

  # Prepare parallel clusters
  iters <- 1:object$nmcmc
  if (cores == 1) {
    chunks <- list(iters)
  } else chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
  if (cores > detectCores()) warning("cores is greater than available nodes")
  cl <- makeCluster(cores)
  registerDoParallel(cl)
    
  thread <- NULL
  alc <- foreach(thread = 1:cores, .combine = "+") %dopar% {
    alc_sum <- rep(0, times = n_new)
      
    for (t in chunks[[thread]]) {
      K <- Exp2Fun(dx, c(1, object$theta[t], object$g[t]))
      Ki <- invdet(K)$Mi
      alc_sum <- alc_sum + alc.C(object$x, Ki, object$theta[t], object$g[t], 
                                 x_new, ref, object$tau2[t])
    } # end of t for loop
    return(alc_sum)
  } # end of foreach statement
    
  stopCluster(cl)
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
  
  # Prepare parallel clusters
  iters <- 1:object$nmcmc
  if (cores == 1) {
    chunks <- list(iters)
  } else chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
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
      
      K <- Exp2Fun(sq_dist(w), c(1, object$theta_y[t], object$g[t]))
      Ki <- invdet(K)$Mi
      alc_sum <- alc_sum + alc.C(w, Ki, object$theta_y[t], object$g[t], w_new, 
                                 ref, object$tau2[t])
    } # end of t for loop
    return(alc_sum)
  } # end of foreach statement

  stopCluster(cl)
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

  # Prepare parallel clusters
  iters <- 1:object$nmcmc
  if (cores == 1) {
    chunks <- list(iters)
  } else chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
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
      
      K <- Exp2Fun(sq_dist(w), c(1, object$theta_y[t], object$g[t]))
      Ki <- invdet(K)$Mi
      alc_sum <- alc_sum + alc.C(w, Ki, object$theta_y[t], object$g[t], w_new, 
                                 ref, object$tau2[t])
    } # end of t for loop
    return(alc_sum)
  } # end of foreach statement
    
  stopCluster(cl)
  toc <- proc.time()[3]
  return(list(value = alc / object$nmcmc, time = toc - tic))
}

# Wij C -----------------------------------------------------------------------

Wij.C <- function(x1, x2, theta, a, b){
  W <- matrix(1, nrow = nrow(x1), ncol = nrow(x2))
  result <- .C("Wij_R",
               X1 = as.double(t(x1)),
               n1 = as.integer(nrow(x1)),
               X2 = as.double(t(x2)),
               n2 = as.integer(nrow(x2)),
               col = as.integer(ncol(x1)),
               theta = as.double(theta),
               a = as.double(a),
               b = as.double(b),
               W = as.double(t(W)))
  return(matrix(result$W, nrow = nrow(x1), ncol = nrow(x2)))
}

# Define IMSE for S3 Objects -------------------------------------------------
#' @title Integrated Mean-Squared (prediction) Error for Sequential Design
#' @description Acts on a \code{gp}, \code{dgp2}, or \code{dgp3} object.
#'     Current version requires squared exponential covariance
#'     (\code{cov = "exp2"}).  Calculates IMSE over the input locations 
#'     \code{x_new}.  Optionally utilizes SNOW parallelization.  User should 
#'     select the point with the lowest IMSE to add to the design.
#'     
#' @details Not yet implemented for Vecchia-approximated fits.
#' 
#'     All iterations in the object are used in the calculation, so samples
#'     should be burned-in.  Thinning the samples using \code{trim} will speed 
#'     up computation.  This function may be used in two ways:
#'     \itemize{
#'         \item Option 1: called on an object with only MCMC iterations, in 
#'         which case \code{x_new} must be specified
#'         \item Option 2: called on an object that has been predicted over, in 
#'         which case the \code{x_new} from \code{predict} is used
#'     }
#'     In Option 2, it is recommended to set \code{store_latent = TRUE} for 
#'     \code{dgp2} and \code{dgp3} objects so latent mappings do not have to 
#'     be re-calculated.  Through \code{predict}, the user may
#'     specify a mean mapping (\code{mean_map = TRUE}) or a full sample from 
#'     the MVN distribution over \code{w_new} (\code{mean_map = FALSE}).  When 
#'     the object has not yet been predicted over (Option 1), the mean mapping 
#'     is used.
#'     
#'     SNOW parallelization reduces computation time but requires more memory storage.
#' 
#' @param object object of class \code{gp}, \code{dgp2}, or \code{dgp3}
#' @param x_new matrix of possible input locations, if object has been run 
#'        through \code{predict} the previously stored \code{x_new} is used
#' @param cores number of cores to utilize in parallel, by default no 
#'        parallelization is used
#' @return list with elements:
#' \itemize{
#'   \item \code{value}: vector of IMSE values, indices correspond to \code{x_new}
#'   \item \code{time}: computation time in seconds
#' }
#' 
#' @references 
#' Sauer, A, RB Gramacy, and D Higdon. 2020. "Active Learning for Deep Gaussian 
#'     Process Surrogates." \emph{Technometrics, to appear:} arXiv:2012.08015. 
#'     \cr\cr
#' Binois, M, J Huang, RB Gramacy, and M Ludkovski. 2019. "Replication or Exploration? 
#'     Sequential Design for Stochastic Simulation Experiments." \emph{Technometrics 
#'     61}, 7-23. Taylor & Francis. doi:10.1080/00401706.2018.1469433.
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
#' # Fit model and calculate IMSE
#' fit <- fit_one_layer(x, y, nmcmc = 100, cov = "exp2")
#' fit <- trim(fit, 50)
#' fit <- predict(fit, x_new, cores = 1, store_latent = TRUE)
#' imse <- IMSE(fit)
#' 
#' \donttest{
#' # --------------------------------------------------------
#' # Example 2: Higdon function
#' # --------------------------------------------------------
#' 
#' f <- function(x) {
#'     i <- which(x <= 0.48)
#'     x[i] <- 2 * sin(pi * x[i] * 4) + 0.4 * cos(pi * x[i] * 16)
#'     x[-i] <- 2 * x[-i] - 1
#'     return(x)
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
#' # Conduct MCMC (can replace fit_three_layer with fit_one_layer/fit_two_layer)
#' fit <- fit_three_layer(x, y, D = 1, nmcmc = 2000, cov = "exp2")
#' plot(fit)
#' fit <- trim(fit, 1000, 2)
#' 
#' # Option 1 - calculate IMSE from only MCMC iterations
#' imse <- IMSE(fit, xx)
#' 
#' # Option 2 - calculate IMSE after predictions
#' fit <- predict(fit, xx, cores = 1, store_latent = TRUE)
#' imse <- IMSE(fit)
#' 
#' # Visualize fit
#' plot(fit)
#' par(new = TRUE) # overlay IMSE
#' plot(xx, imse$value, col = 2, type = 'l', lty = 2, 
#'      axes = FALSE, xlab = '', ylab = '')
#' 
#' # Select next design point
#' x_new <- xx[which.min(imse$value)]
#' }
#' 
#' @rdname IMSE
#' @export

IMSE <- function(object, x_new, cores)
  UseMethod("IMSE", object)

# IMSE One Layer --------------------------------------------------------------
#' @rdname IMSE
#' @export

IMSE.gp <- function(object, x_new = NULL, cores = 1) {
  
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

  n <- nrow(object$x)
  n_new <- nrow(x_new)
  dx <- sq_dist(object$x)
  
  # Define bounds
  a <- apply(x_new, 2, min)
  b <- apply(x_new, 2, max)
  
  Knew_inv <- matrix(nrow = n + 1, ncol = n + 1)
  Wijs <- matrix(nrow = n + 1, ncol = n + 1)

  # Prepare parallel clusters
  iters <- 1:object$nmcmc
  if (cores == 1) {
    chunks <- list(iters)
  } else chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
  if (cores > detectCores()) warning("cores is greater than available nodes")
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  thread <- NULL  
  imse <- foreach(thread = 1:cores, .combine = '+') %dopar% {
    imse_sum <- rep(0, times = n_new)
     
    for (t in chunks[[thread]]) { 
      Kn <- Exp2Fun(dx, c(1, object$theta[t], object$g[t]))
      Kn_inv <- invdet(Kn)$Mi
      kk <- 1 + object$g[t]
      
      # Precalculate all except the last row and last column
      Wijs[1:n, 1:n] <- Wij.C(object$x, object$x, object$theta[t], a, b)
      
      for (i in 1:n_new) {
        x_star <- matrix(x_new[i, ], nrow = 1)
        
        # Calculate new Ki matrix
        k <- Exp2Fun(sq_dist(object$x, x_star), c(1, object$theta[t], eps))
        v <- c(kk - t(k) %*% Kn_inv %*% k)
        g <- (- 1 / v) * Kn_inv %*% k
        Knew_inv[1:n, 1:n] <- Kn_inv + g %*% t(g) * v
        Knew_inv[1:n, n+1] <- g
        Knew_inv[n+1, 1:n] <- g
        Knew_inv[n+1, n+1] <- 1 / v
        
        Wijs[1:n, n+1] <- Wijs[n+1, 1:n] <- Wij.C(object$x, x_star, 
                                                  object$theta[t], a, b)
        Wijs[n+1, n+1] <- Wij.C(x_star, x_star, object$theta[t], a, b)
        imse_sum[i] <- imse_sum[i] + object$tau2[t] * prod(b - a) * 
                                      (1 - sum(Knew_inv * Wijs))
        # Note: sum(Ki * Wijs) == sum(diag(Ki %*% Wijs)) because symmetric
      } # end of i for loop
    } # end of t for loop
    return(imse_sum)
  } # end of foreach statement
    
  stopCluster(cl)
  toc <- proc.time()[3]
  return(list(value = imse / object$nmcmc, time = toc - tic))
}

# IMSE Two Layer --------------------------------------------------------------
#' @rdname IMSE
#' @export

IMSE.dgp2 <- function(object, x_new = NULL, cores = 1) {
  
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
  
  n <- nrow(object$x)
  n_new <- nrow(x_new)
  if (!predicted) {
    D <- ncol(object$w[[1]])
    dx <- sq_dist(object$x)
    d_cross <- sq_dist(x_new, object$x)
  }
  
  Knew_inv <- matrix(nrow = n + 1, ncol = n + 1)
  Wijs <- matrix(nrow = n + 1, ncol = n + 1)
  
  # Prepare parallel clusters
  iters <- 1:object$nmcmc
  if (cores == 1) {
    chunks <- list(iters)
  } else chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
  if (cores > detectCores()) warning("cores is greater than available nodes")
  cl <- makeCluster(cores)
  registerDoParallel(cl)
    
  thread <- NULL
  imse <- foreach(thread = 1:cores, .combine = '+') %dopar% {
    imse_sum <- rep(0, times = n_new)
    
    for (t in chunks[[thread]]) {
      w <- object$w[[t]]
      Kn <- Exp2Fun(sq_dist(w), c(1, object$theta_y[t], object$g[t]))
      Kn_inv <- invdet(Kn)$Mi
      kk <- 1 + object$g[t]
      
      if (predicted) {
        w_new <- object$w_new[[t]]
      } else {
        w_new <- matrix(nrow = n_new, ncol = D)
        for (i in 1:D)
          w_new[, i] <- krig(w[, i], dx, NULL, d_cross, object$theta_w[t, i], 
                             g = eps, v = 999)$mean
      }
      
      # Define bounds
      a <- apply(w_new, 2, min)
      b <- apply(w_new, 2, max)
      
      # Precalculate all except the last row and last column
      Wijs[1:n, 1:n] <- Wij.C(w, w, object$theta_y[t], a, b)
      
      for (i in 1:n_new) {
        w_star <- matrix(w_new[i, ], nrow = 1)
        
        # Calculate new Ki matrix
        k <- Exp2Fun(sq_dist(w, w_star), c(1, object$theta_y[t], eps))
        v <- c(kk - t(k) %*% Kn_inv %*% k)
        g <- (- 1 / v) * Kn_inv %*% k
        Knew_inv[1:n, 1:n] <- Kn_inv + g %*% t(g) * v
        Knew_inv[1:n, n+1] <- g
        Knew_inv[n+1, 1:n] <- g
        Knew_inv[n+1, n+1] <- 1 / v
        
        Wijs[1:n, n+1] <- Wijs[n+1, 1:n] <- Wij.C(w, w_star, object$theta_y[t], a, b)
        Wijs[n+1, n+1] <- Wij.C(w_star, w_star, object$theta_y[t], a, b)
        imse_sum[i] <- imse_sum[i] + object$tau2[t] * prod(b - a) * 
                                        (1 - sum(Knew_inv * Wijs))
        # Note: sum(Ki * Wijs) == sum(diag(Ki %*% Wijs)) because symmetric
      } # end of i for loop
    } # end of t for loop
    return(imse_sum)
  } # end of foreach statement
    
  stopCluster(cl)
  toc <- proc.time()[3]
  return(list(value = imse / object$nmcmc, time = toc - tic))
}

# IMSE Three Layer ------------------------------------------------------------
#' @rdname IMSE
#' @export

IMSE.dgp3 <- function(object, x_new = NULL, cores = 1) {
  
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
  
  n <- nrow(object$x)
  n_new <- nrow(x_new)
  if (!predicted) {
    D <- ncol(object$w[[1]])
    dx <- sq_dist(object$x)
    d_cross <- sq_dist(x_new, object$x)
  }
  
  Knew_inv <- matrix(nrow = n + 1, ncol = n + 1)
  Wijs <- matrix(nrow = n + 1, ncol = n + 1)
  
  # Prepare parallel clusters
  iters <- 1:object$nmcmc
  if (cores == 1) {
    chunks <- list(iters)
  } else chunks <- split(iters, sort(cut(iters, cores, labels = FALSE)))
  if (cores > detectCores()) warning("cores is greater than available nodes")
  cl <- makeCluster(cores)
  registerDoParallel(cl)
    
  thread <- NULL
  imse <- foreach(thread = 1:cores, .combine = '+') %dopar% {
    imse_sum <- rep(0, times = n_new)
    
    for (t in chunks[[thread]]) {
      w <- object$w[[t]]
      Kn <- Exp2Fun(sq_dist(w), c(1, object$theta_y[t], object$g[t]))
      Kn_inv <- invdet(Kn)$Mi
      kk <- 1 + object$g[t]
      
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
      
      # Define bounds
      a <- apply(w_new, 2, min)
      b <- apply(w_new, 2, max)
      
      # Precalculate all except the last row and last column
      Wijs[1:n, 1:n] <- Wij.C(w, w, object$theta_y[t], a, b)
      
      for (i in 1:n_new) {
        w_star <- matrix(w_new[i, ], nrow = 1)
        
        # Calculate new Ki matrix
        k <- Exp2Fun(sq_dist(w, w_star), c(1, object$theta_y[t], eps))
        v <- c(kk - t(k) %*% Kn_inv %*% k)
        g <- (- 1 / v) * Kn_inv %*% k
        Knew_inv[1:n, 1:n] <- Kn_inv + g %*% t(g) * v
        Knew_inv[1:n, n+1] <- g
        Knew_inv[n+1, 1:n] <- g
        Knew_inv[n+1, n+1] <- 1 / v
        
        Wijs[1:n, n+1] <- Wijs[n+1, 1:n] <- Wij.C(w, w_star, 
                                                  object$theta_y[t], a, b)
        Wijs[n+1, n+1] <- Wij.C(w_star, w_star, object$theta_y[t], a, b)
        imse_sum[i] <- imse_sum[i] + object$tau2[t] * prod(b - a) * 
                                        (1 - sum(Knew_inv * Wijs))
        # Note: sum(Ki * Wijs) == sum(diag(Ki %*% Wijs)) because symmetric
      } # end of i for loop
    } # end of t for loop
    return(imse_sum)
  } # end of foreach statement
    
  stopCluster(cl)
  toc <- proc.time()[3]
  return(list(value = imse / object$nmcmc, time = toc - tic))
}

# EI --------------------------------------------------------------------------

exp_improv <- function(mu, sig2, f_min) {
  
  ei_store <- rep(0, times = length(mu))
  i <- which(sig2 > eps)
  mu_i <- mu[i]
  sig_i <- sqrt(sig2[i])
  ei <- (f_min - mu_i) * pnorm(f_min, mean = mu_i, sd = sig_i) + 
    sig_i * dnorm(f_min, mean = mu_i, sd = sig_i)
  ei_store[i] <- ei
  
  return(ei_store)
}
